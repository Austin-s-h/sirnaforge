"""Thermodynamic calculations for siRNA design using ViennaRNA."""

from dataclasses import dataclass

import RNA
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

from sirnaforge.models.sirna import (
    DEFAULT_MIN_ASYMMETRY_SCORE,
    DEFAULT_PLFOLD_MAX_BP_SPAN,
    DEFAULT_PLFOLD_WINDOW,
    SiRNACandidate,
)
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

# Number of terminal base pairs compared when scoring duplex-end asymmetry.
# Applied symmetrically to both ends so the two ΔG values are comparable.
END_WINDOW_NT = 7

# Target-site opening windows, all anchored on the site's 3' end.
#
# Guide and target are antiparallel: for a target site T[1..L] 5'->3' and a guide G[1..L] 5'->3',
# G[i] pairs T[L+1-i]. The guide seed G[2..8] therefore pairs the *3' end* of the target site, so
# the region RISC has to open first is the 3'-most nucleotides -- not the 5'-most. Taking the
# 8-mer at the other end measures accessibility opposite the guide's 3' end and is a measured
# near-null against knockdown (Spearman rho +0.07 vs +0.27; see issue #95).
SEED_END_WINDOW_NT = 8
SEED_ANCHORED_WINDOW_NT = 17

# Conditions for the nearest-neighbour melting temperature: roughly physiological
# ionic strength and a typical transfection-scale duplex concentration.
TM_SODIUM_MM = 100.0
TM_STRAND_CONC_NM = 100.0


@dataclass(frozen=True)
class TargetSiteAccessibility:
    """Local opening probabilities for one siRNA target site on its mRNA.

    Every window is anchored on the site's 3' end, which is the end the guide seed pairs
    (see SEED_END_WINDOW_NT). A value is None when the site is too close to the transcript
    5' end for the window to fit; it is never substituted with a default, because a missing
    input must not score as a good one.

    Attributes:
        seed_end_8mer: P(the 3'-most 8 nt of the site are unpaired). The scored statistic.
        seed_anchored_17mer: P(the 3'-most 17 nt are unpaired). Reported, not scored.
        whole_site: P(the entire site is unpaired). Reported, not scored.
    """

    seed_end_8mer: float | None
    seed_anchored_17mer: float | None
    whole_site: float | None


class TargetAccessibilityProfile:
    """RNAplfold unpaired-probability profile for one transcript, folded once.

    Wraps the matrix returned by ``RNA.pfl_fold_up``, which is 1-based in both axes:
    ``unpaired[i][u]`` is the equilibrium probability that the length-``u`` stretch
    **ending at** 1-based position ``i`` is entirely unpaired. (Confirmed against the
    library rather than its docstring: ``unpaired[1][u]`` is 0 for every u > 1 while
    ``unpaired[N][u]`` is populated for every u, which only holds for the "ending at"
    convention.)

    RNAplfold is used rather than a global partition function because the latter is
    O(n^3) and a transcript is multiple kb; the windowed form is also the quantity
    local-opening models are defined on.
    """

    def __init__(self, unpaired: list[list[float]], length: int, u_max: int) -> None:
        """Wrap an already-computed pfl_fold_up matrix; prefer :meth:`fold` to build one."""
        self._unpaired = unpaired
        self._length = length
        self._u_max = u_max

    @classmethod
    def fold(
        cls,
        sequence: str,
        window_size: int = DEFAULT_PLFOLD_WINDOW,
        max_bp_span: int = DEFAULT_PLFOLD_MAX_BP_SPAN,
        u_max: int = SEED_ANCHORED_WINDOW_NT,
    ) -> "TargetAccessibilityProfile":
        """Fold one transcript and cache its opening probabilities.

        Args:
            sequence: Full transcript sequence; read as RNA, so T and U are equivalent.
            window_size: RNAplfold averaging window (W).
            max_bp_span: RNAplfold maximum base-pair span (L).
            u_max: Longest unpaired stretch to tabulate. Must cover every window queried.
        """
        rna = sequence.upper().replace("T", "U")
        unpaired = RNA.pfl_fold_up(rna, u_max, window_size, max_bp_span)
        return cls(unpaired, len(rna), u_max)

    def site_accessibility(self, start_pos: int, site_length: int) -> TargetSiteAccessibility:
        """Look up the three 3'-anchored windows for a site.

        Args:
            start_pos: Start of the target site (0-based) in the folded transcript.
            site_length: Length of the target site in nucleotides.
        """
        site_end = start_pos + site_length  # 1-based position of the site's 3'-most base
        if start_pos < 0 or site_end > self._length:
            return TargetSiteAccessibility(None, None, None)

        return TargetSiteAccessibility(
            seed_end_8mer=self._opening(site_end, SEED_END_WINDOW_NT),
            seed_anchored_17mer=self._opening(site_end, SEED_ANCHORED_WINDOW_NT),
            whole_site=self._opening(site_end, site_length),
        )

    def _opening(self, end_pos: int, window: int) -> float | None:
        """P(the length-``window`` stretch ending at 1-based ``end_pos`` is unpaired)."""
        if window > self._u_max or end_pos < window or end_pos > self._length:
            return None
        return float(self._unpaired[end_pos][window])


class ThermodynamicCalculator:
    """Calculate thermodynamic properties for siRNA candidates using ViennaRNA."""

    def __init__(self, temperature: float = 37.0):
        """Initialize thermodynamic calculator.

        Args:
            temperature: Temperature in Celsius for calculations
        """
        self.temperature = temperature

        # Set ViennaRNA temperature
        RNA.cvar.temperature = temperature
        self.model_details = RNA.md()
        self.model_details.temperature = temperature

    def calculate_duplex_stability(self, guide: str, passenger: str) -> float:
        """Calculate duplex stability (deltaG) using ViennaRNA.

        Both strands are supplied 5'->3' and are already complementary, so neither
        is reverse-complemented here: ViennaRNA's ``&`` cofold notation pairs the
        two strands antiparallel on its own.
        """
        if len(guide) != len(passenger):
            # TODO: convert to warning and generate opposite sequence from the guide strand
            raise ValueError("Guide and passenger sequences must be same length")
        # TODO we should save this mfe structure to have alongside the other dotplot?
        return self._cofold_mfe(guide, passenger)

    def calculate_asymmetry_score(self, candidate: SiRNACandidate) -> tuple[float, float, float]:
        """Calculate thermodynamic asymmetry score using ViennaRNA.

        Duplex ends are antiparallel, so the guide 5' end pairs with the passenger
        3' end and vice versa. Both windows are ``END_WINDOW_NT`` long so the two
        ΔG values are directly comparable.

        Returns:
            Tuple of (5' end stability, 3' end stability, asymmetry score)
        """
        guide = candidate.guide_sequence
        passenger = candidate.passenger_sequence
        window = min(END_WINDOW_NT, len(guide), len(passenger))

        dg_5p = self._calculate_end_stability(guide[:window], passenger[-window:])
        dg_3p = self._calculate_end_stability(guide[-window:], passenger[:window])

        # Asymmetry score: favor when 5' end is less stable (higher dG)
        # Higher score when dg_5p > dg_3p (5' end less stable)
        asymmetry_raw = dg_5p - dg_3p

        # Normalize to 0-1 scale
        asymmetry_score = max(0.0, min(1.0, (asymmetry_raw + 5.0) / 10.0))

        return dg_5p, dg_3p, asymmetry_score

    def calculate_target_accessibility(
        self,
        target_sequence: str,
        start_pos: int,
        sirna_length: int,
        window_size: int = DEFAULT_PLFOLD_WINDOW,
        max_bp_span: int = DEFAULT_PLFOLD_MAX_BP_SPAN,
    ) -> "TargetSiteAccessibility":
        """Compute local opening probabilities for one target site on the mRNA.

        Convenience single-site entry point. Callers scoring many sites on the same
        transcript should build a :class:`TargetAccessibilityProfile` once instead:
        this folds the whole transcript per call.

        Args:
            target_sequence: Full target mRNA sequence.
            start_pos: Start of the siRNA target site (0-based).
            sirna_length: Length of the target site in nucleotides.
            window_size: RNAplfold averaging window (W).
            max_bp_span: RNAplfold maximum base-pair span (L).

        Returns:
            TargetSiteAccessibility for the site.
        """
        profile = TargetAccessibilityProfile.fold(
            target_sequence,
            window_size=window_size,
            max_bp_span=max_bp_span,
            u_max=max(sirna_length, SEED_ANCHORED_WINDOW_NT),
        )
        return profile.site_accessibility(start_pos, sirna_length)

    def _calculate_end_stability(self, guide_end: str, passenger_end: str) -> float:
        """Calculate stability of duplex end using ViennaRNA.

        Args:
            guide_end: Guide-strand terminal window, 5'->3'
            passenger_end: The passenger window that pairs with it, 5'->3'
        """
        if not guide_end or not passenger_end:
            return 0.0
        return self._cofold_mfe(guide_end, passenger_end)

    def _cofold_mfe(self, strand_a: str, strand_b: str) -> float:
        """Return the cofold MFE of two complementary strands, each given 5'->3'."""
        # Removed RNA.OPTION_EVAL_ONLY to fix segfault
        fc = RNA.fold_compound(f"{strand_a}&{strand_b}", self.model_details)
        _, mfe = fc.mfe()
        return float(mfe)

    def calculate_melting_temperature(self, guide: str, passenger: str) -> float:
        """Calculate duplex melting temperature in °C.

        Uses the RNA nearest-neighbour parameters of Xia et al. (1998) via
        Biopython. Tm needs ΔH and ΔS separately, which ViennaRNA's MFE does not
        provide, so scaling ΔG cannot produce a physical Tm: at the corrected
        duplex ΔG of a 21mer (about -39 kcal/mol) the previous
        ``37 + 2 * -ΔG`` approximation returned 101-124 °C.

        Args:
            guide: Guide strand, 5'->3'
            passenger: Passenger strand, 5'->3' (expected reverse complement of guide)
        """
        if len(guide) != len(passenger):
            raise ValueError("Guide and passenger sequences must be same length")

        guide_rna = guide.upper().replace("T", "U")
        if passenger.upper().replace("T", "U") != str(Seq(guide_rna).reverse_complement_rna()):
            logger.warning(
                "Passenger is not the exact reverse complement of the guide; "
                "melting temperature assumes a perfectly paired duplex."
            )

        tm = MeltingTemp.Tm_NN(
            Seq(guide_rna),
            nn_table=MeltingTemp.RNA_NN2,
            Na=TM_SODIUM_MM,
            dnac1=TM_STRAND_CONC_NM,
            dnac2=TM_STRAND_CONC_NM,
        )
        return float(tm)

    @staticmethod
    def meets_asymmetry_threshold(asymmetry_score: float, threshold: float) -> bool:
        """Check an already-computed asymmetry score against a threshold.

        Callers that have just computed the score use this instead of
        :meth:`is_thermodynamically_favorable` to avoid re-folding both duplex ends.
        """
        return asymmetry_score >= threshold

    def is_thermodynamically_favorable(
        self, candidate: SiRNACandidate, threshold: float = DEFAULT_MIN_ASYMMETRY_SCORE
    ) -> bool:
        """Check if candidate meets thermodynamic asymmetry threshold."""
        _, _, asymmetry_score = self.calculate_asymmetry_score(candidate)
        return self.meets_asymmetry_threshold(asymmetry_score, threshold)

    def calculate_secondary_structure(self, sequence: str) -> tuple[str, float, float]:
        """Calculate secondary structure for a sequence.

        Returns:
            Tuple of (structure, mfe, paired_fraction)
        """
        fc = RNA.fold_compound(sequence, self.model_details)
        structure, mfe = fc.mfe()

        # Calculate paired fraction
        paired_bases = structure.count("(") + structure.count(")")
        paired_fraction = paired_bases / len(structure) if structure else 0.0

        return structure, mfe, paired_fraction

    # # Fallback methods for when ViennaRNA is not available
    # def _fallback_duplex_stability(self, guide: str, passenger: str) -> float:  # noqa: ARG002
    #     """Fallback duplex stability calculation."""
    #     # Simple approximation based on GC content
    #     # Note: passenger parameter kept for interface consistency
    #     gc_content = (guide.count("G") + guide.count("C")) / len(guide)
    #     return -2.0 * gc_content * len(guide)

    # def _fallback_end_stability(self, guide_end: str, passenger_end: str) -> float:  # noqa: ARG002
    #     """Fallback end stability calculation."""
    #     # Note: passenger_end parameter kept for interface consistency
    #     if not guide_end:
    #         return 0.0
    #     gc_content = (guide_end.count("G") + guide_end.count("C")) / len(guide_end)
    #     return -2.0 * gc_content * len(guide_end)

    # def _fallback_accessibility(self, target_sequence: str, start_pos: int, sirna_length: int) -> tuple[float, float]:
    #     """Fallback accessibility calculation."""
    #     target_site = target_sequence[start_pos : start_pos + sirna_length]
    #     at_content = (target_site.count("A") + target_site.count("T") + target_site.count("U")) / len(target_site)
    #     # Higher AT content suggests better accessibility
    #     accessibility = 0.5 + (at_content - 0.5) * 0.5
    #     return max(0.0, min(1.0, accessibility)), -5.0

    # def _fallback_melting_temp(self, guide: str) -> float:
    #     """Fallback melting temperature calculation."""
    #     at_count = guide.count("A") + guide.count("T") + guide.count("U")
    #     gc_count = guide.count("G") + guide.count("C")
    #     return 2 * at_count + 4 * gc_count

    # def _fallback_structure(self, sequence: str) -> tuple[str, float, float]:
    #     """Fallback structure prediction."""
    #     # Return a simple dot-bracket structure
    #     length = len(sequence)
    #     structure = "." * length
    #     mfe = -length * 0.5  # Rough estimate
    #     paired_fraction = 0.0
    #     return structure, mfe, paired_fraction
