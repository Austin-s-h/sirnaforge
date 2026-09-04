"""Core siRNA design algorithms and functionality."""

import hashlib
import logging
import math
import sys
import time

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

from sirnaforge import __version__
from sirnaforge.core.repeat_detection import RepeatObservation, normalize_guide_sequence
from sirnaforge.core.scoring import ScoringError, compute_composite
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.sirna import (
    EMPIRICAL_SCORE_MAX,
    EMPIRICAL_SCORE_MIN,
    DesignParameters,
    DesignResult,
    SiRNACandidate,
)
from sirnaforge.models.sirna import SiRNACandidate as _ModelCandidate

logger = logging.getLogger(__name__)

# Duplex ΔG scales with duplex length, so it is normalised per nucleotide before
# scoring. Window taken from the ViennaRNA range measured for canonical siRNA
# duplexes at 37 °C: strong -> 1.0, weak -> 0.0.
DUPLEX_DG_PER_NT_STRONG = -2.1
DUPLEX_DG_PER_NT_WEAK = -1.4

# Sanitized transcript ids embedded in candidate ids are truncated to this length. Beyond it,
# the last ID_DIGEST_LEN chars are replaced by a deterministic digest of the FULL sanitized id
# (hashlib, not hash(), which is PYTHONHASHSEED-dependent) so two ids that share a long common
# prefix cannot collide after truncation.
TRANSCRIPT_ID_MAX_LEN = 24
ID_DIGEST_LEN = 8


def _as_rna(sequence: str) -> str:
    """Read a stored (DNA) sequence as RNA so T and U compare equal."""
    return sequence.upper().replace("T", "U")


class SiRNADesigner:
    """Main siRNA design engine following the algorithm specification."""

    def __init__(self, parameters: DesignParameters) -> None:
        """Initialize designer with given parameters."""
        self.parameters = parameters
        self.last_guide_to_transcripts: dict[str, set[str]] | None = None

    def design_from_file(self, input_file: str) -> DesignResult:
        """Design siRNAs from input FASTA file."""
        start_time = time.perf_counter()

        # Parse input sequences
        sequences = list(SeqIO.parse(input_file, "fasta"))
        if not sequences:
            raise ValueError(f"No sequences found in {input_file}")

        all_candidates: list[SiRNACandidate] = []
        rejected_pool: list[SiRNACandidate] = []
        # Map guide_sequence -> set of transcript_ids where it appears
        guide_to_transcripts: dict[str, set[str]] = {}

        # Process each sequence
        for seq_record in sequences:
            transcript_id = seq_record.id
            sequence = str(seq_record.seq).upper()

            # Generate candidates for this sequence
            candidates, rejected = self._enumerate_candidates(sequence, transcript_id)
            rejected_pool.extend(rejected)

            # Apply filters
            filtered_candidates = self._apply_filters(candidates)

            # Score candidates
            scored_candidates = self._score_candidates(filtered_candidates)

            # Track which transcripts each guide appears in
            for c in scored_candidates:
                guide_to_transcripts.setdefault(c.guide_sequence, set()).add(c.transcript_id)

            all_candidates.extend(scored_candidates)

        # Sort by composite score (descending)
        all_candidates.sort(key=lambda x: x.composite_score, reverse=True)

        # Get top candidates only from those passing filters; fallback to all if none pass
        passing = [
            c
            for c in all_candidates
            if (c.passes_filters is True)
            or (hasattr(_ModelCandidate, "FilterStatus") and c.passes_filters == _ModelCandidate.FilterStatus.PASS)
        ]
        top_candidates = (passing or all_candidates)[: self.parameters.top_n]

        processing_time = max(0.0, time.perf_counter() - start_time)  # Ensure non-negative
        # Compute transcript hit metrics for each candidate (how many input transcripts contain the guide)
        total_seqs = len(sequences)
        for c in all_candidates:
            hits = len(guide_to_transcripts.get(c.guide_sequence, {c.transcript_id}))
            c.transcript_hit_count = hits
            c.transcript_hit_fraction = hits / total_seqs if total_seqs > 0 else 0.0

        # Expose guide_to_transcripts mapping for workflow's protein-coding coverage computation
        self.last_guide_to_transcripts = guide_to_transcripts

        return DesignResult(
            input_file=input_file,
            parameters=self.parameters,
            candidates=all_candidates,
            top_candidates=top_candidates,
            total_sequences=len(sequences),
            total_candidates=len(all_candidates),
            filtered_candidates=len(
                [
                    c
                    for c in all_candidates
                    if (c.passes_filters is True)
                    or (
                        hasattr(_ModelCandidate, "FilterStatus")
                        and c.passes_filters == _ModelCandidate.FilterStatus.PASS
                    )
                ]
            ),
            processing_time=processing_time,
            tool_versions=self._get_tool_versions(),
            rejected_candidates=rejected_pool,
        )

    def design_from_sequence(self, sequence: str, transcript_id: str = "seq1") -> DesignResult:
        """Design siRNAs from a single sequence."""
        start_time = time.perf_counter()

        sequence = sequence.upper()

        # Generate candidates
        candidates, rejected = self._enumerate_candidates(sequence, transcript_id)

        # Apply filters
        filtered_candidates = self._apply_filters(candidates)

        # Score candidates
        scored_candidates = self._score_candidates(filtered_candidates)

        # Sort by composite score (descending)
        scored_candidates.sort(key=lambda x: x.composite_score, reverse=True)

        # Get top candidates only from those passing filters; fallback to all if none pass
        passing = [
            c
            for c in scored_candidates
            if (c.passes_filters is True)
            or (hasattr(_ModelCandidate, "FilterStatus") and c.passes_filters == _ModelCandidate.FilterStatus.PASS)
        ]
        top_candidates = (passing or scored_candidates)[: self.parameters.top_n]

        processing_time = max(0.0, time.perf_counter() - start_time)
        # For single-sequence runs, transcript hit metrics are trivial (hits=1, fraction=1.0)
        for c in scored_candidates:
            c.transcript_hit_count = 1
            c.transcript_hit_fraction = 1.0
        return DesignResult(
            input_file="<direct_input>",
            parameters=self.parameters,
            candidates=scored_candidates,
            top_candidates=top_candidates,
            total_sequences=1,
            total_candidates=len(scored_candidates),
            filtered_candidates=len(
                [
                    c
                    for c in scored_candidates
                    if (c.passes_filters is True)
                    or (
                        hasattr(_ModelCandidate, "FilterStatus")
                        and c.passes_filters == _ModelCandidate.FilterStatus.PASS
                    )
                ]
            ),
            processing_time=processing_time,
            tool_versions=self._get_tool_versions(),
            rejected_candidates=rejected,
        )

    def _enumerate_candidates(
        self, sequence: str, transcript_id: str
    ) -> tuple[list[SiRNACandidate], list[SiRNACandidate]]:
        """Enumerate all possible siRNA candidates and record those failing early filters."""
        candidates: list[SiRNACandidate] = []
        rejected: list[SiRNACandidate] = []
        sirna_length = self.parameters.sirna_length
        filters = self.parameters.filters

        # Slide window across sequence
        for i in range(len(sequence) - sirna_length + 1):
            target_seq = sequence[i : i + sirna_length]

            # Generate guide (antisense) and passenger (sense) sequences
            guide_seq = str(Seq(target_seq).reverse_complement())
            passenger_seq = target_seq

            # Early filtering for computational efficiency
            gc_content = self._calculate_gc_content(guide_seq)
            fail_reason: SiRNACandidate.FilterStatus | None = None
            if not (filters.gc_min <= gc_content <= filters.gc_max):
                fail_reason = SiRNACandidate.FilterStatus.GC_OUT_OF_RANGE
            elif self._has_poly_runs(guide_seq, filters.max_poly_runs):
                fail_reason = SiRNACandidate.FilterStatus.POLY_RUNS

            # Create candidate ID with project moniker and sanitized transcript id
            # Format: SIRNAF_<TRANSCRIPT>_<start>_<end>
            # Sanitize transcript_id: keep alphanumerics and underscore, replace others with '-'
            safe_tid = "".join([c if (c.isalnum() or c == "_") else "-" for c in transcript_id])
            # Truncate long transcript ids, disambiguating with a digest of the full id so two
            # ids sharing a >=24-char prefix (e.g. de-novo assembly isoforms) cannot collide.
            if len(safe_tid) > TRANSCRIPT_ID_MAX_LEN:
                digest = hashlib.sha256(safe_tid.encode()).hexdigest()[:ID_DIGEST_LEN]
                keep = TRANSCRIPT_ID_MAX_LEN - ID_DIGEST_LEN - 1
                safe_tid = f"{safe_tid[:keep]}-{digest}"
            candidate_id = f"SIRNAF_{safe_tid}_{i + 1}_{i + sirna_length}"

            candidate = SiRNACandidate(
                id=candidate_id,
                transcript_id=transcript_id,
                position=i + 1,  # 1-based
                guide_sequence=guide_seq,
                passenger_sequence=passenger_seq,
                gc_content=gc_content,
                length=sirna_length,
                asymmetry_score=0.0,  # Will be calculated in scoring
                composite_score=0.0,  # Will be calculated in scoring
            )

            if fail_reason is not None:
                candidate.passes_filters = fail_reason
                issues = list(candidate.quality_issues or [])
                label = fail_reason.value if hasattr(fail_reason, "value") else str(fail_reason)
                issues.append(label)
                candidate.quality_issues = issues
                rejected.append(candidate)
                continue

            candidates.append(candidate)

        return candidates, rejected

    def _apply_filters(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Apply remaining filters (early GC and poly-run filtering already done in enumeration)."""
        filtered = []

        for candidate in candidates:
            issues: list[str] = []
            status: bool | _ModelCandidate.FilterStatus = True

            # Note: GC content and poly-run filtering already done in _enumerate_candidates
            # This is mainly for any additional filters or post-processing

            # Update candidate with filter results
            candidate.passes_filters = status
            candidate.quality_issues = issues

            filtered.append(candidate)

        return filtered

    def _score_candidates(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Score candidates using composite scoring algorithm."""
        for candidate in candidates:
            # Calculate component scores
            # Thermodynamic end stabilities and asymmetry
            dg5: float = float("nan")
            dg3: float = float("nan")
            try:
                calc = ThermodynamicCalculator()
                dg5, dg3, asym_score = calc.calculate_asymmetry_score(candidate)
            except Exception:
                # Fallback if ViennaRNA or calculation not available
                asym_score = self._calculate_asymmetry_score(candidate)
                dg5, dg3 = float("nan"), float("nan")
            # Thermodynamic duplex stability (ΔG) and score normalization
            dg_score, duplex_dg = self._calculate_duplex_score(candidate)
            candidate.duplex_stability = duplex_dg
            gc_score = self._calculate_gc_score(candidate.gc_content)
            access_score = self._calculate_accessibility_score(candidate)
            ot_score = self._calculate_off_target_score(candidate)
            empirical_score = self._calculate_empirical_score(candidate)
            self._apply_score_filters(candidate, asym_score, empirical_score)

            # Optional melting temperature estimation (rough)
            tm_c = float("nan")
            try:
                if "calc" in locals():
                    tm_c = calc.calculate_melting_temperature(candidate.guide_sequence, candidate.passenger_sequence)
                else:
                    calc2 = ThermodynamicCalculator()
                    tm_c = calc2.calculate_melting_temperature(candidate.guide_sequence, candidate.passenger_sequence)
            except Exception:
                tm_c = float("nan")

            # Store component scores
            # Combine thermodynamic components: favor asymmetry with contribution from duplex stability
            thermo_combo = 0.7 * asym_score + 0.3 * dg_score

            candidate.component_scores = {
                "asymmetry": asym_score,
                # store as float; if missing, use NaN to satisfy type expectations
                "duplex_stability_dg": float(duplex_dg) if duplex_dg is not None else float("nan"),
                "duplex_stability_score": dg_score,
                "thermo_combo": thermo_combo,
                "dg_5p": float(dg5),
                "dg_3p": float(dg3),
                "delta_dg_end": float(dg5 - dg3) if (not math.isnan(dg5) and not math.isnan(dg3)) else float("nan"),
                "melting_temp_c": float(tm_c),
                "gc_content": gc_score,
                "accessibility": access_score,
                "design_off_target_proxy": ot_score,  # Diagnostic only, not fed into composite
                "empirical": empirical_score,
            }

            # Compute composite score via the shared scorer with design-time terms only.
            # Build features mapping, omitting any NaN values (ViennaRNA failures).
            features: dict[str, float] = {}
            if not math.isnan(asym_score):
                features["asymmetry"] = asym_score
            if not math.isnan(gc_score):
                features["gc_content"] = gc_score
            if not math.isnan(access_score):
                features["accessibility"] = access_score
            if not math.isnan(empirical_score):
                features["empirical"] = empirical_score

            # Compute composite; if all features are NaN, leave composite_score at 0.0.
            if features:
                try:
                    result = compute_composite(features, self.parameters.scoring)
                    candidate.composite_score = result.score
                    candidate.weight_set_version = result.weight_set_version
                    candidate.scored_after_screening = False
                    # Write per-term contributions
                    candidate.score_asymmetry = result.contributions.get("asymmetry")
                    candidate.score_gc_content = result.contributions.get("gc_content")
                    candidate.score_accessibility = result.contributions.get("accessibility")
                    candidate.score_empirical = result.contributions.get("empirical")
                    # Post-screen terms stay None at design time
                    candidate.score_off_target = None
                    candidate.score_isoform_coverage = None
                    candidate.score_conservation = None
                except ScoringError as e:
                    logger.warning(f"Scoring failed for candidate {candidate.id}: {e}. Setting composite_score=0.0.")
                    candidate.composite_score = 0.0
                    candidate.weight_set_version = ""
                    candidate.scored_after_screening = False
            else:
                # All features were NaN; leave composite_score at its default (0.0)
                logger.warning(
                    f"All component scores are NaN for candidate {candidate.id}. Leaving composite_score=0.0."
                )
                candidate.weight_set_version = ""
                candidate.scored_after_screening = False

            candidate.asymmetry_score = asym_score

        return candidates

    def _calculate_duplex_score(self, candidate: SiRNACandidate) -> tuple[float, float | None]:
        """Compute duplex stability ΔG and a normalized score in [0,1].

        ΔG is normalised per nucleotide so 19-23 nt designs stay comparable:
        DUPLEX_DG_PER_NT_STRONG -> 1.0, DUPLEX_DG_PER_NT_WEAK -> 0.0, clamped outside.
        On failure or missing backend, returns (asymmetry_score, None) as a fallback.
        """
        try:
            calc = ThermodynamicCalculator()
            dg = calc.calculate_duplex_stability(candidate.guide_sequence, candidate.passenger_sequence)
            # Normalize per nucleotide: more negative is better
            dg_per_nt = dg / len(candidate.guide_sequence)
            span = DUPLEX_DG_PER_NT_WEAK - DUPLEX_DG_PER_NT_STRONG
            score = (DUPLEX_DG_PER_NT_WEAK - dg_per_nt) / span
            score = max(0.0, min(1.0, score))
            return score, float(dg)
        except Exception:
            # Fallback: use asymmetry as proxy if duplex calc not available
            try:
                asym = self._calculate_asymmetry_score(candidate)
            except Exception:
                asym = 0.5
            return asym, None

    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100

    def _has_poly_runs(self, sequence: str, max_runs: int) -> bool:
        """Check for runs of identical nucleotides exceeding threshold."""
        current_base = sequence[0]
        current_run = 1

        for base in sequence[1:]:
            if base == current_base:
                current_run += 1
                if current_run > max_runs:
                    return True
            else:
                current_base = base
                current_run = 1

        return False

    def _calculate_asymmetry_score(self, candidate: SiRNACandidate) -> float:
        """Calculate thermodynamic asymmetry score using enhanced method."""
        try:
            calc = ThermodynamicCalculator()
            _, _, asymmetry_score = calc.calculate_asymmetry_score(candidate)
            return asymmetry_score
        except ImportError:
            # Fallback to simplified version
            guide = candidate.guide_sequence

            # Calculate stability of 5' end (positions 1-7) vs 3' end (positions 15-21)
            five_prime_end = guide[:7]
            three_prime_end = guide[14:21] if len(guide) >= 21 else guide[14:]

            # Simplified AT/GC ratio as proxy for stability
            five_prime_gc = (five_prime_end.count("G") + five_prime_end.count("C")) / len(five_prime_end)
            three_prime_gc = (three_prime_end.count("G") + three_prime_end.count("C")) / len(three_prime_end)

            # Higher score when 5' end is less stable (lower GC) than 3' end
            asymmetry = three_prime_gc - five_prime_gc

            # Normalize to 0-1 range
            return max(0.0, min(1.0, (asymmetry + 1.0) / 2.0))

    def _calculate_gc_score(self, gc_content: float) -> float:
        """Calculate GC content score with Gaussian penalty around 40%."""
        # GC_score = exp(-((GC-40)/10)^2)
        return math.exp(-(((gc_content - 40) / 10) ** 2))

    def _calculate_accessibility_score(self, candidate: SiRNACandidate) -> float:
        """Calculate target accessibility score using ViennaRNA when available."""
        try:
            calc = ThermodynamicCalculator()

            # For single candidate analysis, we don't have full target sequence context
            # So we'll use the guide sequence as a proxy for structure prediction
            guide = candidate.guide_sequence
            structure, mfe, paired_fraction = calc.calculate_secondary_structure(guide)

            # Store structure info in candidate
            candidate.structure = structure
            candidate.mfe = mfe
            candidate.paired_fraction = paired_fraction

            # Accessibility score: 1 - paired_fraction
            # Flag excessive pairing per filter threshold
            self._flag_excess_pairing(candidate, paired_fraction)
            return 1.0 - paired_fraction

        except ImportError:
            # Fallback to simple heuristic
            guide = candidate.guide_sequence
            at_content = (guide.count("A") + guide.count("T") + guide.count("U")) / len(guide)

            # Moderate AT content suggests better accessibility
            paired_fraction = abs(at_content - 0.5) * 2.0  # heuristic inverse
            candidate.paired_fraction = paired_fraction
            self._flag_excess_pairing(candidate, paired_fraction)
            return 1.0 - paired_fraction

    def _flag_excess_pairing(self, candidate: SiRNACandidate, paired_fraction: float) -> None:
        """Flag a candidate whose guide is too structured to be accessible."""
        if candidate.passes_filters is True and paired_fraction > self.parameters.filters.max_paired_fraction:
            candidate.passes_filters = _ModelCandidate.FilterStatus.EXCESS_PAIRING

    def _calculate_off_target_score(self, candidate: SiRNACandidate) -> float:
        """Score internal sequence repetitiveness as a design-time off-target proxy.

        This looks only at repeated 7-mers *within* the guide. It is not informed by
        transcriptome or miRNA screening: those run after design and are applied as
        a pass/fail gate (see SiRNAWorkflow._integrate_offtarget_results), never fed
        back into composite_score. `off_target_screened` records whether a candidate
        reached that stage.
        """
        # Comprehensive off-target analysis would require external databases and
        # more complex alignment tools
        guide = candidate.guide_sequence

        # Simple penalty for repetitive sequences
        penalty = 0
        for i in range(len(guide) - 6):
            seed = guide[i : i + 7]
            # Count occurrences of this 7-mer in the sequence
            if guide.count(seed) > 1:
                penalty += 10

        # Transform penalty to score: OT_score = exp(-penalty/50)
        candidate.off_target_penalty = penalty
        return math.exp(-penalty / 50)

    def _calculate_empirical_score(self, candidate: SiRNACandidate) -> float:
        """Calculate empirical score using Reynolds et al. rules (simplified).

        Guides are stored as DNA, so the sequence is read as RNA (T is U) before the
        position-19 test; otherwise a T there never earned the A/U bonus. The
        attainable range is EMPIRICAL_SCORE_MIN..EMPIRICAL_SCORE_MAX, not 0..1.
        """
        guide = candidate.guide_sequence.upper().replace("T", "U")
        score = 0.5  # Base score

        # Some simplified Reynolds rules
        # Prefer A/U at position 19 (3' end of guide)
        if len(guide) >= 19 and guide[18] in ("A", "U"):
            score += 0.1

        # Prefer G/C at position 1
        if guide[0] in ("G", "C"):
            score += 0.1

        # Avoid C at position 19
        if len(guide) >= 19 and guide[18] == "C":
            score -= 0.1

        return max(EMPIRICAL_SCORE_MIN, min(EMPIRICAL_SCORE_MAX, score))

    def _apply_score_filters(self, candidate: SiRNACandidate, asymmetry_score: float, empirical_score: float) -> None:
        """Flag candidates failing the asymmetry or empirical-rule thresholds.

        Each threshold gates the quantity it is named after: min_asymmetry_score
        gates the thermodynamic asymmetry score, min_empirical_score gates the
        empirical design-rule score. Only the first failure is recorded, matching
        the earlier GC / poly-run / excess-pairing gates.
        """
        if candidate.passes_filters is not True:
            return

        filters = self.parameters.filters
        if not ThermodynamicCalculator.meets_asymmetry_threshold(asymmetry_score, filters.min_asymmetry_score):
            candidate.passes_filters = _ModelCandidate.FilterStatus.LOW_ASYMMETRY
        elif empirical_score < filters.min_empirical_score:
            candidate.passes_filters = _ModelCandidate.FilterStatus.LOW_EMPIRICAL_SCORE

    @staticmethod
    def stamp_repeat_verdict(candidate: SiRNACandidate, observations: dict[str, RepeatObservation]) -> None:
        """Stamp repeat metadata and verdict on a single candidate if its guide is flagged.

        The REPEAT_ELEMENT verdict is applied only if the candidate is currently passing
        (passes_filters is True or PASS). A candidate that already failed for another
        reason (GC, asymmetry, etc.) retains its earlier verdict — precedence is:
        existing failure > REPEAT_ELEMENT > PASS.

        Args:
            candidate: Candidate to potentially flag.
            observations: Mapping from normalized guide sequence to RepeatObservation.
        """
        norm_guide = normalize_guide_sequence(candidate.guide_sequence)
        obs = observations.get(norm_guide)
        if obs is None:
            return

        # Write repeat metadata regardless of verdict
        candidate.repeat_flagged = obs.is_repeat
        candidate.repeat_transcript_fraction = obs.transcript_fraction

        # Apply REPEAT_ELEMENT verdict only if currently passing
        if obs.is_repeat and (
            candidate.passes_filters is True or candidate.passes_filters == _ModelCandidate.FilterStatus.PASS
        ):
            candidate.passes_filters = _ModelCandidate.FilterStatus.REPEAT_ELEMENT

    def _get_tool_versions(self) -> dict[str, str]:
        """Get versions of tools used in the analysis."""
        try:
            biopython_version = Bio.__version__
        except AttributeError:
            biopython_version = "unknown"

        python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

        return {
            "python": python_version,
            "biopython": biopython_version,
            "sirnaforge": __version__,
        }


# component_scores keys carrying the miRNA biogenesis bonus past the design stage. The bonus is
# folded into composite_score, not into the composite term set, so post-screen rescoring (which
# rebuilds the composite from the term set) has to be able to reapply it from the candidate itself.
MIRNA_BONUS_KEY = "mirna_biogenesis_bonus"
MIRNA_BONUS_MAX_KEY = "mirna_biogenesis_bonus_max"


def apply_mirna_biogenesis_bonus(base_score: float, mirna_bonus: float, max_mirna_bonus: float) -> float:
    """Fold the miRNA biogenesis bonus into a 0-100 composite score.

    The bonuses widen the attainable range, so rescale by the maximum attainable total
    instead of clamping: clamping parked every strong candidate at exactly 100.0 and
    erased the ranking at the top. Order is preserved, since this is monotone in
    (base score + bonus).
    """
    scaled = (base_score + mirna_bonus * 100) / (1.0 + max_mirna_bonus)
    # Guard the model's 0-100 bound for non-default scoring weights
    return max(0.0, min(100.0, scaled))


class MiRNADesigner(SiRNADesigner):
    """miRNA-biogenesis-aware siRNA designer with specialized scoring.

    Extends SiRNADesigner with scoring rules optimized for miRNA-like processing:
    - Argonaute selection preferences (pos1 A/U, mismatch at pos1)
    - 3' supplementary pairing analysis (positions 13-16)
    - Conservative thermodynamic thresholds
    - Seed region quality assessment
    """

    def __init__(self, parameters: DesignParameters) -> None:
        """Initialize miRNA designer with miRNA-specific config validation."""
        super().__init__(parameters)
        # Optionally apply miRNA-specific filter adjustments based on MiRNADesignConfig
        # For now, we rely on the caller to set appropriate filters for miRNA mode

    def _score_candidates(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Score candidates using miRNA-biogenesis-aware composite scoring.

        Adds miRNA-specific scoring components:
        - Argonaute start bonus for A/U at guide position 1
        - Position 1 mismatch/wobble preference
        - 3' supplementary pairing score
        - Enhanced asymmetry requirements
        """
        from sirnaforge.models.sirna import MiRNADesignConfig  # noqa: PLC0415

        mirna_config = MiRNADesignConfig()
        scoring_weights = mirna_config.scoring_weights
        # Every bonus below is capped by its weight (supp_bonus scales a [0,1] score).
        max_mirna_bonus = (
            scoring_weights["ago_start_bonus"]
            + scoring_weights["pos1_mismatch_bonus"]
            + scoring_weights["supp_13_16_bonus"]
        )

        # First, run the standard scoring
        candidates = super()._score_candidates(candidates)

        # Add miRNA-specific scoring enhancements
        for candidate in candidates:
            guide = candidate.guide_sequence
            passenger = candidate.passenger_sequence

            # 1. Argonaute selection: prefer A/U at guide position 1
            # The reported base keeps the stored (DNA) spelling; the test does not.
            guide_pos1_base = guide[0] if guide else ""
            candidate.guide_pos1_base = guide_pos1_base
            ago_start_bonus = scoring_weights["ago_start_bonus"] if _as_rna(guide_pos1_base) in ("A", "U") else 0.0

            # 2. Position 1 pairing state: prefer G:U wobble or mismatch over perfect pair
            pos1_pairing_state = self._classify_pos1_pairing(guide_pos1_base, passenger[-1] if passenger else "")
            candidate.pos1_pairing_state = pos1_pairing_state
            pos1_mismatch_bonus = (
                scoring_weights["pos1_mismatch_bonus"] if pos1_pairing_state in ["wobble", "mismatch"] else 0.0
            )

            # 3. 3' supplementary pairing (positions 13-16)
            supp_score = self._calculate_supplementary_score(guide)
            candidate.supp_13_16_score = supp_score
            supp_bonus = scoring_weights["supp_13_16_bonus"] * supp_score

            # 4. Seed class classification (positions 2-8)
            seed_class = self._classify_seed_region(guide)
            candidate.seed_class = seed_class

            # 5. Apply miRNA-specific bonuses to composite score
            mirna_bonus = ago_start_bonus + pos1_mismatch_bonus + supp_bonus

            # Recorded on the candidate so off-target screening, which recomputes the composite
            # from the term set afterwards, can reapply the same bonus instead of dropping it.
            candidate.component_scores[MIRNA_BONUS_KEY] = mirna_bonus
            candidate.component_scores[MIRNA_BONUS_MAX_KEY] = max_mirna_bonus

            candidate.composite_score = apply_mirna_biogenesis_bonus(
                candidate.composite_score, mirna_bonus, max_mirna_bonus
            )

        return candidates

    def _classify_pos1_pairing(self, guide_base: str, passenger_base: str) -> str:
        """Classify pairing state at guide position 1.

        Args:
            guide_base: Guide strand base at position 1 (5' end)
            passenger_base: Passenger strand base at position 21 (pairs with guide pos1)

        Returns:
            Pairing classification: "perfect", "wobble", or "mismatch"
        """
        # Watson-Crick pairs
        perfect_pairs = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
        # G:U wobble pair
        wobble_pairs = {("G", "U"), ("U", "G")}

        # Bases are stored as DNA, so read them as RNA before lookup: otherwise A:T
        # is not found in perfect_pairs and every A:U pair is called a mismatch.
        pair = (_as_rna(guide_base), _as_rna(passenger_base))
        if pair in perfect_pairs:
            return "perfect"
        if pair in wobble_pairs:
            return "wobble"
        return "mismatch"

    def _calculate_supplementary_score(self, guide: str) -> float:
        """Calculate 3' supplementary pairing potential (positions 13-16).

        Lower score is better (less 3' pairing = better specificity).

        Args:
            guide: Guide strand sequence

        Returns:
            Normalized score [0-1] where 0 = high pairing, 1 = low pairing
        """
        if len(guide) < 16:
            return 0.5  # Default for short sequences

        # Extract positions 13-16 (0-indexed: 12-15)
        supp_region = _as_rna(guide[12:16])

        # Simple heuristic: count A/U content (lower stability)
        au_count = supp_region.count("A") + supp_region.count("U")
        # High A/U = low stability = high score (good for avoiding 3' pairing)
        return au_count / len(supp_region) if supp_region else 0.5

    def _classify_seed_region(self, guide: str) -> str:
        """Classify seed match class based on guide positions 2-8.

        Args:
            guide: Guide strand sequence

        Returns:
            Seed class: "6mer", "7mer-m8", "7mer-a1", or "8mer"
        """
        # This is a simplified classification based on seed length
        # In practice, seed class depends on target matching, which happens during off-target analysis
        # For now, we just categorize based on sequence properties
        if len(guide) < 8:
            return "6mer"
        # This is a placeholder - actual seed class determined during off-target matching
        return "8mer"
