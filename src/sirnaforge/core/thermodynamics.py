"""Thermodynamic calculations for siRNA design using ViennaRNA."""

import RNA
from Bio.Seq import Seq

from sirnaforge.models.sirna import SiRNACandidate


class ThermodynamicCalculator:
    """Calculate thermodynamic properties for siRNA candidates using ViennaRNA."""

    def __init__(self, temperature: float = 37.0):
        """
        Initialize thermodynamic calculator.

        Args:
            temperature: Temperature in Celsius for calculations
        """
        self.temperature = temperature

        # Set ViennaRNA temperature
        RNA.cvar.temperature = temperature
        self.model_details = RNA.md()
        self.model_details.temperature = temperature

    def calculate_duplex_stability(self, guide: str, passenger: str) -> float:
        """Calculate duplex stability (deltaG) using ViennaRNA."""
        if len(guide) != len(passenger):
            # TODO: convert to warning and generate opposite sequence from the guide strand
            raise ValueError("Guide and passenger sequences must be same length")
        # Create duplex string for ViennaRNA
        duplex_seq = guide + "&" + str(Seq(passenger).reverse_complement())

        # Calculate duplex MFE (removed RNA.OPTION_EVAL_ONLY to fix segfault)
        fc = RNA.fold_compound(duplex_seq, self.model_details)
        mfe_structure, mfe = fc.mfe()

        return float(mfe)

    def calculate_asymmetry_score(self, candidate: SiRNACandidate) -> tuple[float, float, float]:
        """
        Calculate thermodynamic asymmetry score using ViennaRNA.

        Returns:
            Tuple of (5' end stability, 3' end stability, asymmetry score)
        """
        guide = candidate.guide_sequence
        passenger = candidate.passenger_sequence

        # Calculate stability of 5' end (positions 1-7)
        dg_5p = self._calculate_end_stability(guide[:7], passenger[:7])

        # Calculate stability of 3' end (positions 15-21)
        guide_3p = guide[14:21] if len(guide) >= 21 else guide[14:]
        passenger_3p = passenger[14:21] if len(passenger) >= 21 else passenger[14:]
        dg_3p = self._calculate_end_stability(guide_3p, passenger_3p)

        # Asymmetry score: favor when 5' end is less stable (higher dG)
        # Higher score when dg_5p > dg_3p (5' end less stable)
        asymmetry_raw = dg_5p - dg_3p

        # Normalize to 0-1 scale
        asymmetry_score = max(0.0, min(1.0, (asymmetry_raw + 5.0) / 10.0))

        return dg_5p, dg_3p, asymmetry_score

    def calculate_target_accessibility(
        self, target_sequence: str, start_pos: int, sirna_length: int
    ) -> tuple[float, float]:
        """
        Calculate target site accessibility using ViennaRNA.

        Args:
            target_sequence: Full target mRNA sequence
            start_pos: Start position of siRNA target site (0-based)
            sirna_length: Length of siRNA

        Returns:
            Tuple of (average_unpaired_probability, mfe)
        """
        # Create fold compound for target sequence
        fc = RNA.fold_compound(target_sequence, self.model_details)

        # Calculate MFE structure
        mfe_structure, mfe = fc.mfe()

        # Calculate partition function and base pair probabilities
        fc.pf()

        # Get unpaired probabilities for target site
        unpaired_probs = []

        for i in range(start_pos, min(start_pos + sirna_length, len(target_sequence))):
            # Get probability that position i is unpaired
            prob = fc.pr_unpaired(i + 1)  # ViennaRNA uses 1-based indexing
            unpaired_probs.append(prob)

        avg_unpaired = sum(unpaired_probs) / len(unpaired_probs) if unpaired_probs else 0.0

        return avg_unpaired, mfe

    def _calculate_end_stability(self, guide_end: str, passenger_end: str) -> float:
        """Calculate stability of duplex end using ViennaRNA."""
        if not guide_end or not passenger_end:
            return 0.0
        # Create duplex for the end region
        duplex_seq = guide_end + "&" + str(Seq(passenger_end).reverse_complement())
        fc = RNA.fold_compound(duplex_seq, self.model_details)
        _, mfe = fc.mfe()

        return float(mfe)

    def calculate_melting_temperature(self, guide: str, passenger: str) -> float:
        """Calculate melting temperature using ViennaRNA thermodynamics."""
        # Calculate duplex stability
        dg = self.calculate_duplex_stability(guide, passenger)

        # Simplified Tm estimation: Tm = dH/dS (assumes dH ≈ -dG for rough estimate)
        # More accurate would require separate enthalpy calculation
        # This is a rough approximation
        if dg >= 0:
            return 25.0  # Low melting temp for unstable duplexes

        # Rough conversion: more negative dG → higher Tm
        tm = 37.0 + (-dg * 2.0)  # Empirical scaling factor

        return max(0.0, tm)

    def is_thermodynamically_favorable(self, candidate: SiRNACandidate, threshold: float = 0.5) -> bool:
        """Check if candidate meets thermodynamic asymmetry threshold."""
        _, _, asymmetry_score = self.calculate_asymmetry_score(candidate)
        return asymmetry_score >= threshold

    def calculate_secondary_structure(self, sequence: str) -> tuple[str, float, float]:
        """
        Calculate secondary structure for a sequence.

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
