"""Core siRNA design algorithms and functionality."""

import math
import time

from Bio import SeqIO
from Bio.Seq import Seq

from sirnaforge.core.off_target import OffTargetAnalyzer
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate


class SiRNADesigner:
    """Main siRNA design engine following the algorithm specification."""

    def __init__(self, parameters: DesignParameters) -> None:
        """Initialize designer with given parameters."""
        self.parameters = parameters

    def design_from_file(self, input_file: str) -> DesignResult:
        """Design siRNAs from input FASTA file."""
        start_time = time.time()

        # Parse input sequences
        sequences = list(SeqIO.parse(input_file, "fasta"))
        if not sequences:
            raise ValueError(f"No sequences found in {input_file}")

        all_candidates = []

        # Process each sequence
        for seq_record in sequences:
            transcript_id = seq_record.id
            sequence = str(seq_record.seq).upper()

            # Generate candidates for this sequence
            candidates = self._enumerate_candidates(sequence, transcript_id)

            # Apply filters
            filtered_candidates = self._apply_filters(candidates)

            # Score candidates
            scored_candidates = self._score_candidates(filtered_candidates)

            all_candidates.extend(scored_candidates)

        # Sort by composite score (descending)
        all_candidates.sort(key=lambda x: x.composite_score, reverse=True)

        # Get top candidates
        top_candidates = all_candidates[: self.parameters.top_n]

        processing_time = time.time() - start_time

        return DesignResult(
            input_file=input_file,
            parameters=self.parameters,
            candidates=all_candidates,
            top_candidates=top_candidates,
            total_sequences=len(sequences),
            total_candidates=len(all_candidates),
            filtered_candidates=len([c for c in all_candidates if c.passes_filters]),
            processing_time=processing_time,
            tool_versions=self._get_tool_versions(),
        )

    def design_from_sequence(self, sequence: str, transcript_id: str = "seq1") -> DesignResult:
        """Design siRNAs from a single sequence."""
        start_time = time.time()

        sequence = sequence.upper()

        # Generate candidates
        candidates = self._enumerate_candidates(sequence, transcript_id)

        # Apply filters
        filtered_candidates = self._apply_filters(candidates)

        # Score candidates
        scored_candidates = self._score_candidates(filtered_candidates)

        # Sort by composite score (descending)
        scored_candidates.sort(key=lambda x: x.composite_score, reverse=True)

        # Get top candidates
        top_candidates = scored_candidates[: self.parameters.top_n]

        processing_time = time.time() - start_time

        return DesignResult(
            input_file="<direct_input>",
            parameters=self.parameters,
            candidates=scored_candidates,
            top_candidates=top_candidates,
            total_sequences=1,
            total_candidates=len(scored_candidates),
            filtered_candidates=len([c for c in scored_candidates if c.passes_filters]),
            processing_time=processing_time,
            tool_versions=self._get_tool_versions(),
        )

    def _enumerate_candidates(self, sequence: str, transcript_id: str) -> list[SiRNACandidate]:
        """Enumerate all possible siRNA candidates using sliding window."""
        candidates = []
        sirna_length = self.parameters.sirna_length

        # Slide window across sequence
        for i in range(len(sequence) - sirna_length + 1):
            target_seq = sequence[i : i + sirna_length]

            # Generate guide (antisense) and passenger (sense) sequences
            guide_seq = str(Seq(target_seq).reverse_complement())
            passenger_seq = target_seq

            # Create candidate ID
            candidate_id = f"{transcript_id}_{i + 1}_{i + sirna_length}"

            # Calculate basic properties
            gc_content = self._calculate_gc_content(guide_seq)

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

            candidates.append(candidate)

        return candidates

    def _apply_filters(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Apply hard filters as specified in the algorithm."""
        filters = self.parameters.filters
        filtered = []

        for candidate in candidates:
            issues = []
            passes = True

            # GC content filter
            if not (filters.gc_min <= candidate.gc_content <= filters.gc_max):
                issues.append(
                    f"GC content {candidate.gc_content:.1f}% outside range {filters.gc_min}-{filters.gc_max}%"
                )
                passes = False

            # Poly run filter (no runs of 4+ identical nucleotides)
            if self._has_poly_runs(candidate.guide_sequence, filters.max_poly_runs):
                issues.append(f"Contains runs of >{filters.max_poly_runs} identical nucleotides")
                passes = False

            # Update candidate with filter results
            candidate.passes_filters = passes
            candidate.quality_issues = issues

            filtered.append(candidate)

        return filtered

    def _score_candidates(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Score candidates using composite scoring algorithm."""
        weights = self.parameters.scoring

        for candidate in candidates:
            # Calculate component scores
            asym_score = self._calculate_asymmetry_score(candidate)
            gc_score = self._calculate_gc_score(candidate.gc_content)
            access_score = self._calculate_accessibility_score(candidate)
            ot_score = self._calculate_off_target_score(candidate)
            empirical_score = self._calculate_empirical_score(candidate)

            # Store component scores
            candidate.component_scores = {
                "asymmetry": asym_score,
                "gc_content": gc_score,
                "accessibility": access_score,
                "off_target": ot_score,
                "empirical": empirical_score,
            }

            # Calculate composite score
            composite = (
                weights.asymmetry * asym_score
                + weights.gc_content * gc_score
                + weights.accessibility * access_score
                + weights.off_target * ot_score
                + weights.empirical * empirical_score
            )

            # Normalize to 0-100 scale
            candidate.composite_score = composite * 100
            candidate.asymmetry_score = asym_score

        return candidates

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
            return 1.0 - paired_fraction

        except ImportError:
            # Fallback to simple heuristic
            guide = candidate.guide_sequence
            at_content = (guide.count("A") + guide.count("T") + guide.count("U")) / len(guide)

            # Moderate AT content suggests better accessibility
            return 1.0 - abs(at_content - 0.5) * 2.0

    def _calculate_off_target_score(self, candidate: SiRNACandidate) -> float:
        """Calculate off-target score using enhanced analysis."""
        try:
            analyzer = OffTargetAnalyzer()
            return analyzer.calculate_off_target_score(candidate)
        except ImportError:
            # Fallback to simplified version
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
        """Calculate empirical score using Reynolds et al. rules (simplified)."""
        guide = candidate.guide_sequence
        score = 0.5  # Base score

        # Some simplified Reynolds rules
        # Prefer A/U at position 19 (3' end of guide)
        if len(guide) >= 19 and guide[18] in ["A", "U"]:
            score += 0.1

        # Prefer G/C at position 1
        if guide[0] in ["G", "C"]:
            score += 0.1

        # Avoid C at position 19
        if len(guide) >= 19 and guide[18] == "C":
            score -= 0.1

        return max(0.0, min(1.0, score))

    def _get_tool_versions(self) -> dict[str, str]:
        """Get versions of tools used (placeholder)."""
        return {
            "python": "3.11.5",
            "biopython": "1.81",
            "sirnaforge": "0.1.0",
        }
