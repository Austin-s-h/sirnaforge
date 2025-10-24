"""Core siRNA design algorithms and functionality."""

import math
import sys
import time
from typing import Optional, Union

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

from sirnaforge import __version__
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.models.sirna import SiRNACandidate as _ModelCandidate


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

        all_candidates: list[SiRNACandidate] = []
        # Map guide_sequence -> set of transcript_ids where it appears
        guide_to_transcripts: dict[str, set[str]] = {}

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

        processing_time = time.time() - start_time
        # Compute transcript hit metrics for each candidate (how many input transcripts contain the guide)
        total_seqs = len(sequences)
        for c in all_candidates:
            hits = len(guide_to_transcripts.get(c.guide_sequence, {c.transcript_id}))
            c.transcript_hit_count = hits
            c.transcript_hit_fraction = hits / total_seqs if total_seqs > 0 else 0.0

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

        # Get top candidates only from those passing filters; fallback to all if none pass
        passing = [
            c
            for c in scored_candidates
            if (c.passes_filters is True)
            or (hasattr(_ModelCandidate, "FilterStatus") and c.passes_filters == _ModelCandidate.FilterStatus.PASS)
        ]
        top_candidates = (passing or scored_candidates)[: self.parameters.top_n]

        processing_time = time.time() - start_time
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
        )

    def _enumerate_candidates(self, sequence: str, transcript_id: str) -> list[SiRNACandidate]:
        """Enumerate all possible siRNA candidates using sliding window with early filtering."""
        candidates = []
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
            if not (filters.gc_min <= gc_content <= filters.gc_max):
                continue

            # Early filtering: check for poly runs before creating candidate object
            if self._has_poly_runs(guide_seq, filters.max_poly_runs):
                continue

            # Create candidate ID with project moniker and sanitized transcript id
            # Format: SIRNAF_<TRANSCRIPT>_<start>_<end>
            # Sanitize transcript_id: keep alphanumerics and underscore, replace others with '-'
            safe_tid = "".join([c if (c.isalnum() or c == "_") else "-" for c in transcript_id])
            # Truncate long transcript ids to keep IDs short while retaining uniqueness
            if len(safe_tid) > 24:
                safe_tid = safe_tid[:24]
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

            candidates.append(candidate)

        return candidates

    def _apply_filters(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Apply remaining filters (early GC and poly-run filtering already done in enumeration)."""
        filtered = []

        for candidate in candidates:
            issues: list[str] = []
            status: Union[bool, _ModelCandidate.FilterStatus] = True

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
            # TODO: The composite_score needs to have basic weighting applied with truth data
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
                "off_target": ot_score,
                "empirical": empirical_score,
            }

            # Calculate composite score using configurable weights from parameters
            # Access the configured scoring weights
            weights = self.parameters.scoring

            # Apply the configured weights to each component
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

    def _calculate_duplex_score(self, candidate: SiRNACandidate) -> tuple[float, Optional[float]]:
        """Compute duplex stability ΔG and a normalized score in [0,1].

        Mapping: dg in [-40, -5] kcal/mol -> score in [1, 0]. Clamp outside this range.
        On failure or missing backend, returns (asymmetry_score, None) as a fallback.
        """
        try:
            calc = ThermodynamicCalculator()
            dg = calc.calculate_duplex_stability(candidate.guide_sequence, candidate.passenger_sequence)
            # Normalize: more negative is better
            # Clamp dg to [-40, -5]
            lo, hi = -40.0, -5.0
            dg_clamped = max(lo, min(hi, dg))
            score = (-(dg_clamped) - 5.0) / (40.0 - 5.0)
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
            try:
                if (
                    hasattr(_ModelCandidate, "FilterStatus")
                    and paired_fraction is not None
                    and paired_fraction > self.parameters.filters.max_paired_fraction
                    and candidate.passes_filters is True
                ):
                    candidate.passes_filters = _ModelCandidate.FilterStatus.EXCESS_PAIRING
            except (AttributeError, ValueError, TypeError):  # nosec B110 acceptable narrow handling
                # Ignore unexpected attribute/value issues; filtering status remains unchanged.
                pass
            return 1.0 - paired_fraction

        except ImportError:
            # Fallback to simple heuristic
            guide = candidate.guide_sequence
            at_content = (guide.count("A") + guide.count("T") + guide.count("U")) / len(guide)

            # Moderate AT content suggests better accessibility
            paired_fraction = abs(at_content - 0.5) * 2.0  # heuristic inverse
            candidate.paired_fraction = paired_fraction
            try:
                if (
                    hasattr(_ModelCandidate, "FilterStatus")
                    and paired_fraction > self.parameters.filters.max_paired_fraction
                    and candidate.passes_filters is True
                ):
                    candidate.passes_filters = _ModelCandidate.FilterStatus.EXCESS_PAIRING
            except (AttributeError, ValueError, TypeError):  # nosec B110
                pass
            return 1.0 - paired_fraction

    def _calculate_off_target_score(self, candidate: SiRNACandidate) -> float:
        """Calculate off-target score using simplified analysis."""
        # Simplified version - comprehensive off-target analysis would require
        # external databases and more complex alignment tools
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

        result = max(0.0, min(1.0, score))
        # Enforce minimal asymmetry threshold as a filter
        try:
            if (
                hasattr(_ModelCandidate, "FilterStatus")
                and result < self.parameters.filters.min_asymmetry_score
                and candidate.passes_filters is True
            ):
                candidate.passes_filters = _ModelCandidate.FilterStatus.LOW_ASYMMETRY
        except (AttributeError, ValueError, TypeError):  # nosec B110
            pass
        return result

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
