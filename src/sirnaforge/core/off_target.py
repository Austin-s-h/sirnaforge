"""Off-target analysis for siRNA design."""

import math
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from sirnaforge.models.sirna import SiRNACandidate


class OffTargetAnalyzer:
    """Analyze off-target potential for siRNA candidates."""

    def __init__(
        self,
        transcriptome_file: Optional[str] = None,
        bowtie_index: Optional[str] = None,
        bwa_index: Optional[str] = None,
        use_external_tools: bool = False,
    ):
        """
        Initialize off-target analyzer.

        Args:
            transcriptome_file: Path to transcriptome FASTA file for off-target search
            bowtie_index: Path to Bowtie index for seed search
            bwa_index: Path to BWA index for full-length search
            use_external_tools: Whether to use external alignment tools
        """
        self.transcriptome_file = transcriptome_file
        self.bowtie_index = bowtie_index
        self.bwa_index = bwa_index
        self.use_external_tools = use_external_tools
        self.transcriptome_seqs: dict[str, str] = {}

        if transcriptome_file and Path(transcriptome_file).exists():
            self._load_transcriptome()

    def _load_transcriptome(self) -> None:
        """Load transcriptome sequences for off-target analysis."""
        # Would use Bio.SeqIO in production
        self.transcriptome_seqs = {}

    def analyze_off_targets(self, candidate: SiRNACandidate) -> tuple[int, float]:
        """
        Analyze off-target potential for a siRNA candidate.

        Returns:
            Tuple of (off_target_count, penalty_score)
        """
        if self.use_external_tools and (self.bowtie_index or self.bwa_index):
            return self._analyze_with_external_tools(candidate)
        if self.transcriptome_seqs:
            return self._search_transcriptome(candidate)
        return self._analyze_sequence_features(candidate)

    def _analyze_with_external_tools(self, candidate: SiRNACandidate) -> tuple[int, float]:
        """Analyze off-targets using Bowtie/BWA like the Nextflow pipeline."""
        guide = candidate.guide_sequence
        off_target_count = 0
        penalty = 0.0

        try:
            # Create temporary FASTA file with the guide sequence
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp_file:
                tmp_file.write(f">{candidate.id}\n{guide}\n")
                tmp_fasta = tmp_file.name

            # Run Bowtie for seed search if available
            if self.bowtie_index:
                bowtie_hits = self._run_bowtie_search(tmp_fasta, guide)
                off_target_count += bowtie_hits
                penalty += bowtie_hits * 50  # High penalty for seed matches

            # Run BWA for full-length search if available
            if self.bwa_index:
                bwa_hits = self._run_bwa_search(tmp_fasta, guide)
                off_target_count += bwa_hits
                penalty += bwa_hits * 20  # Moderate penalty for full-length matches

            # Clean up
            Path(tmp_fasta).unlink(missing_ok=True)

        except Exception:
            # Fall back to sequence-based analysis
            off_target_count, penalty = self._analyze_sequence_features(candidate)

        return off_target_count, penalty

    def _run_bowtie_search(self, query_fasta: str, guide: str) -> int:  # noqa: ARG002
        """Run Bowtie search for seed matches."""
        # Note: query_fasta parameter kept for interface consistency but not currently used
        # as we generate the seed sequence directly from the guide
        try:
            # Extract seed region (positions 2-8)
            seed = guide[1:8] if len(guide) >= 8 else guide[1:]

            # Create seed FASTA
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp_file:
                tmp_file.write(f">seed\n{seed}\n")
                seed_fasta = tmp_file.name

            # Run Bowtie with parameters for seed search
            assert self.bowtie_index is not None, "Bowtie index must be set"
            cmd = [
                "bowtie",
                "-v",
                "0",  # No mismatches in seed
                "-a",  # Report all alignments
                "--sam",  # SAM output
                self.bowtie_index,
                seed_fasta,
            ]

            result = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=30)

            # Count alignments (non-header lines)
            alignments = [line for line in result.stdout.split("\n") if line and not line.startswith("@")]

            # Clean up
            Path(seed_fasta).unlink(missing_ok=True)

            return len(alignments)

        except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
            return 0

    def _run_bwa_search(self, query_fasta: str, guide: str) -> int:  # noqa: ARG002
        """Run BWA search for full-length matches."""
        # Note: guide parameter kept for interface consistency but not currently used
        # as the query_fasta already contains the sequence
        try:
            # Run BWA-MEM with parameters for short sequences
            assert self.bwa_index is not None, "BWA index must be set"
            cmd = [
                "bwa-mem2",
                "mem",
                "-a",  # Report all alignments
                "-k",
                "12",  # Seed length
                "-T",
                "15",  # Minimum score threshold
                self.bwa_index,
                query_fasta,
            ]

            result = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=60)

            # Count alignments (non-header lines)
            alignments = [line for line in result.stdout.split("\n") if line and not line.startswith("@")]

            return len(alignments)

        except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
            return 0

    def _analyze_sequence_features(self, candidate: SiRNACandidate) -> tuple[int, float]:
        """Analyze sequence features that correlate with off-target risk."""
        guide = candidate.guide_sequence
        penalty = 0.0
        off_target_count = 0

        # Check for repetitive elements
        penalty += self._check_low_complexity(guide)
        penalty += self._check_seed_uniqueness(guide)
        penalty += self._check_gc_skew(guide)

        # Estimate off-target count from penalty
        if penalty > 100:
            off_target_count = int(penalty / 20)

        return off_target_count, penalty

    def _search_transcriptome(self, candidate: SiRNACandidate) -> tuple[int, float]:
        """Search transcriptome for potential off-targets."""
        # This would implement transcriptome search using candidate.guide_sequence
        # For now, return placeholder values until transcriptome search is implemented
        _ = candidate  # Acknowledge the parameter is intentionally unused for now
        return 0, 0.0

    def _check_low_complexity(self, sequence: str) -> float:
        """Check for low complexity regions that increase off-target risk."""
        penalty = 0.0

        # Check for homopolymer runs
        for base in ["A", "T", "G", "C", "U"]:
            max_run = self._longest_run(sequence, base)
            if max_run >= 4:  # TODO parameterize threshold
                penalty += max_run * 5

        # Check for dinucleotide repeats
        for i in range(len(sequence) - 3):
            dinuc = sequence[i : i + 2]
            next_dinuc = sequence[i + 2 : i + 4]
            if dinuc == next_dinuc:
                penalty += 10

        return penalty

    def _check_seed_uniqueness(self, guide: str) -> float:
        """Check seed region (positions 2-8) for uniqueness."""
        if len(guide) < 8:
            return 0.0

        seed = guide[1:8]  # Positions 2-8 (0-based indexing)
        penalty = 0.0

        # Check for internal seed matches within the guide itself
        guide_without_seed = guide[:1] + guide[8:]
        seed_matches = len(re.findall(seed, guide_without_seed))
        penalty += seed_matches * 50

        # Simple complexity check for seed
        unique_bases = len(set(seed))
        if unique_bases < 3:
            penalty += 20

        return penalty

    def _check_gc_skew(self, sequence: str) -> float:
        """Check for GC skew that might affect off-target binding."""
        penalty = 0.0

        # Check 5' end GC content (positions 1-10)
        five_prime = sequence[:10]
        five_gc = (five_prime.count("G") + five_prime.count("C")) / len(five_prime)

        # Check 3' end GC content (positions 12-21)
        three_prime = sequence[11:21] if len(sequence) >= 21 else sequence[11:]
        if three_prime:
            three_gc = (three_prime.count("G") + three_prime.count("C")) / len(three_prime)

            # Extreme skew increases off-target risk
            skew = abs(five_gc - three_gc)
            if skew > 0.6:
                penalty += skew * 30

        return penalty

    def _longest_run(self, sequence: str, base: str) -> int:
        """Find the longest run of a specific base in sequence."""
        max_run = 0
        current_run = 0

        for char in sequence:
            if char == base:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0

        return max_run

    def calculate_off_target_score(self, candidate: SiRNACandidate) -> float:
        """
        Calculate off-target score (0-1, higher is better).

        Transforms penalty to score: OT_score = exp(-penalty/50)
        """
        _, penalty = self.analyze_off_targets(candidate)

        # Transform penalty to score using exponential decay
        return math.exp(-penalty / 50.0)

    def get_seed_sequence(self, guide: str) -> str:
        """Extract seed sequence (positions 2-8) from guide."""
        if len(guide) < 8:
            return guide[1:] if len(guide) > 1 else guide
        return guide[1:8]

    def is_seed_unique(self, candidate: SiRNACandidate, max_matches: int = 0) -> bool:
        """Check if seed sequence is sufficiently unique."""
        penalty = self._check_seed_uniqueness(candidate.guide_sequence)
        # Convert penalty to estimated matches
        estimated_matches = int(penalty / 50)
        return estimated_matches <= max_matches
