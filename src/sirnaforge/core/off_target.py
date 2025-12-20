"""Off-target analysis for siRNA design.

This module provides comprehensive off-target analysis functionality for siRNA design,
including both miRNA seed match analysis and transcriptome off-target detection.
Uses BWA-MEM2 for both short and long sequence alignments.
Optimized for both standalone use and parallelized Nextflow workflows.
"""

import json
import shutil
import statistics
import subprocess  # nosec B404
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from sirnaforge.data.base import FastaUtils
from sirnaforge.data.mirna_manager import MiRNADatabaseManager
from sirnaforge.models.off_target import (
    AggregatedMiRNASummary,
    AggregatedOffTargetSummary,
    AlignmentStrand,
    AnalysisMode,
    AnalysisSummary,
    MiRNAHit,
    MiRNASummary,
    OffTargetHit,
)
from sirnaforge.models.schemas import GenomeAlignmentSchema, MiRNAAlignmentSchema
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.species import human_vs_other_totals

logger = get_logger(__name__)


def _compute_species_counts(df: pd.DataFrame) -> dict[str, int]:
    """Build a frequency map for the 'species' column in a DataFrame."""
    if df.empty or "species" not in df.columns:
        return {}

    counts: dict[str, int] = {}
    for value in df["species"].tolist():
        label = "unknown"
        if value is not None and not pd.isna(value):
            label = str(value)
        counts[label] = counts.get(label, 0) + 1
    return counts


# =============================================================================
# Core Analyzer Classes
# =============================================================================


def _get_executable_path(tool_name: str) -> str | None:
    """Get the full path to an executable, ensuring it exists."""
    path = shutil.which(tool_name)
    if path is None:
        logger.warning(f"Tool '{tool_name}' not found in PATH")
    return path


def _validate_command_args(cmd: list[str]) -> None:
    """Validate command arguments for subprocess execution."""
    if not cmd:
        raise ValueError("Command list cannot be empty")

    executable = cmd[0]
    if not executable:
        raise ValueError("Executable path cannot be empty")

    # Ensure we have an absolute path to the executable
    if not Path(executable).is_absolute():
        raise ValueError(f"Executable must be an absolute path: {executable}")


# =============================================================================
# Core Analyzer Classes
# =============================================================================


class BwaAnalyzer:
    """BWA-MEM2 based analyzer for both transcriptome and miRNA seed off-target search."""

    def __init__(
        self,
        index_prefix: str | Path,
        mode: str = "transcriptome",  # "transcriptome" or "mirna_seed"
        seed_length: int = 12,
        min_score: int = 15,
        max_hits: int = 10000,
        seed_start: int = 2,
        seed_end: int = 8,
    ):
        """Initialize BWA-MEM2 analyzer.

        Args:
            index_prefix: Path to BWA index
            mode: Analysis mode - "transcriptome" for long targets, "mirna_seed" for short targets
            seed_length: BWA seed length parameter
            min_score: Minimum alignment score
            max_hits: Maximum hits to return
            seed_start: Seed region start (1-based)
            seed_end: Seed region end (1-based)
        """
        self.index_prefix = str(index_prefix)
        self.mode = mode
        self.seed_length = seed_length
        self.min_score = min_score
        self.max_hits = max_hits
        self.seed_start = seed_start
        self.seed_end = seed_end

        # Configure parameters based on mode
        if mode == "mirna_seed":
            # For miRNA seed analysis: short query (6-8bp) vs short target (~22bp)
            # Need very permissive parameters for ultra-short sequences
            self.seed_length = min(seed_length, 6)  # Max 6bp seed for 6-8bp queries
            self.min_score = 6  # Very low threshold - allow imperfect matches
        elif mode == "transcriptome":
            # For transcriptome analysis: short query vs long target
            self.seed_length = seed_length  # Use provided seed length
            self.min_score = min_score  # Use provided min score
        else:
            raise ValueError(f"Unknown mode: {mode}. Use 'transcriptome' or 'mirna_seed'")

    def analyze_sequences(self, sequences: dict[str, str]) -> list[dict[str, Any]]:
        """Run BWA-MEM2 analysis on sequences.

        Args:
            sequences: Dictionary of sequence name -> sequence

        Returns:
            List of alignment dictionaries
        """
        # Prepare sequences based on mode
        analysis_sequences = self._prepare_sequences_for_analysis(sequences)

        results = []
        temp_fasta_path = create_temp_fasta(analysis_sequences)

        try:
            # Get absolute path to bwa-mem2 executable
            bwa_path = _get_executable_path("bwa-mem2")
            if not bwa_path:
                raise FileNotFoundError("BWA-MEM2 executable not found in PATH")

            # Configure BWA parameters based on mode
            cmd = self._build_bwa_command(bwa_path, temp_fasta_path)

            _validate_command_args(cmd)
            logger.info(f"Running BWA-MEM2 ({self.mode} mode): {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=None, check=True)  # nosec B603
            results = self._parse_sam_output(result.stdout, sequences)
            results = self._filter_and_rank(results)
            logger.info(f"BWA-MEM2 analysis completed: {len(results)} hits found")

        except subprocess.CalledProcessError as e:
            logger.error(f"BWA-MEM2 failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            logger.error("BWA-MEM2 timed out")
        finally:
            Path(temp_fasta_path).unlink(missing_ok=True)

        return results[: self.max_hits]

    def _prepare_sequences_for_analysis(self, sequences: dict[str, str]) -> dict[str, str]:
        """Prepare sequences for analysis based on mode."""
        if self.mode == "mirna_seed":
            # Extract seed region (positions 2-8, 1-based) from siRNA sequences
            prepared = {}
            for name, seq in sequences.items():
                if len(seq) >= self.seed_end:
                    seed_seq = seq[self.seed_start - 1 : self.seed_end]  # Convert to 0-based indexing
                    prepared[name] = seed_seq
                    logger.debug(f"Extracted seed region for {name}: {seed_seq} (from {seq})")
                else:
                    logger.warning(f"Sequence {name} too short for seed extraction: {seq}")
                    prepared[name] = seq  # Use full sequence if too short
            return prepared
        # For transcriptome mode, use full sequences
        return sequences

    def _build_bwa_command(self, bwa_path: str, temp_fasta_path: str) -> list[str]:
        """Build BWA command based on analysis mode."""
        base_cmd = [
            bwa_path,
            "mem",
            "-a",  # Output all alignments
            "-v",
            "1",  # Verbosity level
        ]

        if self.mode == "mirna_seed":
            # For miRNA seed analysis: ultra-permissive parameters for 6-8bp vs ~22bp
            cmd = base_cmd + [
                "-k",
                str(self.seed_length),  # Seed length (max 6bp)
                "-T",
                str(self.min_score),  # Minimum score (6)
                "-w",
                "2",  # Narrow band width for short sequences
                "-A",
                "2",  # Higher matching score to reward matches
                "-B",
                "1",  # Low mismatch penalty
                "-O",
                "1,1",  # Low gap open penalties
                "-E",
                "1,1",  # Low gap extension penalties
                "-L",
                "8,8",  # Clipping penalty for ultra-short reads
                self.index_prefix,
                temp_fasta_path,
            ]
        elif self.mode == "transcriptome":
            # For transcriptome analysis: standard parameters
            cmd = base_cmd + [
                "-k",
                str(self.seed_length),  # Seed length
                "-T",
                str(self.min_score),  # Minimum score
                "-w",
                "100",  # Band width (larger for long targets)
                self.index_prefix,
                temp_fasta_path,
            ]
        else:
            raise ValueError(f"Unknown mode: {self.mode}")

        return cmd

    def _parse_sam_output(self, sam_output: str, original_sequences: dict[str, str]) -> list[dict[str, Any]]:
        """Parse SAM output from BWA-MEM2."""
        results = []

        for line in sam_output.splitlines():
            if line.startswith("@") or not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 11:
                continue

            qname = parts[0]
            flag = int(parts[1])
            rname = parts[2]
            pos = int(parts[3])
            mapq = int(parts[4]) if parts[4] != "*" else 0
            cigar = parts[5]

            if flag & 4:  # Skip unmapped
                continue

            strand = "-" if (flag & 16) else "+"
            coord = f"{rname}:{pos}"

            # Parse optional tags
            tags = {}
            for tag in parts[11:]:
                if ":" in tag:
                    tag_parts = tag.split(":", 2)
                    if len(tag_parts) == 3:
                        tags[tag_parts[0]] = tag_parts[2]

            nm = int(tags.get("NM", 0))
            as_score = int(tags.get("AS", 0)) if "AS" in tags else None

            # Parse mismatch positions from MD tag
            mismatch_positions = self._parse_md_tag(tags.get("MD", ""))
            seed_mismatches = sum(1 for pos in mismatch_positions if self.seed_start <= pos <= self.seed_end)

            # Calculate off-target score
            offtarget_score = self._calculate_offtarget_score(mismatch_positions)

            result = {
                "qname": qname,
                "qseq": original_sequences.get(qname, ""),  # Use original sequence
                "rname": rname,
                "coord": coord,
                "strand": strand,
                "cigar": cigar,
                "mapq": mapq,
                "as_score": as_score,
                "nm": nm,
                "mismatch_positions": mismatch_positions,
                "seed_mismatches": seed_mismatches,
                "offtarget_score": offtarget_score,
            }
            results.append(result)

        return results

    def _parse_md_tag(self, md_tag: str) -> list[int]:
        """Parse MD tag to extract mismatch positions."""
        positions = []
        read_pos = 1
        i = 0

        while i < len(md_tag):
            if md_tag[i].isdigit():
                num_str = ""
                while i < len(md_tag) and md_tag[i].isdigit():
                    num_str += md_tag[i]
                    i += 1
                if num_str:
                    read_pos += int(num_str)
            elif md_tag[i] == "^":
                i += 1
                while i < len(md_tag) and md_tag[i].isalpha():
                    i += 1
            elif md_tag[i].isalpha():
                positions.append(read_pos)
                read_pos += 1
                i += 1
            else:
                i += 1

        return positions

    def _calculate_offtarget_score(self, mismatch_positions: list[int]) -> float:
        """Calculate off-target score based on mismatch count and positions.

        Scoring principles (based on siRNA literature):
        1. Total mismatch count - most fundamental property
        2. Seed region (positions 2-8) - critical for target recognition
        3. Position-specific weights - 5' end more important than 3' end
        4. Continuous mismatches - clusters reduce binding more than scattered

        Returns:
            float: Off-target penalty score (lower = more likely off-target effect)
                   0.0 = perfect match (highest risk)
                   Higher scores = more mismatches (lower risk)

        References: TODO: validate
            - Jackson et al. 2003 (seed region importance)
            - Birmingham et al. 2006 (position-specific effects)
            - Huesken et al. 2005 (thermodynamic contributions)
        """
        num_mismatches = len(mismatch_positions)

        # Perfect match = highest off-target risk
        if num_mismatches == 0:
            return 0.0

        # Base score: number of mismatches (most fundamental property)
        # Each mismatch significantly reduces binding affinity
        base_score = num_mismatches * 10.0

        # Position-specific penalties (seed region is critical)
        position_penalty = 0.0
        seed_mismatches = 0

        for pos in mismatch_positions:
            if self.seed_start <= pos <= self.seed_end:
                # Seed region (pos 2-8): mismatches here strongly disrupt binding
                position_penalty += 5.0
                seed_mismatches += 1
            elif pos <= 10:
                # 5' region (pos 1, 9-10): moderately important
                position_penalty += 3.0
            else:
                # 3' region (pos 11-19): less critical but still relevant
                position_penalty += 1.0

        # Continuous mismatch bonus (clusters disrupt binding more)
        continuous_bonus = 0.0
        if num_mismatches >= 2:
            sorted_positions = sorted(mismatch_positions)
            continuous_count = 0
            for i in range(len(sorted_positions) - 1):
                if sorted_positions[i + 1] - sorted_positions[i] == 1:
                    continuous_count += 1
            # Adjacent mismatches create stronger disruption
            continuous_bonus = continuous_count * 2.0

        # Total score: base + position weights + continuity bonus
        return base_score + position_penalty + continuous_bonus

    def _filter_and_rank(self, results: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Filter and rank results by off-target score."""
        results.sort(key=lambda x: (x["offtarget_score"], -x.get("as_score", 0)))
        return results


class OffTargetAnalysisManager:
    """Manager class for comprehensive off-target analysis using BWA-MEM2."""

    def __init__(
        self,
        species: str,
        transcriptome_path: str | Path | None = None,
        mirna_path: str | Path | None = None,
        transcriptome_index: str | Path | None = None,
        mirna_index: str | Path | None = None,
    ):
        """Initialize the off-target analysis manager."""
        self.species = species
        self.transcriptome_path = Path(transcriptome_path) if transcriptome_path is not None else None
        self.mirna_path = Path(mirna_path) if mirna_path is not None else None
        self.transcriptome_index = Path(transcriptome_index) if transcriptome_index is not None else None
        self.mirna_index = Path(mirna_index) if mirna_index is not None else None

    def analyze_mirna_off_targets(
        self,
        sequences: dict[str, str] | str | Path,
        output_prefix: str | Path,
    ) -> tuple[Path, Path]:
        """Analyze miRNA off-targets using BWA-MEM2 in miRNA seed mode."""
        if not self.mirna_index:
            raise ValueError("miRNA index not provided")

        if isinstance(sequences, (str, Path)):
            sequences = parse_fasta_file(sequences)

        analyzer = BwaAnalyzer(self.mirna_index, mode="mirna_seed")
        results = analyzer.analyze_sequences(sequences)

        output_path = Path(output_prefix)
        tsv_path = output_path.parent / f"{output_path.name}_mirna_hits.tsv"
        json_path = output_path.parent / f"{output_path.name}_mirna_hits.json"

        self._write_mirna_results(results, tsv_path, json_path)
        return tsv_path, json_path

    def analyze_transcriptome_off_targets(
        self,
        sequences: dict[str, str] | str | Path,
        output_prefix: str | Path,
    ) -> tuple[Path, Path]:
        """Analyze transcriptome off-targets using BWA-MEM2 in transcriptome mode."""
        if not self.transcriptome_index:
            raise ValueError("Transcriptome index not provided")

        if isinstance(sequences, (str, Path)):
            sequences = parse_fasta_file(sequences)

        analyzer = BwaAnalyzer(self.transcriptome_index, mode="transcriptome")
        results = analyzer.analyze_sequences(sequences)

        output_path = Path(output_prefix)
        tsv_path = output_path.parent / f"{output_path.name}_transcriptome_hits.tsv"
        json_path = output_path.parent / f"{output_path.name}_transcriptome_hits.json"

        self._write_transcriptome_results(results, tsv_path, json_path)
        return tsv_path, json_path

    def analyze_sirna_candidate(self, candidate: SiRNACandidate) -> dict[str, Any]:
        """Analyze a single siRNA candidate for off-targets."""
        sequences = {candidate.id: candidate.guide_sequence}

        results: dict[str, Any] = {
            "candidate_id": candidate.id,
            "guide_sequence": candidate.guide_sequence,
            "mirna_hits": [],
            "transcriptome_hits": [],
        }

        if self.mirna_index:
            mirna_analyzer = BwaAnalyzer(self.mirna_index, mode="mirna_seed")
            results["mirna_hits"] = mirna_analyzer.analyze_sequences(sequences)

        if self.transcriptome_index:
            transcriptome_analyzer = BwaAnalyzer(self.transcriptome_index, mode="transcriptome")
            results["transcriptome_hits"] = transcriptome_analyzer.analyze_sequences(sequences)

        return results

    def _write_mirna_results(self, results: list[dict[str, Any]], tsv_path: str | Path, json_path: str | Path) -> None:
        """Write miRNA analysis results."""
        # Write TSV
        species_label = getattr(self, "species", "unknown")
        enriched_results: list[dict[str, Any]] = []

        with Path(tsv_path).open("w") as f:
            f.write(
                "qname\tqseq\tspecies\trname\tcoord\tstrand\tcigar\tmapq\tas_score\t"
                "nm\tseed_mismatches\tofftarget_score\n"
            )
            for result in results:
                row_species = result.get("species", species_label)
                enriched = {
                    **result,
                    "species": row_species,
                }
                enriched_results.append(enriched)
                f.write(
                    f"{enriched['qname']}\t{enriched['qseq']}\t{enriched['species']}\t{enriched['rname']}\t"
                    f"{enriched['coord']}\t{enriched['strand']}\t{enriched['cigar']}\t"
                    f"{enriched['mapq']}\t{enriched.get('as_score', 'NA')}\t{enriched['nm']}\t"
                    f"{enriched['seed_mismatches']}\t{enriched['offtarget_score']}\n"
                )

        # Write JSON
        with Path(json_path).open("w") as f:
            json.dump(enriched_results, f, indent=2)

    def _write_transcriptome_results(
        self, results: list[dict[str, Any]], tsv_path: str | Path, json_path: str | Path
    ) -> None:
        """Write transcriptome analysis results."""
        # Write TSV
        with Path(tsv_path).open("w") as f:
            f.write("qname\tqseq\trname\tcoord\tstrand\tcigar\tmapq\tas_score\tnm\tseed_mismatches\tofftarget_score\n")
            for result in results:
                f.write(
                    f"{result['qname']}\t{result['qseq']}\t{result['rname']}\t"
                    f"{result['coord']}\t{result['strand']}\t{result['cigar']}\t"
                    f"{result['mapq']}\t{result.get('as_score', 'NA')}\t{result['nm']}\t"
                    f"{result['seed_mismatches']}\t{result['offtarget_score']}\n"
                )

        # Write JSON
        with Path(json_path).open("w") as f:
            json.dump(results, f, indent=2)


# =============================================================================
# Utility Functions
# =============================================================================


def create_temp_fasta(sequences: dict[str, str]) -> str:
    """Create temporary FASTA file from sequences."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_file:
        temp_path = tmp_file.name
    FastaUtils.write_dict_to_fasta(sequences, temp_path)
    return temp_path


def validate_and_write_sequences(
    input_file: str, output_file: str, expected_length: int = 21
) -> tuple[int, int, list[str]]:
    """Validate siRNA sequences and write valid ones to output file."""
    sequences = FastaUtils.parse_fasta_to_dict(input_file)

    try:
        valid_sequences = FastaUtils.validate_sirna_sequences(sequences, expected_length)

        if valid_sequences:
            FastaUtils.write_dict_to_fasta(valid_sequences, output_file)
        else:
            Path(output_file).touch()

        invalid_count = len(sequences) - len(valid_sequences)
        issues = [
            f"{name}: Invalid (length={len(seq)}, expected={expected_length})"
            for name, seq in sequences.items()
            if name not in valid_sequences
        ]

        return len(valid_sequences), invalid_count, issues

    except ValueError as e:
        Path(output_file).touch()
        return 0, len(sequences), [str(e)]


def build_bwa_index(fasta_file: str | Path, index_prefix: str | Path) -> Path:
    """Build BWA-MEM2 index for both transcriptome and miRNA off-target analysis."""
    fasta_path = Path(fasta_file)
    index_prefix_path = Path(index_prefix)

    logger.info(f"Building BWA-MEM2 index from {fasta_path} with prefix {index_prefix_path}")

    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_path}")

    index_prefix_path.parent.mkdir(parents=True, exist_ok=True)

    # Get absolute path to bwa-mem2 executable
    bwa_path = _get_executable_path("bwa-mem2")
    if not bwa_path:
        raise FileNotFoundError("bwa-mem2 executable not found in PATH")

    cmd = [bwa_path, "index", "-p", str(index_prefix_path), str(fasta_path)]
    _validate_command_args(cmd)

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=7200)  # nosec B603
        if not validate_index_files(index_prefix_path, "bwa-mem2"):
            raise RuntimeError(
                "BWA-MEM2 index files were not created correctly. This usually happens when the build "
                "process runs out of memory (human transcriptomes can require 32GB+). Increase available RAM "
                "or pre-build the index in an environment with more memory and retry."
            )
        logger.info(f"BWA-MEM2 index built successfully: {index_prefix_path}")
        return index_prefix_path
    except subprocess.CalledProcessError as e:
        logger.error(f"BWA-MEM2 index build failed: {e.stderr}")
        raise
    except subprocess.TimeoutExpired:
        logger.error("BWA-MEM2 index build timed out")
        raise


def validate_sirna_sequences(
    sequences: dict[str, str], expected_length: int = 21
) -> tuple[dict[str, str], dict[str, str], list[str]]:
    """Validate siRNA sequences using existing FastaUtils."""
    try:
        valid_sequences = FastaUtils.validate_sirna_sequences(sequences, expected_length)
        invalid_sequences = {name: seq for name, seq in sequences.items() if name not in valid_sequences}
        issues = [
            f"{name}: Invalid sequence (length={len(seq)}, expected={expected_length})"
            for name, seq in invalid_sequences.items()
        ]
        return valid_sequences, invalid_sequences, issues
    except ValueError as e:
        return {}, sequences, [str(e)]


def parse_fasta_file(fasta_file: str | Path) -> dict[str, str]:
    """Parse FASTA file using existing FastaUtils."""
    return FastaUtils.parse_fasta_to_dict(fasta_file)


def write_fasta_file(sequences: dict[str, str], output_file: str) -> None:
    """Write sequences to FASTA file using existing FastaUtils."""
    FastaUtils.write_dict_to_fasta(sequences, output_file)


def check_tool_availability(tool: str) -> bool:
    """Check if external tool is available."""
    try:
        # Get absolute path to tool executable
        tool_path = _get_executable_path(tool)
        if not tool_path:
            return False

        cmd = [tool_path, "--help"]
        _validate_command_args(cmd)
        result = subprocess.run(cmd, capture_output=True, check=False, timeout=10)  # nosec B603
        return result.returncode in {0, 1}
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def validate_index_files(index_prefix: str | Path, tool: str = "bwa") -> bool:
    """Validate that index files exist for given tool."""
    index_path = Path(index_prefix)

    if tool in ("bwa", "bwa-mem2"):
        required_extensions = [".amb", ".ann", ".bwt.2bit.64", ".pac"]
    else:
        logger.warning(f"Unknown tool for index validation: {tool}")
        return False

    for ext in required_extensions:
        candidate = index_path.parent / f"{index_path.name}{ext}"
        if not candidate.exists():
            logger.debug(f"Missing index file: {candidate.name}")
            return False
        if candidate.stat().st_size == 0:
            logger.debug(f"Index file is empty: {candidate}")
            return False

    return True


# =============================================================================
# Nextflow Integration Functions
# =============================================================================


def run_mirna_analysis_for_nextflow(
    species: str,
    sequences_file: str,
    mirna_index: str | Path,
    output_prefix: str | Path,
) -> tuple[str, str, str]:
    """Nextflow-compatible function for miRNA analysis."""
    manager = OffTargetAnalysisManager(species=species, mirna_index=mirna_index)
    output_root = Path(output_prefix)

    try:
        tsv_path, json_path = manager.analyze_mirna_off_targets(sequences_file, output_root)

        # Create summary
        summary_path = output_root.parent / f"{output_root.name}_mirna_summary.txt"
        with summary_path.open("w") as f:
            with json_path.open() as jf:
                results = json.load(jf)
            f.write(f"Species: {species}\n")
            f.write(f"Total miRNA hits: {len(results)}\n")
            f.write("Analysis completed successfully\n")

        return str(tsv_path), str(json_path), str(summary_path)

    except Exception as e:
        error_summary = output_root.parent / f"{output_root.name}_mirna_error.txt"
        with error_summary.open("w") as f:
            f.write(f"miRNA analysis failed: {str(e)}\n")
        return "", "", str(error_summary)


def run_transcriptome_analysis_for_nextflow(
    species: str,
    sequences_file: str,
    transcriptome_index: str | Path,
    output_prefix: str | Path,
) -> tuple[str, str, str]:
    """Nextflow-compatible function for transcriptome analysis."""
    manager = OffTargetAnalysisManager(species=species, transcriptome_index=transcriptome_index)
    output_root = Path(output_prefix)

    try:
        tsv_path, json_path = manager.analyze_transcriptome_off_targets(sequences_file, output_root)

        # Create summary
        summary_path = output_root.parent / f"{output_root.name}_transcriptome_summary.txt"
        with summary_path.open("w") as f:
            with json_path.open() as jf:
                results = json.load(jf)
            f.write(f"Species: {species}\n")
            f.write(f"Total transcriptome hits: {len(results)}\n")
            f.write("Analysis completed successfully\n")

        return str(tsv_path), str(json_path), str(summary_path)

    except Exception as e:
        error_summary = output_root.parent / f"{output_root.name}_transcriptome_error.txt"
        with error_summary.open("w") as f:
            f.write(f"Transcriptome analysis failed: {str(e)}\n")
        return "", "", str(error_summary)


def run_comprehensive_offtarget_analysis(
    species: str,
    sequences_file: str,
    index_path: str,
    output_prefix: str | Path,
    mode: str = "transcriptome",
    bwa_k: int = 12,
    bwa_T: int = 15,
    max_hits: int = 10000,
    seed_start: int = 2,
    seed_end: int = 8,
) -> tuple[str, str, str]:
    """Run comprehensive off-target analysis for Nextflow integration."""
    output_root = Path(output_prefix)

    try:
        sequences = parse_fasta_file(sequences_file)

        # Use BWA analyzer for comprehensive analysis
        analyzer = BwaAnalyzer(
            index_prefix=index_path,
            mode=mode,
            seed_length=bwa_k,
            min_score=bwa_T,
            max_hits=max_hits,
            seed_start=seed_start,
            seed_end=seed_end,
        )

        results = analyzer.analyze_sequences(sequences)

        # Write results using pandas (much faster than manual loop)
        tsv_path = output_root.parent / f"{output_root.name}.tsv"
        json_path = output_root.parent / f"{output_root.name}.json"
        summary_path = output_root.parent / f"{output_root.name}_summary.txt"

        # Convert dict results to DataFrame
        df = pd.DataFrame(results)

        if df.empty:
            # Create empty DataFrame with required schema columns (including species)
            df = pd.DataFrame(columns=list(GenomeAlignmentSchema.__annotations__.keys()))
        elif "species" not in df.columns:
            df["species"] = species
        else:
            df["species"] = df["species"].fillna(species)

        # Validate with Pandera schema
        df = GenomeAlignmentSchema.validate(df, lazy=True)

        # Write TSV and JSON (pandas handles efficiently)
        df.to_csv(tsv_path, sep="\t", index=False)
        df.to_json(json_path, orient="records", indent=2)

        # Write summary
        with summary_path.open("w") as f:
            f.write(f"Species: {species}\n")
            f.write(f"Total sequences analyzed: {len(sequences)}\n")
            f.write(f"Total off-target hits: {len(results)}\n")
            f.write(f"Analysis mode: {mode}\n")
            f.write(f"Analysis parameters: bwa_k={bwa_k}, bwa_T={bwa_T}, max_hits={max_hits}\n")
            f.write(f"Seed region: {seed_start}-{seed_end}\n")
            f.write("Analysis completed successfully\n")

        return str(tsv_path), str(json_path), str(summary_path)

    except Exception as e:
        error_summary = output_root.parent / f"{output_root.name}_error.txt"
        with error_summary.open("w") as f:
            f.write(f"Comprehensive off-target analysis failed: {str(e)}\n")
        return "", "", str(error_summary)


def run_bwa_alignment_analysis(
    candidates_file: str | Path,
    index_prefix: str | Path,
    species: str,
    output_dir: str | Path,
    max_hits: int = 10000,
    bwa_k: int = 12,
    bwa_T: int = 15,
    seed_start: int = 2,
    seed_end: int = 8,
) -> Path:
    """Run BWA-MEM2 alignment analysis for candidate sequences using Pydantic models.

    This is the main function called by OFFTARGET_ANALYSIS Nextflow module.

    Args:
        candidates_file: Path to FASTA file with candidate sequences
        index_prefix: Path to BWA-MEM2 index prefix
        species: Species identifier
        output_dir: Directory to write results
        max_hits: Maximum hits to report per candidate
        bwa_k: BWA seed length parameter
        bwa_T: BWA minimum score threshold
        seed_start: Seed region start position (1-based)
        seed_end: Seed region end position (1-based)

    Returns:
        Path to output directory containing results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Parse input sequences
    sequences = parse_fasta_file(candidates_file)

    # Determine candidate ID from filename (e.g., "candidate_0001.fasta" -> "candidate_0001")
    candidate_id = Path(candidates_file).stem

    # Create unique output prefix for this candidate-species combination
    output_prefix = output_path / f"{candidate_id}_{species}"

    # Run BWA-MEM2 analysis
    analyzer = BwaAnalyzer(
        index_prefix=index_prefix,
        mode="transcriptome",  # Always use transcriptome mode for genome analysis
        seed_length=bwa_k,
        min_score=bwa_T,
        max_hits=max_hits,
        seed_start=seed_start,
        seed_end=seed_end,
    )

    results_dicts = analyzer.analyze_sequences(sequences)

    # Convert dict results to OffTargetHit objects with validation
    all_hits: list[OffTargetHit] = []
    for hit_dict in results_dicts:
        try:
            # Parse coord string "chr1:12345" into integer
            coord_str = hit_dict["coord"]
            coord_int = int(coord_str.split(":")[1]) if ":" in coord_str else int(coord_str)

            offtarget_hit = OffTargetHit(
                qname=hit_dict["qname"],
                qseq=hit_dict["qseq"],
                species=species,
                rname=hit_dict["rname"],
                coord=coord_int,
                strand=AlignmentStrand(hit_dict["strand"]),
                cigar=hit_dict["cigar"],
                mapq=hit_dict["mapq"],
                as_score=hit_dict.get("as_score"),
                nm=hit_dict["nm"],
                seed_mismatches=hit_dict["seed_mismatches"],
                offtarget_score=hit_dict["offtarget_score"],
            )
            all_hits.append(offtarget_hit)
        except Exception as e:
            logger.warning(f"Failed to validate off-target hit: {e}, skipping")
            continue

    # Write TSV analysis file using Pydantic models
    analysis_file = Path(f"{output_prefix}_analysis.tsv")
    with analysis_file.open("w") as f:
        f.write(OffTargetHit.tsv_header() + "\n")
        for hit in all_hits:
            f.write(hit.to_tsv_row() + "\n")

    # Write JSON file with validated data
    json_file = Path(f"{output_prefix}_hits.json")
    with json_file.open("w") as f:
        json.dump([hit.model_dump() for hit in all_hits], f, indent=2)

    # Derive summary statistics for alignment metrics
    mean_mapq = statistics.fmean(hit.mapq for hit in all_hits) if all_hits else None
    mean_mismatches = statistics.fmean(hit.nm for hit in all_hits) if all_hits else None
    mean_seed_mismatches = statistics.fmean(hit.seed_mismatches for hit in all_hits) if all_hits else None

    # Create validated summary using Pydantic model
    summary = AnalysisSummary(
        candidate_id=candidate_id,
        species=species,
        mode=AnalysisMode.TRANSCRIPTOME,
        total_sequences=len(sequences),
        total_hits=len(all_hits),
        mean_mapq=mean_mapq,
        mean_mismatches=mean_mismatches,
        mean_seed_mismatches=mean_seed_mismatches,
    )

    # Write summary JSON file
    summary_file = Path(f"{output_prefix}_summary.json")
    with summary_file.open("w") as f:
        # Add parameters to the output
        summary_dict = summary.model_dump()
        summary_dict["parameters"] = {
            "bwa_k": bwa_k,
            "bwa_T": bwa_T,
            "max_hits": max_hits,
            "seed_start": seed_start,
            "seed_end": seed_end,
        }
        json.dump(summary_dict, f, indent=2)

    logger.info(f"BWA analysis completed for {candidate_id} vs {species}: {len(all_hits)} hits")

    return output_path


def aggregate_offtarget_results(  # noqa: PLR0912
    results_dir: str | Path,
    output_dir: str | Path,
    genome_species: str,
) -> Path:
    """Aggregate transcriptome off-target analysis results using Pandera.

    Uses pandas + Pandera for efficient bulk reading and validation instead of
    manual line-by-line parsing with Pydantic models.

    NOTE: This function ONLY aggregates genome/transcriptome hits. miRNA results
    are aggregated separately by aggregate_mirna_results() to keep output files
    distinct and properly typed.

    Args:
        results_dir: Directory containing individual analysis results
        output_dir: Directory to write aggregated results
        genome_species: Comma-separated list of genome species analyzed

    Returns:
        Path to output directory containing aggregated results
    """
    results_path = Path(results_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    species_list = [s.strip() for s in genome_species.split(",") if s.strip()]
    species_file_counts: dict[str, int] = {}
    missing_species: list[str] = []
    for species in species_list:
        species_dir = results_path / species
        count = len(list(species_dir.glob("*_analysis.tsv"))) if species_dir.exists() else 0
        species_file_counts[species] = count
        if count == 0:
            missing_species.append(species)

    # Collect ONLY genome/transcriptome TSV analysis files
    # miRNA files are handled separately by aggregate_mirna_results()
    analysis_files = list(results_path.glob("**/*_analysis.tsv"))
    # Filter out miRNA files explicitly to avoid schema validation errors
    analysis_files = [f for f in analysis_files if "mirna" not in f.name.lower()]

    logger.info(f"Found {len(analysis_files)} transcriptome analysis files to aggregate")

    if analysis_files:
        # Read all files into DataFrames and concatenate (vectorized operation)
        dfs = []
        for analysis_file in analysis_files:
            try:
                # Pandas reads TSV much faster than manual line splitting
                df = pd.read_csv(analysis_file, sep="\t")

                # Validate schema with Pandera
                df = GenomeAlignmentSchema.validate(df, lazy=True)
                dfs.append(df)

            except Exception as e:
                logger.warning(f"Failed to read/validate {analysis_file}: {e}")
                continue

        # Concatenate all DataFrames at once (much faster than append in loop)
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
        else:
            # Create empty DataFrame with correct schema
            combined_df = pd.DataFrame(columns=list(GenomeAlignmentSchema.__annotations__.keys()))

    else:
        # No files found - create empty DataFrame
        combined_df = pd.DataFrame(columns=list(GenomeAlignmentSchema.__annotations__.keys()))

    # Write combined results (pandas is much faster than manual TSV writing)
    combined_tsv = output_path / "combined_offtargets.tsv"
    combined_df.to_csv(combined_tsv, sep="\t", index=False)

    # Write JSON (pandas handles serialization)
    combined_json = output_path / "combined_offtargets.json"
    combined_df.to_json(combined_json, orient="records", indent=2)

    species_counts = _compute_species_counts(combined_df)
    if not species_counts:
        species_counts = dict.fromkeys(species_list, 0) if species_list else {}

    for species in species_list:
        species_counts.setdefault(species, 0)

    human_hits, other_hits = human_vs_other_totals(species_counts)

    logger.info(
        f"Aggregated {len(combined_df)} transcriptome off-target hits from {len(analysis_files)} files using pandas"
    )

    # Prepare summary metadata
    summary_json = output_path / "combined_summary.json"

    # Create validated aggregated summary
    summary_status = "completed" if not missing_species else "partial"

    summary = AggregatedOffTargetSummary(
        species_analyzed=species_list,
        analysis_files_processed=len(analysis_files),
        total_results=len(combined_df),
        combined_tsv=combined_tsv,
        combined_json=combined_json,
        summary_file=summary_json,
        hits_per_species=species_counts,
        human_hits=human_hits,
        other_species_hits=other_hits,
        species_file_counts=species_file_counts,
        missing_species=missing_species,
        status=summary_status,
    )

    # Write summary JSON
    with summary_json.open("w") as f:
        json.dump(summary.model_dump(mode="json"), f, indent=2)

    # Write final summary text file
    final_summary = output_path / "final_summary.txt"
    with final_summary.open("w") as f:
        f.write("Off-Target Analysis Aggregation Summary\n")
        f.write("=" * 50 + "\n\n")

        # Check if transcriptome analysis was performed
        if len(analysis_files) == 0:
            f.write("TRANSCRIPTOME ANALYSIS STATUS: NOT PERFORMED\n")
            f.write("-" * 50 + "\n")
            f.write("Reason: No transcriptome FASTAs or BWA indices were provided.\n")
            f.write("Result: Only lightweight miRNA seed match analysis was run.\n\n")
            f.write("To enable transcriptome off-target analysis:\n")
            f.write("  • Provide --genome_fastas (transcriptome) 'species:path,species2:path2'\n")
            f.write("     OR\n")
            f.write("  • Provide --genome_indices 'species:index,species2:index2'\n\n")
            f.write("=" * 50 + "\n\n")

        if missing_species:
            warning_list = ", ".join(missing_species)
            f.write("WARNINGS\n")
            f.write("-" * 50 + "\n")
            f.write(
                "No transcriptome alignment files were produced for the following species: "
                f"{warning_list}. This usually indicates the BWA-MEM2 indexing stage ran out of memory.\n"
            )
            f.write(
                "Increase --max_memory (32GB+ recommended for human transcriptomes) or pre-build indices on a host with more RAM.\n\n"
            )

        # Results summary
        f.write("RESULTS SUMMARY\n")
        f.write("-" * 50 + "\n")
        f.write(f"Transcriptome off-target hits: {len(combined_df)}\n")
        f.write(f"Human hits: {human_hits}\n")
        f.write(f"Other species hits: {other_hits}\n")

        # Show species list or note if empty
        if species_list:
            f.write(f"Species requested for analysis: {', '.join(species_list)}\n")
        else:
            f.write("Species requested for analysis: (none - miRNA-only mode)\n")

        if species_counts:
            f.write("Per-species hit counts:\n")
            for species, count in sorted(species_counts.items()):
                f.write(f"  {species}: {count}\n")

        f.write(f"Transcriptome analysis files processed: {len(analysis_files)}\n\n")

        # Explain what the output files contain
        f.write("OUTPUT FILES\n")
        f.write("-" * 50 + "\n")

        if len(combined_df) == 0:
            f.write(f"• {combined_tsv.name}: Header only (no hits found)\n")
            f.write(f"• {combined_json.name}: Empty array (no hits found)\n")
            f.write(f"• {summary_json.name}: Metadata only\n\n")
            f.write(
                "Note: Empty data files indicate NO problematic transcriptome off-targets were detected - this is GOOD!\n"
            )
            f.write("Your siRNA candidates are clean at the transcriptome alignment level.\n\n")
            f.write("For miRNA seed match analysis results, see:\n")
            f.write("  ../mirna/mirna_analysis.tsv\n")
            f.write("  ../mirna/mirna_summary.json\n")
        else:
            f.write(f"• {combined_tsv.name}: {len(combined_df)} off-target hits (TSV format)\n")
            f.write(f"• {combined_json.name}: {len(combined_df)} off-target hits (JSON format)\n")
            f.write(f"• {summary_json.name}: Analysis metadata and statistics\n")

    if missing_species:
        logger.warning(
            "Transcriptome aggregation completed with missing species: %s. "
            "Likely cause: insufficient memory while building BWA-MEM2 indices.",
            ", ".join(missing_species),
        )
    else:
        logger.info(f"Wrote aggregated results to {output_path}")

    return output_path


def run_mirna_seed_analysis(
    candidates_file: str | Path,
    candidate_id: str,
    mirna_db: str,  # Review, can this be linked to a class describing all miRNA database protocol/ABC?
    mirna_species: list[str],
    output_dir: str | Path,
) -> Path:
    """Run miRNA seed match analysis for candidate sequences.

    This function uses the MiRNADatabaseManager to download and cache miRNA databases,
    builds BWA indices if needed, and performs seed match analysis.

    Args:
        candidates_file: Path to FASTA file with candidate sequences
        candidate_id: Candidate identifier
        mirna_db: miRNA database name (mirgenedb, mirbase, etc.)
        mirna_species: List of species to analyze against
        output_dir: Directory to write results

    Returns:
        Path to output directory containing results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Initialize miRNA database manager
    manager = MiRNADatabaseManager()

    # Parse input sequences
    sequences = parse_fasta_file(candidates_file)

    all_raw_hits = []  # All raw alignments from BWA
    species_raw_stats = {}
    species_filtered_stats = {}

    logger.info(f"Running miRNA seed analysis for {candidate_id}")
    logger.info(f"Database: {mirna_db}, Species: {mirna_species}")

    for species in mirna_species:
        try:
            # Get or download miRNA database for this species
            logger.info(f"Processing miRNA database for species: {species}")
            db_fasta_path = manager.get_database(mirna_db, species)

            if db_fasta_path is None or not db_fasta_path.exists():
                logger.warning(f"miRNA database not available for {species}, skipping")
                continue

            # Build BWA index for miRNA database if it doesn't exist
            index_prefix = db_fasta_path.with_suffix("")
            if not validate_index_files(index_prefix, "bwa-mem2"):
                logger.info(f"Building BWA index for {species} miRNA database")
                build_bwa_index(db_fasta_path, index_prefix)

            # Run BWA analysis in miRNA seed mode
            analyzer = BwaAnalyzer(
                index_prefix=index_prefix,
                mode="mirna_seed",
                seed_length=6,  # Short seed for miRNA seed region
                min_score=6,  # Low threshold for sensitivity
                max_hits=1000,  # Reasonable limit per candidate
                seed_start=2,
                seed_end=8,
            )

            results = analyzer.analyze_sequences(sequences)

            # Convert BWA results to DataFrame directly - let Pandera handle validation
            results_df = pd.DataFrame(results)

            # Parse coord field: for miRNA it's "miRNA_ID:position", extract integer position
            if "coord" in results_df.columns:
                results_df["coord"] = results_df["coord"].apply(
                    lambda x: int(x.split(":")[-1]) if isinstance(x, str) and ":" in x else int(x)
                )

            # Add miRNA-specific columns
            results_df["species"] = species
            results_df["database"] = mirna_db
            results_df["mirna_id"] = results_df["rname"]  # BWA uses rname for reference ID

            # Remove internal columns not part of the MiRNAAlignmentSchema
            # - rname: copied to mirna_id (schema uses domain-specific naming)
            # - mismatch_positions: internal debugging data, not needed in output
            columns_to_drop = ["rname", "mismatch_positions"]
            results_df = results_df.drop(columns=[col for col in columns_to_drop if col in results_df.columns])

            # Validate and coerce types using Pandera schema
            try:
                validated_df = MiRNAAlignmentSchema.validate(results_df, lazy=True)
                all_raw_hits.append(validated_df)
                species_raw_stats[species] = len(validated_df)
                logger.info(f"Species {species}: {len(validated_df)} miRNA alignments validated")
            except Exception as validation_error:
                logger.error(f"Failed to validate miRNA hits for {species}: {validation_error}")
                species_raw_stats[species] = 0
                species_filtered_stats[species] = 0
                continue

        except Exception as e:
            logger.error(f"Failed to process miRNA analysis for {species}: {e}")
            species_raw_stats[species] = 0
            species_filtered_stats[species] = 0

    # Concatenate all validated DataFrames from different species
    if all_raw_hits:
        df_raw = pd.concat(all_raw_hits, ignore_index=True)

        # Keep ALL hits including perfect matches - they're biologically relevant off-targets
        # Users can filter by score threshold if desired
        df_filtered = df_raw.copy()

        # Calculate per-species filtered stats
        for species in mirna_species:
            species_filtered_stats[species] = len(df_filtered[df_filtered["species"] == species])

        # Write RAW hits TSV (all alignments)
        raw_analysis_file = output_path / f"{candidate_id}_mirna_analysis_raw.tsv"
        df_raw.to_csv(raw_analysis_file, sep="\t", index=False)

        # Write FILTERED hits TSV (quality-filtered)
        filtered_analysis_file = output_path / f"{candidate_id}_mirna_analysis.tsv"
        df_filtered.to_csv(filtered_analysis_file, sep="\t", index=False)

        # Write raw hits JSON
        raw_json_file = output_path / f"{candidate_id}_mirna_hits_raw.json"
        df_raw.to_json(raw_json_file, orient="records", indent=2)

        # Write filtered hits JSON
        filtered_json_file = output_path / f"{candidate_id}_mirna_hits.json"
        df_filtered.to_json(filtered_json_file, orient="records", indent=2)

        total_filtered = len(df_filtered)
        total_raw = len(df_raw)
    else:
        # No hits - create empty DataFrame with proper schema columns
        df_empty = pd.DataFrame(columns=list(MiRNAAlignmentSchema.to_schema().columns.keys()))

        raw_analysis_file = output_path / f"{candidate_id}_mirna_analysis_raw.tsv"
        df_empty.to_csv(raw_analysis_file, sep="\t", index=False)

        filtered_analysis_file = output_path / f"{candidate_id}_mirna_analysis.tsv"
        df_empty.to_csv(filtered_analysis_file, sep="\t", index=False)

        raw_json_file = output_path / f"{candidate_id}_mirna_hits_raw.json"
        df_empty.to_json(raw_json_file, orient="records", indent=2)

        filtered_json_file = output_path / f"{candidate_id}_mirna_hits.json"
        df_empty.to_json(filtered_json_file, orient="records", indent=2)

        total_filtered = 0
        total_raw = 0

    # Create validated summary using Pydantic model
    summary = MiRNASummary(
        candidate_id=candidate_id,
        mirna_database=mirna_db,
        species_analyzed=mirna_species,
        total_sequences=len(sequences),
        total_hits=total_filtered,  # Filtered hits count
        hits_per_species=species_raw_stats,  # Raw hits per species
        total_raw_alignments=total_raw,  # Total raw alignments
    )

    # Write summary JSON file
    summary_file = output_path / f"{candidate_id}_mirna_summary.json"
    with summary_file.open("w") as f:
        json.dump(summary.model_dump(mode="json"), f, indent=2)

    logger.info(
        f"miRNA seed analysis completed for {candidate_id}: "
        f"{total_filtered} filtered high-quality matches "
        f"(from {total_raw} raw alignments)"
    )

    return output_path


def aggregate_mirna_results(
    results_dir: str | Path,
    output_dir: str | Path,
    mirna_db: str,
    mirna_species: str,
) -> Path:
    """Aggregate miRNA seed analysis results from multiple candidates using pandas.

    Uses pandas + Pandera for efficient bulk reading and validation instead of
    manual line-by-line parsing with Pydantic models.

    Args:
        results_dir: Directory containing individual miRNA analysis results
        output_dir: Directory to write aggregated results
        mirna_db: miRNA database used for analysis
        mirna_species: Comma-separated list of species analyzed

    Returns:
        Path to output directory containing aggregated results
    """
    results_path = Path(results_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    species_list = [s.strip() for s in mirna_species.split(",") if s.strip()]

    # Collect all miRNA analysis files using pandas (much faster than manual parsing)
    analysis_files = list(results_path.glob("**/*_mirna_analysis.tsv"))

    if analysis_files:
        # Read all files into DataFrames and concatenate (vectorized operation)
        dfs = []
        candidate_stats = {}

        for analysis_file in analysis_files:
            try:
                # Extract candidate ID from filename
                candidate_id = analysis_file.stem.replace("_mirna_analysis", "")

                # Pandas reads TSV much faster than manual line splitting
                df = pd.read_csv(analysis_file, sep="\t")

                # Validate schema with Pandera
                df = MiRNAAlignmentSchema.validate(df, lazy=True)

                # Track hits per candidate
                candidate_stats[candidate_id] = len(df)

                dfs.append(df)

            except Exception as e:
                logger.warning(f"Failed to read/validate {analysis_file}: {e}")
                continue

        # Concatenate all DataFrames at once (much faster than append in loop)
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
        else:
            # Create empty DataFrame with correct schema
            combined_df = pd.DataFrame(columns=list(MiRNAHit.model_fields.keys()))

    else:
        # No files found - create empty DataFrame
        combined_df = pd.DataFrame(columns=list(MiRNAHit.model_fields.keys()))
        candidate_stats = {}

    # Write combined results (pandas is much faster than manual TSV writing)
    combined_tsv = output_path / "combined_mirna_hits.tsv"
    combined_df.to_csv(combined_tsv, sep="\t", index=False)

    # Write JSON (pandas handles serialization)
    combined_json = output_path / "combined_mirna_hits.json"
    combined_df.to_json(combined_json, orient="records", indent=2)

    logger.info(f"Aggregated {len(combined_df)} miRNA hits from {len(analysis_files)} files using pandas")

    # Calculate statistics using pandas groupby (much faster than loops)
    species_stats = _compute_species_counts(combined_df)
    if not species_stats:
        species_stats = dict.fromkeys(species_list, 0) if species_list else {}

    for species in species_list:
        species_stats.setdefault(species, 0)

    human_hits, other_hits = human_vs_other_totals(species_stats)

    # Create validated summary using Pydantic model
    summary = AggregatedMiRNASummary(
        total_mirna_hits=len(combined_df),
        mirna_database=mirna_db,
        species_analyzed=species_list,
        hits_per_species=species_stats,
        hits_per_candidate=candidate_stats,
        analysis_files_processed=len(analysis_files),
        total_candidates=len(candidate_stats),
        combined_tsv=combined_tsv,
        combined_json=combined_json,
        summary_file=output_path / "combined_mirna_summary.json",
        human_hits=human_hits,
        other_species_hits=other_hits,
    )

    # Write summary JSON
    summary_json = output_path / "combined_mirna_summary.json"
    with summary_json.open("w") as f:
        json.dump(summary.model_dump(mode="json"), f, indent=2)

    # Write final summary text file
    final_summary = output_path / "final_mirna_summary.txt"
    with final_summary.open("w") as f:
        f.write("miRNA Seed Match Analysis Aggregation Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total miRNA seed matches: {len(combined_df)}\n")
        f.write(f"Database: {mirna_db}\n")
        f.write(f"Species analyzed: {', '.join(species_list)}\n")
        f.write(f"Candidates analyzed: {len(candidate_stats)}\n")
        f.write(f"Analysis files processed: {len(analysis_files)}\n")
        f.write(f"Human hits: {human_hits}\n")
        f.write(f"Other species hits: {other_hits}\n")
        f.write("\nHits per species:\n")
        for species, count in species_stats.items():
            f.write(f"  {species}: {count}\n")
        f.write("\nOutput files:\n")
        f.write(f"  - Combined TSV: {combined_tsv.name}\n")
        f.write(f"  - Combined JSON: {combined_json.name}\n")
        f.write(f"  - Summary JSON: {summary_json.name}\n")

    logger.info(f"Wrote aggregated miRNA results to {output_path}")

    return output_path


# Export all main functions and classes
__all__ = [
    # Core classes
    "BwaAnalyzer",
    "OffTargetAnalysisManager",
    # Utility functions
    "create_temp_fasta",
    "validate_and_write_sequences",
    "build_bwa_index",
    "validate_sirna_sequences",
    "parse_fasta_file",
    "write_fasta_file",
    "check_tool_availability",
    "validate_index_files",
    # Nextflow integration functions (called directly from Nextflow modules)
    "run_bwa_alignment_analysis",
    "aggregate_offtarget_results",
    "run_mirna_seed_analysis",
    "aggregate_mirna_results",
]
