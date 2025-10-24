"""Off-target analysis for siRNA design.

This module provides comprehensive off-target analysis functionality for siRNA design,
including both miRNA seed match analysis and transcriptome off-target detection.
Uses BWA-MEM2 for both short and long sequence alignments.
Optimized for both standalone use and parallelized Nextflow workflows.
"""

import json
import shutil
import subprocess  # nosec B404
import tempfile
from pathlib import Path
from typing import Any, Optional, Union

from sirnaforge.data.base import FastaUtils
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


# =============================================================================
# Core Analyzer Classes
# =============================================================================


def _get_executable_path(tool_name: str) -> Optional[str]:
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
        index_prefix: Union[str, Path],
        mode: str = "transcriptome",  # "transcriptome" or "mirna_seed"
        seed_length: int = 12,
        min_score: int = 15,
        max_hits: int = 10000,
        seed_start: int = 2,
        seed_end: int = 8,
    ):
        """
        Initialize BWA-MEM2 analyzer.

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
        """
        Run BWA-MEM2 analysis on sequences.

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
        """Calculate off-target score based on mismatch positions."""
        score = 0.0

        for pos in mismatch_positions:
            if self.seed_start <= pos <= self.seed_end:
                score += 5.0  # Seed region mismatches are more critical
            else:
                score += 1.0  # Non-seed mismatches

        if len(mismatch_positions) == 0:
            score += 10.0  # Perfect matches get high penalty (bad for off-target)

        return score

    def _filter_and_rank(self, results: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Filter and rank results by off-target score."""
        results.sort(key=lambda x: (x["offtarget_score"], -x.get("as_score", 0)))
        return results


class OffTargetAnalysisManager:
    """Manager class for comprehensive off-target analysis using BWA-MEM2."""

    def __init__(
        self,
        species: str,
        transcriptome_path: Optional[Union[str, Path]] = None,
        mirna_path: Optional[Union[str, Path]] = None,
        transcriptome_index: Optional[Union[str, Path]] = None,
        mirna_index: Optional[Union[str, Path]] = None,
    ):
        self.species = species
        self.transcriptome_path = Path(transcriptome_path) if transcriptome_path is not None else None
        self.mirna_path = Path(mirna_path) if mirna_path is not None else None
        self.transcriptome_index = Path(transcriptome_index) if transcriptome_index is not None else None
        self.mirna_index = Path(mirna_index) if mirna_index is not None else None

    def analyze_mirna_off_targets(
        self,
        sequences: Union[dict[str, str], str, Path],
        output_prefix: Union[str, Path],
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
        sequences: Union[dict[str, str], str, Path],
        output_prefix: Union[str, Path],
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

    def _write_mirna_results(
        self, results: list[dict[str, Any]], tsv_path: Union[str, Path], json_path: Union[str, Path]
    ) -> None:
        """Write miRNA analysis results."""
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

    def _write_transcriptome_results(
        self, results: list[dict[str, Any]], tsv_path: Union[str, Path], json_path: Union[str, Path]
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


def build_bwa_index(fasta_file: Union[str, Path], index_prefix: Union[str, Path]) -> Path:
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


def parse_fasta_file(fasta_file: Union[str, Path]) -> dict[str, str]:
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


def validate_index_files(index_prefix: Union[str, Path], tool: str = "bwa") -> bool:
    """Validate that index files exist for given tool."""
    index_path = Path(index_prefix)

    if tool in ("bwa", "bwa-mem2"):
        required_extensions = [".amb", ".ann", ".bwt.2bit.64", ".pac"]
    else:
        logger.warning(f"Unknown tool for index validation: {tool}")
        return False

    for ext in required_extensions:
        if not (index_path.parent / f"{index_path.name}{ext}").exists():
            logger.debug(f"Missing index file: {index_path.name}{ext}")
            return False

    return True


# =============================================================================
# Nextflow Integration Functions
# =============================================================================


def run_mirna_analysis_for_nextflow(
    species: str,
    sequences_file: str,
    mirna_index: Union[str, Path],
    output_prefix: Union[str, Path],
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
    transcriptome_index: Union[str, Path],
    output_prefix: Union[str, Path],
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
    output_prefix: Union[str, Path],
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

        # Write results
        tsv_path = output_root.parent / f"{output_root.name}.tsv"
        json_path = output_root.parent / f"{output_root.name}.json"
        summary_path = output_root.parent / f"{output_root.name}_summary.txt"

        # Write TSV
        with tsv_path.open("w") as f:
            f.write("qname\tqseq\trname\tcoord\tstrand\tcigar\tmapq\tas_score\tnm\tseed_mismatches\tofftarget_score\n")
            for result in results:
                f.write(
                    f"{result['qname']}\t{result['qseq']}\t{result['rname']}\t"
                    f"{result['coord']}\t{result['strand']}\t{result['cigar']}\t"
                    f"{result['mapq']}\t{result.get('as_score', 'NA')}\t{result['nm']}\t"
                    f"{result['seed_mismatches']}\t{result['offtarget_score']}\n"
                )

        # Write JSON
        with json_path.open("w") as f:
            json.dump(results, f, indent=2)

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
    # Nextflow functions
    "run_mirna_analysis_for_nextflow",
    "run_transcriptome_analysis_for_nextflow",
    "run_comprehensive_offtarget_analysis",
]
