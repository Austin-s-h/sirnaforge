"""Off-target analysis for siRNA design.

This module provides comprehensive off-target analysis functionality for siRNA design,
including both miRNA seed match analysis and transcriptome off-target detection.
Optimized for both standalone use and parallelized Nextflow workflows.
"""

import json
import subprocess
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


class BowtieAnalyzer:
    """Bowtie-based analyzer for miRNA seed matches."""

    def __init__(
        self,
        index_prefix: str,
        mismatches: int = 0,
        seed_start: int = 2,
        seed_end: int = 8,
    ):
        """
        Initialize Bowtie analyzer for miRNA seed search.

        Args:
            index_prefix: Path to Bowtie index
            mismatches: Number of allowed mismatches
            seed_start: Seed region start (1-based)
            seed_end: Seed region end (1-based)
        """
        self.index_prefix = index_prefix
        self.mismatches = mismatches
        self.seed_start = seed_start
        self.seed_end = seed_end

    def analyze_sequences(self, sequences: dict[str, str]) -> list[dict[str, Any]]:
        """
        Run Bowtie analysis on sequences.

        Args:
            sequences: Dictionary of sequence name -> sequence

        Returns:
            List of alignment dictionaries
        """
        results = []
        temp_fasta_path = create_temp_fasta(sequences)

        try:
            cmd = ["bowtie", "-v", str(self.mismatches), "-a", "--sam", "--quiet", self.index_prefix, temp_fasta_path]
            logger.info(f"Running Bowtie: {' '.join(cmd)}")

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, check=True)
            results = self._parse_sam_output(result.stdout, sequences)
            logger.info(f"Bowtie analysis completed: {len(results)} hits found")

        except subprocess.CalledProcessError as e:
            logger.error(f"Bowtie failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            logger.error("Bowtie timed out")
        finally:
            Path(temp_fasta_path).unlink(missing_ok=True)

        return results

    def _parse_sam_output(self, sam_output: str, sequences: dict[str, str]) -> list[dict[str, Any]]:
        """Parse SAM output from Bowtie."""
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

            if flag & 4:  # Skip unmapped
                continue

            strand = "-" if (flag & 16) else "+"
            coord = f"{rname}:{pos}"

            # Count mismatches (for exact Bowtie -v mode)
            total_mismatches = self.mismatches
            seed_mismatches = 0  # Calculate based on mismatch positions if available

            # Calculate score (lower is better for miRNA seed matches)
            score = total_mismatches + (seed_mismatches * 3)

            result = {
                "qname": qname,
                "qseq": sequences.get(qname, ""),
                "rname": rname,
                "coord": coord,
                "strand": strand,
                "mismatches": total_mismatches,
                "seed_mismatches": seed_mismatches,
                "score": score,
            }
            results.append(result)

        return results


class BwaAnalyzer:
    """BWA-MEM2 based analyzer for transcriptome off-target search."""

    def __init__(
        self,
        index_prefix: str,
        seed_length: int = 12,
        min_score: int = 15,
        max_hits: int = 10000,
        seed_start: int = 2,
        seed_end: int = 8,
    ):
        """
        Initialize BWA-MEM2 analyzer for transcriptome search.

        Args:
            index_prefix: Path to BWA index
            seed_length: BWA seed length parameter
            min_score: Minimum alignment score
            max_hits: Maximum hits to return
            seed_start: Seed region start (1-based)
            seed_end: Seed region end (1-based)
        """
        self.index_prefix = index_prefix
        self.seed_length = seed_length
        self.min_score = min_score
        self.max_hits = max_hits
        self.seed_start = seed_start
        self.seed_end = seed_end

    def analyze_sequences(self, sequences: dict[str, str]) -> list[dict[str, Any]]:
        """
        Run BWA-MEM2 analysis on sequences.

        Args:
            sequences: Dictionary of sequence name -> sequence

        Returns:
            List of alignment dictionaries
        """
        results = []
        temp_fasta_path = create_temp_fasta(sequences)

        try:
            cmd = [
                "bwa-mem2",
                "mem",
                "-a",
                "-k",
                str(self.seed_length),
                "-T",
                str(self.min_score),
                "-v",
                "1",
                self.index_prefix,
                temp_fasta_path,
            ]

            logger.info(f"Running BWA-MEM2: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=None, check=True)
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

    def _parse_sam_output(self, sam_output: str, sequences: dict[str, str]) -> list[dict[str, Any]]:
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
                "qseq": sequences.get(qname, ""),
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
    """Manager class for comprehensive off-target analysis."""

    def __init__(
        self,
        species: str,
        transcriptome_path: Optional[Union[str, Path]] = None,
        mirna_path: Optional[Union[str, Path]] = None,
        transcriptome_index: Optional[str] = None,
        mirna_index: Optional[str] = None,
    ):
        self.species = species
        self.transcriptome_path = transcriptome_path
        self.mirna_path = mirna_path
        self.transcriptome_index = transcriptome_index
        self.mirna_index = mirna_index

    def analyze_mirna_off_targets(
        self,
        sequences: Union[dict[str, str], str, Path],
        output_prefix: str,
    ) -> tuple[str, str]:
        """Analyze miRNA off-targets using Bowtie."""
        if not self.mirna_index:
            raise ValueError("miRNA index not provided")

        if isinstance(sequences, (str, Path)):
            sequences = parse_fasta_file(str(sequences))

        analyzer = BowtieAnalyzer(self.mirna_index)
        results = analyzer.analyze_sequences(sequences)

        tsv_path = f"{output_prefix}_mirna_hits.tsv"
        json_path = f"{output_prefix}_mirna_hits.json"

        self._write_mirna_results(results, tsv_path, json_path)
        return tsv_path, json_path

    def analyze_transcriptome_off_targets(
        self,
        sequences: Union[dict[str, str], str, Path],
        output_prefix: str,
    ) -> tuple[str, str]:
        """Analyze transcriptome off-targets using BWA-MEM2."""
        if not self.transcriptome_index:
            raise ValueError("Transcriptome index not provided")

        if isinstance(sequences, (str, Path)):
            sequences = parse_fasta_file(str(sequences))

        analyzer = BwaAnalyzer(self.transcriptome_index)
        results = analyzer.analyze_sequences(sequences)

        tsv_path = f"{output_prefix}_transcriptome_hits.tsv"
        json_path = f"{output_prefix}_transcriptome_hits.json"

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
            mirna_analyzer: BowtieAnalyzer = BowtieAnalyzer(self.mirna_index)
            results["mirna_hits"] = mirna_analyzer.analyze_sequences(sequences)

        if self.transcriptome_index:
            transcriptome_analyzer = BwaAnalyzer(self.transcriptome_index)
            results["transcriptome_hits"] = transcriptome_analyzer.analyze_sequences(sequences)

        return results

    def _write_mirna_results(self, results: list[dict[str, Any]], tsv_path: str, json_path: str) -> None:
        """Write miRNA analysis results."""
        # Write TSV
        with Path(tsv_path).open("w") as f:
            f.write("qname\tqseq\trname\tcoord\tstrand\tmismatches\tseed_mismatches\tscore\n")
            for result in results:
                f.write(
                    f"{result['qname']}\t{result['qseq']}\t{result['rname']}\t"
                    f"{result['coord']}\t{result['strand']}\t{result['mismatches']}\t"
                    f"{result['seed_mismatches']}\t{result['score']}\n"
                )

        # Write JSON
        with Path(json_path).open("w") as f:
            json.dump(results, f, indent=2)

    def _write_transcriptome_results(self, results: list[dict[str, Any]], tsv_path: str, json_path: str) -> None:
        """Write transcriptome analysis results."""
        # Write TSV
        with Path(tsv_path).open("w") as f:
            f.write(
                "qname\tqseq\trname\tcoord\tstrand\tcigar\tmapq\tas_score\tnm\t" "seed_mismatches\tofftarget_score\n"
            )
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


def build_bowtie_index(fasta_file: str, index_prefix: str) -> str:
    """Build Bowtie index for miRNA analysis."""
    logger.info(f"Building Bowtie index from {fasta_file} with prefix {index_prefix}")

    if not Path(fasta_file).exists():
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_file}")

    Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)

    cmd = ["bowtie-build", fasta_file, index_prefix]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=3600)
        logger.info(f"Bowtie index built successfully: {index_prefix}")
        return index_prefix
    except subprocess.CalledProcessError as e:
        logger.error(f"Bowtie index build failed: {e.stderr}")
        raise
    except subprocess.TimeoutExpired:
        logger.error("Bowtie index build timed out")
        raise


def build_bwa_index(fasta_file: str, index_prefix: str) -> str:
    """Build BWA-MEM2 index for transcriptome analysis."""
    logger.info(f"Building BWA-MEM2 index from {fasta_file} with prefix {index_prefix}")

    if not Path(fasta_file).exists():
        raise FileNotFoundError(f"Input FASTA file not found: {fasta_file}")

    Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)

    cmd = ["bwa-mem2", "index", "-p", index_prefix, fasta_file]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=7200)
        logger.info(f"BWA-MEM2 index built successfully: {index_prefix}")
        return index_prefix
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


def parse_fasta_file(fasta_file: str) -> dict[str, str]:
    """Parse FASTA file using existing FastaUtils."""
    return FastaUtils.parse_fasta_to_dict(fasta_file)


def write_fasta_file(sequences: dict[str, str], output_file: str) -> None:
    """Write sequences to FASTA file using existing FastaUtils."""
    FastaUtils.write_dict_to_fasta(sequences, output_file)


def check_tool_availability(tool: str) -> bool:
    """Check if external tool is available."""
    try:
        result = subprocess.run([tool, "--help"], capture_output=True, check=False, timeout=10)
        return result.returncode in {0, 1}
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def validate_index_files(index_prefix: str, tool: str) -> bool:
    """Validate that index files exist for given tool."""
    index_path = Path(index_prefix)

    if tool == "bowtie":
        required_extensions = [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]
    elif tool == "bwa":
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
    mirna_index: str,
    output_prefix: str,
) -> tuple[str, str, str]:
    """Nextflow-compatible function for miRNA analysis."""
    manager = OffTargetAnalysisManager(species=species, mirna_index=mirna_index)

    try:
        tsv_path, json_path = manager.analyze_mirna_off_targets(sequences_file, output_prefix)

        # Create summary
        summary_path = f"{output_prefix}_mirna_summary.txt"
        with Path(summary_path).open("w") as f:
            with Path(json_path).open() as jf:
                results = json.load(jf)
            f.write(f"Species: {species}\n")
            f.write(f"Total miRNA hits: {len(results)}\n")
            f.write("Analysis completed successfully\n")

        return tsv_path, json_path, summary_path

    except Exception as e:
        error_summary = f"{output_prefix}_mirna_error.txt"
        with Path(error_summary).open("w") as f:
            f.write(f"miRNA analysis failed: {str(e)}\n")
        return "", "", error_summary


def run_transcriptome_analysis_for_nextflow(
    species: str,
    sequences_file: str,
    transcriptome_index: str,
    output_prefix: str,
) -> tuple[str, str, str]:
    """Nextflow-compatible function for transcriptome analysis."""
    manager = OffTargetAnalysisManager(species=species, transcriptome_index=transcriptome_index)

    try:
        tsv_path, json_path = manager.analyze_transcriptome_off_targets(sequences_file, output_prefix)

        # Create summary
        summary_path = f"{output_prefix}_transcriptome_summary.txt"
        with Path(summary_path).open("w") as f:
            with Path(json_path).open() as jf:
                results = json.load(jf)
            f.write(f"Species: {species}\n")
            f.write(f"Total transcriptome hits: {len(results)}\n")
            f.write("Analysis completed successfully\n")

        return tsv_path, json_path, summary_path

    except Exception as e:
        error_summary = f"{output_prefix}_transcriptome_error.txt"
        with Path(error_summary).open("w") as f:
            f.write(f"Transcriptome analysis failed: {str(e)}\n")
        return "", "", error_summary


def run_comprehensive_offtarget_analysis(
    species: str,
    sequences_file: str,
    index_path: str,
    output_prefix: str,
    bwa_k: int = 12,
    bwa_T: int = 15,
    max_hits: int = 10000,
    seed_start: int = 2,
    seed_end: int = 8,
) -> tuple[str, str, str]:
    """Run comprehensive off-target analysis for Nextflow integration."""
    try:
        sequences = parse_fasta_file(sequences_file)

        # Use BWA analyzer for comprehensive analysis
        analyzer = BwaAnalyzer(
            index_prefix=index_path,
            seed_length=bwa_k,
            min_score=bwa_T,
            max_hits=max_hits,
            seed_start=seed_start,
            seed_end=seed_end,
        )

        results = analyzer.analyze_sequences(sequences)

        # Write results
        tsv_path = f"{output_prefix}.tsv"
        json_path = f"{output_prefix}.json"
        summary_path = f"{output_prefix}_summary.txt"

        # Write TSV
        with Path(tsv_path).open("w") as f:
            f.write(
                "qname\tqseq\trname\tcoord\tstrand\tcigar\tmapq\tas_score\tnm\t" "seed_mismatches\tofftarget_score\n"
            )
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

        # Write summary
        with Path(summary_path).open("w") as f:
            f.write(f"Species: {species}\n")
            f.write(f"Total sequences analyzed: {len(sequences)}\n")
            f.write(f"Total off-target hits: {len(results)}\n")
            f.write(f"Analysis parameters: bwa_k={bwa_k}, bwa_T={bwa_T}, max_hits={max_hits}\n")
            f.write(f"Seed region: {seed_start}-{seed_end}\n")
            f.write("Analysis completed successfully\n")

        return tsv_path, json_path, summary_path

    except Exception as e:
        error_summary = f"{output_prefix}_error.txt"
        with Path(error_summary).open("w") as f:
            f.write(f"Comprehensive off-target analysis failed: {str(e)}\n")
        return "", "", error_summary


# Export all main functions and classes
__all__ = [
    # Core classes
    "BowtieAnalyzer",
    "BwaAnalyzer",
    "OffTargetAnalysisManager",
    # Utility functions
    "create_temp_fasta",
    "validate_and_write_sequences",
    "build_bowtie_index",
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
