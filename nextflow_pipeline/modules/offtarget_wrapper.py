#!/usr/bin/env python3
"""
Enhanced off-target wrapper supporting multiple alignment tools and comprehensive analysis.

This module provides comprehensive off-target analysis for siRNA candidates using:
- Bowtie for fast seed-based alignment
- BWA-MEM2 for sensitive full-length alignment
- Combined analysis with scoring and filtering

Usage:
    python offtarget_wrapper.py --mode combined --index-prefix /path/to/genome --queries input.fasta --out-prefix results

Features:
- Seed region-aware mismatch scoring
- JSON and TSV output formats
- Configurable alignment parameters
- Support for multiple genome species
- Comprehensive hit filtering and ranking

Author: siRNA Design Toolkit
Version: 1.0.0
"""

import argparse
import json
import logging
import subprocess
import sys
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


class OffTargetAnalyzer:
    """Main class for off-target analysis of siRNA candidates."""

    def __init__(
        self,
        index_prefix: str,
        mode: str = "combined",
        bwa_k: int = 12,
        bwa_T: int = 15,
        max_hits: int = 10000,
        seed_start: int = 2,
        seed_end: int = 8,
    ):
        """
        Initialize the off-target analyzer.

        Args:
            index_prefix: Path to genome index
            mode: Analysis mode ('bowtie-only', 'bwa-mem2-only', 'combined')
            bwa_k: BWA seed length
            bwa_T: BWA minimum score threshold
            max_hits: Maximum hits per query
            seed_start: Seed region start (1-based)
            seed_end: Seed region end (1-based)
        """
        self.index_prefix = index_prefix
        self.mode = mode
        self.bwa_k = bwa_k
        self.bwa_T = bwa_T
        self.max_hits = max_hits
        self.seed_start = seed_start
        self.seed_end = seed_end

        # Validate inputs
        self._validate_parameters()

    def _validate_parameters(self):
        """Validate analysis parameters."""
        if not Path(self.index_prefix).exists():
            logger.warning(f"Index path {self.index_prefix} does not exist")

        if self.mode not in ["bowtie-only", "bwa-mem2-only", "combined"]:
            raise ValueError(f"Invalid mode: {self.mode}")

        if self.seed_start >= self.seed_end:
            raise ValueError("seed_start must be less than seed_end")

    def run_bowtie(self, query_fasta: str, v_mismatches: int = 0) -> str:
        """
        Run Bowtie alignment for seed-based search.

        Args:
            query_fasta: Path to query sequences
            v_mismatches: Number of allowed mismatches

        Returns:
            SAM format alignment output
        """
        cmd = ["bowtie", "-v", str(v_mismatches), "-a", "--sam", "--quiet", self.index_prefix, query_fasta]

        logger.info(f"Running Bowtie: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour timeout
                check=True,
            )
            logger.info("Bowtie completed successfully")
            return result.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"Bowtie failed with return code {e.returncode}")
            logger.error(f"STDERR: {e.stderr}")
            return ""

        except subprocess.TimeoutExpired:
            logger.error("Bowtie timed out after 1 hour")
            return ""

    def run_bwa_mem2(self, query_fasta: str) -> str:
        """
        Run BWA-MEM2 alignment for sensitive search.

        Args:
            query_fasta: Path to query sequences

        Returns:
            SAM format alignment output
        """
        cmd = [
            "bwa-mem2",
            "mem",
            "-a",  # Report all alignments
            "-k",
            str(self.bwa_k),
            "-T",
            str(self.bwa_T),
            "-v",
            "1",  # Verbose level 1
            self.index_prefix,
            query_fasta,
        ]

        logger.info(f"Running BWA-MEM2: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200,  # 2 hour timeout
                check=True,
            )
            logger.info("BWA-MEM2 completed successfully")
            return result.stdout

        except subprocess.CalledProcessError as e:
            logger.error(f"BWA-MEM2 failed with return code {e.returncode}")
            logger.error(f"STDERR: {e.stderr}")
            return ""

        except subprocess.TimeoutExpired:
            logger.error("BWA-MEM2 timed out after 2 hours")
            return ""

    def parse_sam_alignments(self, sam_text: str, query_seqs: dict[str, str]) -> list[dict[str, Any]]:
        """
        Parse SAM format alignments.

        Args:
            sam_text: SAM format text
            query_seqs: Dictionary of query sequences by name

        Returns:
            List of alignment dictionaries
        """
        alignments = []

        for line in sam_text.splitlines():
            if not line or line.startswith("@"):
                continue

            try:
                parts = line.split("\t")
                if len(parts) < 11:
                    continue

                qname = parts[0]
                flag = int(parts[1])
                rname = parts[2]
                pos = int(parts[3])
                mapq = int(parts[4]) if parts[4] != "*" else 0
                cigar = parts[5]
                seq = parts[9]

                # Skip unmapped reads
                if flag & 4:
                    continue

                # Parse optional tags
                tags = {}
                for tag in parts[11:]:
                    if ":" in tag:
                        tag_parts = tag.split(":", 2)
                        if len(tag_parts) == 3:
                            tags[tag_parts[0]] = tag_parts[2]

                # Extract alignment information
                strand = "-" if (flag & 16) else "+"
                nm = int(tags.get("NM", 0))
                as_score = int(tags.get("AS", 0)) if "AS" in tags else None

                # Parse mismatch positions from MD tag
                mismatch_positions = []
                if "MD" in tags:
                    mismatch_positions = self._parse_md_tag(tags["MD"])

                alignment = {
                    "qname": qname,
                    "qseq": query_seqs.get(qname, seq),
                    "rname": rname,
                    "pos": pos,
                    "strand": strand,
                    "cigar": cigar,
                    "mapq": mapq,
                    "NM": nm,
                    "AS": as_score,
                    "mismatch_positions": mismatch_positions,
                    "raw_tags": tags,
                }

                alignments.append(alignment)

            except (ValueError, IndexError) as e:
                logger.warning(f"Error parsing SAM line: {e}")
                continue

        return alignments

    def _parse_md_tag(self, md_tag: str) -> list[int]:
        """
        Parse MD tag to extract mismatch positions.

        Args:
            md_tag: MD tag value

        Returns:
            List of 1-based mismatch positions
        """
        positions = []
        read_pos = 1
        i = 0

        while i < len(md_tag):
            if md_tag[i].isdigit():
                # Count matching bases
                num_str = ""
                while i < len(md_tag) and md_tag[i].isdigit():
                    num_str += md_tag[i]
                    i += 1
                read_pos += int(num_str)

            elif md_tag[i] == "^":
                # Deletion in reference - skip deletion bases
                i += 1
                while i < len(md_tag) and md_tag[i].isalpha():
                    i += 1

            elif md_tag[i].isalpha():
                # Mismatch
                positions.append(read_pos)
                read_pos += 1
                i += 1

            else:
                i += 1

        return positions

    def load_query_sequences(self, fasta_path: str) -> dict[str, str]:
        """
        Load query sequences from FASTA file.

        Args:
            fasta_path: Path to FASTA file

        Returns:
            Dictionary mapping sequence names to sequences
        """
        sequences = {}
        current_name = None
        current_seq = []

        try:
            with Path(fasta_path).open() as f:
                for file_line in f:
                    line = file_line.strip()
                    if not line:
                        continue

                    if line.startswith(">"):
                        if current_name:
                            sequences[current_name] = "".join(current_seq).upper().replace("U", "T")
                        current_name = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)

                # Add last sequence
                if current_name:
                    sequences[current_name] = "".join(current_seq).upper().replace("U", "T")

        except OSError as e:
            logger.error(f"Error reading FASTA file {fasta_path}: {e}")
            return {}

        logger.info(f"Loaded {len(sequences)} query sequences")
        return sequences

    def score_alignment(self, alignment: dict[str, Any]) -> float:
        """
        Calculate off-target score for an alignment.

        Args:
            alignment: Alignment dictionary

        Returns:
            Off-target score (higher = more problematic)
        """
        mismatch_positions = alignment.get("mismatch_positions", [])

        score = 0.0
        for pos in mismatch_positions:
            if self.seed_start <= pos <= self.seed_end:
                # Seed region mismatches are more critical
                score += 5.0
            else:
                # Non-seed mismatches
                score += 1.0

        # Bonus penalty for perfect matches
        if len(mismatch_positions) == 0:
            score += 10.0

        return score

    def filter_alignments(self, alignments: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """
        Filter and rank alignments by off-target score.

        Args:
            alignments: List of alignment dictionaries

        Returns:
            Filtered and sorted alignments
        """
        # Add scores to alignments
        for alignment in alignments:
            alignment["offtarget_score"] = self.score_alignment(alignment)

        # Group by query
        by_query = defaultdict(list)
        for alignment in alignments:
            by_query[alignment["qname"]].append(alignment)

        # Filter top hits per query
        filtered = []
        for _qname, query_alignments in by_query.items():
            # Sort by score (ascending - lower is better for matches)
            sorted_alignments = sorted(query_alignments, key=lambda x: (x["offtarget_score"], -x.get("AS", 0)))

            # Take top hits
            filtered.extend(sorted_alignments[: self.max_hits])

        logger.info(f"Filtered to {len(filtered)} alignments")
        return filtered

    def analyze(self, query_fasta: str, output_prefix: str) -> tuple[str, str]:
        """
        Run complete off-target analysis.

        Args:
            query_fasta: Path to query FASTA file
            output_prefix: Output file prefix

        Returns:
            Tuple of (TSV path, JSON path)
        """
        logger.info(f"Starting off-target analysis with mode: {self.mode}")

        # Load query sequences
        query_seqs = self.load_query_sequences(query_fasta)
        if not query_seqs:
            raise ValueError("No query sequences loaded")

        all_alignments = []

        # Run alignments based on mode
        if self.mode in ("bowtie-only", "combined"):
            logger.info("Running Bowtie analysis...")
            bowtie_out = self.run_bowtie(query_fasta, v_mismatches=0)
            if bowtie_out:
                bowtie_alignments = self.parse_sam_alignments(bowtie_out, query_seqs)
                for alignment in bowtie_alignments:
                    alignment["source"] = "bowtie"
                all_alignments.extend(bowtie_alignments)
                logger.info(f"Bowtie found {len(bowtie_alignments)} alignments")

        if self.mode in ("bwa-mem2-only", "combined"):
            logger.info("Running BWA-MEM2 analysis...")
            bwa_out = self.run_bwa_mem2(query_fasta)
            if bwa_out:
                bwa_alignments = self.parse_sam_alignments(bwa_out, query_seqs)
                for alignment in bwa_alignments:
                    alignment["source"] = "bwa-mem2"

                # For combined mode, merge unique alignments
                if self.mode == "combined":
                    existing_coords = {(a["rname"], a["pos"], a["strand"], a["cigar"]) for a in all_alignments}
                    new_alignments = [
                        a
                        for a in bwa_alignments
                        if (a["rname"], a["pos"], a["strand"], a["cigar"]) not in existing_coords
                    ]
                    all_alignments.extend(new_alignments)
                    logger.info(f"BWA-MEM2 added {len(new_alignments)} unique alignments")
                else:
                    all_alignments.extend(bwa_alignments)
                    logger.info(f"BWA-MEM2 found {len(bwa_alignments)} alignments")

        # Filter and score alignments
        filtered_alignments = self.filter_alignments(all_alignments)

        # Write output files
        tsv_path, json_path = self._write_outputs(filtered_alignments, output_prefix)

        logger.info(f"Analysis complete. Results written to {tsv_path} and {json_path}")
        return tsv_path, json_path

    def _write_outputs(self, alignments: list[dict[str, Any]], output_prefix: str) -> tuple[str, str]:
        """Write alignments to TSV and JSON files."""

        # Ensure output directory exists
        prefix_path = Path(output_prefix)
        output_dir = prefix_path.parent if prefix_path.parent != prefix_path else Path()
        output_dir.mkdir(parents=True, exist_ok=True)

        tsv_path = f"{output_prefix}.tsv"
        json_path = f"{output_prefix}.json"

        # Write TSV
        headers = [
            "qname",
            "qseq",
            "rname",
            "coord",
            "strand",
            "cigar",
            "mapq",
            "AS",
            "NM",
            "mismatch_count",
            "mismatch_positions",
            "offtarget_score",
            "source",
        ]

        with Path(tsv_path).open("w") as f:
            f.write("\t".join(headers) + "\n")

            for alignment in alignments:
                coord = f"{alignment['rname']}:{alignment['pos']}"
                mm_positions = alignment.get("mismatch_positions", [])
                mm_pos_str = ";".join(map(str, mm_positions))

                row = [
                    alignment.get("qname", ""),
                    alignment.get("qseq", ""),
                    alignment.get("rname", ""),
                    coord,
                    alignment.get("strand", ""),
                    alignment.get("cigar", ""),
                    str(alignment.get("mapq", "")),
                    str(alignment.get("AS", "")),
                    str(alignment.get("NM", "")),
                    str(len(mm_positions)),
                    mm_pos_str,
                    f"{alignment.get('offtarget_score', 0):.2f}",
                    alignment.get("source", ""),
                ]

                f.write("\t".join(row) + "\n")

        # Write JSON
        json_output = []
        for alignment in alignments:
            coord = f"{alignment['rname']}:{alignment['pos']}"
            mm_positions = alignment.get("mismatch_positions", [])

            # Flag seed vs non-seed mismatches
            flagged_positions = []
            for pos in mm_positions:
                flagged_positions.append({"pos": pos, "is_seed": self.seed_start <= pos <= self.seed_end})

            json_item = {
                "qname": alignment.get("qname"),
                "qseq": alignment.get("qseq"),
                "rname": alignment.get("rname"),
                "coord": coord,
                "strand": alignment.get("strand"),
                "cigar": alignment.get("cigar"),
                "mapq": alignment.get("mapq"),
                "AS": alignment.get("AS"),
                "NM": alignment.get("NM"),
                "mismatch_count": len(mm_positions),
                "mismatch_positions": flagged_positions,
                "offtarget_score": alignment.get("offtarget_score", 0),
                "source": alignment.get("source"),
                "raw_tags": alignment.get("raw_tags", {}),
            }

            json_output.append(json_item)

        with Path(json_path).open("w") as f:
            json.dump(json_output, f, indent=2)

        return tsv_path, json_path


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Comprehensive off-target analysis for siRNA candidates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic combined analysis
  %(prog)s --mode combined --index-prefix /data/genomes/human/bwa_index --queries candidates.fasta --out-prefix results/human

  # BWA-only with custom parameters
  %(prog)s --mode bwa-mem2-only --index-prefix /data/genomes/human/bwa_index --queries candidates.fasta --out-prefix results/human --bwa-k 10 --bwa-T 20

  # Custom seed region
  %(prog)s --mode combined --index-prefix /data/genomes/human/bwa_index --queries candidates.fasta --out-prefix results/human --seed-start 1 --seed-end 7
        """,
    )

    parser.add_argument(
        "--mode",
        choices=["bowtie-only", "bwa-mem2-only", "combined"],
        default="combined",
        help="Analysis mode (default: combined)",
    )

    parser.add_argument("--index-prefix", required=True, help="Genome index prefix path")

    parser.add_argument("--queries", required=True, help="Query sequences FASTA file")

    parser.add_argument("--out-prefix", required=True, help="Output files prefix")

    parser.add_argument("--bwa-k", type=int, default=12, help="BWA-MEM2 seed length (default: 12)")

    parser.add_argument("--bwa-T", type=int, default=15, help="BWA-MEM2 minimum score threshold (default: 15)")

    parser.add_argument("--max-hits", type=int, default=10000, help="Maximum hits per query (default: 10000)")

    parser.add_argument("--seed-start", type=int, default=2, help="Seed region start position, 1-based (default: 2)")

    parser.add_argument("--seed-end", type=int, default=8, help="Seed region end position, 1-based (default: 8)")

    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Initialize analyzer
        analyzer = OffTargetAnalyzer(
            index_prefix=args.index_prefix,
            mode=args.mode,
            bwa_k=args.bwa_k,
            bwa_T=args.bwa_T,
            max_hits=args.max_hits,
            seed_start=args.seed_start,
            seed_end=args.seed_end,
        )

        # Run analysis
        tsv_path, json_path = analyzer.analyze(args.queries, args.out_prefix)

        print(f"TSV: {tsv_path}")
        print(f"JSON: {json_path}")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
