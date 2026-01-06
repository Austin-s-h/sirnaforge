#!/usr/bin/env python3
"""Transcriptome filtering utilities for reducing reference size.

Provides functions to filter transcriptome FASTA files based on:
- Biotype (protein_coding, etc.)
- Canonical transcript status
- Specific regions (e.g., 3' UTR)

This helps reduce memory requirements for BWA-MEM2 indexing on low-RAM machines.
"""

import logging
import re
from collections.abc import Callable
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class TranscriptFilter:
    """Filter transcripts based on header metadata and biotype."""

    # Ensembl cDNA header format examples:
    # >ENST00000456328.2 cdna chromosome:GRCh38:1:11869:14409:1 gene:ENSG00000223972.5 gene_biotype:transcribed_unprocessed_pseudogene transcript_biotype:processed_transcript gene_symbol:DDX11L1 description:DEAD/H-box helicase 11 like 1 [Source:HGNC Symbol;Acc:HGNC:37102]
    # >ENST00000450305.2 cdna chromosome:GRCh38:1:12010:13670:1 gene:ENSG00000223972.5 gene_biotype:transcribed_unprocessed_pseudogene transcript_biotype:transcribed_unprocessed_pseudogene gene_symbol:DDX11L1 description:DEAD/H-box helicase 11 like 1 [Source:HGNC Symbol;Acc:HGNC:37102]

    @staticmethod
    def parse_ensembl_header(description: str) -> dict[str, str]:
        """Parse Ensembl FASTA header to extract metadata.

        Args:
            description: Full FASTA description line (without '>')

        Returns:
            Dictionary of parsed fields (gene_biotype, transcript_biotype, gene_symbol, etc.)
        """
        metadata: dict[str, str] = {}

        # Extract fields in format "key:value"
        # Pattern matches: word_key:value_with_spaces (stops at next key: or end)
        # Example: "gene_biotype:protein_coding" or "description:text with spaces"
        for match in re.finditer(r"(\w+):([^\s]+(?:\s+[^\s:]+)*?)(?=\s+\w+:|$)", description):
            key = match.group(1)
            value = match.group(2).strip()
            metadata[key] = value

        return metadata

    @staticmethod
    def is_protein_coding(record: SeqRecord) -> bool:
        """Check if transcript is protein coding.

        Args:
            record: SeqRecord from FASTA

        Returns:
            True if protein_coding, False otherwise
        """
        metadata = TranscriptFilter.parse_ensembl_header(record.description)

        # Check both gene_biotype and transcript_biotype
        gene_biotype = metadata.get("gene_biotype", "")
        transcript_biotype = metadata.get("transcript_biotype", "")

        return "protein_coding" in gene_biotype or "protein_coding" in transcript_biotype

    @staticmethod
    def is_canonical(record: SeqRecord) -> bool:
        """Check if transcript is canonical (MANE Select or Ensembl canonical).

        Args:
            record: SeqRecord from FASTA

        Returns:
            True if canonical, False otherwise
        """
        description = record.description.lower()

        # Look for canonical markers in description
        # Ensembl marks canonical transcripts with "canonical:1" or in tags
        # MANE Select transcripts are also considered canonical
        canonical_markers = [
            "canonical:1",
            "mane_select",
            "mane select",
            "appris_principal",
            "tag:basic",  # Basic tag indicates high-quality transcript
        ]

        return any(marker in description for marker in canonical_markers)

    @staticmethod
    def filter_fasta(
        input_fasta: Path,
        output_fasta: Path,
        filter_func: Callable[[SeqRecord], bool],
        filter_name: str = "custom",
    ) -> int:
        """Filter FASTA file using a custom filter function.

        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output FASTA file
            filter_func: Function that returns True for records to keep
            filter_name: Name of filter for logging

        Returns:
            Number of sequences kept
        """
        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

        kept_count = 0
        total_count = 0

        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        with output_fasta.open("w") as out_handle:
            for record in SeqIO.parse(input_fasta, "fasta"):
                total_count += 1
                if filter_func(record):
                    SeqIO.write(record, out_handle, "fasta")
                    kept_count += 1

        logger.info(
            f"Filtered {input_fasta.name} with '{filter_name}': "
            f"kept {kept_count}/{total_count} sequences ({kept_count / total_count * 100:.1f}%)"
        )

        return kept_count

    @staticmethod
    def apply_protein_coding_filter(input_fasta: Path, output_fasta: Path) -> int:
        """Filter to protein-coding transcripts only.

        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output FASTA file

        Returns:
            Number of sequences kept
        """
        return TranscriptFilter.filter_fasta(
            input_fasta, output_fasta, TranscriptFilter.is_protein_coding, "protein_coding"
        )

    @staticmethod
    def apply_canonical_filter(input_fasta: Path, output_fasta: Path) -> int:
        """Filter to canonical transcripts only.

        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output FASTA file

        Returns:
            Number of sequences kept
        """
        return TranscriptFilter.filter_fasta(input_fasta, output_fasta, TranscriptFilter.is_canonical, "canonical_only")

    @staticmethod
    def apply_combined_filter(input_fasta: Path, output_fasta: Path, filters: list[str]) -> int:
        """Apply multiple filters in sequence.

        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output FASTA file
            filters: List of filter names (e.g., ['protein_coding', 'canonical_only'])

        Returns:
            Number of sequences kept
        """

        def combined_filter_func(record: SeqRecord) -> bool:
            for filter_name in filters:
                if filter_name == "protein_coding":
                    if not TranscriptFilter.is_protein_coding(record):
                        return False
                elif filter_name == "canonical_only":
                    if not TranscriptFilter.is_canonical(record):
                        return False
                else:
                    logger.warning(f"Unknown filter: {filter_name}")
            return True

        filter_desc = "+".join(filters)
        return TranscriptFilter.filter_fasta(input_fasta, output_fasta, combined_filter_func, filter_desc)


def get_filter_spec(filter_string: str | None) -> list[str]:
    """Parse filter specification string into list of filter names.

    Args:
        filter_string: Comma-separated filter names (e.g., "protein_coding,canonical_only")

    Returns:
        List of filter names, or empty list if None
    """
    if not filter_string:
        return []

    filters = [f.strip() for f in filter_string.split(",") if f.strip()]
    valid_filters = {"protein_coding", "canonical_only"}

    invalid = [f for f in filters if f not in valid_filters]
    if invalid:
        raise ValueError(f"Invalid filter(s): {', '.join(invalid)}. Valid: {', '.join(valid_filters)}")

    return filters
