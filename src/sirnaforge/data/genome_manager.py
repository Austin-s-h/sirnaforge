#!/usr/bin/env python3
"""Genome reference manager built on the composable transcriptome manager core."""

from __future__ import annotations

from pathlib import Path

from .transcriptome_manager import TranscriptomeManager, TranscriptomeSource


class GenomeManager(TranscriptomeManager):
    """Genome FASTA manager with caching and optional BWA-MEM2 index building."""

    SOURCE_LABEL = "genome"

    SOURCES = {
        "ensembl_human_hg38_primary": TranscriptomeSource(
            name="ensembl_human_hg38_primary",
            url=(
                "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
                "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            ),
            species="human",
            format="fasta",
            compressed=True,
            description="Ensembl human GRCh38 primary assembly genomic reference",
        ),
        # TODO: Add other genomes...
    }

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        cache_ttl_days: int = 30,
        auto_build_indices: bool = True,
    ):
        """Initialize genome manager.

        Args:
            cache_dir: Directory for genome cache
            cache_ttl_days: Cache time-to-live in days
            auto_build_indices: Automatically build BWA-MEM2 indices when missing
        """
        super().__init__(
            cache_dir=cache_dir,
            cache_ttl_days=cache_ttl_days,
            auto_build_indices=auto_build_indices,
            cache_subdir="genomes",
            sources=self.SOURCES,
            source_label=self.SOURCE_LABEL,
        )

    def get_genome(
        self,
        source_name: str,
        force_refresh: bool = False,
        build_index: bool = True,
    ) -> dict[str, Path] | None:
        """Get genome reference, downloading and indexing as needed."""
        return self.get_transcriptome(source_name=source_name, force_refresh=force_refresh, build_index=build_index)

    def get_custom_genome(
        self,
        fasta_path: str | Path,
        build_index: bool = True,
        cache_name: str | None = None,
    ) -> dict[str, Path] | None:
        """Get custom genome FASTA through shared caching/indexing pipeline."""
        return self.get_custom_transcriptome(fasta_path=fasta_path, build_index=build_index, cache_name=cache_name)
