"""Data handling and retrieval modules."""

from .gene_search import (
    DatabaseType,
    GeneInfo,
    GeneSearcher,
    GeneSearchResult,
    TranscriptInfo,
    search_gene_sync,
    search_multiple_databases_sync,
)
from .genome_manager import GenomeManager
from .species_registry import normalize_species_name
from .transcriptome_manager import TranscriptomeManager

__all__ = [
    "DatabaseType",
    "GeneInfo",
    "GeneSearcher",
    "GeneSearchResult",
    "GenomeManager",
    "TranscriptInfo",
    "TranscriptomeManager",
    "normalize_species_name",
    "search_gene_sync",
    "search_multiple_databases_sync",
]
