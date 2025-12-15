#!/usr/bin/env python3
"""Unified cache manager for all sirnaforge reference databases.

Provides a single interface to manage miRNA databases, transcriptomes, and indices
using composition and protocol-based design.
"""

from dataclasses import dataclass
from typing import Protocol, TypedDict, Union, cast


class CacheStats(TypedDict):
    """Statistics for a cache component."""

    cache_directory: str
    total_files: int
    total_size_mb: float
    cache_ttl_days: int


class ClearResult(TypedDict):
    """Result of a cache clear operation."""

    files_deleted: int
    size_freed_mb: float
    status: str


class CacheComponent(Protocol):
    """Protocol for cache components."""

    def cache_info(self) -> CacheStats:
        """Get cache statistics."""
        ...

    def clear_cache(self, confirm: bool = False) -> ClearResult:
        """Clear the cache."""
        ...


@dataclass
class UnifiedCacheManager:
    """Unified manager for all sirnaforge caches.

    Uses composition to combine miRNA and transcriptome caches
    into a single, easy-to-use interface.
    """

    def __init__(self) -> None:
        """Initialize with all cache components."""
        from sirnaforge.data.mirna_manager import MiRNADatabaseManager  # noqa: PLC0415
        from sirnaforge.data.transcriptome_manager import TranscriptomeManager  # noqa: PLC0415

        self.mirna = MiRNADatabaseManager()
        self.transcriptome = TranscriptomeManager()

    def get_info(self, include_mirna: bool = True, include_transcriptome: bool = True) -> dict[str, CacheStats]:
        """Get cache info for selected components.

        Args:
            include_mirna: Include miRNA cache stats
            include_transcriptome: Include transcriptome cache stats

        Returns:
            Dictionary mapping component name to cache stats
        """
        info: dict[str, CacheStats] = {}
        if include_mirna:
            info["mirna"] = cast(CacheStats, self.mirna.cache_info())
        if include_transcriptome:
            info["transcriptome"] = self._get_transcriptome_stats()
        return info

    def _get_transcriptome_stats(self) -> CacheStats:
        """Get transcriptome cache statistics."""
        base_info = self.transcriptome.cache_info()
        return CacheStats(
            cache_directory=base_info["cache_directory"],
            total_files=base_info["total_fasta_files"] + base_info["index_files"],
            total_size_mb=base_info["total_size_mb"],
            cache_ttl_days=base_info["cache_ttl_days"],
        )

    def clear(
        self, clear_mirna: bool = False, clear_transcriptome: bool = False, dry_run: bool = False
    ) -> dict[str, ClearResult]:
        """Clear selected cache components.

        Args:
            clear_mirna: Clear miRNA databases
            clear_transcriptome: Clear transcriptomes and indices
            dry_run: Show what would be deleted without deleting

        Returns:
            Dictionary mapping component name to clear results
        """
        results: dict[str, ClearResult] = {}

        if clear_mirna:
            results["mirna"] = cast(ClearResult, self.mirna.clear_cache(confirm=not dry_run))

        if clear_transcriptome:
            results["transcriptome"] = self._clear_transcriptome(dry_run)

        return results

    def _clear_transcriptome(self, dry_run: bool) -> ClearResult:
        """Clear transcriptome cache.

        Delegates to TranscriptomeManager.clear_cache so cache deletion
        behavior stays centralized.
        """
        result = self.transcriptome.clear_cache(confirm=not dry_run)
        return ClearResult(
            files_deleted=int(result.get("files_deleted", 0)),
            size_freed_mb=float(result.get("size_freed_mb", 0.0)),
            status="dry run - no files deleted" if dry_run else "cleared",
        )

    def get_total_stats(self) -> dict[str, Union[float, int]]:
        """Get combined statistics across all caches."""
        info = self.get_info()
        total_files = sum(stats["total_files"] for stats in info.values())
        total_size = sum(stats["total_size_mb"] for stats in info.values())

        return {"total_files": total_files, "total_size_mb": total_size}
