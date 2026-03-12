#!/usr/bin/env python3
"""Genomic annotation manager built on the transcriptome/reference cache pipeline."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from .transcriptome_manager import TranscriptomeManager, TranscriptomeSource

logger = logging.getLogger(__name__)


@dataclass
class AnnotationSource(TranscriptomeSource):
    """Annotation-specific source model."""

    pass


class AnnotationManager(TranscriptomeManager):
    """Annotation manager with caching and no index-building side effects."""

    SOURCE_LABEL = "annotation"

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        cache_ttl_days: int = 30,
    ):
        """Initialize annotation manager.

        Args:
            cache_dir: Optional cache directory override
            cache_ttl_days: Cache TTL in days
        """
        super().__init__(
            cache_dir=cache_dir,
            cache_ttl_days=cache_ttl_days,
            auto_build_indices=False,
            cache_subdir="annotations",
            sources={},
            source_label=self.SOURCE_LABEL,
        )

    def get_custom_annotation(self, annotation_path_or_url: str | Path, cache_name: str | None = None) -> Path | None:
        """Resolve and cache a user-supplied annotation resource.

        Supports local paths and remote URLs through the shared reference pipeline.
        """
        annotation_str = str(annotation_path_or_url)
        if self._is_remote_location(annotation_str):
            return self._handle_url_annotation(annotation_str, cache_name)

        result = self.get_custom_transcriptome(
            fasta_path=annotation_path_or_url,
            build_index=False,
            cache_name=cache_name,
        )
        if result is None:
            return None
        return result["fasta"]

    def _handle_url_annotation(self, url: str, cache_name: str | None = None) -> Path | None:
        """Download/cache remote annotation as-is (keep .gz compressed payloads)."""
        inferred_name = cache_name or self._default_cache_name_for_resource(url)
        compressed = url.endswith(".gz")
        suffix = ".gtf.gz" if compressed else ".gtf"

        source = AnnotationSource(
            name=inferred_name,
            url=url,
            species="custom",
            format="gtf",
            compressed=False,
            description=f"Custom annotation from {url}",
        )
        cache_key = source.cache_key()
        cache_file = self.cache_dir / f"{cache_key}{suffix}"

        if self._is_cache_valid(cache_key):
            logger.info("✅ Using cached custom annotation: %s", cache_file)
            return cache_file

        if not self._download_to_path(source, cache_file):
            return None

        self._record_cache_entry(cache_key, source, cache_file, persist=True)
        logger.info("✅ Cached custom annotation: %s (%s bytes)", cache_file, cache_file.stat().st_size)
        return cache_file
