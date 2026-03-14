#!/usr/bin/env python3
"""Transcriptome Database Manager with local caching and automatic index building.

This module provides a clean interface for downloading, caching, and managing
transcriptome FASTA files with automatic BWA-MEM2 index building and cache management.
"""

import hashlib
import logging
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from .reference_manager import CacheMetadata, ReferenceManager, ReferenceSource

logger = logging.getLogger(__name__)


@dataclass
class TranscriptomeSource(ReferenceSource):
    """Transcriptome-specific database source configuration.

    Inherits from ReferenceSource with transcriptome-specific extensions.
    """

    pass


class TranscriptomeManager(ReferenceManager[TranscriptomeSource]):
    """Transcriptome database manager with caching and automatic BWA-MEM2 index building."""

    SOURCE_LABEL = "transcriptome"

    # Common transcriptome sources
    SOURCES = {
        "ensembl_human_cdna": TranscriptomeSource(
            name="ensembl_cdna",
            url="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            species="human",
            format="fasta",
            compressed=True,
            description="Ensembl human cDNA sequences (GRCh38)",
        ),
        "ensembl_mouse_cdna": TranscriptomeSource(
            name="ensembl_cdna",
            url="https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
            species="mouse",
            format="fasta",
            compressed=True,
            description="Ensembl mouse cDNA sequences (GRCm39)",
        ),
        "ensembl_rat_cdna": TranscriptomeSource(
            name="ensembl_cdna",
            url="https://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/cdna/Rattus_norvegicus.GRCr8.cdna.all.fa.gz",
            species="rat",
            format="fasta",
            compressed=True,
            description="Ensembl rat cDNA sequences (GRCr8)",
        ),
        "ensembl_macaque_cdna": TranscriptomeSource(
            name="ensembl_cdna",
            url="https://ftp.ensembl.org/pub/current_fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz",
            species="macaque",
            format="fasta",
            compressed=True,
            description="Ensembl rhesus macaque cDNA sequences (Mmul_10)",
        ),
    }

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        cache_ttl_days: int = 90,
        auto_build_indices: bool = True,
        local_content_dedupe: bool = True,
        cache_subdir: str = "transcriptomes",
        sources: dict[str, TranscriptomeSource] | None = None,
        source_label: str | None = None,
    ):
        """Initialize the transcriptome database manager.

        Args:
            cache_dir: Directory for caching transcriptomes (default: ~/.cache/sirnaforge/transcriptomes)
            cache_ttl_days: Cache time-to-live in days (default: 90 days for large files)
            auto_build_indices: Automatically build BWA-MEM2 indices when missing
            local_content_dedupe: Reuse cached local files by content hash across different input paths
            cache_subdir: Cache subdirectory override for composable manager reuse
            sources: Optional source dictionary override
            source_label: Optional source label used in status/error reporting
        """
        super().__init__(cache_subdir=cache_subdir, cache_dir=cache_dir, cache_ttl_days=cache_ttl_days)
        self.auto_build_indices = auto_build_indices
        self.local_content_dedupe = local_content_dedupe
        self.sources: dict[str, TranscriptomeSource] = dict(sources or self.SOURCES)
        self.source_label = source_label or self.SOURCE_LABEL
        self.local_content_index: dict[str, str] = {}
        self._rebuild_local_content_index()

    def _rebuild_local_content_index(self) -> None:
        """Rebuild local content-hash index from loaded metadata."""
        self.local_content_index.clear()
        for cache_key, meta in self.metadata.items():
            extra = meta.extra or {}
            content_hash = extra.get("local_content_hash")
            if isinstance(content_hash, str) and content_hash:
                self.local_content_index[content_hash] = cache_key

    def describe_source_status(self, source_name: str) -> dict[str, Any]:
        """Return cache metadata for a configured source."""
        source = self.sources.get(source_name)
        if source is None:
            return {
                "namespace": self.source_label,
                "identifier": source_name,
                "species": "unknown",
                "description": f"unknown {self.source_label} source",
                "cached": False,
                "indexed": False,
                "cache_path": None,
                "index_path": None,
                "size_mb": None,
                "last_updated": None,
            }

        cache_key = source.cache_key()
        meta = self.metadata.get(cache_key)
        cache_path: Path | None = None
        if meta:
            cache_path = Path(meta.file_path)
            if not cache_path.exists():
                cache_path = None
        cached = self._is_cache_valid(cache_key)
        index_path = None
        indexed = False

        if meta:
            candidate_index = self._get_index_path(meta)
            if candidate_index and candidate_index.with_suffix(".amb").exists():
                index_path = candidate_index
                indexed = True
            elif candidate_index:
                index_path = candidate_index

        size_mb = (meta.file_size / (1024 * 1024)) if meta else None
        last_updated = meta.downloaded_at if meta else None

        return {
            "namespace": self.source_label,
            "identifier": source_name,
            "species": source.species,
            "description": source.description,
            "cached": cached,
            "indexed": indexed,
            "cache_path": str(cache_path) if cache_path else None,
            "index_path": str(index_path) if index_path else None,
            "size_mb": size_mb,
            "last_updated": last_updated,
        }

    def describe_sources_status(self, source_names: list[str] | tuple[str, ...] | None = None) -> list[dict[str, Any]]:
        """Return cache statuses for multiple transcriptome sources."""
        targets = source_names or list(self.sources.keys())
        return [self.describe_source_status(name) for name in targets]

    def _metadata_from_dict(self, data: dict[str, Any]) -> CacheMetadata:
        """Create CacheMetadata with TranscriptomeSource."""
        return CacheMetadata.from_dict(data, source_class=TranscriptomeSource)

    def _get_index_path(self, meta: CacheMetadata) -> Path | None:
        """Get index path from metadata's extra field."""
        if meta.extra and "index_path" in meta.extra:
            return Path(meta.extra["index_path"])
        return None

    def _set_index_path(self, meta: CacheMetadata, index_path: Path) -> None:
        """Set index path in metadata's extra field."""
        if meta.extra is None:
            meta.extra = {}
        meta.extra["index_path"] = str(index_path)
        meta.extra["index_built_at"] = datetime.now().isoformat()

    def _is_index_complete(self, index_prefix: Path) -> bool:
        """Check if all required BWA-MEM2 index files exist and are non-empty.

        Args:
            index_prefix: Path prefix for index files

        Returns:
            True if all index files are complete, False otherwise
        """
        try:
            # Late import to avoid circular dependency
            from sirnaforge.core.off_target import validate_index_files  # noqa: PLC0415

            return validate_index_files(index_prefix, tool="bwa-mem2")
        except Exception as e:
            logger.debug(f"Index validation failed: {e}")
            return False

    def _build_index(self, fasta_path: Path, index_prefix: Path) -> bool:
        """Build BWA-MEM2 index for transcriptome FASTA.

        Args:
            fasta_path: Path to FASTA file
            index_prefix: Path prefix for index files

        Returns:
            True if successful, False otherwise
        """
        try:
            # Late import to avoid circular dependency
            from sirnaforge.core.off_target import build_bwa_index as _build_bwa_index  # noqa: PLC0415

            logger.info(f"🔨 Building BWA-MEM2 index for {fasta_path.name}...")
            _build_bwa_index(fasta_path, index_prefix)
            return True
        except Exception as e:
            logger.error(f"❌ Failed to build BWA-MEM2 index: {e}")
            return False

    def _ensure_index_marker(self, index_prefix: Path) -> None:
        """Ensure the index prefix path exists as a filesystem entry.

        BWA(-MEM2) produces multiple files that share a prefix (e.g. <prefix>.amb,
        <prefix>.ann, ...). The prefix itself is not a file created by the tool.

        Some higher-level code/tests treat the prefix as a `Path` and call
        `.exists()`. Creating a tiny marker file at the prefix path makes that
        check meaningful without changing the prefix semantics.
        """
        try:
            index_prefix.parent.mkdir(parents=True, exist_ok=True)
            index_prefix.touch(exist_ok=True)
        except Exception as e:
            logger.debug(f"Could not create index marker for {index_prefix}: {e}")

    def get_transcriptome(  # noqa: PLR0911
        self, source_name: str, force_refresh: bool = False, build_index: bool = True
    ) -> dict[str, Path] | None:
        """Get transcriptome database, downloading and building index if needed.

        Args:
            source_name: Pre-configured source name (e.g., "ensembl_human_cdna")
            force_refresh: Force re-download even if cached
            build_index: Build BWA-MEM2 index if missing

        Returns:
            Dictionary with 'fasta' and optionally 'index' paths, or None if failed
        """
        if source_name not in self.sources:
            available = ", ".join(self.sources.keys())
            logger.error(f"Unknown {self.source_label} source: {source_name}. Available: {available}")
            return None

        source = self.sources[source_name]
        cache_key = self._cache_key_for_remote_uri(source.url) or source.cache_key()
        cache_file = (
            Path(self.metadata[cache_key].file_path)
            if cache_key in self.metadata
            else self.cache_dir / f"{cache_key}.fa"
        )
        index_prefix = self.cache_dir / f"{cache_key}_index"

        # Use cached version if valid
        if not force_refresh and self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached {source.name} ({source.species}): {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        # Recover from interrupted prior runs where file exists but metadata is missing.
        if not force_refresh and self._recover_remote_cache_entry(
            source=source, cache_key=cache_key, cache_file=cache_file
        ):
            logger.info(f"✅ Using recovered cached {source.name} ({source.species}): {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        # Download transcriptome
        logger.info(f"🔄 Downloading {source.name} ({source.species})...")
        if not self._download_to_path(source, cache_file):
            return None

        # Update metadata
        self._record_cache_entry(cache_key, source, cache_file)
        logger.info(f"✅ Cached {source.name}: {cache_file} ({cache_file.stat().st_size:,} bytes)")

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def get_custom_transcriptome(
        self, fasta_path: str | Path, build_index: bool = True, cache_name: str | None = None
    ) -> dict[str, Path] | None:
        """Process a custom transcriptome FASTA with caching and index building.

        Args:
            fasta_path: Path or URL to transcriptome FASTA file
            build_index: Build BWA-MEM2 index if missing
            cache_name: Custom cache name (default: derived from filename)

        Returns:
            Dictionary with 'fasta' and optionally 'index' paths, or None if failed
        """
        fasta_str = str(fasta_path)

        # Handle URL downloads
        if self._is_remote_location(fasta_str):
            return self._handle_url_transcriptome(fasta_str, cache_name, build_index)

        # Handle local files
        input_path = Path(fasta_path)
        if not input_path.exists():
            logger.error(f"FASTA file not found: {input_path}")
            return None

        # File already in cache dir? Use it directly
        if input_path.parent == self.cache_dir:
            return self._handle_cached_file(input_path, cache_name or input_path.stem, build_index)

        # Copy file to cache
        return self._cache_local_file(input_path, cache_name or input_path.stem, build_index)

    def _handle_url_transcriptome(self, url: str, cache_name: str | None, build_index: bool) -> dict[str, Path] | None:
        """Download and cache transcriptome from URL."""
        cache_name = cache_name or self._default_cache_name_for_resource(url)
        source = TranscriptomeSource(
            name=cache_name,
            url=url,
            species="custom",
            compressed=url.endswith(".gz"),
            description=f"Custom transcriptome from {url}",
        )

        cache_key = self._cache_key_for_remote_uri(url) or source.cache_key()
        cache_file = (
            Path(self.metadata[cache_key].file_path)
            if cache_key in self.metadata
            else self.cache_dir / f"{cache_key}.fa"
        )
        index_prefix = self.cache_dir / f"{cache_key}_index"

        # Check cache
        if self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached custom transcriptome: {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        if self._recover_remote_cache_entry(source=source, cache_key=cache_key, cache_file=cache_file):
            logger.info(f"✅ Using recovered cached custom transcriptome: {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        # Download
        if not self._download_to_path(source, cache_file):
            return None

        # Save metadata
        self._record_cache_entry(cache_key, source, cache_file)
        logger.info(f"✅ Cached custom transcriptome: {cache_file} ({cache_file.stat().st_size:,} bytes)")

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def _handle_cached_file(self, file_path: Path, cache_name: str, build_index: bool) -> dict[str, Path]:
        """Handle transcriptome file already in cache directory."""
        extra: dict[str, Any] | None = None
        if self.local_content_dedupe:
            content_hash = self._compute_file_checksum(file_path)
            cache_key = hashlib.md5(f"local_content_{content_hash}".encode()).hexdigest()[:12]
            extra = {"local_content_hash": content_hash}
            self.local_content_index[content_hash] = cache_key
        else:
            cache_key = hashlib.md5(f"local_{cache_name}_{file_path}".encode()).hexdigest()[:12]
        index_prefix = file_path.parent / f"{file_path.stem}_index"

        # Ensure metadata exists
        if cache_key not in self.metadata:
            self._record_cache_entry(
                cache_key,
                TranscriptomeSource(
                    name=cache_name, url=str(file_path), species="custom", description=f"Local: {file_path}"
                ),
                file_path,
                extra=extra,
            )

        return self._prepare_result_with_index(file_path, index_prefix, cache_key, build_index)

    def _cache_local_file(self, input_path: Path, cache_name: str, build_index: bool) -> dict[str, Path]:
        """Copy local file to cache and prepare it."""
        extra: dict[str, Any] | None = None
        if self.local_content_dedupe:
            content_hash = self._compute_file_checksum(input_path)
            known_cache_key = self.local_content_index.get(content_hash)
            if known_cache_key and known_cache_key in self.metadata:
                known_meta = self.metadata[known_cache_key]
                known_path = Path(known_meta.file_path)
                if known_path.exists():
                    logger.info("✅ Reusing local cached transcriptome by content hash: %s", known_path)
                    index_prefix = self._get_index_path(known_meta) or (known_path.parent / f"{known_path.stem}_index")
                    return self._prepare_result_with_index(known_path, index_prefix, known_cache_key, build_index)

            cache_key = hashlib.md5(f"local_content_{content_hash}".encode()).hexdigest()[:12]
            cache_file = self.cache_dir / f"{cache_key}.fa"
            index_prefix = self.cache_dir / f"{cache_key}_index"
            extra = {"local_content_hash": content_hash}
            self.local_content_index[content_hash] = cache_key
        else:
            cache_key = hashlib.md5(f"local_{cache_name}_{input_path}".encode()).hexdigest()[:12]
            cache_file = self.cache_dir / f"{cache_name}_{cache_key}.fa"
            index_prefix = self.cache_dir / f"{cache_name}_{cache_key}_index"

        # Check if already cached
        if not cache_file.exists():
            logger.info(f"🔄 Copying {input_path.name} to cache...")
            shutil.copyfile(input_path, cache_file)

        # Save metadata
        self._record_cache_entry(
            cache_key,
            TranscriptomeSource(
                name=cache_name, url=str(input_path), species="custom", description=f"Local: {input_path}"
            ),
            cache_file,
            extra=extra,
        )

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def get_filtered_transcriptome(
        self,
        source_name: str,
        filters: list[str],
        force_refresh: bool = False,
        build_index: bool = True,
    ) -> dict[str, Path] | None:
        """Get a filtered transcriptome with caching.

        Args:
            source_name: Pre-configured source name (e.g., "ensembl_human_cdna")
            filters: List of filter names (e.g., ['protein_coding', 'canonical_only'])
            force_refresh: Force re-download and re-filter
            build_index: Build BWA-MEM2 index if missing

        Returns:
            Dictionary with 'fasta' and optionally 'index' paths, or None if failed
        """
        if not filters:
            return self.get_transcriptome(source_name, force_refresh, build_index)

        # Ensure base transcriptome is cached (without building index)
        base_result = self.get_transcriptome(source_name, force_refresh, build_index=False)
        if not base_result:
            return None

        from .transcriptome_filter import TranscriptFilter  # noqa: PLC0415

        source = self.sources[source_name]
        filter_spec = "+".join(sorted(filters))
        filtered_cache_key = f"{source.cache_key()}_{filter_spec}"
        filtered_fasta = self.cache_dir / f"{filtered_cache_key}.fa"
        filtered_index = self.cache_dir / f"{filtered_cache_key}_index"

        # Check cache
        if not force_refresh and filtered_cache_key in self.metadata:
            cached_path = Path(self.metadata[filtered_cache_key].file_path)
            if cached_path.exists():
                return self._prepare_result_with_index(cached_path, filtered_index, filtered_cache_key, build_index)

        # Apply filters
        logger.info(f"🔍 Applying filters to {source_name}: {', '.join(filters)}")
        kept = TranscriptFilter.apply_combined_filter(base_result["fasta"], filtered_fasta, filters)
        if kept == 0:
            filtered_fasta.unlink(missing_ok=True)
            return None

        # Cache metadata
        self._record_cache_entry(
            filtered_cache_key,
            TranscriptomeSource(
                name=f"{source.name}_filtered",
                url=source.url,
                species=source.species,
                description=f"{source.description} [filtered: {filter_spec}]",
            ),
            filtered_fasta,
            extra={"filters": filters, "kept_count": kept},
        )

        return self._prepare_result_with_index(filtered_fasta, filtered_index, filtered_cache_key, build_index)

    def _prepare_result_with_index(
        self, fasta: Path, index_prefix: Path, cache_key: str, build_index: bool
    ) -> dict[str, Path]:
        """Helper to prepare result dict with optional index building."""
        if not (build_index and self.auto_build_indices):
            self._save_metadata()
            return {"fasta": fasta}

        meta = self.metadata[cache_key]
        index_path = self._get_index_path(meta) or index_prefix

        if self._is_index_complete(index_path):
            self._ensure_index_marker(index_path)
            logger.info(f"✅ Using cached BWA-MEM2 index: {index_path}")
            return {"fasta": fasta, "index": index_path}

        logger.info(f"⚠️  Building index: {index_prefix}")
        if self._build_index(fasta, index_prefix):
            self._ensure_index_marker(index_prefix)
            self._set_index_path(meta, index_prefix)
            self._save_metadata()
            return {"fasta": fasta, "index": index_prefix}

        logger.warning("Index build failed, returning FASTA without index")
        self._save_metadata()
        return {"fasta": fasta}

    def list_available_sources(self) -> dict[str, TranscriptomeSource]:
        """List all pre-configured transcriptome sources."""
        return self.sources

    def _cache_info_files(self) -> list[Path]:
        """Return transcriptome FASTA files to include in total cache stats."""
        return list(self.cache_dir.glob("*.fa"))

    def _cache_info_extra(self, total_files: int, _total_size_bytes: int) -> dict[str, Any]:
        """Augment base cache stats with transcriptome/index-specific fields."""
        index_files = len(list(self.cache_dir.glob("*_index.amb")))

        return {
            "total_fasta_files": total_files,
            "index_files": index_files,
            "auto_build_indices": self.auto_build_indices,
            "local_content_dedupe": self.local_content_dedupe,
            "cached_references": list(self.metadata.keys()),
            "cached_transcriptomes": list(self.metadata.keys()),
        }

    def _clean_extra_files(self, meta: CacheMetadata) -> None:
        """Remove BWA index files associated with a cached transcriptome entry."""
        index_path = self._get_index_path(meta)
        if index_path is None:
            return

        # Remove marker file for the index prefix (if present)
        index_path.unlink(missing_ok=True)
        for ext in [".amb", ".ann", ".bwt.2bit.64", ".pac"]:
            index_path.with_suffix(ext).unlink(missing_ok=True)
