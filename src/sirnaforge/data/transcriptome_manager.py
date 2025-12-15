#!/usr/bin/env python3
"""Transcriptome Database Manager with local caching and automatic index building.

This module provides a clean interface for downloading, caching, and managing
transcriptome FASTA files with automatic BWA-MEM2 index building and cache management.
"""

import hashlib
import logging
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional, Union

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
        self, cache_dir: Optional[Union[str, Path]] = None, cache_ttl_days: int = 90, auto_build_indices: bool = True
    ):
        """Initialize the transcriptome database manager.

        Args:
            cache_dir: Directory for caching transcriptomes (default: ~/.cache/sirnaforge/transcriptomes)
            cache_ttl_days: Cache time-to-live in days (default: 90 days for large files)
            auto_build_indices: Automatically build BWA-MEM2 indices when missing
        """
        super().__init__(cache_subdir="transcriptomes", cache_dir=cache_dir, cache_ttl_days=cache_ttl_days)
        self.auto_build_indices = auto_build_indices

    def describe_source_status(self, source_name: str) -> dict[str, Any]:
        """Return cache metadata for a configured source."""
        source = self.SOURCES.get(source_name)
        if source is None:
            return {
                "namespace": "transcriptome",
                "identifier": source_name,
                "species": "unknown",
                "description": "unknown transcriptome source",
                "cached": False,
                "indexed": False,
                "cache_path": None,
                "index_path": None,
                "size_mb": None,
                "last_updated": None,
            }

        cache_key = source.cache_key()
        meta = self.metadata.get(cache_key)
        cache_path: Optional[Path] = None
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
            "namespace": "transcriptome",
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

    def describe_sources_status(
        self, source_names: Optional[Union[list[str], tuple[str, ...]]] = None
    ) -> list[dict[str, Any]]:
        """Return cache statuses for multiple transcriptome sources."""
        targets = source_names or list(self.SOURCES.keys())
        return [self.describe_source_status(name) for name in targets]

    def _metadata_from_dict(self, data: dict[str, Any]) -> CacheMetadata:
        """Create CacheMetadata with TranscriptomeSource."""
        return CacheMetadata.from_dict(data, source_class=TranscriptomeSource)

    def _get_index_path(self, meta: CacheMetadata) -> Optional[Path]:
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
    ) -> Optional[dict[str, Path]]:
        """Get transcriptome database, downloading and building index if needed.

        Args:
            source_name: Pre-configured source name (e.g., "ensembl_human_cdna")
            force_refresh: Force re-download even if cached
            build_index: Build BWA-MEM2 index if missing

        Returns:
            Dictionary with 'fasta' and optionally 'index' paths, or None if failed
        """
        if source_name not in self.SOURCES:
            available = ", ".join(self.SOURCES.keys())
            logger.error(f"Unknown transcriptome source: {source_name}. Available: {available}")
            return None

        source = self.SOURCES[source_name]
        cache_key = source.cache_key()
        cache_file = self.cache_dir / f"{cache_key}.fa"
        index_prefix = self.cache_dir / f"{cache_key}_index"

        # Use cached version if valid
        if not force_refresh and self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached {source.name} ({source.species}): {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        # Download transcriptome
        logger.info(f"🔄 Downloading {source.name} ({source.species})...")
        content = self._download_file(source)
        if not content or not content.strip():
            return None

        cache_file.write_text(content, encoding="utf-8")
        if cache_file.stat().st_size == 0:
            cache_file.unlink()
            return None

        # Update metadata
        self.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=self._compute_file_checksum(cache_file),
            file_path=str(cache_file),
        )
        logger.info(f"✅ Cached {source.name}: {cache_file} ({cache_file.stat().st_size:,} bytes)")

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def get_custom_transcriptome(
        self, fasta_path: Union[str, Path], build_index: bool = True, cache_name: Optional[str] = None
    ) -> Optional[dict[str, Path]]:
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
        if fasta_str.startswith(("http://", "https://", "ftp://")):
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

    def _handle_url_transcriptome(
        self, url: str, cache_name: Optional[str], build_index: bool
    ) -> Optional[dict[str, Path]]:
        """Download and cache transcriptome from URL."""
        cache_name = cache_name or url.split("/")[-1].replace(".gz", "").replace(".fa", "").replace(".fasta", "")
        source = TranscriptomeSource(
            name=cache_name,
            url=url,
            species="custom",
            compressed=url.endswith(".gz"),
            description=f"Custom transcriptome from {url}",
        )

        cache_key = source.cache_key()
        cache_file = self.cache_dir / f"{cache_key}.fa"
        index_prefix = self.cache_dir / f"{cache_key}_index"

        # Check cache
        if self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached custom transcriptome: {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        # Download
        content = self._download_file(source)
        if not content or not content.strip():
            return None

        cache_file.write_text(content, encoding="utf-8")
        if cache_file.stat().st_size == 0:
            cache_file.unlink()
            return None

        # Save metadata
        self.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=self._compute_file_checksum(cache_file),
            file_path=str(cache_file),
        )
        logger.info(f"✅ Cached custom transcriptome: {cache_file} ({cache_file.stat().st_size:,} bytes)")

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def _handle_cached_file(self, file_path: Path, cache_name: str, build_index: bool) -> dict[str, Path]:
        """Handle transcriptome file already in cache directory."""
        cache_key = hashlib.md5(f"local_{cache_name}_{file_path}".encode()).hexdigest()[:12]
        index_prefix = file_path.parent / f"{file_path.stem}_index"

        # Ensure metadata exists
        if cache_key not in self.metadata:
            self.metadata[cache_key] = CacheMetadata(
                source=TranscriptomeSource(
                    name=cache_name, url=str(file_path), species="custom", description=f"Local: {file_path}"
                ),
                downloaded_at=datetime.now().isoformat(),
                file_size=file_path.stat().st_size,
                checksum=self._compute_file_checksum(file_path),
                file_path=str(file_path),
            )

        return self._prepare_result_with_index(file_path, index_prefix, cache_key, build_index)

    def _cache_local_file(self, input_path: Path, cache_name: str, build_index: bool) -> dict[str, Path]:
        """Copy local file to cache and prepare it."""
        cache_key = hashlib.md5(f"local_{cache_name}_{input_path}".encode()).hexdigest()[:12]
        cache_file = self.cache_dir / f"{cache_name}_{cache_key}.fa"
        index_prefix = self.cache_dir / f"{cache_name}_{cache_key}_index"

        # Check if already cached
        if not cache_file.exists():
            logger.info(f"🔄 Copying {input_path.name} to cache...")
            cache_file.write_text(input_path.read_text())

        # Save metadata
        self.metadata[cache_key] = CacheMetadata(
            source=TranscriptomeSource(
                name=cache_name, url=str(input_path), species="custom", description=f"Local: {input_path}"
            ),
            downloaded_at=datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=self._compute_file_checksum(cache_file),
            file_path=str(cache_file),
        )

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def get_filtered_transcriptome(
        self,
        source_name: str,
        filters: list[str],
        force_refresh: bool = False,
        build_index: bool = True,
    ) -> Optional[dict[str, Path]]:
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

        source = self.SOURCES[source_name]
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
        self.metadata[filtered_cache_key] = CacheMetadata(
            source=TranscriptomeSource(
                name=f"{source.name}_filtered",
                url=source.url,
                species=source.species,
                description=f"{source.description} [filtered: {filter_spec}]",
            ),
            downloaded_at=datetime.now().isoformat(),
            file_size=filtered_fasta.stat().st_size,
            checksum=self._compute_file_checksum(filtered_fasta),
            file_path=str(filtered_fasta),
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
        return self.SOURCES

    def cache_info(self) -> dict[str, Any]:
        """Get information about the current cache state."""
        total_files = len(list(self.cache_dir.glob("*.fa")))
        total_size = sum(f.stat().st_size for f in self.cache_dir.glob("*.fa"))
        index_files = len(list(self.cache_dir.glob("*_index.amb")))

        return {
            "cache_directory": str(self.cache_dir),
            "total_fasta_files": total_files,
            "total_size_mb": total_size / (1024 * 1024),
            "index_files": index_files,
            "cache_ttl_days": self.cache_ttl.days,
            "auto_build_indices": self.auto_build_indices,
            "cached_transcriptomes": list(self.metadata.keys()),
        }

    def clean_cache(self, older_than_days: Optional[int] = None) -> None:
        """Clean old cache files.

        Args:
            older_than_days: Remove files older than this (default: use TTL)
        """
        if older_than_days is None:
            older_than_days = self.cache_ttl.days

        cutoff = datetime.now() - timedelta(days=older_than_days)
        removed_count = 0

        for cache_key in list(self.metadata.keys()):
            meta = self.metadata[cache_key]
            downloaded_at = datetime.fromisoformat(meta.downloaded_at)

            if downloaded_at < cutoff:
                # Remove FASTA file
                fasta_path = Path(meta.file_path)
                if fasta_path.exists():
                    fasta_path.unlink()
                    removed_count += 1

                # Remove index files if they exist
                index_path = self._get_index_path(meta)
                if index_path:
                    # Remove marker file for the index prefix (if present)
                    index_path.unlink(missing_ok=True)
                    for ext in [".amb", ".ann", ".bwt.2bit.64", ".pac"]:
                        index_file = index_path.with_suffix(ext)
                        if index_file.exists():
                            index_file.unlink()

                del self.metadata[cache_key]

        if removed_count > 0:
            self._save_metadata()
            logger.info(f"🧹 Cleaned {removed_count} old cache files")
        else:
            logger.info("🧹 No old cache files to clean")
