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
            url="https://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz",
            species="rat",
            format="fasta",
            compressed=True,
            description="Ensembl rat cDNA sequences (mRatBN7.2)",
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
            logger.error(f"Unknown transcriptome source: {source_name}")
            logger.info(f"Available sources: {', '.join(self.SOURCES.keys())}")
            return None

        source = self.SOURCES[source_name]
        cache_key = source.cache_key()
        cache_file = self.cache_dir / f"{cache_key}.fa"
        index_prefix = self.cache_dir / f"{cache_key}_index"

        # Check if we can use cached version
        if not force_refresh and self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached {source.name} ({source.species}): {cache_file}")
            meta = self.metadata[cache_key]

            # Check if index exists and is needed
            if build_index and self.auto_build_indices:
                index_path = self._get_index_path(meta)
                if index_path and index_path.with_suffix(".amb").exists():
                    return {"fasta": cache_file, "index": index_path}
                # Index missing, build it
                if self._build_index(cache_file, index_prefix):
                    self._set_index_path(meta, index_prefix)
                    self._save_metadata()
                    return {"fasta": cache_file, "index": index_prefix}
                return {"fasta": cache_file}
            return {"fasta": cache_file}

        # Download transcriptome
        logger.info(f"🔄 Downloading {source.name} ({source.species})...")
        content = self._download_file(source)

        if content is None or not content.strip():
            logger.error(f"Failed to download transcriptome: {source_name}")
            return None

        # Save to cache
        with cache_file.open("w", encoding="utf-8") as f:
            f.write(content)

        if cache_file.stat().st_size == 0:
            logger.error(f"Downloaded content for {source_name} is empty")
            cache_file.unlink(missing_ok=True)
            return None

        # Update metadata
        checksum = self._compute_file_checksum(cache_file)
        self.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=checksum,
            file_path=str(cache_file),
        )

        logger.info(f"✅ Cached {source.name} ({source.species}): {cache_file} ({cache_file.stat().st_size:,} bytes)")

        # Build index if requested
        if build_index and self.auto_build_indices:
            if self._build_index(cache_file, index_prefix):
                self._set_index_path(self.metadata[cache_key], index_prefix)
                self._save_metadata()
                return {"fasta": cache_file, "index": index_prefix}
            self._save_metadata()
            return {"fasta": cache_file}
        self._save_metadata()
        return {"fasta": cache_file}

    def get_custom_transcriptome(  # noqa: PLR0911, PLR0912
        self, fasta_path: Union[str, Path], build_index: bool = True, cache_name: Optional[str] = None
    ) -> Optional[dict[str, Path]]:
        """Process a custom transcriptome FASTA file with caching and index building.

        Args:
            fasta_path: Path or URL to transcriptome FASTA file
            build_index: Build BWA-MEM2 index if missing
            cache_name: Custom cache name (default: derived from filename)

        Returns:
            Dictionary with 'fasta' and optionally 'index' paths, or None if failed
        """
        fasta_str = str(fasta_path)

        # Check if it's a URL
        if fasta_str.startswith(("http://", "https://", "ftp://")):
            # Create temporary source for URL
            if cache_name is None:
                cache_name = fasta_str.split("/")[-1].replace(".gz", "").replace(".fa", "").replace(".fasta", "")

            source = TranscriptomeSource(
                name=cache_name,
                url=fasta_str,
                species="custom",
                compressed=fasta_str.endswith(".gz"),
                description=f"Custom transcriptome from {fasta_str}",
            )

            cache_key = source.cache_key()
            cache_file = self.cache_dir / f"{cache_key}.fa"
            index_prefix = self.cache_dir / f"{cache_key}_index"

            # Check cache
            if self._is_cache_valid(cache_key):
                logger.info(f"✅ Using cached custom transcriptome: {cache_file}")
                meta = self.metadata[cache_key]

                if build_index and self.auto_build_indices:
                    index_path = self._get_index_path(meta)
                    if index_path and index_path.with_suffix(".amb").exists():
                        return {"fasta": cache_file, "index": index_path}
                    if self._build_index(cache_file, index_prefix):
                        self._set_index_path(meta, index_prefix)
                        self._save_metadata()
                        return {"fasta": cache_file, "index": index_prefix}
                    return {"fasta": cache_file}
                return {"fasta": cache_file}

            # Download and cache
            content = self._download_file(source)
            if content is None or not content.strip():
                return None

            with cache_file.open("w", encoding="utf-8") as f:
                f.write(content)

            if cache_file.stat().st_size == 0:
                cache_file.unlink(missing_ok=True)
                return None

            # Update metadata
            checksum = self._compute_file_checksum(cache_file)
            self.metadata[cache_key] = CacheMetadata(
                source=source,
                downloaded_at=datetime.now().isoformat(),
                file_size=cache_file.stat().st_size,
                checksum=checksum,
                file_path=str(cache_file),
            )

            logger.info(f"✅ Cached custom transcriptome: {cache_file} ({cache_file.stat().st_size:,} bytes)")

            # Build index if requested
            if build_index and self.auto_build_indices:
                if self._build_index(cache_file, index_prefix):
                    self._set_index_path(self.metadata[cache_key], index_prefix)
                    self._save_metadata()
                    return {"fasta": cache_file, "index": index_prefix}
                self._save_metadata()
                return {"fasta": cache_file}
            self._save_metadata()
            return {"fasta": cache_file}

        # Local file path
        input_path = Path(fasta_path)
        if not input_path.exists():
            logger.error(f"FASTA file not found: {input_path}")
            return None

        # Determine cache name and key
        if cache_name is None:
            cache_name = input_path.stem

        # Check if already in cache directory
        if input_path.parent == self.cache_dir:
            # Already cached - create metadata if missing
            cache_key = hashlib.md5(f"local_{cache_name}_{input_path}".encode()).hexdigest()[:12]
            cache_file = input_path
            index_prefix = input_path.parent / f"{input_path.stem}_index"

            # Create or update metadata
            if cache_key not in self.metadata:
                checksum = self._compute_file_checksum(cache_file)
                source = TranscriptomeSource(
                    name=cache_name,
                    url=str(input_path),
                    species="custom",
                    description=f"Local transcriptome from {input_path}",
                )

                self.metadata[cache_key] = CacheMetadata(
                    source=source,
                    downloaded_at=datetime.now().isoformat(),
                    file_size=cache_file.stat().st_size,
                    checksum=checksum,
                    file_path=str(cache_file),
                )
                self._save_metadata()

            # Build index if needed
            if build_index and self.auto_build_indices:
                if not index_prefix.with_suffix(".amb").exists():
                    if self._build_index(input_path, index_prefix):
                        self._set_index_path(self.metadata[cache_key], index_prefix)
                        self._save_metadata()
                        return {"fasta": input_path, "index": index_prefix}
                    return {"fasta": input_path}
                # Index exists - update metadata
                self._set_index_path(self.metadata[cache_key], index_prefix)
                self._save_metadata()
                return {"fasta": input_path, "index": index_prefix}
            return {"fasta": input_path}

        cache_key = hashlib.md5(f"local_{cache_name}_{input_path}".encode()).hexdigest()[:12]
        cache_file = self.cache_dir / f"{cache_name}_{cache_key}.fa"
        index_prefix = self.cache_dir / f"{cache_name}_{cache_key}_index"

        # Check if already cached
        if cache_file.exists():
            logger.info(f"✅ Using cached copy: {cache_file}")
        else:
            logger.info(f"🔄 Copying {input_path.name} to cache...")
            with input_path.open("r") as src, cache_file.open("w") as dst:
                dst.write(src.read())

        # Update metadata
        checksum = self._compute_file_checksum(cache_file)
        source = TranscriptomeSource(
            name=cache_name,
            url=str(input_path),
            species="custom",
            description=f"Local transcriptome from {input_path}",
        )

        self.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=checksum,
            file_path=str(cache_file),
        )

        # Build index if requested
        if build_index and self.auto_build_indices:
            if not index_prefix.with_suffix(".amb").exists():
                if self._build_index(cache_file, index_prefix):
                    self._set_index_path(self.metadata[cache_key], index_prefix)
                    self._save_metadata()
                    return {"fasta": cache_file, "index": index_prefix}
                self._save_metadata()
                return {"fasta": cache_file}
            self._set_index_path(self.metadata[cache_key], index_prefix)
            self._save_metadata()
            return {"fasta": cache_file, "index": index_prefix}
        self._save_metadata()
        return {"fasta": cache_file}

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
