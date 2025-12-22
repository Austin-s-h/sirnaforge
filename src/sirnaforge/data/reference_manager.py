#!/usr/bin/env python3
"""Base Reference Database Manager for shared caching functionality.

This module provides the foundation for managing reference databases (miRNA, transcriptome, genome)
with automatic caching, download management, and integrity validation.
"""

import gzip
import hashlib
import html
import json
import logging
import urllib.request
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Generic, TypeVar

from sirnaforge.utils.cache_utils import resolve_cache_subdir

logger = logging.getLogger(__name__)


@dataclass
class ReferenceSource:
    """Base configuration for a reference database source."""

    name: str
    url: str
    species: str
    format: str = "fasta"
    compressed: bool = False
    description: str = ""

    def cache_key(self) -> str:
        """Generate a unique cache key for this source."""
        content = f"{self.name}_{self.species}_{self.url}"
        return hashlib.md5(content.encode()).hexdigest()[:12]


@dataclass
class CacheMetadata:
    """Metadata for cached reference files."""

    source: ReferenceSource
    downloaded_at: str
    file_size: int
    checksum: str
    file_path: str
    version: str = "1.0"
    extra: dict[str, Any] | None = None  # For subclass-specific metadata

    @classmethod
    def from_dict(cls, data: dict, source_class: type = ReferenceSource) -> "CacheMetadata":
        """Create CacheMetadata from dictionary.

        Args:
            data: Dictionary containing metadata
            source_class: Class to use for source (allows subclasses)
        """
        source = source_class(**data["source"])
        return cls(
            source=source,
            downloaded_at=data["downloaded_at"],
            file_size=data["file_size"],
            checksum=data["checksum"],
            file_path=data["file_path"],
            version=data.get("version", "1.0"),
            extra=data.get("extra"),
        )

    def to_dict(self) -> dict:
        """Convert CacheMetadata to dictionary."""
        return asdict(self)


SourceT = TypeVar("SourceT", bound=ReferenceSource)


class ReferenceManager(ABC, Generic[SourceT]):
    """Base class for reference database managers with caching support.

    This provides common functionality for:
    - Cache directory management with multiple fallback locations
    - File downloading with retry and progress
    - Checksum validation
    - TTL-based cache invalidation
    - Metadata persistence

    Subclasses should implement source-specific operations like
    filtering, indexing, or multi-source aggregation.
    """

    def __init__(
        self,
        cache_subdir: str,
        cache_dir: str | Path | None = None,
        cache_ttl_days: int = 30,
    ):
        """Initialize the reference manager.

        Args:
            cache_subdir: Subdirectory name under cache root (e.g., 'mirna', 'transcriptomes')
            cache_dir: Directory for caching databases (default: ~/.cache/sirnaforge/{cache_subdir})
            cache_ttl_days: Cache time-to-live in days
        """
        self.cache_dir = self._resolve_cache_directory(cache_dir, cache_subdir)
        self.cache_ttl = timedelta(days=cache_ttl_days)
        self.metadata_file = self.cache_dir / "cache_metadata.json"
        self._load_metadata()

    def _resolve_cache_directory(self, cache_dir: str | Path | None, cache_subdir: str) -> Path:
        """Resolve cache directory with fallback locations.

        Single source of truth lives in `sirnaforge.utils.cache_utils.resolve_cache_subdir`.
        This method delegates to it for consistency across subsystems.
        """
        if cache_dir is not None:
            return resolve_cache_subdir(cache_subdir, override=cache_dir)
        return resolve_cache_subdir(cache_subdir)

    def _load_metadata(self) -> None:
        """Load cache metadata from disk."""
        self.metadata: dict[str, CacheMetadata] = {}

        if self.metadata_file.exists():
            try:
                with self.metadata_file.open("r") as f:
                    data = json.load(f)
                    for key, meta_dict in data.items():
                        self.metadata[key] = self._metadata_from_dict(meta_dict)
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Failed to load cache metadata: {e}")

    @abstractmethod
    def _metadata_from_dict(self, data: dict) -> CacheMetadata:
        """Create CacheMetadata from dictionary.

        Subclasses should override to use their specific source class.
        """
        pass

    def _save_metadata(self) -> None:
        """Save cache metadata to disk."""
        try:
            data = {key: meta.to_dict() for key, meta in self.metadata.items()}
            # Compute from cache_dir so callers/tests that monkeypatch cache_dir stay consistent.
            metadata_path = self.cache_dir / "cache_metadata.json"
            self.metadata_file = metadata_path
            with metadata_path.open("w") as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save cache metadata: {e}")

    def clear_cache(self, confirm: bool = False) -> dict[str, Any]:
        """Clear all cached files for this manager.

        By default this only reports what would be removed. Set `confirm=True`
        to actually delete files.
        """
        if not self.cache_dir.exists():
            return {
                "cache_directory": str(self.cache_dir),
                "files_deleted": 0,
                "size_freed_mb": 0.0,
                "status": "Cache directory does not exist",
            }

        files = [p for p in self.cache_dir.glob("*") if p.is_file()]
        total_size = sum(p.stat().st_size for p in files)
        if not confirm:
            return {
                "cache_directory": str(self.cache_dir),
                "files_deleted": len(files),
                "size_freed_mb": total_size / (1024 * 1024),
                "status": f"Would delete {len(files)} files ({total_size / (1024 * 1024):.2f} MB)",
            }

        deleted = 0
        for file_path in files:
            try:
                file_path.unlink()
                deleted += 1
            except FileNotFoundError:
                continue

        self.metadata.clear()
        return {
            "cache_directory": str(self.cache_dir),
            "files_deleted": deleted,
            "size_freed_mb": total_size / (1024 * 1024),
            "status": "Cache cleared successfully",
        }

    def _compute_file_checksum(self, file_path: Path) -> str:
        """Compute MD5 checksum of a file."""
        hash_md5 = hashlib.md5()
        with file_path.open("rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def _is_cache_valid(self, cache_key: str) -> bool:
        """Check if cached data is still valid.

        Args:
            cache_key: Cache key to validate

        Returns:
            True if cache is valid, False otherwise
        """
        if cache_key not in self.metadata:
            return False

        meta = self.metadata[cache_key]
        cache_file = Path(meta.file_path)

        # Check if file exists
        if not cache_file.exists():
            return False

        # Reject zero-byte cache entries
        if cache_file.stat().st_size == 0:
            logger.warning("Cache file %s is empty; marking as invalid", cache_file)
            return False

        # Check TTL
        downloaded_at = datetime.fromisoformat(meta.downloaded_at)
        if datetime.now() - downloaded_at > self.cache_ttl:
            return False

        # Check file integrity
        if self._compute_file_checksum(cache_file) != meta.checksum:
            logger.warning(f"Cache file {cache_file} corrupted, will re-download")
            return False

        return True

    def _download_file(self, source: ReferenceSource, timeout: int = 600) -> str | None:
        """Download file from source URL and return as text.

        Args:
            source: ReferenceSource configuration
            timeout: Download timeout in seconds

        Returns:
            File content as string, or None if download failed
        """
        try:
            logger.info(f"📥 Downloading {source.name} ({source.species}): {source.url}")

            request = urllib.request.Request(
                source.url,
                headers={
                    "User-Agent": "sirnaforge/1.0 (+https://github.com/austin-s-h/sirnaforge)",
                    "Accept": "text/plain,application/octet-stream",
                },
            )

            with urllib.request.urlopen(request, timeout=timeout) as response:
                data = response.read()

            if source.compressed and source.url.endswith(".gz"):
                logger.info("🔄 Decompressing gzipped file...")
                data = gzip.decompress(data)

            # Decode as text
            try:
                content: str = data.decode("utf-8")
            except UnicodeDecodeError:
                try:
                    content = data.decode("latin-1")
                except UnicodeDecodeError:
                    logger.error(f"❌ Cannot decode {source.url} as text")
                    return None

            # Fix HTML entities
            content = html.unescape(content)
            content = content.replace("<br>", "\n").replace("<BR>", "\n")

            if not content.strip():
                logger.error("Received empty response from %s", source.url)
                return None

            logger.info(f"✅ Downloaded {len(content):,} characters")
            return content

        except Exception as e:
            logger.error(f"❌ Failed to download {source.url}: {e}")
            return None

    def cache_info(self) -> dict[str, Any]:
        """Get information about the current cache state.

        Returns:
            Dictionary containing cache statistics
        """
        files = list(self.cache_dir.glob("*.fa")) + list(self.cache_dir.glob("*.fasta"))
        total_files = len(files)
        total_size = sum(f.stat().st_size for f in files if f.exists())

        return {
            "cache_directory": str(self.cache_dir),
            "total_files": total_files,
            "total_size_mb": total_size / (1024 * 1024),
            "cache_ttl_days": self.cache_ttl.days,
            "cached_items": list(self.metadata.keys()),
        }

    def clean_cache(self, older_than_days: int | None = None) -> None:
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
                cache_file = Path(meta.file_path)
                if cache_file.exists():
                    cache_file.unlink()
                    removed_count += 1

                # Allow subclasses to clean extra files
                self._clean_extra_files(meta)

                del self.metadata[cache_key]

        if removed_count > 0:
            self._save_metadata()
            logger.info(f"🧹 Cleaned {removed_count} old cache files")
        else:
            logger.info("🧹 No old cache files to clean")

    def _clean_extra_files(self, meta: CacheMetadata) -> None:
        """Clean extra files associated with a cache entry.

        Subclasses can override to clean index files, etc.

        Args:
            meta: Metadata for the cache entry being cleaned
        """
        pass
