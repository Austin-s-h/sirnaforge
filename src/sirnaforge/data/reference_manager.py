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
import shutil
import urllib.request
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Generic, TypeVar

from sirnaforge.utils.cache_utils import (
    LEGACY_PRODUCER_VERSION,
    current_producer_version,
    is_producer_version_current,
    log_discarded_artifact,
    normalize_producer_version,
    resolve_cache_subdir,
)

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
    # Version of the code that *produced* this artifact, compared on load by
    # `ReferenceManager._is_cache_valid`. See `cache_utils.PRODUCER_VERSIONS`.
    version: str = LEGACY_PRODUCER_VERSION
    extra: dict[str, Any] | None = None  # For subclass-specific metadata

    @classmethod
    def from_dict(cls, data: dict[str, Any], source_class: type[ReferenceSource] = ReferenceSource) -> "CacheMetadata":
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
            version=normalize_producer_version(data.get("version")),
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
        artifact_class: str | None = None,
    ):
        """Initialize the reference manager.

        Args:
            cache_subdir: Subdirectory name under cache root (e.g., 'mirna', 'transcriptomes')
            cache_dir: Directory for caching databases (default: ~/.cache/sirnaforge/{cache_subdir})
            cache_ttl_days: Cache time-to-live in days
            artifact_class: Producer-version scope for this cache (default: `cache_subdir`).
                One subdirectory holds one kind of artifact, so the subdirectory name is
                also the natural invalidation scope; override only when two managers share
                a directory but not a producer.
        """
        self.artifact_class = artifact_class or cache_subdir
        self.producer_version = current_producer_version(self.artifact_class)
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
        self.uri_index: dict[str, str] = {}

        if self.metadata_file.exists():
            try:
                with self.metadata_file.open("r") as f:
                    data = json.load(f)
                    if isinstance(data, dict) and "entries" in data:
                        entries = data.get("entries", {})
                        for key, meta_dict in entries.items():
                            self.metadata[key] = self._metadata_from_dict(meta_dict)
                        raw_uri_index = data.get("uri_index", {})
                        if isinstance(raw_uri_index, dict):
                            self.uri_index = {str(uri): str(key) for uri, key in raw_uri_index.items()}
                    else:
                        for key, meta_dict in data.items():
                            self.metadata[key] = self._metadata_from_dict(meta_dict)
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Failed to load cache metadata: {e}")

        # Backfill URI index from metadata so older cache files upgrade in-place.
        for cache_key, meta in self.metadata.items():
            source_url = meta.source.url
            if self._is_remote_location(source_url):
                self.uri_index[source_url] = cache_key

    @abstractmethod
    def _metadata_from_dict(self, data: dict) -> CacheMetadata:
        """Create CacheMetadata from dictionary.

        Subclasses should override to use their specific source class.
        """
        pass

    def _save_metadata(self) -> None:
        """Save cache metadata to disk."""
        try:
            data = {
                "version": "2.0",
                "entries": {key: meta.to_dict() for key, meta in self.metadata.items()},
                "uri_index": self.uri_index,
            }
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
        self.uri_index.clear()
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

    def _record_cache_entry(
        self,
        cache_key: str,
        source: ReferenceSource,
        cache_file: Path,
        *,
        extra: dict[str, Any] | None = None,
        downloaded_at: str | None = None,
        persist: bool = False,
        producer_version: str | None = None,
    ) -> CacheMetadata:
        """Create/replace metadata for a cache entry from an on-disk file."""
        metadata = CacheMetadata(
            source=source,
            downloaded_at=downloaded_at or datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=self._compute_file_checksum(cache_file),
            file_path=str(cache_file),
            version=producer_version or self.producer_version,
            extra=extra,
        )
        self.metadata[cache_key] = metadata
        if self._is_remote_location(source.url):
            self.uri_index[source.url] = cache_key
        if persist:
            self._save_metadata()
        return metadata

    def _cache_key_for_remote_uri(self, uri: str) -> str | None:
        """Return an indexed cache key for a remote URI when available."""
        cache_key = self.uri_index.get(uri)
        if cache_key is None:
            return None
        if cache_key not in self.metadata:
            # Remove stale mappings so future lookups do not repeatedly miss.
            self.uri_index.pop(uri, None)
            return None
        return cache_key

    def _recover_remote_cache_entry(
        self,
        *,
        source: ReferenceSource,
        cache_key: str,
        cache_file: Path,
        persist: bool = True,
    ) -> bool:
        """Recover metadata for an already-downloaded remote artifact.

        This handles interrupted runs where the file exists on disk but metadata
        was never flushed, preventing unnecessary re-downloads.
        """
        if cache_key in self.metadata:
            return self._is_cache_valid(cache_key)

        if not cache_file.exists() or cache_file.stat().st_size == 0:
            return False

        # Recovery only ever adopts a raw download whose bytes came straight from
        # upstream, so no producer of ours shaped its content and stamping it with
        # the current version is honest. Derived artifacts must never be recovered
        # this way; they carry their own stamps (see cache_utils artifact stamps).
        self._record_cache_entry(cache_key, source, cache_file, persist=persist)
        logger.info(
            "Recovered remote cache metadata for %s from existing file %s",
            source.url,
            cache_file,
        )
        return self._is_cache_valid(cache_key)

    def is_producer_version_current(self, meta: CacheMetadata) -> bool:
        """Whether `meta` was written by the producer version this build ships.

        Public because it is the one question every cache consumer needs to be able
        to ask (and test) directly, independently of TTL and checksums.
        """
        return is_producer_version_current(self.artifact_class, meta.version)

    def stale_producer_entries(self) -> list[str]:
        """Cache keys whose recorded producer version is no longer current."""
        return [key for key, meta in self.metadata.items() if not self.is_producer_version_current(meta)]

    def _is_cache_valid(self, cache_key: str) -> bool:  # noqa: PLR0911
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

        # Reject artifacts written by a producer we have since fixed. This has to
        # come before the checksum check, which compares the file against the md5
        # recorded when that same file was written and therefore can only detect
        # corruption after the fact, never wrong content the writer itself put there.
        if not self.is_producer_version_current(meta):
            log_discarded_artifact(
                self.artifact_class,
                cache_file,
                f"producer version {normalize_producer_version(meta.version)} != {self.producer_version}",
            )
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

    def _download_to_path(self, source: ReferenceSource, destination: Path, timeout: int = 600) -> bool:
        """Stream a source URL directly to destination, optionally decompressing gzip content."""
        try:
            logger.info(f"📥 Downloading {source.name} ({source.species}): {source.url}")
            request = urllib.request.Request(
                source.url,
                headers={
                    "User-Agent": "sirnaforge/1.0 (+https://github.com/austin-s-h/sirnaforge)",
                    "Accept": "text/plain,application/octet-stream",
                },
            )

            destination.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                if source.compressed and source.url.endswith(".gz"):
                    logger.info("🔄 Decompressing gzipped file...")
                    with gzip.GzipFile(fileobj=response) as decompressed, destination.open("wb") as handle:
                        shutil.copyfileobj(decompressed, handle, length=1024 * 1024)
                else:
                    with destination.open("wb") as handle:
                        shutil.copyfileobj(response, handle, length=1024 * 1024)

            if destination.stat().st_size == 0:
                logger.error("Received empty response from %s", source.url)
                destination.unlink(missing_ok=True)
                return False

            logger.info("✅ Downloaded %s bytes to %s", destination.stat().st_size, destination)
            return True
        except Exception as e:
            logger.error(f"❌ Failed to download {source.url}: {e}")
            destination.unlink(missing_ok=True)
            return False

    def _is_remote_location(self, resource_location: str) -> bool:
        """Return whether a user resource string points to a remote URL."""
        return resource_location.startswith(("http://", "https://", "ftp://"))

    def _default_cache_name_for_resource(self, resource_location: str) -> str:
        """Derive a stable cache name from a user-supplied local path or URL."""
        name = resource_location.rstrip("/").split("/")[-1] if resource_location else "resource"
        normalized = name.replace(".gz", "").replace(".fa", "").replace(".fasta", "")
        return normalized or "resource"

    def _cache_info_files(self) -> list[Path]:
        """Return cache files included in aggregate size/count stats."""
        return list(self.cache_dir.glob("*.fa")) + list(self.cache_dir.glob("*.fasta"))

    def _cache_info_extra(self, _total_files: int, _total_size_bytes: int) -> dict[str, Any]:
        """Return subclass-specific fields for cache_info."""
        return {}

    def cache_info(self) -> dict[str, Any]:
        """Get information about the current cache state.

        Returns:
            Dictionary containing cache statistics
        """
        files = self._cache_info_files()
        total_files = len(files)
        total_size = sum(f.stat().st_size for f in files if f.exists())
        info: dict[str, Any] = {
            "cache_directory": str(self.cache_dir),
            "total_files": total_files,
            "total_size_mb": total_size / (1024 * 1024),
            "cache_ttl_days": self.cache_ttl.days,
            "cached_items": list(self.metadata.keys()),
            "artifact_class": self.artifact_class,
            "producer_version": self.producer_version,
            # Surfaced so `sirnaforge cache --info` can explain an upcoming refresh.
            "stale_producer_entries": self.stale_producer_entries(),
        }
        info.update(self._cache_info_extra(total_files, total_size))
        return info

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
