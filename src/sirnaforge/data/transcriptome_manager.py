#!/usr/bin/env python3
"""Transcriptome Database Manager with local caching and automatic index building.

This module provides a clean interface for downloading, caching, and managing
transcriptome FASTA files with automatic BWA-MEM2 index building and cache management.
"""

import hashlib
import logging
import shutil
import uuid
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from sirnaforge.utils.cache_utils import (
    ARTIFACT_STAMP_SUFFIX,
    ARTIFACT_TRANSCRIPTOME_INDEX,
    discard_artifact_stamp,
    is_artifact_stamp_current,
    log_discarded_artifact,
    read_artifact_stamp,
    write_artifact_stamp,
)

from .ensembl_references import build_transcriptome_sources
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

    # File suffixes written by `bwa-mem2 index` for a given prefix. Deliberately a
    # superset of what off_target.validate_index_files() requires: `.0123` is not needed
    # to *use* an index but bwa-mem2 does write it (see modules/local/build_bwa_index.nf),
    # and it is one of the largest files, so cleanup must cover it.
    INDEX_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")

    # Appended to the base download URL to give a filtered artifact its own cache
    # identity. Filtering is a subset operation, so a filtered FASTA must never be
    # reachable through the URI of the unfiltered download it was cut from.
    FILTER_URI_FRAGMENT = "#filters="

    # `extra` key naming an index prefix whose files outlived the content they were
    # built from and could not be deleted. Persisted, so a later process that finds the
    # cache entry valid (and therefore never re-reconciles) still refuses the leftovers.
    UNREMOVABLE_INDEX_KEY = "unremovable_stale_index"

    # Common transcriptome sources, generated from the shared Ensembl assembly table
    # (sirnaforge.data.ensembl_references) so cDNA and genome references stay in lockstep
    # and adding a species is a single table entry. Keys/URLs are unchanged from the
    # previously hand-written entries (ensembl_human_cdna, ensembl_mouse_cdna, ...).
    SOURCES = build_transcriptome_sources()

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

    def _cache_key_for_remote_uri(self, uri: str) -> str | None:
        """Resolve a remote URI to a cache key, never to a derived (filtered) artifact.

        Deciding this at resolution time rather than by rewriting `uri_index` on load is
        deliberate: `ReferenceManager._load_metadata` backfills the index from each
        entry's own `source.url`, so metadata written before filtered entries got their
        own filter-qualified URI re-adds the bad mapping on *every* load and a one-off
        repair would only hold for one process. A filtered FASTA is a subset of the base
        download, so only a URI carrying the filter fragment may resolve to it - an
        unfiltered request resolving here would screen off-targets against a subset and
        report false negatives.
        """
        cache_key = super()._cache_key_for_remote_uri(uri)
        if cache_key is None:
            return None

        meta = self.metadata.get(cache_key)
        if meta is not None and (meta.extra or {}).get("filters") and self.FILTER_URI_FRAGMENT not in uri:
            return None
        return cache_key

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

    @staticmethod
    def _cache_dir_is_writable(cache_dir: Path) -> bool:
        """Return whether a cache directory can be created and written to."""
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            probe = cache_dir / f".sirnaforge_write_probe_{uuid.uuid4().hex}"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink(missing_ok=True)
            return True
        except OSError:
            return False

    def _resolve_remote_cache_file(self, cache_key: str, suffix: str = ".fa") -> Path:
        """Resolve a writable remote-cache destination for a cache key.

        Stale metadata can reference old host paths (for example from a prior
        non-container run). When that parent directory is no longer writable,
        re-home new downloads to the active manager cache directory.
        """
        default_cache_file = self.cache_dir / f"{cache_key}{suffix}"
        if cache_key not in self.metadata:
            return default_cache_file

        metadata_cache_file = Path(self.metadata[cache_key].file_path)
        if self._is_cache_valid(cache_key):
            return metadata_cache_file

        if self._cache_dir_is_writable(metadata_cache_file.parent):
            return metadata_cache_file

        logger.warning(
            "Transcriptome metadata path '%s' is not writable; re-homing download to '%s'",
            metadata_cache_file,
            default_cache_file,
        )
        return default_cache_file

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

    @staticmethod
    def _index_outputs(index_prefix: Path) -> list[Path]:
        """The files a BWA-MEM2 index actually consists of, for output fingerprinting.

        The prefix itself is not a file (only a marker `_ensure_index_marker` touches),
        so the members are discovered by glob rather than from a hard-coded extension
        list: bwa-mem2 releases differ in which suffixes they emit, and a stamp that
        only covered the four we happen to name would let the others rot unnoticed.
        The sidecar stamp is excluded - it cannot fingerprint itself.
        """
        return sorted(
            path
            for path in index_prefix.parent.glob(f"{index_prefix.name}.*")
            if path.is_file() and not path.name.endswith(ARTIFACT_STAMP_SUFFIX)
        )

    @staticmethod
    def _index_inputs(meta: CacheMetadata) -> dict[str, str]:
        """Fingerprint of the content an index must have been built from.

        Uses the checksum already recorded for the cached FASTA rather than hashing
        it again: a transcriptome is tens of GB and the digest is computed anyway
        when the file is cached.
        """
        return {"fasta": meta.checksum}

    def _index_is_trustworthy(self, meta: CacheMetadata, index_path: Path) -> bool:
        """Whether an existing index may be reused for the currently cached content.

        bwa-mem2 aligns against the index files alone and never re-reads the FASTA,
        so an index that outlived the content it was built from silently reports hits
        with the previous release's transcript IDs and coordinates.
        """
        inputs = self._index_inputs(meta)
        if read_artifact_stamp(index_path) is not None:
            # Stamped: the stamp is the whole answer (and logs the mismatch itself).
            return is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, index_path, inputs=inputs)

        # Unstamped index from an older siRNAforge. Rebuilding a human transcriptome
        # index costs an hour of CPU, so do not force that on everyone; our own
        # metadata already answers the question. An index recorded as built *after*
        # the cached content was written must have been built from that content
        # (_record_cache_entry rewrites downloaded_at whenever content is replaced),
        # so adopt it and stamp it. Anything we cannot date is discarded.
        built_at = (meta.extra or {}).get("index_built_at")
        if isinstance(built_at, str) and built_at >= meta.downloaded_at:
            write_artifact_stamp(
                ARTIFACT_TRANSCRIPTOME_INDEX,
                index_path,
                inputs=inputs,
                outputs=self._index_outputs(index_path),
            )
            logger.info("Adopted pre-existing BWA-MEM2 index %s (built after the cached content)", index_path)
            return True

        log_discarded_artifact(
            ARTIFACT_TRANSCRIPTOME_INDEX,
            index_path,
            "it cannot be shown to have been built from the cached reference",
        )
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

    def _remove_index_files(self, index_prefix: Path) -> None:
        """Delete the marker and all bwa-mem2 files sharing an index prefix.

        This now runs on the hot path of every re-download, where the cache parent is
        not guaranteed writable (the reason `_resolve_remote_cache_file` exists). A
        stale index we cannot delete must not abort the run that replaced its FASTA, so
        failures are reported per file and the caller re-checks the prefix.
        """
        candidates = [index_prefix, *(index_prefix.parent / f"{index_prefix.name}{ext}" for ext in self.INDEX_SUFFIXES)]
        for candidate in candidates:
            try:
                candidate.unlink(missing_ok=True)
            except OSError as e:
                logger.warning("Could not remove stale index file %s: %s", candidate, e)

    def _reconcile_index_for_new_content(
        self, previous: CacheMetadata | None, current: CacheMetadata, index_prefix: Path
    ) -> None:
        """Carry over or destroy the BWA index after cached content is rewritten.

        Index files are named from the cache key, which hashes a release-agnostic URL
        (ENSEMBL_FTP_BASE points at `current_fasta`), so a refreshed reference lands on
        the same prefix. bwa-mem2 aligns against the index alone and never reads the
        FASTA, so an index left over from the previous content would silently report
        hits with the previous release's transcript IDs and coordinates. Identical bytes
        are the only case where the existing index remains trustworthy.
        """
        if previous is not None and previous.checksum == current.checksum:
            # `_record_cache_entry` replaces metadata wholesale; keep index bookkeeping.
            merged = {**(previous.extra or {}), **(current.extra or {})}
            current.extra = merged or None
            return

        stale_prefixes = {index_prefix}
        previous_index = self._get_index_path(previous) if previous is not None else None
        if previous_index is not None:
            stale_prefixes.add(previous_index)

        for prefix in stale_prefixes:
            if self._is_index_complete(prefix):
                logger.info("🗑️  Discarding BWA-MEM2 index built from replaced content: %s", prefix)
            self._remove_index_files(prefix)

        # An index we failed to delete still looks complete to _is_index_complete, so
        # record it as unusable. Escalating and then letting _prepare_result_with_index
        # hand the same files back would be the corruption the escalation warns about.
        if self._is_index_complete(index_prefix):
            logger.error(
                "❌ Stale BWA-MEM2 index at %s could not be removed and no longer matches the cached "
                "reference; it will not be used. Delete %s.* or run `sirnaforge cache "
                "--clear-transcriptome` to get a cached index back.",
                index_prefix,
                index_prefix,
            )
            if current.extra is None:
                current.extra = {}
            current.extra[self.UNREMOVABLE_INDEX_KEY] = str(index_prefix)
            # Belt and braces: an index_path carried in from anywhere must not outlive
            # the content it described.
            current.extra.pop("index_path", None)

    def _unremovable_index_blocks_reuse(self, meta: CacheMetadata) -> bool:
        """Whether this entry's index prefix is quarantined leftovers we must not serve.

        Set by `_reconcile_index_for_new_content` when a stale index outlived the content
        it was built from and the filesystem refused to delete it. Rebuilding onto the
        same prefix is not an option either: a directory that rejects `unlink` rejects
        the writes too, and a partly-overwritten prefix would still pass
        `_is_index_complete` while mixing old and new files. So the run proceeds with no
        index until the leftovers are gone, which is the one thing that clears this.
        """
        quarantined = (meta.extra or {}).get(self.UNREMOVABLE_INDEX_KEY)
        if not quarantined:
            return False

        if self._is_index_complete(Path(quarantined)):
            logger.error(
                "🚫 Not using BWA-MEM2 index %s: it was built from reference content that has since "
                "been replaced and could not be deleted. Proceeding without a cached index.",
                quarantined,
            )
            return True

        logger.info("♻️  Stale BWA-MEM2 index %s is gone; the prefix can be rebuilt.", quarantined)
        if meta.extra is not None:
            meta.extra.pop(self.UNREMOVABLE_INDEX_KEY, None)
        return False

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
        cache_file = self._resolve_remote_cache_file(cache_key)
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
        previous = self.metadata.get(cache_key)
        if not self._download_to_path(source, cache_file):
            return None

        # Update metadata
        current = self._record_cache_entry(cache_key, source, cache_file)
        self._reconcile_index_for_new_content(previous, current, index_prefix)
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
        cache_file = self._resolve_remote_cache_file(cache_key)
        index_prefix = self.cache_dir / f"{cache_key}_index"

        # Check cache
        if self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached custom transcriptome: {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        if self._recover_remote_cache_entry(source=source, cache_key=cache_key, cache_file=cache_file):
            logger.info(f"✅ Using recovered cached custom transcriptome: {cache_file}")
            return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

        # Download
        previous = self.metadata.get(cache_key)
        if not self._download_to_path(source, cache_file):
            return None

        # Save metadata
        current = self._record_cache_entry(cache_key, source, cache_file)
        self._reconcile_index_for_new_content(previous, current, index_prefix)
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
        # Name the index from the cache key, not the file stem: the stem is stable while
        # the file's contents are not, so a rewritten FASTA sitting in the cache dir used
        # to be handed the previous content's index. bwa-mem2 never reads the FASTA back,
        # so that mismatch surfaces only as hits with the previous content's IDs.
        index_prefix = file_path.parent / f"{cache_key}_index"

        # Re-record on every call. The recorded checksum is what tells us whether the
        # index still matches these bytes, and with dedupe disabled the cache key is
        # path-derived, so an existing entry can describe an older revision of the file.
        previous = self.metadata.get(cache_key)
        current = self._record_cache_entry(
            cache_key,
            TranscriptomeSource(
                name=cache_name, url=str(file_path), species="custom", description=f"Local: {file_path}"
            ),
            file_path,
            extra=extra,
        )
        self._reconcile_index_for_new_content(previous, current, index_prefix)

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
            # The cache file is content-addressed, so an existing copy holds these bytes.
            copy_needed = not cache_file.exists()
        else:
            cache_key = hashlib.md5(f"local_{cache_name}_{input_path}".encode()).hexdigest()[:12]
            cache_file = self.cache_dir / f"{cache_name}_{cache_key}.fa"
            index_prefix = self.cache_dir / f"{cache_name}_{cache_key}_index"
            # Here the cache key is derived from the *path*, so the cached copy can be an
            # older revision of it. Existence alone would pin the first version we ever
            # saw and ignore every later edit to the user's FASTA.
            copy_needed = not cache_file.exists() or self._compute_file_checksum(cache_file) != (
                self._compute_file_checksum(input_path)
            )

        if copy_needed:
            logger.info(f"🔄 Copying {input_path.name} to cache...")
            shutil.copyfile(input_path, cache_file)

        # Save metadata, then invalidate an index built from the copy we just replaced.
        previous = self.metadata.get(cache_key)
        current = self._record_cache_entry(
            cache_key,
            TranscriptomeSource(
                name=cache_name, url=str(input_path), species="custom", description=f"Local: {input_path}"
            ),
            cache_file,
            extra=extra,
        )
        self._reconcile_index_for_new_content(previous, current, index_prefix)

        return self._prepare_result_with_index(cache_file, index_prefix, cache_key, build_index)

    def _base_content_identity(self, source: TranscriptomeSource, base_fasta: Path) -> tuple[str, str]:
        """Return the (cache key, checksum) of the base download a filter derives from."""
        base_cache_key = self._cache_key_for_remote_uri(source.url) or source.cache_key()
        base_meta = self.metadata.get(base_cache_key)
        # get_transcriptome() has just validated (checksum-verified) or rewritten this
        # entry, so its recorded checksum matches the bytes on disk; only hash when the
        # entry is missing or points somewhere else.
        if base_meta is not None and Path(base_meta.file_path) == base_fasta:
            return base_cache_key, base_meta.checksum
        return base_cache_key, self._compute_file_checksum(base_fasta)

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

        # A filtered FASTA is derived, so its identity is the base bytes it was cut from.
        # The Ensembl URL is release-agnostic, so a new release reuses the base cache key
        # and file name; recording the base checksum is what lets an existing filtered
        # entry (which has no download of its own to TTL out) notice that.
        base_cache_key, base_checksum = self._base_content_identity(source, base_result["fasta"])

        # Check cache: valid on its own terms (present, non-empty, in-TTL, uncorrupted)
        # *and* still derived from the base content we just resolved.
        if not force_refresh and self._is_cache_valid(filtered_cache_key):
            cached_meta = self.metadata[filtered_cache_key]
            if (cached_meta.extra or {}).get("base_checksum") == base_checksum:
                cached_path = Path(cached_meta.file_path)
                return self._prepare_result_with_index(cached_path, filtered_index, filtered_cache_key, build_index)
            logger.info(
                "♻️  %s changed since it was filtered (%s); re-filtering and discarding the derived index",
                source_name,
                filter_spec,
            )

        # Apply filters
        logger.info(f"🔍 Applying filters to {source_name}: {', '.join(filters)}")
        kept = TranscriptFilter.apply_combined_filter(base_result["fasta"], filtered_fasta, filters)
        if kept == 0:
            filtered_fasta.unlink(missing_ok=True)
            return None

        # Cache metadata. The filtered artifact must not be recorded under the bare
        # source URL: that URL indexes the *unfiltered* download, and a later
        # get_transcriptome() would resolve to this filtered FASTA and screen
        # off-targets against a subset of the transcriptome. `base_checksum` records
        # which base bytes these are a subset of.
        previous = self.metadata.get(filtered_cache_key)
        current = self._record_cache_entry(
            filtered_cache_key,
            TranscriptomeSource(
                name=f"{source.name}_filtered",
                url=f"{source.url}{self.FILTER_URI_FRAGMENT}{filter_spec}",
                species=source.species,
                description=f"{source.description} [filtered: {filter_spec}]",
            ),
            filtered_fasta,
            extra={
                "filters": filters,
                "kept_count": kept,
                "base_cache_key": base_cache_key,
                "base_checksum": base_checksum,
            },
        )
        self._reconcile_index_for_new_content(previous, current, filtered_index)

        return self._prepare_result_with_index(filtered_fasta, filtered_index, filtered_cache_key, build_index)

    def _prepare_result_with_index(
        self, fasta: Path, index_prefix: Path, cache_key: str, build_index: bool
    ) -> dict[str, Path]:
        """Helper to prepare result dict with optional index building."""
        if not (build_index and self.auto_build_indices):
            self._save_metadata()
            return {"fasta": fasta}

        meta = self.metadata[cache_key]

        # Checked before the completeness test below, which cannot tell a fresh index
        # from leftovers of the content this entry replaced.
        if self._unremovable_index_blocks_reuse(meta):
            self._save_metadata()
            return {"fasta": fasta}

        index_path = self._get_index_path(meta) or index_prefix

        if self._is_index_complete(index_path) and self._index_is_trustworthy(meta, index_path):
            self._ensure_index_marker(index_path)
            logger.info(f"✅ Using cached BWA-MEM2 index: {index_path}")
            return {"fasta": fasta, "index": index_path}

        logger.info(f"⚠️  Building index: {index_prefix}")
        # Drop any stamp first so a build that dies half-way cannot leave a stamp
        # vouching for a partial index.
        discard_artifact_stamp(index_prefix)
        if self._build_index(fasta, index_prefix):
            self._ensure_index_marker(index_prefix)
            self._set_index_path(meta, index_prefix)
            write_artifact_stamp(
                ARTIFACT_TRANSCRIPTOME_INDEX,
                index_prefix,
                inputs=self._index_inputs(meta),
                outputs=self._index_outputs(index_prefix),
            )
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

        # `_remove_index_files` supersedes the inline extension list it replaces: it covers `.0123`
        # (which bwa-mem2 writes), builds paths by concatenation so dotted prefixes work where
        # `with_suffix` silently did not, and guards each unlink because this runs on the hot path
        # of every re-download.
        self._remove_index_files(index_path)
        # Leave no stamp behind to vouch for an index that is gone.
        discard_artifact_stamp(index_path)
