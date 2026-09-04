"""Improved variant caching using Parquet for efficient storage and retrieval."""

import json
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import pandas as pd

from sirnaforge.models.variant import VariantRecord, VariantSource
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class VariantParquetCache:
    """Efficient variant cache using Parquet files for better performance than JSON.

    Benefits over JSON:
    - Columnar storage format is much more efficient for variant data
    - Built-in compression reduces disk usage
    - Fast filtering and querying with pandas
    - Batch operations instead of individual file I/O
    """

    # Canonical on-disk schema: every column this version reads and writes, one per
    # persisted ``VariantRecord`` field plus the two bookkeeping columns. ``put``
    # writes all of them rather than whatever the existing file happens to have, so
    # a cache written by an older version gains new columns instead of silently
    # dropping them from new rows; columns it does not recognise are carried
    # through untouched (see ``put``). Keep in step with ``VariantRecord`` --
    # ``test_columns_match_variant_record_fields`` fails if the two drift.
    _COLUMNS = (
        "cache_key",
        "id",
        "chr",
        "pos",
        "ref",
        "alt",
        "assembly",
        "sources",
        "clinvar_significance",
        "af",
        "population_afs",
        "annotations",
        "provenance",
        "cached_at",
    )

    def __init__(self, cache_dir: Path, ttl_days: int = 90):
        """Initialize the Parquet-based variant cache.

        Args:
            cache_dir: Directory for cache storage
            ttl_days: Time-to-live for cached entries in days (default: 90)
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.ttl_days = ttl_days
        self.cache_file = self.cache_dir / "variants.parquet"
        self._write_failure_reported = False

        # Initialize empty cache if it doesn't exist
        if not self.cache_file.exists():
            self._init_empty_cache()

    @staticmethod
    def _serialize_value(value: Any) -> str:
        """Serialize Python values to JSON for stable, safe storage."""
        return json.dumps(value, sort_keys=True)

    @staticmethod
    def _deserialize_value(value: Any, default: Any) -> Any:
        """Deserialize cached JSON payloads."""
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return default
        if isinstance(value, str):
            return json.loads(value)
        return value

    @staticmethod
    def _normalize_timestamp(value: Any) -> Any:
        """Coerce a timestamp value to the ISO string form used on disk."""
        # ``pd.NaT`` subclasses ``datetime`` and its isoformat() is the literal
        # string "NaT", which round-trips back to NaT and made the TTL readers
        # disagree about the row, so it has to be caught before the datetime
        # branch. Persist a real null instead: every reader treats an unknown
        # timestamp the same way (stale).
        if value is pd.NaT:
            return None
        if isinstance(value, datetime):
            return value.isoformat()
        return value

    @staticmethod
    def _parse_cached_at(values: Any) -> Any:
        """Parse ``cached_at`` values into datetimes, tolerating mixed ISO layouts.

        ``datetime.now().isoformat()`` omits the ``.%f`` suffix when microsecond is
        exactly 0, so one column can legitimately hold two ISO layouts; pandas 2.x
        infers a single format from the first element and raises on the rest, which
        silently disabled both TTL readers. ``format="ISO8601"`` accepts either
        layout and ``errors="coerce"`` maps anything still unreadable (including a
        null written by ``_normalize_timestamp``) to NaT.
        """
        return pd.to_datetime(values, format="ISO8601", errors="coerce")

    def _fresh_mask(self, df: pd.DataFrame) -> pd.Series:
        """Return a boolean mask of entries still inside the TTL.

        ``cleanup_stale_entries`` and ``get_stats`` share this so they can never
        disagree about a row: NaT comparisons are always False, so a row with an
        unreadable timestamp is stale for both.
        """
        cutoff_date = datetime.now() - timedelta(days=self.ttl_days)
        fresh: pd.Series = self._parse_cached_at(df["cached_at"]) > cutoff_date
        return fresh

    def _is_servable(self, frame: pd.DataFrame, row: "pd.Series[Any]", cache_key: str) -> bool:
        """Decide whether a matched row may be handed back to the caller."""
        # An unreadable timestamp parses to NaT and counts as stale here too, so
        # get, cleanup_stale_entries and get_stats agree about the same row.
        cached_at = self._parse_cached_at(row["cached_at"])
        if pd.isna(cached_at):
            logger.debug(f"Cache entry for {cache_key} has an unreadable cached_at; treating as stale")
            return False

        age = datetime.now() - cached_at.to_pydatetime()
        if age > timedelta(days=self.ttl_days):
            logger.debug(f"Cache entry for {cache_key} is stale (age: {age.days} days)")
            # Don't delete here, let cleanup handle it
            return False

        # ``population_afs`` was not persisted at all before this fix, and the
        # current writer always stores at least "{}". A null therefore means
        # "never written", not "this variant has no population AFs" -- and
        # avoid-mode filtering keys off exactly those AFs, so serving the row
        # would silently reinstate the bug for the 90-day TTL of every cache that
        # already exists. Report a miss so the caller re-fetches.
        population_afs = row["population_afs"] if "population_afs" in frame.columns else None
        if population_afs is None or (isinstance(population_afs, float) and pd.isna(population_afs)):
            logger.info(
                f"Discarding pre-0.6.0 cache entry for {cache_key}: population allele frequencies "
                "were never stored, so the entry cannot be trusted for avoid-mode filtering"
            )
            return False

        return True

    def _init_empty_cache(self) -> None:
        """Initialize an empty cache file."""
        empty_df = pd.DataFrame(columns=list(self._COLUMNS))
        empty_df.to_parquet(self.cache_file, index=False, engine="pyarrow", compression="snappy")
        logger.info(f"Initialized empty variant cache at {self.cache_file}")

    def get(self, cache_key: str) -> VariantRecord | None:
        """Retrieve a variant from cache by key.

        Args:
            cache_key: Cache key for the variant

        Returns:
            VariantRecord if found and not stale, None otherwise
        """
        try:
            df = pd.read_parquet(self.cache_file, engine="pyarrow")

            if df.empty:
                return None

            # Filter by cache key
            matches = df[df["cache_key"] == cache_key]

            if matches.empty:
                return None

            # Check TTL and that the row was written by a version that stored
            # everything the caller will read back.
            row = matches.iloc[0]
            if not self._is_servable(matches, row, cache_key):
                return None

            # Reconstruct VariantRecord
            sources = self._deserialize_value(row["sources"], [])
            population_afs = self._deserialize_value(row["population_afs"], {})
            annotations = self._deserialize_value(row["annotations"], {})
            provenance = self._deserialize_value(row["provenance"], {})

            variant = VariantRecord(
                id=row["id"] if pd.notna(row["id"]) else None,
                chr=str(row["chr"]),
                pos=int(row["pos"]),
                ref=str(row["ref"]),
                alt=str(row["alt"]),
                assembly=str(row["assembly"]),
                sources=[VariantSource(s) for s in sources],
                clinvar_significance=row["clinvar_significance"] if pd.notna(row["clinvar_significance"]) else None,
                af=float(row["af"]) if pd.notna(row["af"]) else None,
                population_afs=population_afs,
                annotations=annotations,
                provenance=provenance,
            )

            logger.debug(f"Cache hit for {cache_key}")
            return variant

        except Exception as e:
            logger.warning(f"Error reading from cache: {e}")
            return None

    def put(self, cache_key: str, variant: VariantRecord) -> None:
        """Store a variant in the cache.

        Args:
            cache_key: Cache key for storage
            variant: VariantRecord to cache
        """
        try:
            # Read existing cache
            df = pd.read_parquet(self.cache_file, engine="pyarrow")

            # Remove existing entry with same key if present
            kept = df[df["cache_key"] != cache_key]

            # Caches written by an older cleanup_stale_entries hold datetime64 here;
            # normalise so the rebuilt column is one parquet-writable type.
            if "cached_at" in kept.columns:
                kept = kept.assign(cached_at=kept["cached_at"].map(self._normalize_timestamp))

            # Create new row
            new_row: dict[str, Any] = {
                "cache_key": cache_key,
                "id": variant.id,
                "chr": variant.chr,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "assembly": variant.assembly,
                "sources": self._serialize_value([s.value for s in variant.sources]),
                "clinvar_significance": variant.clinvar_significance.value if variant.clinvar_significance else None,
                "af": variant.af,
                "population_afs": self._serialize_value(variant.population_afs),
                "annotations": self._serialize_value(variant.annotations),
                "provenance": self._serialize_value(variant.provenance),
                "cached_at": datetime.now().isoformat(),
            }

            # Rebuild the frame column-wise rather than appending in place. An
            # in-place append had to pick an index label (masking above preserves
            # the original labels, so the next free position is already taken and
            # the write clobbered an unrelated entry) and had to fit the existing
            # column dtypes (a datetime64 ``cached_at`` rejected the ISO string and
            # left the cache permanently unwritable). Rebuilding sidesteps both,
            # and also avoids the pandas FutureWarning that concat emits on
            # empty/all-NA frames.
            #
            # Columns the canonical schema does not know about are carried through
            # rather than dropped: with a 90-day TTL one cache directory is easily
            # shared with a newer sirnaforge, and dropping its columns would destroy
            # data the other version still needs. The new row simply has no value
            # for them.
            unknown_columns = [column for column in kept.columns if column not in self._COLUMNS]
            df = pd.DataFrame(
                {
                    column: [
                        *(kept[column].tolist() if column in kept.columns else [None] * len(kept)),
                        new_row.get(column),
                    ]
                    for column in (*self._COLUMNS, *unknown_columns)
                }
            )

            # Write back to file
            df.to_parquet(self.cache_file, index=False, engine="pyarrow", compression="snappy")

            logger.debug(f"Cached variant with key {cache_key}")

        except Exception as e:
            # A cache that cannot be written stays broken for the rest of the run,
            # so surface the first failure loudly and only trace the repeats.
            if self._write_failure_reported:
                logger.debug(f"Error writing to cache: {e}")
            else:
                self._write_failure_reported = True
                logger.warning(f"Error writing to cache, variant cache will not be updated: {e}", exc_info=True)

    def cleanup_stale_entries(self) -> int:
        """Remove entries older than TTL.

        Returns:
            Number of entries removed
        """
        try:
            df = pd.read_parquet(self.cache_file, engine="pyarrow")

            if df.empty:
                return 0

            original_count = len(df)

            # Compare on a temporary mask and leave ``cached_at`` untouched: writing
            # datetime64 back into it would break every later put, which appends an
            # ISO string. Whatever type the file already holds is written back
            # unchanged here; the next ``put`` is what normalises it to ISO strings.
            df = df[self._fresh_mask(df)].reset_index(drop=True)

            # Write back
            df.to_parquet(self.cache_file, index=False, engine="pyarrow", compression="snappy")

            removed = original_count - len(df)
            if removed > 0:
                logger.info(f"Cleaned up {removed} stale cache entries")

            return removed

        except Exception as e:
            logger.warning(f"Error cleaning up cache: {e}")
            return 0

    def get_stats(self) -> dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        try:
            df = pd.read_parquet(self.cache_file, engine="pyarrow")

            if df.empty:
                return {"total_entries": 0, "stale_entries": 0}

            total = len(df)
            # Complement of the mask cleanup uses, so the two always agree.
            stale = int((~self._fresh_mask(df)).sum())

            return {
                "total_entries": total,
                "stale_entries": stale,
                "cache_file": str(self.cache_file),
                "cache_size_mb": self.cache_file.stat().st_size / (1024 * 1024) if self.cache_file.exists() else 0,
            }

        except Exception as e:
            logger.warning(f"Error getting cache stats: {e}")
            return {"total_entries": 0, "stale_entries": 0, "error": str(e)}

    def clear(self) -> None:
        """Clear all cache entries."""
        self._init_empty_cache()
        logger.info("Cleared variant cache")
