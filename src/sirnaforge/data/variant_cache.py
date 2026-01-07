"""Improved variant caching using Parquet for efficient storage and retrieval."""

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

        # Initialize empty cache if it doesn't exist
        if not self.cache_file.exists():
            self._init_empty_cache()

    def _init_empty_cache(self) -> None:
        """Initialize an empty cache file."""
        empty_df = pd.DataFrame(
            columns=[
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
                "annotations",
                "provenance",
                "cached_at",
            ]
        )
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

            # Check TTL
            row = matches.iloc[0]
            cached_at = pd.to_datetime(row["cached_at"])
            age = datetime.now() - cached_at.to_pydatetime()

            if age > timedelta(days=self.ttl_days):
                logger.debug(f"Cache entry for {cache_key} is stale (age: {age.days} days)")
                # Don't delete here, let cleanup handle it
                return None

            # Reconstruct VariantRecord
            sources_str = str(row["sources"])
            annotations_str = str(row["annotations"])
            provenance_str = str(row["provenance"])

            sources = eval(sources_str) if sources_str != "nan" else []
            annotations = eval(annotations_str) if annotations_str != "nan" else {}
            provenance = eval(provenance_str) if provenance_str != "nan" else {}

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
            df = df[df["cache_key"] != cache_key]

            # Create new row
            new_row: dict[str, Any] = {
                "cache_key": cache_key,
                "id": variant.id,
                "chr": variant.chr,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "assembly": variant.assembly,
                "sources": str([s.value for s in variant.sources]),
                "clinvar_significance": variant.clinvar_significance.value if variant.clinvar_significance else None,
                "af": variant.af,
                "annotations": str(variant.annotations),
                "provenance": str(variant.provenance),
                "cached_at": datetime.now().isoformat(),
            }

            # Append new row
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

            # Write back to file
            df.to_parquet(self.cache_file, index=False, engine="pyarrow", compression="snappy")

            logger.debug(f"Cached variant with key {cache_key}")

        except Exception as e:
            logger.warning(f"Error writing to cache: {e}")

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

            # Convert cached_at to datetime
            df["cached_at"] = pd.to_datetime(df["cached_at"])
            cutoff_date = datetime.now() - timedelta(days=self.ttl_days)

            # Filter out stale entries
            df = df[df["cached_at"] > cutoff_date]

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

            df["cached_at"] = pd.to_datetime(df["cached_at"])
            cutoff_date = datetime.now() - timedelta(days=self.ttl_days)

            total = len(df)
            stale = len(df[df["cached_at"] <= cutoff_date])

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
