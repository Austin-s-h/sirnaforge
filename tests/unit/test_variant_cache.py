"""Unit tests for Parquet-based variant cache."""

import logging
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.data.variant_cache import VariantParquetCache
from sirnaforge.models.variant import ClinVarSignificance, VariantMode, VariantRecord, VariantSource


class TestVariantParquetCache:
    """Tests for VariantParquetCache."""

    def test_init_creates_empty_cache(self, tmp_path: Path):
        """Test that initialization creates an empty cache file."""
        cache = VariantParquetCache(tmp_path)

        assert cache.cache_file.exists()
        assert cache.cache_dir == tmp_path

    def test_put_and_get(self, tmp_path: Path):
        """Test storing and retrieving a variant."""
        cache = VariantParquetCache(tmp_path)

        variant = VariantRecord(
            id="rs1234",
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.CLINVAR],
            clinvar_significance=ClinVarSignificance.PATHOGENIC,
            af=0.05,
        )

        cache_key = "test_key_1"
        cache.put(cache_key, variant)

        retrieved = cache.get(cache_key)

        assert retrieved is not None
        assert retrieved.id == "rs1234"
        assert retrieved.chr == "chr17"
        assert retrieved.pos == 7577121
        assert retrieved.ref == "G"
        assert retrieved.alt == "A"
        assert retrieved.af == 0.05

    def test_get_nonexistent_key(self, tmp_path: Path):
        """Test that getting a nonexistent key returns None."""
        cache = VariantParquetCache(tmp_path)

        result = cache.get("nonexistent_key")

        assert result is None

    def test_update_existing_entry(self, tmp_path: Path):
        """Test that putting the same key updates the entry."""
        cache = VariantParquetCache(tmp_path)

        variant1 = VariantRecord(
            id="rs1234",
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
            af=0.05,
        )

        variant2 = VariantRecord(
            id="rs1234",
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
            af=0.10,  # Different AF
        )

        cache_key = "test_key"
        cache.put(cache_key, variant1)
        cache.put(cache_key, variant2)

        retrieved = cache.get(cache_key)

        assert retrieved is not None
        assert retrieved.af == 0.10  # Should have updated value

        # Check that we don't have duplicates
        stats = cache.get_stats()
        assert stats["total_entries"] == 1

    def test_multiple_variants(self, tmp_path: Path):
        """Test storing and retrieving multiple variants."""
        cache = VariantParquetCache(tmp_path)

        variants = [VariantRecord(id=f"rs{i}", chr="chr1", pos=100 + i, ref="A", alt="T") for i in range(5)]

        for i, variant in enumerate(variants):
            cache.put(f"key_{i}", variant)

        # Retrieve all
        for i, _variant in enumerate(variants):
            retrieved = cache.get(f"key_{i}")
            assert retrieved is not None
            assert retrieved.id == f"rs{i}"

    def test_get_stats(self, tmp_path: Path):
        """Test getting cache statistics."""
        cache = VariantParquetCache(tmp_path)

        # Empty cache
        stats = cache.get_stats()
        assert stats["total_entries"] == 0
        assert stats["stale_entries"] == 0

        # Add some entries
        for i in range(3):
            variant = VariantRecord(
                id=f"rs{i}",
                chr="chr1",
                pos=100 + i,
                ref="A",
                alt="T",
            )
            cache.put(f"key_{i}", variant)

        stats = cache.get_stats()
        assert stats["total_entries"] == 3
        assert stats["stale_entries"] == 0
        assert "cache_size_mb" in stats

    def test_clear_cache(self, tmp_path: Path):
        """Test clearing the cache."""
        cache = VariantParquetCache(tmp_path)

        # Add entries
        for i in range(3):
            variant = VariantRecord(
                id=f"rs{i}",
                chr="chr1",
                pos=100 + i,
                ref="A",
                alt="T",
            )
            cache.put(f"key_{i}", variant)

        assert cache.get_stats()["total_entries"] == 3

        # Clear
        cache.clear()

        assert cache.get_stats()["total_entries"] == 0

        # Entries should not be retrievable
        assert cache.get("key_0") is None

    def test_cache_with_complex_annotations(self, tmp_path: Path):
        """Test caching variants with complex annotation dictionaries."""
        cache = VariantParquetCache(tmp_path)

        variant = VariantRecord(
            id="rs1234",
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            annotations={
                "gene": "TP53",
                "consequence": "missense_variant",
                "impact": "MODERATE",
                "nested": {"key1": "value1", "key2": [1, 2, 3]},
            },
            provenance={
                "source": "test",
                "timestamp": "2025-12-19T00:00:00",
            },
        )

        cache.put("test_key", variant)
        retrieved = cache.get("test_key")

        assert retrieved is not None
        assert retrieved.annotations["gene"] == "TP53"
        assert retrieved.provenance["source"] == "test"

    def test_ttl_not_enforced_on_get(self, tmp_path: Path):
        """Test that TTL is checked but old entries are not deleted on get."""
        # Use short TTL for testing
        cache = VariantParquetCache(tmp_path, ttl_days=0)

        variant = VariantRecord(
            id="rs1234",
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
        )

        cache.put("test_key", variant)

        # Entry is immediately stale due to TTL=0
        # Should return None but not delete
        result = cache.get("test_key")
        assert result is None

        # Entry should still be in cache (cleanup must be called explicitly)
        stats = cache.get_stats()
        assert stats["total_entries"] == 1
        assert stats["stale_entries"] == 1


@pytest.mark.unit
class TestVariantParquetCacheIntegrity:
    """Regression tests for cache writes corrupting or dropping unrelated rows."""

    @staticmethod
    def _variant(pos: int, **kwargs) -> VariantRecord:
        return VariantRecord(chr="chr1", pos=pos, ref="A", alt="T", **kwargs)

    def test_overwrite_keeps_other_entries(self, tmp_path: Path):
        """Re-putting a key must not evict entries stored after it."""
        cache = VariantParquetCache(tmp_path)

        for i in range(3):
            cache.put(f"key_{i}", self._variant(100 + i, id=f"rs{i}"))

        # Overwrites the middle entry; with label-preserving indexing the append
        # landed on label 2 and silently destroyed key_2.
        cache.put("key_1", self._variant(201, id="rs1_updated"))

        assert cache.get("key_0") is not None
        updated = cache.get("key_1")
        assert updated is not None
        assert updated.id == "rs1_updated"
        assert cache.get("key_2") is not None, "unrelated entry evicted by overwrite"
        assert cache.get_stats()["total_entries"] == 3

    def test_put_still_works_after_cleanup(self, tmp_path: Path):
        """cleanup_stale_entries must not leave the cache unwritable."""
        cache = VariantParquetCache(tmp_path)
        cache.put("first", self._variant(100, id="rs1"))

        assert cache.cleanup_stale_entries() == 0

        cache.put("second", self._variant(101, id="rs2"))

        assert cache.get("second") is not None, "cache went write-dead after cleanup"
        assert cache.get("first") is not None
        assert cache.get_stats()["total_entries"] == 2

    def test_put_recovers_cache_with_datetime_cached_at(self, tmp_path: Path):
        """A cache file left with datetime64 cached_at by an older build stays writable."""
        cache = VariantParquetCache(tmp_path)
        cache.put("first", self._variant(100, id="rs1"))

        corrupted = pd.read_parquet(cache.cache_file, engine="pyarrow")
        corrupted["cached_at"] = pd.to_datetime(corrupted["cached_at"])
        corrupted.to_parquet(cache.cache_file, index=False, engine="pyarrow", compression="snappy")

        cache.put("second", self._variant(101, id="rs2"))

        assert cache.get("second") is not None
        assert cache.get("first") is not None

    def test_population_afs_survive_round_trip(self, tmp_path: Path):
        """Population-specific AFs must not be dropped by the cache round trip."""
        cache = VariantParquetCache(tmp_path)
        variant = self._variant(100, id="rs1", af=0.005, population_afs={"AFR": 0.15, "EUR": 0.02})

        cache.put("key", variant)
        retrieved = cache.get("key")

        assert retrieved is not None
        assert retrieved.population_afs == {"AFR": 0.15, "EUR": 0.02}
        # Drives avoid-mode filtering, so a dropped dict changes which variants pass.
        assert retrieved.get_effective_af_for_mode(VariantMode.AVOID) == 0.15

    def test_reads_legacy_cache_without_population_afs(self, tmp_path: Path):
        """A cache file written before the population_afs column stays readable."""
        cache = VariantParquetCache(tmp_path)
        cache.put("key", self._variant(100, id="rs1", af=0.005))

        legacy = pd.read_parquet(cache.cache_file, engine="pyarrow").drop(columns=["population_afs"])
        legacy.to_parquet(cache.cache_file, index=False, engine="pyarrow", compression="snappy")

        retrieved = cache.get("key")
        assert retrieved is not None
        assert retrieved.population_afs == {}

        # And the missing column must not be dropped from newly written rows.
        cache.put("key2", self._variant(101, id="rs2", population_afs={"AFR": 0.2}))
        retrieved2 = cache.get("key2")
        assert retrieved2 is not None
        assert retrieved2.population_afs == {"AFR": 0.2}

    def test_write_failure_warns_once(self, tmp_path: Path, caplog):
        """A cache that cannot be written warns on the first failure, not on every one."""
        cache = VariantParquetCache(tmp_path)
        cache.cache_file.unlink()
        cache.cache_file.mkdir()  # makes read_parquet/to_parquet fail

        with caplog.at_level(logging.DEBUG, logger="sirnaforge.data.variant_cache"):
            cache.put("key_1", self._variant(100))
            cache.put("key_2", self._variant(101))

        warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "will not be updated" in warnings[0].message
