"""Unit tests for Parquet-based variant cache."""

import pytest

from sirnaforge.data.variant_cache import VariantParquetCache
from sirnaforge.models.variant import ClinVarSignificance, VariantRecord, VariantSource


class TestVariantParquetCache:
    """Tests for VariantParquetCache."""

    def test_init_creates_empty_cache(self, tmp_path):
        """Test that initialization creates an empty cache file."""
        cache = VariantParquetCache(tmp_path)

        assert cache.cache_file.exists()
        assert cache.cache_dir == tmp_path

    def test_put_and_get(self, tmp_path):
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

    def test_get_nonexistent_key(self, tmp_path):
        """Test that getting a nonexistent key returns None."""
        cache = VariantParquetCache(tmp_path)

        result = cache.get("nonexistent_key")

        assert result is None

    def test_update_existing_entry(self, tmp_path):
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

    def test_multiple_variants(self, tmp_path):
        """Test storing and retrieving multiple variants."""
        cache = VariantParquetCache(tmp_path)

        variants = [VariantRecord(id=f"rs{i}", chr="chr1", pos=100 + i, ref="A", alt="T") for i in range(5)]

        for i, variant in enumerate(variants):
            cache.put(f"key_{i}", variant)

        # Retrieve all
        for i, variant in enumerate(variants):
            retrieved = cache.get(f"key_{i}")
            assert retrieved is not None
            assert retrieved.id == f"rs{i}"

    def test_get_stats(self, tmp_path):
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

    def test_clear_cache(self, tmp_path):
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

    def test_cache_with_complex_annotations(self, tmp_path):
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

    def test_ttl_not_enforced_on_get(self, tmp_path):
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
