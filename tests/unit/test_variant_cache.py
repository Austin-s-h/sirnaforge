"""Unit tests for Parquet-based variant cache."""

import logging
from datetime import datetime
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
        # Avoid-mode filtering reads these, so dropping them made the record the
        # resolver hands back on a warm run differ from the cold-run record.
        assert retrieved.get_effective_af_for_mode(VariantMode.AVOID) == 0.15

    @pytest.mark.parametrize("legacy_style", ["column_absent", "column_null"])
    def test_legacy_row_without_population_afs_is_a_miss(self, tmp_path: Path, legacy_style: str, caplog):
        """A row that predates the population_afs column must not be served.

        The current writer always stores at least ``"{}"``, so a null is a reliable
        "never written" marker rather than "no population AFs". Reading it as an
        empty dict would keep serving the original defect out of every cache that
        already exists, for the full 90-day TTL.
        """
        cache = VariantParquetCache(tmp_path)
        cache.put("key", self._variant(100, id="rs1", af=0.005))

        legacy = pd.read_parquet(cache.cache_file, engine="pyarrow")
        if legacy_style == "column_absent":
            legacy = legacy.drop(columns=["population_afs"])
        else:
            legacy["population_afs"] = None
        legacy.to_parquet(cache.cache_file, index=False, engine="pyarrow", compression="snappy")

        with caplog.at_level(logging.INFO, logger="sirnaforge.data.variant_cache"):
            assert cache.get("key") is None, "pre-fix row served as if it had no population AFs"

        # The user has to be told, otherwise the extra re-fetch looks like a bug.
        assert any("pre-0.6.0 cache entry" in r.message for r in caplog.records)

        # And the missing column must not be dropped from newly written rows.
        cache.put("key2", self._variant(101, id="rs2", population_afs={"AFR": 0.2}))
        retrieved = cache.get("key2")
        assert retrieved is not None
        assert retrieved.population_afs == {"AFR": 0.2}

    def test_columns_match_variant_record_fields(self):
        """The hand-maintained schema must cover exactly the model's fields.

        The root cause of the dropped population AFs was schema drift: ``_COLUMNS``
        and the field-by-field serialisation in ``put``/``get`` are maintained by
        hand against a model with ``extra="forbid"``, so a field added to
        ``VariantRecord`` is silently absent from the cache instead of failing.
        """
        bookkeeping = {"cache_key", "cached_at"}
        persisted = set(VariantParquetCache._COLUMNS) - bookkeeping

        assert persisted == set(VariantRecord.model_fields), (
            "VariantParquetCache._COLUMNS has drifted from VariantRecord: "
            f"missing {sorted(set(VariantRecord.model_fields) - persisted)}, "
            f"unexpected {sorted(persisted - set(VariantRecord.model_fields))}. "
            "Add the field to _COLUMNS and to both put() and get()."
        )
        assert bookkeeping <= set(VariantParquetCache._COLUMNS)

    def test_every_field_survives_round_trip(self, tmp_path: Path):
        """Covering _COLUMNS is not enough; put/get must actually carry each field."""
        cache = VariantParquetCache(tmp_path)
        variant = VariantRecord(
            id="rs1",
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            assembly="GRCh38",
            sources=[VariantSource.CLINVAR, VariantSource.ENSEMBL],
            clinvar_significance=ClinVarSignificance.PATHOGENIC,
            af=0.005,
            population_afs={"AFR": 0.15},
            annotations={"gene": "TP53"},
            provenance={"queried_at": "2026-01-01"},
        )

        cache.put("key", variant)

        assert cache.get("key") == variant

    def test_unknown_on_disk_columns_are_preserved(self, tmp_path: Path):
        """A column written by a newer sirnaforge must survive an older put.

        With a 90-day TTL one cache directory is easily shared between versions, so
        dropping everything outside ``_COLUMNS`` would destroy the newer version's
        data on the first write from the older one.
        """
        cache = VariantParquetCache(tmp_path)
        cache.put("key", self._variant(100, id="rs1"))

        with_future = pd.read_parquet(cache.cache_file, engine="pyarrow")
        with_future["future_col"] = ["from a newer version"]
        with_future.to_parquet(cache.cache_file, index=False, engine="pyarrow", compression="snappy")

        cache.put("key2", self._variant(101, id="rs2"))

        written = pd.read_parquet(cache.cache_file, engine="pyarrow")
        assert "future_col" in written.columns, "unknown column dropped by put"
        assert written.loc[written["cache_key"] == "key", "future_col"].iloc[0] == "from a newer version"
        # The new row has no value for a column this version knows nothing about.
        assert written.loc[written["cache_key"] == "key2", "future_col"].iloc[0] is None

    def test_ttl_readers_tolerate_mixed_iso_timestamps(self, tmp_path: Path):
        """Mixed ISO layouts in cached_at must not disable the TTL readers.

        ``datetime.now().isoformat()`` omits ``.%f`` when microsecond is exactly 0,
        so the column heterogeneity arises on its own; pandas then infers one format
        and raises, which left cleanup and get_stats permanently no-ops.
        """
        cache = VariantParquetCache(tmp_path)
        cache.put("second_precision", self._variant(100, id="rs1"))
        cache.put("microsecond_precision", self._variant(101, id="rs2"))

        mixed = pd.read_parquet(cache.cache_file, engine="pyarrow")
        mixed.loc[0, "cached_at"] = datetime.now().replace(microsecond=0).isoformat()
        mixed.to_parquet(cache.cache_file, index=False, engine="pyarrow", compression="snappy")

        stats = cache.get_stats()
        assert "error" not in stats
        assert stats["total_entries"] == 2
        assert stats["stale_entries"] == 0
        assert cache.cleanup_stale_entries() == 0
        assert cache.get("second_precision") is not None

    def test_unreadable_timestamp_is_stale_everywhere(self, tmp_path: Path):
        """cleanup, get_stats and get must agree about a row with no usable timestamp.

        ``pd.NaT`` subclasses ``datetime``, so normalising it produced the literal
        string "NaT": cleanup removed the row while get_stats called it fresh and
        get kept serving it.
        """
        cache = VariantParquetCache(tmp_path)
        cache.put("bad", self._variant(100, id="rs1"))
        cache.put("good", self._variant(101, id="rs2"))

        with_nat = pd.read_parquet(cache.cache_file, engine="pyarrow")
        with_nat["cached_at"] = pd.to_datetime(with_nat["cached_at"], format="ISO8601")
        with_nat.loc[with_nat["cache_key"] == "bad", "cached_at"] = pd.NaT
        with_nat.to_parquet(cache.cache_file, index=False, engine="pyarrow", compression="snappy")

        # A put rewrites the column; NaT must become a null, never the string "NaT".
        cache.put("third", self._variant(102, id="rs3"))
        assert "NaT" not in pd.read_parquet(cache.cache_file, engine="pyarrow")["cached_at"].tolist()

        assert cache.get("bad") is None, "row with no usable timestamp still served"
        assert cache.get_stats()["stale_entries"] == 1
        assert cache.cleanup_stale_entries() == 1
        assert cache.get("good") is not None
        assert cache.get("third") is not None

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
