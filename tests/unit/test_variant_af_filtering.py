"""Test variant AF filtering behavior.

Tests that variants without allele frequency data are included in 'all' results
but excluded from 'passed' results when AF filtering is applied.
"""

from pathlib import Path

import pytest

from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import (
    ClinVarSignificance,
    VariantRecord,
    VariantSource,
)


class TestVariantAFFiltering:
    """Test allele frequency filtering behavior."""

    def test_variant_with_af_passes_filter(self):
        """Test that variants with sufficient AF pass the filter."""
        resolver = VariantResolver(min_af=0.01)

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.05,  # 5% AF - above threshold
        )

        assert resolver._passes_filters(variant) is True

    def test_variant_with_low_af_fails_filter(self):
        """Test that variants with insufficient AF fail the filter."""
        resolver = VariantResolver(min_af=0.01)

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.0005,  # 0.05% AF - below threshold
        )

        assert resolver._passes_filters(variant) is False

    def test_variant_without_af_data_passes_filter(self):
        """Test that variants without AF data pass the filter (with warning)."""
        resolver = VariantResolver(min_af=0.01)

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=None,  # No AF data
        )

        # Should pass filter since no AF data is available
        assert resolver._passes_filters(variant) is True

    def test_variant_with_population_af_avoid_mode(self):
        """Test that avoid mode uses max population AF for filtering."""
        resolver = VariantResolver(min_af=0.01, variant_mode="avoid")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.005,  # 0.5% global AF - below threshold
            population_afs={"AFR": 0.15, "EUR": 0.02, "EAS": 0.08},  # Max is 15%
        )

        # In avoid mode, should use max population AF (15%) which is above threshold
        assert resolver._passes_filters(variant) is True

    def test_variant_with_population_af_target_mode(self):
        """Test that target mode uses global AF for filtering."""
        resolver = VariantResolver(min_af=0.01, variant_mode="target")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.005,  # 0.5% global AF - below threshold
            population_afs={"AFR": 0.15, "EUR": 0.02, "EAS": 0.08},
        )

        # In target mode, should use global AF (0.5%) which is below threshold
        assert resolver._passes_filters(variant) is False

    def test_variant_with_population_af_both_mode(self):
        """Test that both mode uses global AF for filtering."""
        resolver = VariantResolver(min_af=0.01, variant_mode="both")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.005,  # 0.5% global AF - below threshold
            population_afs={"AFR": 0.15, "EUR": 0.02, "EAS": 0.08},
        )

        # In both mode, should use global AF (0.5%) which is below threshold
        assert resolver._passes_filters(variant) is False

    def test_variant_without_af_data_avoid_mode(self):
        """Test that variants without AF data pass avoid mode filter."""
        resolver = VariantResolver(min_af=0.01, variant_mode="avoid")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=None,  # No AF data
            population_afs={},  # No population data
        )

        # Should pass filter since no AF data is available
        assert resolver._passes_filters(variant) is True

    def test_avoid_mode_uses_population_af_when_global_missing(self):
        """Avoid mode should still enforce AF threshold with population-only AFs."""
        resolver = VariantResolver(min_af=0.05, variant_mode="avoid")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=None,
            population_afs={"AFR": 0.08, "EUR": 0.01},
        )

        assert resolver._passes_filters(variant) is True

    def test_avoid_mode_population_af_below_threshold_fails(self):
        """Avoid mode should reject variants if max population AF is below threshold."""
        resolver = VariantResolver(min_af=0.05, variant_mode="avoid")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=None,
            population_afs={"AFR": 0.03, "EUR": 0.01},
        )

        assert resolver._passes_filters(variant) is False

    def test_target_mode_with_only_population_af_skips_af_filter(self):
        """Target mode uses global AF only and therefore skips AF filtering when global AF is missing."""
        resolver = VariantResolver(min_af=0.05, variant_mode="target")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=None,
            population_afs={"AFR": 0.20, "EUR": 0.18},
        )

        assert resolver._passes_filters(variant) is True

    def test_both_mode_with_only_population_af_skips_af_filter(self):
        """Both mode follows target-mode AF semantics when global AF is missing."""
        resolver = VariantResolver(min_af=0.05, variant_mode="both")

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=None,
            population_afs={"AFR": 0.20, "EUR": 0.18},
        )

        assert resolver._passes_filters(variant) is True

    def test_clinvar_significance_filter(self):
        """Test that ClinVar significance filtering works correctly."""
        resolver = VariantResolver(min_af=0.01, clinvar_filters=[ClinVarSignificance.PATHOGENIC])

        # Pathogenic variant should pass
        pathogenic_variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.CLINVAR],
            af=0.05,
            clinvar_significance=ClinVarSignificance.PATHOGENIC,
        )

        assert resolver._passes_filters(pathogenic_variant) is True

        # Benign variant should fail
        benign_variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.CLINVAR],
            af=0.05,
            clinvar_significance=ClinVarSignificance.BENIGN,
        )

        assert resolver._passes_filters(benign_variant) is False

        # Non-ClinVar variant should pass (no ClinVar significance filter applied)
        non_clinvar_variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.05,
            clinvar_significance=None,
        )

        assert resolver._passes_filters(non_clinvar_variant) is True


@pytest.mark.unit
class TestCachedVariantAFFiltering:
    """AF filtering must give the same answer on a warm cache as on a cold one."""

    def test_avoid_mode_verdict_survives_cache_round_trip(self, tmp_path: Path):
        """A cached variant must keep the population AFs that avoid mode filters on."""
        resolver = VariantResolver(min_af=0.01, variant_mode="avoid", cache_dir=tmp_path)

        variant = VariantRecord(
            id="rs28934576",
            chr="chr17",
            pos=7673802,
            ref="C",
            alt="A",
            assembly="GRCh38",
            sources=[VariantSource.ENSEMBL],
            af=0.005,  # 0.5% global AF - below threshold
            population_afs={"AFR": 0.15, "EUR": 0.02, "EAS": 0.08},  # Max is 15%
        )

        # Cold cache: passes on max population AF.
        assert resolver._passes_filters(variant) is True

        resolver._put_to_cache("cache_key", variant)
        cached = resolver._get_from_cache("cache_key")

        assert cached is not None
        assert cached.population_afs == variant.population_afs
        assert resolver._passes_filters(cached) is True
