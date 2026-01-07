"""Tests for population-specific allele frequency handling."""

from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import VariantMode, VariantRecord, VariantSource


class TestPopulationAF:
    """Tests for population-specific allele frequency methods."""

    def test_get_max_population_af(self):
        """Test getting maximum population AF."""
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.02,
            population_afs={"AFR": 0.15, "EUR": 0.01, "EAS": 0.03},
        )

        max_af = variant.get_max_population_af()
        assert max_af == 0.15

    def test_get_max_population_af_no_data(self):
        """Test max population AF when no population data."""
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.02,
        )

        max_af = variant.get_max_population_af()
        assert max_af is None

    def test_get_effective_af_avoid_mode_with_population_data(self):
        """Test effective AF in avoid mode with population data."""
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.02,
            population_afs={"AFR": 0.15, "EUR": 0.01},
        )

        effective_af = variant.get_effective_af_for_mode(VariantMode.AVOID)
        assert effective_af == 0.15

    def test_get_effective_af_avoid_mode_without_population_data(self):
        """Test effective AF in avoid mode without population data."""
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.02,
        )

        effective_af = variant.get_effective_af_for_mode(VariantMode.AVOID)
        assert effective_af == 0.02

    def test_get_effective_af_target_mode(self):
        """Test effective AF in target mode uses global AF."""
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.02,
            population_afs={"AFR": 0.15, "EUR": 0.01},
        )

        effective_af = variant.get_effective_af_for_mode(VariantMode.TARGET)
        assert effective_af == 0.02

    def test_get_effective_af_both_mode(self):
        """Test effective AF in both mode uses global AF."""
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.02,
            population_afs={"AFR": 0.15, "EUR": 0.01},
        )

        effective_af = variant.get_effective_af_for_mode(VariantMode.BOTH)
        assert effective_af == 0.02


class TestAvoidModeFiltering:
    """Tests for avoid mode filtering with population AFs."""

    def test_avoid_mode_filters_high_population_af(self):
        """Test that avoid mode filters variants with high population AF."""
        resolver = VariantResolver(
            min_af=0.05,
            variant_mode="avoid",
        )

        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.ENSEMBL],
            af=0.02,
            population_afs={"AFR": 0.12, "EUR": 0.01},
        )

        assert resolver._passes_filters(variant) is True

        resolver.min_af = 0.15
        assert resolver._passes_filters(variant) is False

    def test_avoid_mode_without_population_data_uses_global(self):
        """Test that avoid mode falls back to global AF without population data."""
        resolver = VariantResolver(
            min_af=0.05,
            variant_mode="avoid",
        )

        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.ENSEMBL],
            af=0.08,
        )

        assert resolver._passes_filters(variant) is True


class TestTargetModeFiltering:
    """Tests for target mode filtering."""

    def test_target_mode_uses_global_af(self):
        """Test that target mode uses global AF regardless of population data."""
        resolver = VariantResolver(
            min_af=0.05,
            variant_mode="target",
        )

        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.ENSEMBL],
            af=0.03,
            population_afs={"AFR": 0.15},
        )

        assert resolver._passes_filters(variant) is False

    def test_target_mode_targets_common_variants(self):
        """Test that target mode correctly filters for common variants."""
        resolver = VariantResolver(
            min_af=0.05,
            variant_mode="target",
        )

        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.ENSEMBL],
            af=0.10,
            population_afs={"AFR": 0.08, "EUR": 0.12},
        )

        assert resolver._passes_filters(variant) is True
