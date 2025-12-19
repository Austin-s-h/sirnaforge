"""Tests for variant workflow integration."""

import pytest

from sirnaforge.models.variant import ClinVarSignificance, VariantMode, VariantRecord, VariantSource
from sirnaforge.workflow_variant import (
    VariantWorkflowConfig,
    parse_clinvar_filter_string,
)


class TestVariantWorkflowConfig:
    """Tests for VariantWorkflowConfig."""

    def test_default_config(self):
        """Test default configuration."""
        config = VariantWorkflowConfig()

        assert config.variant_ids == []
        assert config.vcf_file is None
        assert config.variant_mode == VariantMode.AVOID
        assert config.min_af == 0.01
        assert ClinVarSignificance.PATHOGENIC in config.clinvar_filter_levels
        assert ClinVarSignificance.LIKELY_PATHOGENIC in config.clinvar_filter_levels
        assert config.assembly == "GRCh38"
        assert config.has_variants is False

    def test_with_variant_ids(self):
        """Test configuration with variant IDs."""
        config = VariantWorkflowConfig(variant_ids=["rs1234", "rs5678"])

        assert len(config.variant_ids) == 2
        assert config.has_variants is True

    def test_with_vcf_file(self, tmp_path):
        """Test configuration with VCF file."""
        vcf_path = tmp_path / "variants.vcf"
        vcf_path.touch()

        config = VariantWorkflowConfig(vcf_file=vcf_path)

        assert config.vcf_file == vcf_path
        assert config.has_variants is True

    def test_custom_mode(self):
        """Test custom variant mode."""
        config = VariantWorkflowConfig(variant_mode=VariantMode.TARGET)

        assert config.variant_mode == VariantMode.TARGET


class TestParseClinVarFilterString:
    """Tests for parse_clinvar_filter_string."""

    def test_parse_single_level(self):
        """Test parsing a single significance level."""
        result = parse_clinvar_filter_string("Pathogenic")

        assert len(result) == 1
        assert ClinVarSignificance.PATHOGENIC in result

    def test_parse_multiple_levels(self):
        """Test parsing multiple significance levels."""
        result = parse_clinvar_filter_string("Pathogenic,Likely pathogenic")

        assert len(result) == 2
        assert ClinVarSignificance.PATHOGENIC in result
        assert ClinVarSignificance.LIKELY_PATHOGENIC in result

    def test_parse_with_whitespace(self):
        """Test parsing with whitespace."""
        result = parse_clinvar_filter_string("Pathogenic , Likely pathogenic , Benign")

        assert len(result) == 3
        assert ClinVarSignificance.PATHOGENIC in result
        assert ClinVarSignificance.LIKELY_PATHOGENIC in result
        assert ClinVarSignificance.BENIGN in result

    def test_parse_case_insensitive(self):
        """Test case-insensitive parsing."""
        result = parse_clinvar_filter_string("pathogenic,LIKELY PATHOGENIC")

        assert len(result) == 2
        assert ClinVarSignificance.PATHOGENIC in result
        assert ClinVarSignificance.LIKELY_PATHOGENIC in result

    def test_parse_invalid_level(self):
        """Test that invalid levels raise ValueError."""
        with pytest.raises(ValueError, match="Invalid ClinVar significance level"):
            parse_clinvar_filter_string("Invalid,Pathogenic")


class TestVariantDeduplication:
    """Tests for variant deduplication logic."""

    def test_deduplicate_variants(self):
        """Test deduplication of variants."""
        from sirnaforge.workflow_variant import _deduplicate_variants

        variants = [
            VariantRecord(
                id="rs1234",
                chr="chr1",
                pos=100,
                ref="A",
                alt="T",
                sources=[VariantSource.ENSEMBL],
            ),
            VariantRecord(
                id="rs1234",
                chr="chr1",
                pos=100,
                ref="A",
                alt="T",
                sources=[VariantSource.CLINVAR],  # Higher priority
            ),
            VariantRecord(
                id="rs5678",
                chr="chr2",
                pos=200,
                ref="G",
                alt="C",
                sources=[VariantSource.DBSNP],
            ),
        ]

        result = _deduplicate_variants(variants)

        # Should keep only 2 unique variants
        assert len(result) == 2

        # First variant should be the ClinVar one (higher priority)
        chr1_variants = [v for v in result if v.chr == "chr1"]
        assert len(chr1_variants) == 1
        assert VariantSource.CLINVAR in chr1_variants[0].sources


class TestVariantCounting:
    """Tests for variant counting utilities."""

    def test_count_by_source(self):
        """Test counting variants by source."""
        from sirnaforge.workflow_variant import _count_by_source

        variants = [
            VariantRecord(chr="chr1", pos=100, ref="A", alt="T", sources=[VariantSource.CLINVAR]),
            VariantRecord(chr="chr1", pos=200, ref="G", alt="C", sources=[VariantSource.CLINVAR]),
            VariantRecord(chr="chr2", pos=300, ref="T", alt="G", sources=[VariantSource.ENSEMBL]),
        ]

        counts = _count_by_source(variants)

        assert counts["clinvar"] == 2
        assert counts["ensembl"] == 1

    def test_count_by_chromosome(self):
        """Test counting variants by chromosome."""
        from sirnaforge.workflow_variant import _count_by_chromosome

        variants = [
            VariantRecord(chr="chr1", pos=100, ref="A", alt="T"),
            VariantRecord(chr="chr1", pos=200, ref="G", alt="C"),
            VariantRecord(chr="chr2", pos=300, ref="T", alt="G"),
            VariantRecord(chr="chrX", pos=400, ref="A", alt="C"),
        ]

        counts = _count_by_chromosome(variants)

        assert counts["chr1"] == 2
        assert counts["chr2"] == 1
        assert counts["chrX"] == 1
