"""Unit tests for variant models and parsing."""

import pytest

from sirnaforge.models.variant import (
    ClinVarSignificance,
    VariantMode,
    VariantQuery,
    VariantQueryType,
    VariantRecord,
    VariantSource,
)


class TestVariantMode:
    """Tests for VariantMode enum."""

    def test_variant_mode_values(self):
        """Test VariantMode enum values."""
        assert VariantMode.TARGET.value == "target"
        assert VariantMode.AVOID.value == "avoid"
        assert VariantMode.BOTH.value == "both"


class TestVariantSource:
    """Tests for VariantSource enum."""

    def test_variant_source_priority(self):
        """Test VariantSource priority ordering."""
        assert VariantSource.CLINVAR.value == "clinvar"
        assert VariantSource.ENSEMBL.value == "ensembl"
        assert VariantSource.DBSNP.value == "dbsnp"
        assert VariantSource.LOCAL_VCF.value == "local-vcf"


class TestClinVarSignificance:
    """Tests for ClinVarSignificance enum."""

    def test_clinvar_significance_values(self):
        """Test ClinVarSignificance enum values."""
        assert ClinVarSignificance.PATHOGENIC.value == "Pathogenic"
        assert ClinVarSignificance.LIKELY_PATHOGENIC.value == "Likely pathogenic"
        assert ClinVarSignificance.BENIGN.value == "Benign"


class TestVariantRecord:
    """Tests for VariantRecord model."""

    def test_variant_record_creation(self):
        """Test creating a valid VariantRecord."""
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

        assert variant.id == "rs1234"
        assert variant.chr == "chr17"
        assert variant.pos == 7577121
        assert variant.ref == "G"
        assert variant.alt == "A"
        assert variant.assembly == "GRCh38"
        assert VariantSource.CLINVAR in variant.sources
        assert variant.clinvar_significance == ClinVarSignificance.PATHOGENIC
        assert variant.af == 0.05

    def test_variant_record_to_vcf_style(self):
        """Test converting VariantRecord to VCF-style string."""
        variant = VariantRecord(
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
        )

        assert variant.to_vcf_style() == "chr17:7577121:G:A"

    def test_variant_record_get_primary_source(self):
        """Test getting primary source from VariantRecord."""
        # Single source
        variant = VariantRecord(
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
            sources=[VariantSource.ENSEMBL],
        )
        assert variant.get_primary_source() == VariantSource.ENSEMBL

        # Multiple sources - should return highest priority
        variant_multi = VariantRecord(
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
            sources=[VariantSource.DBSNP, VariantSource.CLINVAR, VariantSource.ENSEMBL],
        )
        assert variant_multi.get_primary_source() == VariantSource.CLINVAR

    def test_variant_record_validation(self):
        """Test VariantRecord validation."""
        # Test minimum position constraint
        with pytest.raises(ValueError):
            VariantRecord(
                chr="chr17",
                pos=0,  # Invalid: must be >= 1
                ref="G",
                alt="A",
            )

        # Test empty ref/alt
        with pytest.raises(ValueError):
            VariantRecord(
                chr="chr17",
                pos=100,
                ref="",  # Invalid: must have at least 1 character
                alt="A",
            )

    def test_variant_record_af_validation(self):
        """Test allele frequency validation."""
        # Valid AF
        variant = VariantRecord(chr="chr1", pos=100, ref="A", alt="T", af=0.5)
        assert variant.af == 0.5

        # AF out of range
        with pytest.raises(ValueError):
            VariantRecord(chr="chr1", pos=100, ref="A", alt="T", af=1.5)

        with pytest.raises(ValueError):
            VariantRecord(chr="chr1", pos=100, ref="A", alt="T", af=-0.1)


class TestVariantQuery:
    """Tests for VariantQuery model."""

    def test_variant_query_rsid(self):
        """Test VariantQuery for rsID."""
        query = VariantQuery(
            raw_input="rs1234",
            query_type=VariantQueryType.RSID,
            rsid="rs1234",
        )

        assert query.raw_input == "rs1234"
        assert query.query_type == VariantQueryType.RSID
        assert query.rsid == "rs1234"
        assert query.chr is None
        assert query.pos is None

    def test_variant_query_coordinate(self):
        """Test VariantQuery for coordinate."""
        query = VariantQuery(
            raw_input="chr17:7577121:G:A",
            query_type=VariantQueryType.COORDINATE,
            chr="chr17",
            pos=7577121,
            ref="G",
            alt="A",
        )

        assert query.raw_input == "chr17:7577121:G:A"
        assert query.query_type == VariantQueryType.COORDINATE
        assert query.chr == "chr17"
        assert query.pos == 7577121
        assert query.ref == "G"
        assert query.alt == "A"

    def test_variant_query_hgvs(self):
        """Test VariantQuery for HGVS."""
        query = VariantQuery(
            raw_input="NM_000546.6:c.215C>G",
            query_type=VariantQueryType.HGVS,
            hgvs="NM_000546.6:c.215C>G",
        )

        assert query.raw_input == "NM_000546.6:c.215C>G"
        assert query.query_type == VariantQueryType.HGVS
        assert query.hgvs == "NM_000546.6:c.215C>G"
        assert query.chr is None

    def test_variant_query_assembly_default(self):
        """Test VariantQuery assembly defaults to GRCh38."""
        query = VariantQuery(
            raw_input="rs1234",
            query_type=VariantQueryType.RSID,
            rsid="rs1234",
        )

        assert query.assembly == "GRCh38"
