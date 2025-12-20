"""Unit tests for VariantResolver."""

import pytest

from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import (
    ClinVarSignificance,
    VariantQuery,
    VariantQueryType,
    VariantRecord,
    VariantSource,
)


class TestVariantResolverParsing:
    """Tests for variant identifier parsing."""

    def test_parse_rsid(self):
        """Test parsing rsID format."""
        resolver = VariantResolver()

        # Standard rsID
        query = resolver.parse_identifier("rs1234")
        assert query.query_type == VariantQueryType.RSID
        assert query.rsid == "rs1234"
        assert query.raw_input == "rs1234"

        # Case insensitive
        query2 = resolver.parse_identifier("RS9876")
        assert query2.query_type == VariantQueryType.RSID
        assert query2.rsid == "RS9876"

    def test_parse_coordinate_with_chr_prefix(self):
        """Test parsing coordinate format with chr prefix."""
        resolver = VariantResolver()

        query = resolver.parse_identifier("chr17:7577121:G:A")
        assert query.query_type == VariantQueryType.COORDINATE
        assert query.chr == "chr17"
        assert query.pos == 7577121
        assert query.ref == "G"
        assert query.alt == "A"

    def test_parse_coordinate_without_chr_prefix(self):
        """Test parsing coordinate format without chr prefix."""
        resolver = VariantResolver()

        query = resolver.parse_identifier("17:7577121:G:A")
        assert query.query_type == VariantQueryType.COORDINATE
        assert query.chr == "chr17"  # Should be normalized to chr17
        assert query.pos == 7577121
        assert query.ref == "G"
        assert query.alt == "A"

    def test_parse_coordinate_sex_chromosomes(self):
        """Test parsing coordinates for X, Y, MT chromosomes."""
        resolver = VariantResolver()

        # X chromosome
        query_x = resolver.parse_identifier("chrX:12345:A:T")
        assert query_x.chr == "chrX"

        # Y chromosome
        query_y = resolver.parse_identifier("Y:67890:C:G")
        assert query_y.chr == "chrY"

        # Mitochondrial
        query_mt = resolver.parse_identifier("MT:100:A:G")
        assert query_mt.chr == "chrMT"

    def test_parse_hgvs_refseq(self):
        """Test parsing HGVS RefSeq format."""
        resolver = VariantResolver()

        query = resolver.parse_identifier("NM_000546.6:c.215C>G")
        assert query.query_type == VariantQueryType.HGVS
        assert query.hgvs == "NM_000546.6:c.215C>G"

    def test_parse_hgvs_ensembl(self):
        """Test parsing HGVS Ensembl format."""
        resolver = VariantResolver()

        query = resolver.parse_identifier("ENST00000269305.9:c.215C>G")
        assert query.query_type == VariantQueryType.HGVS
        assert query.hgvs == "ENST00000269305.9:c.215C>G"

    def test_parse_invalid_format(self):
        """Test parsing invalid format raises ValueError."""
        resolver = VariantResolver()

        with pytest.raises(ValueError, match="Unrecognized variant identifier format"):
            resolver.parse_identifier("invalid_format")

        with pytest.raises(ValueError):
            resolver.parse_identifier("chr17:not_a_number:G:A")

        with pytest.raises(ValueError):
            resolver.parse_identifier("rs")  # rsID without number

    def test_parse_coordinate_multichar_alleles(self):
        """Test parsing coordinates with multi-character alleles."""
        resolver = VariantResolver()

        # Deletion
        query = resolver.parse_identifier("chr1:12345:ATG:A")
        assert query.ref == "ATG"
        assert query.alt == "A"

        # Insertion
        query2 = resolver.parse_identifier("chr1:12345:A:ATGC")
        assert query2.ref == "A"
        assert query2.alt == "ATGC"


class TestVariantResolverFiltering:
    """Tests for variant filtering logic."""

    def test_af_filter_pass(self):
        """Test variants passing AF filter."""
        resolver = VariantResolver(min_af=0.01)
        # AF above threshold
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.05,
        )
        assert resolver._passes_filters(variant) is True

        # AF exactly at threshold
        variant2 = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.01,
        )
        assert resolver._passes_filters(variant2) is True

        # No AF info - should pass
        variant3 = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
        )
        assert resolver._passes_filters(variant3) is True

    def test_af_filter_fail(self):
        """Test variants failing AF filter."""
        resolver = VariantResolver(min_af=0.01)
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.005,  # Below threshold
        )
        assert resolver._passes_filters(variant) is False

    def test_clinvar_filter_pass(self):
        """Test variants passing ClinVar filter."""
        resolver = VariantResolver(
            clinvar_filters=[ClinVarSignificance.PATHOGENIC, ClinVarSignificance.LIKELY_PATHOGENIC]
        )

        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.CLINVAR],
            clinvar_significance=ClinVarSignificance.PATHOGENIC,
        )
        assert resolver._passes_filters(variant) is True

    def test_clinvar_filter_fail(self):
        """Test variants failing ClinVar filter."""
        resolver = VariantResolver(
            clinvar_filters=[ClinVarSignificance.PATHOGENIC, ClinVarSignificance.LIKELY_PATHOGENIC]
        )
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.CLINVAR],
            clinvar_significance=ClinVarSignificance.BENIGN,  # Not in allowed list
        )
        assert resolver._passes_filters(variant) is False

    def test_clinvar_filter_non_clinvar_source(self):
        """Test that ClinVar filter only applies to ClinVar sources."""
        resolver = VariantResolver(
            clinvar_filters=[ClinVarSignificance.PATHOGENIC, ClinVarSignificance.LIKELY_PATHOGENIC]
        )
        # Variant from Ensembl, not ClinVar - should pass regardless of significance
        variant = VariantRecord(
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            sources=[VariantSource.ENSEMBL],
        )
        assert resolver._passes_filters(variant) is True


class TestVariantResolverCache:
    """Tests for variant caching."""

    def test_cache_key_generation(self):
        """Test cache key generation."""
        resolver = VariantResolver(min_af=0.01)

        query = VariantQuery(
            raw_input="rs1234",
            query_type=VariantQueryType.RSID,
            rsid="rs1234",
        )

        cache_key = resolver._get_cache_key(query)
        assert isinstance(cache_key, str)
        assert len(cache_key) > 0

        # Same query should produce same cache key
        cache_key2 = resolver._get_cache_key(query)
        assert cache_key == cache_key2

    def test_cache_put_and_get(self, tmp_path):
        """Test caching and retrieving variants."""
        resolver = VariantResolver(cache_dir=tmp_path)

        variant = VariantRecord(
            id="rs1234",
            chr="chr1",
            pos=100,
            ref="A",
            alt="T",
            af=0.05,
        )

        # Cache the variant
        cache_key = "test_key"
        resolver._put_to_cache(cache_key, variant)

        # Retrieve from cache
        cached_variant = resolver._get_from_cache(cache_key)
        assert cached_variant is not None
        assert cached_variant.id == "rs1234"
        assert cached_variant.chr == "chr1"
        assert cached_variant.pos == 100
        assert cached_variant.af == 0.05


class TestVariantResolverInitialization:
    """Tests for VariantResolver initialization."""

    def test_init_defaults(self):
        """Test default initialization."""
        resolver = VariantResolver()

        assert resolver.min_af == 0.01
        assert ClinVarSignificance.PATHOGENIC in resolver.clinvar_filters
        assert ClinVarSignificance.LIKELY_PATHOGENIC in resolver.clinvar_filters
        assert resolver.assembly == "GRCh38"
        assert resolver.source_priority == [VariantSource.CLINVAR, VariantSource.ENSEMBL, VariantSource.DBSNP]

    def test_init_custom_values(self):
        """Test custom initialization."""
        resolver = VariantResolver(
            min_af=0.05,
            clinvar_filters=[ClinVarSignificance.PATHOGENIC],
            timeout=60,
        )

        assert resolver.min_af == 0.05
        assert resolver.clinvar_filters == [ClinVarSignificance.PATHOGENIC]
        assert resolver.timeout == 60

    def test_init_invalid_assembly(self):
        """Test that non-GRCh38 assembly raises error."""
        with pytest.raises(ValueError, match="Only GRCh38 assembly is supported"):
            VariantResolver(assembly="GRCh37")
