"""Unit tests for variant models and parsing."""

import pytest

from sirnaforge.models.variant import (
    ClinVarSignificance,
    ClinVarVariationResponse,
    EnsemblMapping,
    EnsemblPopulationFrequency,
    EnsemblVariationResponse,
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


class TestEnsemblMapping:
    """Tests for EnsemblMapping model."""

    def test_ensembl_mapping_creation(self):
        """Test creating a valid EnsemblMapping."""
        mapping = EnsemblMapping(
            location="17:7673802-7673802",
            allele_string="C/A",
            assembly_name="GRCh38",
            seq_region_name="17",
            strand=1,
            start=7673802,
            end=7673802,
            coord_system="chromosome",
        )

        assert mapping.location == "17:7673802-7673802"
        assert mapping.allele_string == "C/A"
        assert mapping.assembly_name == "GRCh38"
        assert mapping.seq_region_name == "17"
        assert mapping.strand == 1
        assert mapping.start == 7673802
        assert mapping.end == 7673802

    def test_ensembl_mapping_extra_fields(self):
        """Test EnsemblMapping allows extra fields."""
        mapping = EnsemblMapping(
            location="17:7673802-7673802",
            allele_string="C/A",
            assembly_name="GRCh38",
            seq_region_name="17",
            strand=1,
            start=7673802,
            end=7673802,
            coord_system="chromosome",
            extra_field="extra_value",  # Should be allowed
        )

        assert mapping.extra_field == "extra_value"


class TestEnsemblPopulationFrequency:
    """Tests for EnsemblPopulationFrequency model."""

    def test_population_frequency_creation(self):
        """Test creating a valid EnsemblPopulationFrequency."""
        freq = EnsemblPopulationFrequency(
            population="AFR",
            frequency=0.15,
            allele_count=30,
            allele_number=200,
        )

        assert freq.population == "AFR"
        assert freq.frequency == 0.15
        assert freq.allele_count == 30
        assert freq.allele_number == 200

    def test_population_frequency_optional_fields(self):
        """Test EnsemblPopulationFrequency with optional fields."""
        freq = EnsemblPopulationFrequency(
            population="EUR",
            frequency=None,
            allele_count=None,
            allele_number=None,
        )

        assert freq.population == "EUR"
        assert freq.frequency is None
        assert freq.allele_count is None
        assert freq.allele_number is None

    def test_population_frequency_validation(self):
        """Test EnsemblPopulationFrequency field validation."""
        with pytest.raises(ValueError):
            EnsemblPopulationFrequency(
                population="AFR",
                frequency=1.5,  # Invalid: > 1.0
            )

        with pytest.raises(ValueError):
            EnsemblPopulationFrequency(
                population="AFR",
                frequency=-0.1,  # Invalid: < 0.0
            )


class TestEnsemblVariationResponse:
    """Tests for EnsemblVariationResponse model."""

    def test_ensembl_response_creation(self):
        """Test creating a valid EnsemblVariationResponse."""
        response = EnsemblVariationResponse(
            name="rs28934576",
            var_class="SNP",
            source="Variants (including SNPs and indels) imported from dbSNP",
            most_severe_consequence="missense_variant",
            MAF=0.001,
            minor_allele="A",
            ambiguity="N",
            mappings=[
                EnsemblMapping(
                    location="17:7673802-7673802",
                    allele_string="C/A",
                    assembly_name="GRCh38",
                    seq_region_name="17",
                    strand=1,
                    start=7673802,
                    end=7673802,
                    coord_system="chromosome",
                )
            ],
            clinical_significance=["pathogenic"],
            synonyms=["RCV000376655"],
            evidence=["Frequency", "1000Genomes"],
            populations=[],
        )

        assert response.name == "rs28934576"
        assert response.var_class == "SNP"
        assert response.MAF == 0.001
        assert len(response.mappings) == 1
        assert response.clinical_significance == ["pathogenic"]

    def test_ensembl_response_extra_fields(self):
        """Test EnsemblVariationResponse allows extra fields."""
        response = EnsemblVariationResponse(
            name="rs28934576",
            var_class="SNP",
            source="dbSNP",
            most_severe_consequence=None,
            MAF=None,
            minor_allele=None,
            ambiguity="N",
            mappings=[],
            clinical_significance=None,
            synonyms=None,
            evidence=None,
            populations=[],
            extra_api_field="extra_value",  # Should be allowed
        )

        assert response.extra_api_field == "extra_value"


class TestClinVarVariationResponse:
    """Tests for ClinVarVariationResponse model."""

    def test_clinvar_response_creation(self):
        """Test creating a valid ClinVarVariationResponse."""
        response = ClinVarVariationResponse(
            uid="376655",
            obj_type="single nucleotide variant",
            accession="VCV000376655",
            title="NM_000546.6(TP53):c.818G>T (p.Arg273Leu)",
            germline_classification={
                "description": "Pathogenic",
                "last_evaluated": "2024/09/29 00:00",
                "review_status": "criteria provided, multiple submitters, no conflicts",
            },
            clinical_impact_classification=None,
            variation_set=[],
            genes=[],
            molecular_consequence_list=["missense variant"],
        )

        assert response.uid == "376655"
        assert response.accession == "VCV000376655"
        assert response.title == "NM_000546.6(TP53):c.818G>T (p.Arg273Leu)"
        assert response.germline_classification["description"] == "Pathogenic"

    def test_clinvar_response_minimal(self):
        """Test ClinVarVariationResponse with minimal required fields."""
        response = ClinVarVariationResponse(
            uid="376655",
            obj_type=None,
            accession=None,
            title=None,
            germline_classification=None,
            clinical_impact_classification=None,
            variation_set=None,
            genes=None,
            molecular_consequence_list=None,
        )

        assert response.uid == "376655"
        assert response.obj_type is None
        assert response.germline_classification is None

    def test_clinvar_response_extra_fields(self):
        """Test ClinVarVariationResponse allows extra fields."""
        response = ClinVarVariationResponse(
            uid="376655",
            obj_type="single nucleotide variant",
            accession="VCV000376655",
            title="Test variant",
            germline_classification=None,
            clinical_impact_classification=None,
            variation_set=None,
            genes=None,
            molecular_consequence_list=None,
            extra_api_field="extra_value",  # Should be allowed
        )

        assert response.extra_api_field == "extra_value"
