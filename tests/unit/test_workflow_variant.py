"""Tests for variant workflow integration."""

import json
from pathlib import Path

import pytest

from sirnaforge.models.variant import ClinVarSignificance, VariantMode, VariantRecord, VariantSource
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig
from sirnaforge.workflow_variant import (
    VariantWorkflowConfig,
    _count_by_chromosome,
    _count_by_source,
    _deduplicate_variants,
    _save_variant_report,
    normalize_variant_mode,
    parse_clinvar_filter_string,
    resolve_workflow_variants,
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


class TestNormalizeVariantMode:
    """Tests for normalize_variant_mode function."""

    def test_enum_passthrough(self):
        """Test that enum values pass through unchanged."""
        assert normalize_variant_mode(VariantMode.AVOID) == VariantMode.AVOID
        assert normalize_variant_mode(VariantMode.TARGET) == VariantMode.TARGET
        assert normalize_variant_mode(VariantMode.BOTH) == VariantMode.BOTH

    def test_string_lowercase(self):
        """Test string values are normalized to enum."""
        assert normalize_variant_mode("avoid") == VariantMode.AVOID
        assert normalize_variant_mode("target") == VariantMode.TARGET
        assert normalize_variant_mode("both") == VariantMode.BOTH

    def test_string_uppercase(self):
        """Test uppercase strings are normalized."""
        assert normalize_variant_mode("AVOID") == VariantMode.AVOID
        assert normalize_variant_mode("TARGET") == VariantMode.TARGET
        assert normalize_variant_mode("BOTH") == VariantMode.BOTH

    def test_string_mixedcase(self):
        """Test mixed case strings are normalized."""
        assert normalize_variant_mode("Avoid") == VariantMode.AVOID
        assert normalize_variant_mode("Target") == VariantMode.TARGET
        assert normalize_variant_mode("Both") == VariantMode.BOTH

    def test_invalid_string_raises(self):
        """Test invalid string values raise ValueError."""
        with pytest.raises(ValueError, match="Invalid variant mode"):
            normalize_variant_mode("invalid")

        with pytest.raises(ValueError, match="Invalid variant mode"):
            normalize_variant_mode("skip")

    def test_none_raises(self):
        """Test None value raises ValueError."""
        with pytest.raises(ValueError, match="Invalid variant mode"):
            normalize_variant_mode(None)  # type: ignore[arg-type]


class TestResolveWorkflowVariants:
    """Tests for the asynchronous variant resolution workflow helper."""

    @pytest.mark.asyncio
    async def test_resolve_variants_from_vcf_writes_report(self, tmp_path):
        """Resolve variants from a local VCF and persist the summary report."""
        vcf_path = Path(__file__).resolve().parents[2] / "examples" / "variant_demo.vcf"
        output_dir = tmp_path / "variant_workflow"
        config = VariantWorkflowConfig(
            vcf_file=vcf_path,
            variant_mode=VariantMode.AVOID,
            min_af=0.01,
            assembly="GRCh38",
        )

        variants = await resolve_workflow_variants(config=config, gene_name="TP53", output_dir=output_dir)

        # Expect the low-AF record in the demo VCF to be filtered out at min_af=0.01
        assert len(variants) == 2

        report_path = output_dir / "logs" / "resolved_variants.json"
        assert report_path.exists()

        with report_path.open() as f:
            report = json.load(f)

        assert report["gene"] == "TP53"
        assert report["variant_mode"] == "avoid"
        assert report["summary"]["total_variants"] == 2
        assert report["summary"]["chromosomes"] == {"chr1": 1, "chr2": 1}

    def test_save_variant_report_allows_empty(self, tmp_path):
        """Ensure variant reports are written even when no variants pass filters."""
        config = VariantWorkflowConfig(variant_ids=["rsNope"], variant_mode=VariantMode.AVOID)
        output_path = tmp_path / "resolved_variants.json"

        _save_variant_report([], output_path, "TP53", config)

        with output_path.open() as f:
            report = json.load(f)

        assert report["gene"] == "TP53"
        assert report["summary"]["total_variants"] == 0
        assert report["variants"] == []

    def test_write_candidate_variant_links(self, tmp_path):
        """Verify candidate-variant links are persisted to JSON."""

        class DummyCandidate:
            def __init__(self, cid: str):
                self.id = cid
                self.transcript_id = "ENST0001"
                self.variant_mode = "target"
                self.allele_specific = True
                self.targeted_alleles = ["alt"]
                self.overlapped_variants = [{"id": "rs123", "chr": "chr1", "pos": 100}]

        candidates = [DummyCandidate("c1"), DummyCandidate("c2")]
        config = WorkflowConfig(output_dir=tmp_path, gene_query="TP53")
        wf = SiRNAWorkflow(config)
        output_path = tmp_path / "logs" / "candidate_variants.json"

        wf._write_candidate_variant_links(candidates, output_path)

        with output_path.open() as fh:
            payload = json.load(fh)

        assert payload["gene"] == "TP53"
        assert payload["variant_annotated_candidates"] == 2
        assert payload["candidates"][0]["id"] == "c1"
        assert payload["candidates"][0]["overlapped_variants"][0]["id"] == "rs123"
