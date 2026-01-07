"""Unit tests for transcript annotation models."""

import pytest
from pydantic import ValidationError

from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.models.transcript_annotation import Interval, TranscriptAnnotation, TranscriptAnnotationBundle


class TestInterval:
    """Tests for Interval model."""

    def test_interval_creation(self):
        """Test basic interval creation."""
        interval = Interval(
            seq_region_name="17",
            start=1000,
            end=2000,
            strand=1,
        )

        assert interval.seq_region_name == "17"
        assert interval.start == 1000
        assert interval.end == 2000
        assert interval.strand == 1

    def test_interval_optional_strand(self):
        """Test interval creation without strand."""
        interval = Interval(
            seq_region_name="X",
            start=500,
            end=1500,
        )

        assert interval.strand is None

    def test_interval_immutable(self):
        """Test that intervals are immutable."""
        interval = Interval(seq_region_name="1", start=100, end=200)

        with pytest.raises(ValidationError):
            interval.start = 300  # type: ignore

    def test_interval_hashable(self):
        """Test that intervals can be used in sets/dicts."""
        interval1 = Interval(seq_region_name="1", start=100, end=200, strand=1)
        interval2 = Interval(seq_region_name="1", start=100, end=200, strand=1)
        interval3 = Interval(seq_region_name="1", start=100, end=300, strand=1)

        # Same intervals should hash to same value
        assert hash(interval1) == hash(interval2)

        # Different intervals should have different hashes (usually)
        assert hash(interval1) != hash(interval3)

        # Can be used in sets
        interval_set = {interval1, interval2, interval3}
        assert len(interval_set) == 2  # interval1 and interval2 are duplicates


class TestTranscriptAnnotation:
    """Tests for TranscriptAnnotation model."""

    def test_minimal_annotation(self):
        """Test creation with minimal required fields."""
        annotation = TranscriptAnnotation(
            transcript_id="ENST00000269305",
            gene_id="ENSG00000141510",
            seq_region_name="17",
            start=7661779,
            end=7687550,
            strand=1,
            provider="ensembl_rest",
        )

        assert annotation.transcript_id == "ENST00000269305"
        assert annotation.gene_id == "ENSG00000141510"
        assert annotation.symbol is None
        assert annotation.biotype is None
        assert len(annotation.exons) == 0
        assert len(annotation.cds) == 0

    def test_complete_annotation(self):
        """Test creation with all fields."""
        exon1 = Interval(seq_region_name="17", start=7661779, end=7661910, strand=1)
        exon2 = Interval(seq_region_name="17", start=7668402, end=7669000, strand=1)
        cds1 = Interval(seq_region_name="17", start=7661779, end=7661910, strand=1)

        annotation = TranscriptAnnotation(
            transcript_id="ENST00000269305",
            gene_id="ENSG00000141510",
            symbol="TP53",
            biotype="protein_coding",
            seq_region_name="17",
            start=7661779,
            end=7687550,
            strand=1,
            exons=[exon1, exon2],
            cds=[cds1],
            provider="ensembl_rest",
            endpoint="https://rest.ensembl.org/lookup/id/ENST00000269305",
            reference_choice="GRCh38",
        )

        assert annotation.symbol == "TP53"
        assert annotation.biotype == "protein_coding"
        assert len(annotation.exons) == 2
        assert len(annotation.cds) == 1
        assert annotation.endpoint == "https://rest.ensembl.org/lookup/id/ENST00000269305"

    def test_strand_values(self):
        """Test that strand accepts valid values."""
        # Forward strand
        annotation = TranscriptAnnotation(
            transcript_id="TEST1",
            gene_id="GENE1",
            seq_region_name="1",
            start=1000,
            end=2000,
            strand=1,
            provider="test",
        )
        assert annotation.strand == 1

        # Reverse strand
        annotation = TranscriptAnnotation(
            transcript_id="TEST2",
            gene_id="GENE2",
            seq_region_name="1",
            start=1000,
            end=2000,
            strand=-1,
            provider="test",
        )
        assert annotation.strand == -1


class TestTranscriptAnnotationBundle:
    """Tests for TranscriptAnnotationBundle model."""

    def test_empty_bundle(self):
        """Test creation of empty bundle."""
        reference = ReferenceChoice.disabled("test disabled")
        bundle = TranscriptAnnotationBundle(reference_choice=reference)

        assert bundle.resolved_count == 0
        assert bundle.unresolved_count == 0
        assert bundle.total_requested == 0

    def test_bundle_with_transcripts(self):
        """Test bundle with resolved transcripts."""
        annotation1 = TranscriptAnnotation(
            transcript_id="ENST00000269305",
            gene_id="ENSG00000141510",
            seq_region_name="17",
            start=7661779,
            end=7687550,
            strand=1,
            provider="ensembl_rest",
        )

        annotation2 = TranscriptAnnotation(
            transcript_id="ENST00000359597",
            gene_id="ENSG00000141510",
            seq_region_name="17",
            start=7661779,
            end=7687550,
            strand=1,
            provider="ensembl_rest",
        )

        reference = ReferenceChoice.explicit("GRCh38", reason="test")
        bundle = TranscriptAnnotationBundle(
            transcripts={
                "ENST00000269305": annotation1,
                "ENST00000359597": annotation2,
            },
            unresolved=["ENST99999999"],
            reference_choice=reference,
        )

        assert bundle.resolved_count == 2
        assert bundle.unresolved_count == 1
        assert bundle.total_requested == 3
        assert "ENST00000269305" in bundle.transcripts
        assert "ENST99999999" in bundle.unresolved

    def test_bundle_reference_tracking(self):
        """Test that bundle tracks reference choice properly."""
        # Explicit reference
        ref_explicit = ReferenceChoice.explicit("GRCh38", reason="user provided")
        bundle = TranscriptAnnotationBundle(reference_choice=ref_explicit)
        assert bundle.reference_choice.value == "GRCh38"
        assert bundle.reference_choice.state.value == "explicit"

        # Default reference
        ref_default = ReferenceChoice.default("GRCh38", reason="auto-selected")
        bundle = TranscriptAnnotationBundle(reference_choice=ref_default)
        assert bundle.reference_choice.state.value == "default"

        # Disabled reference
        ref_disabled = ReferenceChoice.disabled("no annotation configured")
        bundle = TranscriptAnnotationBundle(reference_choice=ref_disabled)
        assert bundle.reference_choice.value is None
        assert bundle.reference_choice.state.value == "disabled"
