"""Tests for transcriptome filtering functionality."""

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from sirnaforge.data.transcriptome_filter import TranscriptFilter, get_filter_spec


@pytest.mark.dev
class TestTranscriptFilter:
    """Test suite for transcriptome filtering."""

    def test_parse_ensembl_header(self):
        """Test parsing Ensembl FASTA headers."""
        header = "ENST00000456328.2 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:BRCA1"
        metadata = TranscriptFilter.parse_ensembl_header(header)
        
        assert metadata["gene_biotype"] == "protein_coding"
        assert metadata["transcript_biotype"] == "protein_coding"
        assert metadata["gene_symbol"] == "BRCA1"

    def test_is_protein_coding_positive(self):
        """Test identifying protein-coding transcripts."""
        record = SeqRecord(
            Seq("ATCG"),
            id="ENST001",
            description="gene_biotype:protein_coding"
        )
        assert TranscriptFilter.is_protein_coding(record)

    def test_is_protein_coding_negative(self):
        """Test identifying non-protein-coding transcripts."""
        record = SeqRecord(
            Seq("ATCG"),
            id="ENST002",
            description="gene_biotype:lncRNA"
        )
        assert not TranscriptFilter.is_protein_coding(record)

    def test_is_canonical_positive(self):
        """Test identifying canonical transcripts."""
        test_cases = [
            "ENST001 canonical:1",
            "ENST002 tag:MANE_Select",
            "ENST003 tag:basic",
        ]
        
        for desc in test_cases:
            record = SeqRecord(Seq("ATCG"), id="test", description=desc)
            assert TranscriptFilter.is_canonical(record), f"Failed for: {desc}"

    def test_is_canonical_negative(self):
        """Test identifying non-canonical transcripts."""
        record = SeqRecord(
            Seq("ATCG"),
            id="ENST002",
            description="gene_biotype:protein_coding"
        )
        assert not TranscriptFilter.is_canonical(record)

    def test_filter_fasta_protein_coding(self, tmp_path):
        """Test filtering FASTA for protein-coding transcripts."""
        input_fasta = tmp_path / "input.fa"
        output_fasta = tmp_path / "output.fa"
        
        # Create input with mixed biotypes
        input_fasta.write_text(
            ">ENST001 gene_biotype:protein_coding\nATCG\n"
            ">ENST002 gene_biotype:lncRNA\nGCTA\n"
            ">ENST003 gene_biotype:protein_coding\nTTAA\n"
        )
        
        kept = TranscriptFilter.apply_protein_coding_filter(input_fasta, output_fasta)
        
        assert kept == 2
        output_text = output_fasta.read_text()
        assert "ENST001" in output_text
        assert "ENST003" in output_text
        assert "ENST002" not in output_text

    def test_get_filter_spec_valid(self):
        """Test parsing valid filter specifications."""
        filters = get_filter_spec("protein_coding,canonical_only")
        assert filters == ["protein_coding", "canonical_only"]

    def test_get_filter_spec_invalid(self):
        """Test parsing invalid filter specifications."""
        with pytest.raises(ValueError, match="Invalid filter"):
            get_filter_spec("invalid_filter")

    def test_get_filter_spec_empty(self):
        """Test parsing empty filter specification."""
        filters = get_filter_spec(None)
        assert filters == []
        
        filters = get_filter_spec("")
        assert filters == []
