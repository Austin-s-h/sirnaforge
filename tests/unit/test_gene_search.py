"""Tests for gene search functionality."""

from unittest.mock import patch

import pytest

from sirnaforge.data.gene_search import (
    DatabaseType,
    GeneInfo,
    GeneSearcher,
    GeneSearchResult,
    TranscriptInfo,
    search_gene_sync,
)


class TestGeneSearchModels:
    """Test gene search data models."""

    def test_gene_info_creation(self):
        """Test GeneInfo model creation."""
        gene = GeneInfo(
            gene_id="ENSG00000123456", gene_name="TEST_GENE", gene_type="protein_coding", database=DatabaseType.ENSEMBL
        )
        assert gene.gene_id == "ENSG00000123456"
        assert gene.database == DatabaseType.ENSEMBL

    def test_transcript_info_validation(self):
        """Test TranscriptInfo sequence validation."""
        # Valid RNA sequence
        transcript = TranscriptInfo(
            transcript_id="ENST00000123456",
            gene_id="ENSG00000123456",
            sequence="AUGCUGAUC",
            database=DatabaseType.ENSEMBL,
        )
        assert transcript.sequence == "AUGCUGAUC"

        # Invalid sequence should raise error
        with pytest.raises(ValueError, match="invalid RNA bases"):
            TranscriptInfo(
                transcript_id="ENST00000123456",
                gene_id="ENSG00000123456",
                sequence="ATGCXGAUC",  # X is invalid
                database=DatabaseType.ENSEMBL,
            )

    def test_gene_search_result_success(self):
        """Test GeneSearchResult success property."""
        gene_info = GeneInfo(gene_id="ENSG00000123456", database=DatabaseType.ENSEMBL)
        transcript = TranscriptInfo(
            transcript_id="ENST00000123456", gene_id="ENSG00000123456", database=DatabaseType.ENSEMBL
        )

        # Successful result
        result = GeneSearchResult(
            query="TEST_GENE", database=DatabaseType.ENSEMBL, gene_info=gene_info, transcripts=[transcript]
        )
        assert result.success is True

        # Failed result
        failed_result = GeneSearchResult(
            query="NONEXISTENT_GENE", database=DatabaseType.ENSEMBL, error="Gene not found"
        )
        assert failed_result.success is False


class TestGeneSearcher:
    """Test GeneSearcher functionality."""

    def test_searcher_initialization(self):
        """Test GeneSearcher initialization."""
        searcher = GeneSearcher(preferred_db=DatabaseType.ENSEMBL, timeout=60)
        assert searcher.preferred_db == DatabaseType.ENSEMBL
        assert searcher.timeout == 60

    @pytest.mark.asyncio
    async def test_search_gene_ensembl_mock(self):
        """Test gene search with mocked Ensembl response."""
        searcher = GeneSearcher()

        # Mock the Ensembl-specific methods
        mock_gene_info = GeneInfo(gene_id="ENSG00000141510", gene_name="TP53", database=DatabaseType.ENSEMBL)

        mock_transcript = TranscriptInfo(
            transcript_id="ENST00000269305",
            gene_id="ENSG00000141510",
            gene_name="TP53",
            sequence="AUGCUG" * 100,  # Mock sequence
            database=DatabaseType.ENSEMBL,
            is_canonical=True,
        )

        with patch.object(searcher, "_search_ensembl") as mock_ensembl:
            mock_ensembl.return_value = GeneSearchResult(
                query="TP53", database=DatabaseType.ENSEMBL, gene_info=mock_gene_info, transcripts=[mock_transcript]
            )

            result = await searcher.search_gene("TP53", DatabaseType.ENSEMBL)

            assert result.success is True
            assert result.gene_info.gene_name == "TP53"
            assert len(result.transcripts) == 1
            assert result.transcripts[0].is_canonical is True

    @pytest.mark.asyncio
    async def test_search_multiple_databases(self):
        """Test searching multiple databases."""
        searcher = GeneSearcher()

        with patch.object(searcher, "search_gene") as mock_search:
            mock_search.return_value = GeneSearchResult(
                query="TP53", database=DatabaseType.ENSEMBL, error="Mock result"
            )

            results = await searcher.search_multiple_databases("TP53", [DatabaseType.ENSEMBL, DatabaseType.REFSEQ])

            assert len(results) == 2
            assert all(result.query == "TP53" for result in results)

    def test_save_transcripts_fasta(self, tmp_path):
        """Test saving transcripts to FASTA format."""
        searcher = GeneSearcher()

        transcripts = [
            TranscriptInfo(
                transcript_id="ENST00000269305",
                gene_id="ENSG00000141510",
                gene_name="TP53",
                sequence="AUGCUGAUCCUG",
                database=DatabaseType.ENSEMBL,
                is_canonical=True,
            ),
            TranscriptInfo(
                transcript_id="ENST00000359597",
                gene_id="ENSG00000141510",
                gene_name="TP53",
                sequence="AUGCUGAUCCUGAAA",
                database=DatabaseType.ENSEMBL,
            ),
        ]

        output_file = tmp_path / "test_output.fasta"
        searcher.save_transcripts_fasta(transcripts, output_file)

        # Check file was created and contains expected content
        assert output_file.exists()
        content = output_file.read_text()

        assert ">ENST00000269305" in content
        assert ">ENST00000359597" in content
        assert "AUGCUGAUCCUG" in content
        assert "canonical:true" in content

    def test_synchronous_wrapper(self):
        """Test synchronous wrapper function."""
        with patch("sirnaforge.data.gene_search.asyncio.run") as mock_run:
            mock_result = GeneSearchResult(query="TP53", database=DatabaseType.ENSEMBL)
            mock_run.return_value = mock_result

            result = search_gene_sync("TP53")

            assert result == mock_result
            mock_run.assert_called_once()


class TestDatabaseSpecificMethods:
    """Test database-specific search methods."""

    @pytest.mark.asyncio
    async def test_refseq_placeholder(self):
        """Test RefSeq search placeholder."""
        searcher = GeneSearcher()
        result = await searcher._search_refseq("NM_000546", True)

        assert result.database == DatabaseType.REFSEQ
        assert result.error == "RefSeq search not yet implemented"

    @pytest.mark.asyncio
    async def test_gencode_placeholder(self):
        """Test GENCODE search placeholder."""
        searcher = GeneSearcher()
        result = await searcher._search_gencode("TP53", True)

        assert result.database == DatabaseType.GENCODE
        assert result.error == "GENCODE search not yet implemented"


if __name__ == "__main__":
    pytest.main([__file__])
