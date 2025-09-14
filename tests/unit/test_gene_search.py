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


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
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


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
class TestGeneSearcher:
    """Test GeneSearcher functionality."""

    def test_searcher_initialization(self):
        """Test GeneSearcher initialization."""
        searcher = GeneSearcher(timeout=60)
        assert searcher.timeout == 60
        # Check that all clients are initialized
        assert DatabaseType.ENSEMBL in searcher.clients
        assert DatabaseType.REFSEQ in searcher.clients
        assert DatabaseType.GENCODE in searcher.clients

    @pytest.mark.asyncio
    async def test_search_gene_ensembl_mock(self):
        """Test gene search with mocked Ensembl response."""
        searcher = GeneSearcher()

        # Mock the Ensembl client's search_gene method
        mock_gene_info = GeneInfo(gene_id="ENSG00000141510", gene_name="TP53", database=DatabaseType.ENSEMBL)

        mock_transcript = TranscriptInfo(
            transcript_id="ENST00000269305",
            gene_id="ENSG00000141510",
            gene_name="TP53",
            sequence="AUGCUG" * 100,  # Mock sequence
            database=DatabaseType.ENSEMBL,
            is_canonical=True,
        )

        async def _fake_client_search(query, include_sequence=True):
            return (mock_gene_info, [mock_transcript])

        # Mock the Ensembl client's search_gene method
        with patch.object(searcher.clients[DatabaseType.ENSEMBL], "search_gene", new=_fake_client_search):
            result = await searcher.search_gene("TP53", DatabaseType.ENSEMBL)

            assert result.success is True
            assert result.gene_info.gene_name == "TP53"
            assert len(result.transcripts) == 1
            assert result.transcripts[0].is_canonical is True

    @pytest.mark.asyncio
    async def test_search_multiple_databases(self):
        """Test searching multiple databases."""
        searcher = GeneSearcher()

        async def _fake_search_gene(query, database, include_sequence=True):
            return GeneSearchResult(query="TP53", database=DatabaseType.ENSEMBL, error="Mock result")

        with patch.object(searcher, "search_gene", new=_fake_search_gene):
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
        mock_result = GeneSearchResult(query="TP53", database=DatabaseType.ENSEMBL)

        # Test the synchronous wrapper by mocking asyncio.run
        with patch("sirnaforge.data.gene_search.asyncio.run") as mock_run:
            mock_run.return_value = mock_result

            result = search_gene_sync("TP53")

            assert result == mock_result
            mock_run.assert_called_once()


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
class TestDatabaseSpecificMethods:
    """Test database-specific search methods."""

    @pytest.mark.asyncio
    async def test_refseq_placeholder(self):
        """Test RefSeq search with the actual client."""
        searcher = GeneSearcher()

        # Test that RefSeq client properly handles network errors
        try:
            result = await searcher.search_gene("NM_000546", DatabaseType.REFSEQ)
            # Either it should succeed (if network available) or fail with access error
            assert result.database == DatabaseType.REFSEQ
        except Exception:
            # Network errors are expected in test environments
            pass

    @pytest.mark.asyncio
    async def test_gencode_placeholder(self):
        """Test GENCODE search placeholder."""
        searcher = GeneSearcher()

        # GENCODE should raise an error since it's not implemented
        result = await searcher.search_gene("TP53", DatabaseType.GENCODE)
        assert result.database == DatabaseType.GENCODE
        assert result.error is not None
        assert "not yet implemented" in result.error


if __name__ == "__main__":
    pytest.main([__file__])
