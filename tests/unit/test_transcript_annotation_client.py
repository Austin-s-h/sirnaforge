"""Unit tests for transcript annotation clients."""

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.transcript_annotation import EnsemblTranscriptModelClient, VepConsequenceClient
from sirnaforge.models.transcript_annotation import TranscriptAnnotationBundle


class TestEnsemblTranscriptModelClient:
    """Tests for EnsemblTranscriptModelClient."""

    @pytest.fixture
    def client(self):
        """Create a test client instance."""
        return EnsemblTranscriptModelClient(timeout=10)

    @pytest.fixture
    def mock_lookup_response(self):
        """Mock response from Ensembl lookup/id endpoint."""
        return {
            "id": "ENST00000269305",
            "Parent": "ENSG00000141510",
            "display_name": "TP53-201",
            "external_name": "TP53",
            "biotype": "protein_coding",
            "seq_region_name": "17",
            "start": 7661779,
            "end": 7687550,
            "strand": -1,
            "Exon": [
                {
                    "id": "ENSE00003659301",
                    "seq_region_name": "17",
                    "start": 7661779,
                    "end": 7661910,
                    "strand": -1,
                },
                {
                    "id": "ENSE00003659302",
                    "seq_region_name": "17",
                    "start": 7668402,
                    "end": 7669000,
                    "strand": -1,
                },
            ],
            "Translation": {
                "start": 7661779,
                "end": 7669000,
            },
        }

    @pytest.fixture
    def mock_region_response(self):
        """Mock response from Ensembl overlap/region endpoint."""
        return [
            {
                "feature_type": "transcript",
                "id": "ENST00000269305",
                "Parent": "ENSG00000141510",
                "external_name": "TP53",
                "biotype": "protein_coding",
                "seq_region_name": "17",
                "start": 7661779,
                "end": 7687550,
                "strand": -1,
            },
            {
                "feature_type": "exon",
                "id": "ENSE00003659301",
                "Parent": "ENST00000269305",
                "seq_region_name": "17",
                "start": 7661779,
                "end": 7661910,
                "strand": -1,
            },
            {
                "feature_type": "cds",
                "id": "CDS1",
                "Parent": "ENST00000269305",
                "seq_region_name": "17",
                "start": 7661779,
                "end": 7669000,
                "strand": -1,
            },
        ]

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_by_ids_success(self, client, mock_lookup_response):
        """Test successful fetch by IDs."""
        reference = ReferenceChoice.explicit("GRCh38", reason="test")

        with patch("sirnaforge.data.transcript_annotation.aiohttp.ClientSession") as mock_session_class:
            # Create mock response object
            mock_response = MagicMock()
            mock_response.status = 200
            mock_response.json = AsyncMock(return_value=mock_lookup_response)

            # Create mock session with context manager support
            mock_session = MagicMock()
            mock_session.get = MagicMock(return_value=mock_response)
            mock_session.__aenter__ = AsyncMock(return_value=mock_session)
            mock_session.__aexit__ = AsyncMock()

            # Mock response context manager
            mock_response.__aenter__ = AsyncMock(return_value=mock_response)
            mock_response.__aexit__ = AsyncMock()

            # Setup session class to return our mock
            mock_session_class.return_value = mock_session

            # Test fetch
            bundle = await client.fetch_by_ids(
                ["ENST00000269305"],
                species="human",
                reference=reference,
            )

            assert bundle.resolved_count == 1
            assert bundle.unresolved_count == 0
            assert "ENST00000269305" in bundle.transcripts

            annotation = bundle.transcripts["ENST00000269305"]
            assert annotation.gene_id == "ENSG00000141510"
            # Symbol is extracted from display_name which includes transcript variant
            assert annotation.symbol in ("TP53", "TP53-201")  # Can be either depending on response
            assert annotation.biotype == "protein_coding"
            assert annotation.seq_region_name == "17"
            assert annotation.strand == -1
            assert len(annotation.exons) == 2
            assert len(annotation.cds) == 1

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_by_ids_not_found(self, client):
        """Test fetch by IDs with 404 response."""
        reference = ReferenceChoice.explicit("GRCh38", reason="test")

        with patch("sirnaforge.data.transcript_annotation.aiohttp.ClientSession") as mock_session_class:
            # Setup mock 404 response
            mock_response = MagicMock()
            mock_response.status = 404
            mock_response.__aenter__ = AsyncMock(return_value=mock_response)
            mock_response.__aexit__ = AsyncMock()

            # Create mock session
            mock_session = MagicMock()
            mock_session.get = MagicMock(return_value=mock_response)
            mock_session.__aenter__ = AsyncMock(return_value=mock_session)
            mock_session.__aexit__ = AsyncMock()
            mock_session_class.return_value = mock_session

            # Test fetch
            bundle = await client.fetch_by_ids(
                ["ENST99999999"],
                species="human",
                reference=reference,
            )

            assert bundle.resolved_count == 0
            assert bundle.unresolved_count == 1
            assert "ENST99999999" in bundle.unresolved

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_by_ids_access_error(self, client):
        """Test fetch by IDs with access error."""
        reference = ReferenceChoice.explicit("GRCh38", reason="test")

        with patch("sirnaforge.data.transcript_annotation.aiohttp.ClientSession") as mock_session_class:
            # Setup mock 503 response (service unavailable)
            mock_response = MagicMock()
            mock_response.status = 503
            mock_response.__aenter__ = AsyncMock(return_value=mock_response)
            mock_response.__aexit__ = AsyncMock()

            # Create mock session
            mock_session = MagicMock()
            mock_session.get = MagicMock(return_value=mock_response)
            mock_session.__aenter__ = AsyncMock(return_value=mock_session)
            mock_session.__aexit__ = AsyncMock()
            mock_session_class.return_value = mock_session

            # Test fetch - should add to unresolved rather than raising
            bundle = await client.fetch_by_ids(
                ["ENST00000269305"],
                species="human",
                reference=reference,
            )

            assert bundle.unresolved_count == 1

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_fetch_by_regions_success(self, client, mock_region_response):
        """Test successful fetch by genomic regions."""
        reference = ReferenceChoice.explicit("GRCh38", reason="test")

        with patch("sirnaforge.data.transcript_annotation.aiohttp.ClientSession") as mock_session_class:
            # Setup mock response
            mock_response = MagicMock()
            mock_response.status = 200
            mock_response.json = AsyncMock(return_value=mock_region_response)
            mock_response.__aenter__ = AsyncMock(return_value=mock_response)
            mock_response.__aexit__ = AsyncMock()

            # Create mock session
            mock_session = MagicMock()
            mock_session.get = MagicMock(return_value=mock_response)
            mock_session.__aenter__ = AsyncMock(return_value=mock_session)
            mock_session.__aexit__ = AsyncMock()
            mock_session_class.return_value = mock_session

            # Test fetch
            bundle = await client.fetch_by_regions(
                ["17:7661779-7687550"],
                species="human",
                reference=reference,
            )

            assert bundle.resolved_count == 1
            assert "ENST00000269305" in bundle.transcripts

            annotation = bundle.transcripts["ENST00000269305"]
            assert annotation.gene_id == "ENSG00000141510"
            assert len(annotation.exons) == 1
            assert len(annotation.cds) == 1

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_caching_behavior(self, client, mock_lookup_response):
        """Test that caching works correctly."""
        reference = ReferenceChoice.explicit("GRCh38", reason="test")

        with patch("sirnaforge.data.transcript_annotation.aiohttp.ClientSession") as mock_session_class:
            # Setup mock response
            mock_response = MagicMock()
            mock_response.status = 200
            mock_response.json = AsyncMock(return_value=mock_lookup_response)
            mock_response.__aenter__ = AsyncMock(return_value=mock_response)
            mock_response.__aexit__ = AsyncMock()

            # Create mock session
            mock_session = MagicMock()
            mock_get = MagicMock(return_value=mock_response)
            mock_session.get = mock_get
            mock_session.__aenter__ = AsyncMock(return_value=mock_session)
            mock_session.__aexit__ = AsyncMock()
            mock_session_class.return_value = mock_session

            # First fetch - should hit the API
            bundle1 = await client.fetch_by_ids(
                ["ENST00000269305"],
                species="human",
                reference=reference,
            )
            assert bundle1.resolved_count == 1
            assert mock_get.call_count == 1

            # Second fetch - should use cache
            bundle2 = await client.fetch_by_ids(
                ["ENST00000269305"],
                species="human",
                reference=reference,
            )
            assert bundle2.resolved_count == 1
            # Should not make additional API calls
            assert mock_get.call_count == 1

    @pytest.mark.unit
    def test_species_normalization(self, client):
        """Test species name normalization."""
        assert client._normalize_species("human") == "homo_sapiens"
        assert client._normalize_species("mouse") == "mus_musculus"
        assert client._normalize_species("rat") == "rattus_norvegicus"
        assert client._normalize_species("Homo sapiens") == "homo_sapiens"
        assert client._normalize_species("custom_species") == "custom_species"


class TestVepConsequenceClient:
    """Tests for VepConsequenceClient."""

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_vep_client_initialization(self):
        """Test VEP client can be initialized."""
        client = VepConsequenceClient()
        assert client.timeout == 30
        assert client.base_url == "https://rest.ensembl.org"

    @pytest.mark.unit
    @pytest.mark.asyncio
    async def test_vep_enrich_placeholder(self):
        """Test VEP enrichment (currently a placeholder)."""
        client = VepConsequenceClient()
        reference = ReferenceChoice.explicit("GRCh38", reason="test")
        bundle = TranscriptAnnotationBundle(reference_choice=reference)

        # Should return the same bundle (placeholder implementation)
        enriched = await client.enrich_annotations(bundle)
        assert enriched.resolved_count == bundle.resolved_count
        assert enriched.reference_choice == bundle.reference_choice
