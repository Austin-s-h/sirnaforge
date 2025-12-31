"""Integration tests for transcript annotation with real Ensembl REST API."""

import pytest

from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.transcript_annotation import EnsemblTranscriptModelClient


@pytest.mark.integration
@pytest.mark.requires_network
@pytest.mark.asyncio
async def test_ensembl_fetch_tp53_by_id():
    """Test fetching TP53 transcript annotation from real Ensembl REST API.
    
    This test requires network access and may be skipped in offline environments.
    """
    client = EnsemblTranscriptModelClient(timeout=30)
    reference = ReferenceChoice.explicit("GRCh38", reason="integration test")
    
    # Fetch TP53 canonical transcript
    bundle = await client.fetch_by_ids(
        ids=["ENST00000269305"],  # TP53-201 canonical
        species="human",
        reference=reference,
    )
    
    # Basic validations
    assert bundle.resolved_count >= 1, "Should resolve at least one transcript"
    assert bundle.unresolved_count == 0, "Should not have unresolved IDs"
    assert "ENST00000269305" in bundle.transcripts, "Should contain requested transcript"
    
    # Validate annotation structure
    annotation = bundle.transcripts["ENST00000269305"]
    assert annotation.transcript_id == "ENST00000269305"
    assert annotation.gene_id.startswith("ENSG"), f"Gene ID should be Ensembl format: {annotation.gene_id}"
    assert annotation.symbol == "TP53", "Should have correct gene symbol"
    assert annotation.biotype == "protein_coding", "TP53 should be protein-coding"
    assert annotation.seq_region_name in ("17", "chr17"), "Should be on chromosome 17"
    assert annotation.strand in (1, -1), "Should have valid strand"
    
    # Validate structural features
    assert len(annotation.exons) > 0, "Protein-coding transcript should have exons"
    assert len(annotation.cds) > 0, "Protein-coding transcript should have CDS"
    
    # Validate provenance
    assert annotation.provider == "ensembl_rest"
    assert annotation.endpoint is not None
    assert annotation.reference_choice == "GRCh38"


@pytest.mark.integration
@pytest.mark.requires_network
@pytest.mark.asyncio
async def test_ensembl_fetch_by_region():
    """Test fetching transcripts by genomic region from real Ensembl REST API."""
    client = EnsemblTranscriptModelClient(timeout=30)
    reference = ReferenceChoice.explicit("GRCh38", reason="integration test")
    
    # TP53 region on chr17
    bundle = await client.fetch_by_regions(
        regions=["17:7661779-7687550"],
        species="human",
        reference=reference,
    )
    
    # Should find at least TP53 transcripts
    assert bundle.resolved_count > 0, "Should find transcripts in TP53 region"
    
    # Check that at least one transcript is TP53
    tp53_transcripts = [
        t for t in bundle.transcripts.values()
        if t.symbol == "TP53"
    ]
    assert len(tp53_transcripts) > 0, "Should find at least one TP53 transcript"


@pytest.mark.integration
@pytest.mark.requires_network
@pytest.mark.asyncio
async def test_ensembl_unresolved_id():
    """Test handling of nonexistent transcript ID."""
    client = EnsemblTranscriptModelClient(timeout=30)
    reference = ReferenceChoice.explicit("GRCh38", reason="integration test")
    
    # Use a clearly invalid ID
    bundle = await client.fetch_by_ids(
        ids=["ENST99999999999"],
        species="human",
        reference=reference,
    )
    
    assert bundle.resolved_count == 0, "Invalid ID should not resolve"
    assert bundle.unresolved_count == 1, "Should mark ID as unresolved"
    assert "ENST99999999999" in bundle.unresolved


@pytest.mark.integration
@pytest.mark.requires_network
@pytest.mark.asyncio
async def test_ensembl_mixed_ids():
    """Test fetching mix of valid and invalid transcript IDs."""
    client = EnsemblTranscriptModelClient(timeout=30)
    reference = ReferenceChoice.explicit("GRCh38", reason="integration test")
    
    bundle = await client.fetch_by_ids(
        ids=[
            "ENST00000269305",  # Valid: TP53
            "ENST99999999999",  # Invalid
        ],
        species="human",
        reference=reference,
    )
    
    assert bundle.resolved_count == 1, "Should resolve one valid ID"
    assert bundle.unresolved_count == 1, "Should mark one ID as unresolved"
    assert "ENST00000269305" in bundle.transcripts
    assert "ENST99999999999" in bundle.unresolved
