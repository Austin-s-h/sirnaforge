"""Explicit local-only checks of the remote NCBI RefSeq retrieval path."""

import pytest

from sirnaforge.data.gene_search import DatabaseType, GeneSearcher


@pytest.mark.integration
@pytest.mark.requires_network
@pytest.mark.asyncio
async def test_refseq_gene_search_live():
    """Require real gene/transcript retrieval rather than accepting network errors."""
    result = await GeneSearcher().search_gene("TP53", DatabaseType.REFSEQ)
    assert result.database == DatabaseType.REFSEQ
    assert result.success, result.error
    assert any(transcript.sequence for transcript in result.transcripts)
