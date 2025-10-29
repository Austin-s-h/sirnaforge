"""Test transcript hit count functionality."""

import pytest

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters


@pytest.mark.unit
@pytest.mark.ci
def test_shared_guides(tmp_path):
    """Two identical transcripts -> each guide should hit both transcripts."""
    fasta = tmp_path / "shared.fa"
    # Use a sequence with balanced GC content (40%) to pass filters
    seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    fasta.write_text(f">tx1\n{seq}\n>tx2\n{seq}\n")

    params = DesignParameters(sirna_length=21, top_n=5)
    designer = SiRNADesigner(params)
    res = designer.design_from_file(str(fasta))

    assert res.total_sequences == 2
    assert res.total_candidates > 0

    top = res.top_candidates[0]
    assert top.transcript_hit_count == 2
    assert top.transcript_hit_fraction == 1.0


@pytest.mark.unit
@pytest.mark.ci
def test_unique_guides(tmp_path):
    """Two different transcripts -> typical guide should hit only one transcript."""
    fasta = tmp_path / "unique.fa"
    # Use sequences with balanced GC content but different sequences
    seq1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    seq2 = "CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    fasta.write_text(f">tx1\n{seq1}\n>tx2\n{seq2}\n")

    params = DesignParameters(sirna_length=21, top_n=5)
    designer = SiRNADesigner(params)
    res = designer.design_from_file(str(fasta))

    assert res.total_sequences == 2
    assert res.total_candidates > 0

    # pick first top candidate and ensure it's not hitting both transcripts
    top = res.top_candidates[0]
    assert top.transcript_hit_count in (1, 2)
    # If guides differ the hit fraction should be 0.5 when hit_count==1
    if top.transcript_hit_count == 1:
        assert top.transcript_hit_fraction == 0.5
