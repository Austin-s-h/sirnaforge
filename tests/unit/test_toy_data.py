"""Tests for toy data used in CI/CD smoke tests."""

from pathlib import Path

import pytest


@pytest.mark.unit
@pytest.mark.ci
@pytest.mark.smoke
def test_toy_transcripts_exist():
    """Test that toy transcript data exists and is valid."""
    toy_file = Path(__file__).parent / "data" / "toy_transcripts.fasta"
    assert toy_file.exists(), "Toy transcripts file should exist"
    
    content = toy_file.read_text()
    assert content.startswith(">"), "Should be valid FASTA format"
    assert "toy_transcript" in content, "Should contain toy transcript headers"
    
    # Should be small for fast CI/CD
    assert len(content) < 500, "Toy data should be minimal for fast CI/CD"


@pytest.mark.unit
@pytest.mark.ci  
@pytest.mark.smoke
def test_toy_candidates_exist():
    """Test that toy siRNA candidate data exists and is valid."""
    toy_file = Path(__file__).parent / "data" / "toy_candidates.fasta"
    assert toy_file.exists(), "Toy candidates file should exist"
    
    content = toy_file.read_text()
    assert content.startswith(">"), "Should be valid FASTA format"
    assert "toy_sirna" in content, "Should contain toy siRNA headers"
    
    # Should be minimal for smoke tests
    assert len(content) < 200, "Toy data should be ultra-minimal"


@pytest.mark.unit
@pytest.mark.ci
def test_toy_data_sequences_valid():
    """Test that toy sequences contain valid RNA bases."""
    toy_file = Path(__file__).parent / "data" / "toy_transcripts.fasta"
    content = toy_file.read_text()
    
    lines = content.strip().split('\n')
    sequences = [line for line in lines if not line.startswith('>')]
    
    for seq in sequences:
        # Should only contain valid RNA bases
        valid_bases = set('AUGC')
        seq_bases = set(seq.upper())
        assert seq_bases.issubset(valid_bases), f"Invalid bases in sequence: {seq_bases - valid_bases}"
        
        # Should be reasonable length for siRNA design
        assert 20 <= len(seq) <= 100, f"Sequence length {len(seq)} not suitable for siRNA design"