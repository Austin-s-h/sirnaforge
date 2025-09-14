"""Tests for toy database functionality used in fast CI/CD workflows."""

import tempfile
from pathlib import Path

import pytest

from sirnaforge.core.off_target import build_bwa_index
from sirnaforge.data.base import FastaUtils


@pytest.mark.smoke
def test_toy_transcriptome_database_format():
    """Test that toy transcriptome database has valid FASTA format."""
    toy_db_path = Path(__file__).parent / "data" / "toy_transcriptome_db.fasta"
    assert toy_db_path.exists(), "Toy transcriptome database file should exist"

    # Parse sequences
    sequences = FastaUtils.parse_fasta_to_dict(str(toy_db_path))

    # Should have multiple transcripts
    assert len(sequences) >= 3, "Should have at least 3 toy transcripts"

    # Check sequence format
    for seq_id, sequence in sequences.items():
        assert len(sequence) > 200, f"Transcript {seq_id} should be substantial (>200 bp)"
        assert all(base in "ATGC" for base in sequence), f"Transcript {seq_id} should contain only DNA bases"
        assert sequence.startswith("ATG"), f"Transcript {seq_id} should start with start codon"


@pytest.mark.smoke
def test_toy_mirna_database_format():
    """Test that toy miRNA database has valid FASTA format."""
    toy_db_path = Path(__file__).parent / "data" / "toy_mirna_db.fasta"
    assert toy_db_path.exists(), "Toy miRNA database file should exist"

    # Parse sequences
    sequences = FastaUtils.parse_fasta_to_dict(str(toy_db_path))

    # Should have multiple miRNAs
    assert len(sequences) >= 10, "Should have at least 10 toy miRNAs"

    # Check sequence format
    for seq_id, sequence in sequences.items():
        assert 18 <= len(sequence) <= 25, f"miRNA {seq_id} should be typical miRNA length (18-25 nt)"
        assert all(base in "ATGC" for base in sequence), f"miRNA {seq_id} should contain only DNA bases"
        assert "hsa-" in seq_id, f"miRNA {seq_id} should follow human miRNA naming convention"


@pytest.mark.smoke
@pytest.mark.skipif(
    not Path("/usr/bin/bwa-mem2").exists() and not Path("/usr/local/bin/bwa-mem2").exists(),
    reason="BWA-MEM2 not available",
)
def test_toy_transcriptome_bwa_index_building():
    """Test that BWA index can be built from toy transcriptome database."""
    toy_db_path = Path(__file__).parent / "data" / "toy_transcriptome_db.fasta"

    with tempfile.TemporaryDirectory() as temp_dir:
        index_prefix = str(Path(temp_dir) / "toy_transcriptome")

        # Should not crash
        result_prefix = build_bwa_index(str(toy_db_path), index_prefix)
        assert result_prefix == index_prefix

        # Check that index files were created
        index_files = [f"{index_prefix}.{ext}" for ext in ["0123", "amb", "ann", "bwt.2bit.64", "pac"]]
        for index_file in index_files:
            assert Path(index_file).exists(), f"BWA index file should exist: {index_file}"


@pytest.mark.smoke
@pytest.mark.skipif(
    not Path("/usr/bin/bwa-mem2").exists() and not Path("/usr/local/bin/bwa-mem2").exists(),
    reason="BWA-MEM2 not available",
)
def test_toy_mirna_bwa_index_building():
    """Test that BWA-MEM2 index can be built from toy miRNA database."""
    toy_db_path = Path(__file__).parent / "data" / "toy_mirna_db.fasta"

    with tempfile.TemporaryDirectory() as temp_dir:
        index_prefix = str(Path(temp_dir) / "toy_mirna")

        # Should not crash
        result_prefix = build_bwa_index(str(toy_db_path), index_prefix)
        assert result_prefix == index_prefix

        # Check that BWA-MEM2 index files were created
        index_extensions = ["0123", "amb", "ann", "bwt.2bit.64", "pac"]
        for ext in index_extensions:
            index_file = f"{index_prefix}.{ext}"
            assert Path(index_file).exists(), f"BWA-MEM2 index file should exist: {index_file}"


@pytest.mark.smoke
def test_toy_databases_size_optimization():
    """Test that toy databases are optimized for fast CI/CD."""
    transcriptome_path = Path(__file__).parent / "data" / "toy_transcriptome_db.fasta"
    mirna_path = Path(__file__).parent / "data" / "toy_mirna_db.fasta"

    # Check file sizes are reasonable for fast CI/CD
    transcriptome_size = transcriptome_path.stat().st_size
    mirna_size = mirna_path.stat().st_size

    # Should be small enough for fast processing but large enough for meaningful testing
    assert transcriptome_size < 10000, f"Transcriptome DB too large: {transcriptome_size} bytes"
    assert transcriptome_size > 2000, f"Transcriptome DB too small: {transcriptome_size} bytes"

    assert mirna_size < 5000, f"miRNA DB too large: {mirna_size} bytes"
    assert mirna_size > 1000, f"miRNA DB too small: {mirna_size} bytes"
