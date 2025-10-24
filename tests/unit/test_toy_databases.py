"""Tests for toy database functionality used in fast CI/CD workflows."""

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
    sequences = FastaUtils.parse_fasta_to_dict(toy_db_path)

    # Should have multiple transcripts
    assert len(sequences) >= 3, "Should have at least 3 toy transcripts"

    # Check sequence format
    for seq_id, sequence in sequences.items():
        assert len(sequence) > 200, f"Transcript {seq_id} should be substantial (>200 bp)"
        assert all(base in "ATGC" for base in sequence), f"Transcript {seq_id} should contain only DNA bases"


@pytest.mark.smoke
def test_toy_mirna_database_format():
    """Test that toy miRNA database has valid FASTA format."""
    toy_db_path = Path(__file__).parent / "data" / "toy_mirna_db.fasta"
    assert toy_db_path.exists(), "Toy miRNA database file should exist"

    # Parse sequences
    sequences = FastaUtils.parse_fasta_to_dict(toy_db_path)

    # Should have multiple miRNAs
    assert len(sequences) >= 10, "Should have at least 10 toy miRNAs"

    # Check sequence format
    for seq_id, sequence in sequences.items():
        assert 18 <= len(sequence) <= 25, f"miRNA {seq_id} should be typical miRNA length (18-25 nt)"
        assert all(base in "ATGC" for base in sequence), f"miRNA {seq_id} should contain only DNA bases"
        assert "hsa-" in seq_id, f"miRNA {seq_id} should follow human miRNA naming convention"


@pytest.mark.smoke
@pytest.mark.skipif(not Path("bwa-mem2").exists(), reason="BWA-MEM2 not available")
def test_toy_transcriptome_bwa_index_building(tmp_path: Path):
    """Test that BWA index can be built from toy transcriptome database."""
    toy_db_path = Path(__file__).parent / "data" / "toy_transcriptome_db.fasta"

    index_prefix = build_bwa_index(toy_db_path, tmp_path / "toy_transcriptome")

    # Should not crash (implicit by reaching here)
    assert index_prefix == tmp_path / "toy_transcriptome"

    # Check that index files were created
    for ext in ["0123", "amb", "ann", "bwt.2bit.64", "pac"]:
        index_file = index_prefix.parent / f"{index_prefix.name}.{ext}"
        assert index_file.exists(), f"BWA index file should exist: {index_file}"


@pytest.mark.smoke
@pytest.mark.skipif(not Path("bwa-mem2").exists(), reason="BWA-MEM2 not available")
def test_toy_mirna_bwa_index_building(tmp_path: Path):
    """Test that BWA-MEM2 index can be built from toy miRNA database."""
    toy_db_path = Path(__file__).parent / "data" / "toy_mirna_db.fasta"

    index_prefix = build_bwa_index(toy_db_path, tmp_path / "toy_mirna")

    assert index_prefix == tmp_path / "toy_mirna"

    # Check that BWA-MEM2 index files were created
    for ext in ["0123", "amb", "ann", "bwt.2bit.64", "pac"]:
        index_file = index_prefix.parent / f"{index_prefix.name}.{ext}"
        assert index_file.exists(), f"BWA-MEM2 index file should exist: {index_file}"
