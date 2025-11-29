"""Integration tests for toy databases with offtarget workflow."""

from pathlib import Path

import pytest

from sirnaforge.core.off_target import OffTargetAnalysisManager, build_bwa_index
from sirnaforge.data.base import FastaUtils


@pytest.mark.smoke
@pytest.mark.requires_docker
@pytest.mark.integration
def test_toy_transcriptome_analysis_workflow(tmp_path: Path):
    """Test complete transcriptome analysis workflow with toy database."""
    # Get toy database path
    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    transcriptome_db = test_data_dir / "toy_transcriptome_db.fasta"
    candidate_file = test_data_dir / "toy_candidates.fasta"

    assert transcriptome_db.exists(), "Toy transcriptome database should exist"
    assert candidate_file.exists(), "Toy candidate sequences should exist"

    # Build BWA index for transcriptome
    index_prefix = build_bwa_index(transcriptome_db, tmp_path / "transcriptome_index")

    # Run transcriptome analysis
    manager = OffTargetAnalysisManager(species="toy_human", transcriptome_index=index_prefix)

    output_prefix = tmp_path / "test_output"
    tsv_path, json_path = manager.analyze_transcriptome_off_targets(candidate_file, output_prefix)

    # Verify outputs exist
    assert tsv_path.exists(), "TSV output should be created"
    assert json_path.exists(), "JSON output should be created"

    # Verify outputs have content
    assert tsv_path.stat().st_size > 0, "TSV should not be empty"
    assert json_path.stat().st_size > 10, "JSON should have meaningful content"


@pytest.mark.smoke
@pytest.mark.requires_docker
@pytest.mark.integration
def test_toy_mirna_analysis_workflow(tmp_path: Path):
    """Test complete miRNA analysis workflow with toy database."""
    # Get toy database path
    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    mirna_db = test_data_dir / "toy_mirna_db.fasta"
    candidate_file = test_data_dir / "toy_candidates.fasta"

    assert mirna_db.exists(), "Toy miRNA database should exist"
    assert candidate_file.exists(), "Toy candidate sequences should exist"

    # Build BWA index for miRNA (can use BWA for miRNA analysis too)
    index_prefix = build_bwa_index(mirna_db, tmp_path / "mirna_index")

    # Run miRNA analysis
    manager = OffTargetAnalysisManager(species="toy_human", mirna_index=index_prefix)

    output_prefix = tmp_path / "test_output"
    tsv_path, json_path = manager.analyze_mirna_off_targets(candidate_file, output_prefix)

    # Verify outputs exist
    assert tsv_path.exists(), "TSV output should be created"
    assert json_path.exists(), "JSON output should be created"

    # Verify outputs have content
    assert tsv_path.stat().st_size > 0, "TSV should not be empty"
    assert json_path.stat().st_size > 10, "JSON should have meaningful content"


@pytest.mark.requires_docker
@pytest.mark.integration
def test_combined_offtarget_analysis(tmp_path: Path):
    """Test combined transcriptome and miRNA analysis with toy databases."""
    # Get toy database paths
    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    transcriptome_db = test_data_dir / "toy_transcriptome_db.fasta"
    mirna_db = test_data_dir / "toy_mirna_db.fasta"
    candidate_file = test_data_dir / "toy_candidates.fasta"

    # Build indices
    transcriptome_index = build_bwa_index(transcriptome_db, tmp_path / "transcriptome_index")
    mirna_index = build_bwa_index(mirna_db, tmp_path / "mirna_index")

    # Run combined analysis
    manager = OffTargetAnalysisManager(
        species="toy_human", transcriptome_index=transcriptome_index, mirna_index=mirna_index
    )

    # Test combined analysis
    sequences = FastaUtils.parse_fasta_to_dict(candidate_file)

    # Test both transcriptome and miRNA analysis
    transcriptome_tsv, transcriptome_json = manager.analyze_transcriptome_off_targets(
        sequences, tmp_path / "combined_test"
    )
    mirna_tsv, mirna_json = manager.analyze_mirna_off_targets(sequences, tmp_path / "combined_test")

    # Verify output files were created
    assert transcriptome_tsv.exists()
    assert transcriptome_json.exists()
    assert mirna_tsv.exists()
    assert mirna_json.exists()
