"""Integration tests for toy databases with offtarget workflow."""

import shutil
import tempfile
from pathlib import Path

import pytest

from sirnaforge.core.off_target import OffTargetAnalysisManager, build_bwa_index
from sirnaforge.data.base import FastaUtils


@pytest.mark.smoke
@pytest.mark.docker
@pytest.mark.skipif(not shutil.which("bwa-mem2"), reason="bwa-mem2 executable not available - run in Docker container")
def test_toy_transcriptome_analysis_workflow():
    """Test complete transcriptome analysis workflow with toy database."""
    # Get toy database path
    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    transcriptome_db = test_data_dir / "toy_transcriptome_db.fasta"
    candidate_file = test_data_dir / "toy_candidates.fasta"

    assert transcriptome_db.exists(), "Toy transcriptome database should exist"
    assert candidate_file.exists(), "Toy candidate sequences should exist"

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Build BWA index for transcriptome
        index_prefix = str(temp_path / "transcriptome_index")
        build_bwa_index(str(transcriptome_db), index_prefix)

        # Run transcriptome analysis
        manager = OffTargetAnalysisManager(species="toy_human", transcriptome_index=index_prefix)

        output_prefix = str(temp_path / "test_output")
        tsv_path, json_path = manager.analyze_transcriptome_off_targets(str(candidate_file), output_prefix)

        # Verify outputs exist
        assert Path(tsv_path).exists(), "TSV output should be created"
        assert Path(json_path).exists(), "JSON output should be created"

        # Verify outputs have content
        assert Path(tsv_path).stat().st_size > 0, "TSV should not be empty"
        assert Path(json_path).stat().st_size > 10, "JSON should have meaningful content"


@pytest.mark.smoke
@pytest.mark.docker
@pytest.mark.skipif(not shutil.which("bwa-mem2"), reason="bwa-mem2 executable not available - run in Docker container")
def test_toy_mirna_analysis_workflow():
    """Test complete miRNA analysis workflow with toy database."""
    # Get toy database path
    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    mirna_db = test_data_dir / "toy_mirna_db.fasta"
    candidate_file = test_data_dir / "toy_candidates.fasta"

    assert mirna_db.exists(), "Toy miRNA database should exist"
    assert candidate_file.exists(), "Toy candidate sequences should exist"

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Build BWA index for miRNA (can use BWA for miRNA analysis too)
        index_prefix = str(temp_path / "mirna_index")
        build_bwa_index(str(mirna_db), index_prefix)

        # Run miRNA analysis
        manager = OffTargetAnalysisManager(species="toy_human", mirna_index=index_prefix)

        output_prefix = str(temp_path / "test_output")
        tsv_path, json_path = manager.analyze_mirna_off_targets(str(candidate_file), output_prefix)

        # Verify outputs exist
        assert Path(tsv_path).exists(), "TSV output should be created"
        assert Path(json_path).exists(), "JSON output should be created"

        # Verify outputs have content
        assert Path(tsv_path).stat().st_size > 0, "TSV should not be empty"
        assert Path(json_path).stat().st_size > 10, "JSON should have meaningful content"


@pytest.mark.docker_integration
@pytest.mark.docker
@pytest.mark.skipif(not shutil.which("bwa-mem2"), reason="bwa-mem2 executable not available - run in Docker container")
def test_combined_offtarget_analysis():
    """Test combined transcriptome and miRNA analysis with toy databases."""
    # Get toy database paths
    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    transcriptome_db = test_data_dir / "toy_transcriptome_db.fasta"
    mirna_db = test_data_dir / "toy_mirna_db.fasta"
    candidate_file = test_data_dir / "toy_candidates.fasta"

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Build indices
        transcriptome_index = str(temp_path / "transcriptome_index")
        mirna_index = str(temp_path / "mirna_index")

        build_bwa_index(str(transcriptome_db), transcriptome_index)
        build_bwa_index(str(mirna_db), mirna_index)

        # Run combined analysis
        manager = OffTargetAnalysisManager(
            species="toy_human", transcriptome_index=transcriptome_index, mirna_index=mirna_index
        )

        # Test combined analysis
        sequences = FastaUtils.parse_fasta_to_dict(str(candidate_file))

        results = manager.analyze_combined_off_targets(sequences)

        # Should have both transcriptome and miRNA results
        assert "transcriptome_hits" in results
        assert "mirna_hits" in results
        assert isinstance(results["transcriptome_hits"], list)
        assert isinstance(results["mirna_hits"], list)
