"""Integration tests for toy databases with offtarget workflow.

These tests require BWA-MEM2 which is only available inside the Docker container.
They are marked with `runs_in_container` to run via `make docker-test`.
"""

import json
import shutil
from pathlib import Path

import pytest

from sirnaforge.core.off_target import (
    MiRNASeedBackend,
    OffTargetAnalysisManager,
    build_bwa_index,
    mirna_seed_hit_identity,
    run_mirna_seed_analysis,
)
from sirnaforge.data.base import FastaUtils


def _load_semantic_mirna_hit_identities(
    output_dir: Path, candidate_id: str
) -> set[tuple[str, str, int, int, int, float]]:
    """Load normalized miRNA hits from a backend run and project them onto the parity contract."""
    hits_path = output_dir / f"{candidate_id}_mirna_hits_raw.json"
    hits = json.loads(hits_path.read_text())
    return {mirna_seed_hit_identity(hit) for hit in hits}


def _assert_semantic_hit_parity(
    *,
    baseline_hits: set[tuple[str, str, int, int, int, float]],
    candidate_hits: set[tuple[str, str, int, int, int, float]],
    baseline_label: str,
    candidate_label: str,
) -> None:
    """Assert semantic parity with a diff that makes backend drift obvious to maintainers."""
    missing_hits = sorted(baseline_hits - candidate_hits)
    extra_hits = sorted(candidate_hits - baseline_hits)
    assert candidate_hits == baseline_hits, (
        f"{candidate_label} semantic miRNA seed hits diverged from {baseline_label}: "
        f"missing={missing_hits[:10]}, extra={extra_hits[:10]}, "
        f"baseline_count={len(baseline_hits)}, candidate_count={len(candidate_hits)}"
    )


@pytest.mark.runs_in_container
@pytest.mark.smoke
@pytest.mark.integration
def test_toy_transcriptome_analysis_workflow(tmp_path: Path):
    """Test complete transcriptome analysis workflow with toy database."""
    if shutil.which("bwa-mem2") is None:
        pytest.skip("bwa-mem2 not available - run this test in the Docker container")

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


@pytest.mark.runs_in_container
@pytest.mark.smoke
@pytest.mark.integration
def test_toy_mirna_analysis_workflow(tmp_path: Path):
    """Test complete miRNA analysis workflow with toy database."""
    if shutil.which("bwa-mem2") is None:
        pytest.skip("bwa-mem2 not available - run this test in the Docker container")

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


@pytest.mark.runs_in_container
@pytest.mark.integration
def test_combined_offtarget_analysis(tmp_path: Path):
    """Test combined transcriptome and miRNA analysis with toy databases."""
    if shutil.which("bwa-mem2") is None:
        pytest.skip("bwa-mem2 not available - run this test in the Docker container")

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


@pytest.mark.runs_in_container
@pytest.mark.integration
def test_toy_mirna_seed_backend_matches_bwa_semantic_hits(tmp_path: Path):
    """The optimized miRNA backend should match BWA on toy-data semantic hit identities."""
    if shutil.which("bwa-mem2") is None:
        pytest.skip("bwa-mem2 not available - run this test in the Docker container")

    test_data_dir = Path(__file__).parent.parent / "unit" / "data"
    mirna_db = test_data_dir / "toy_mirna_bwa_parity_db.fasta"
    candidate_file = test_data_dir / "toy_mirna_bwa_parity_candidates.fasta"
    candidate_id = "toy"

    assert mirna_db.exists(), "Semantic-parity miRNA toy database should exist"
    assert candidate_file.exists(), "Toy candidate sequences should exist"

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(
        "sirnaforge.data.mirna_manager.MiRNADatabaseManager.get_database",
        lambda *_args, **_kwargs: mirna_db,
    )
    try:
        bwa_output_dir = run_mirna_seed_analysis(
            candidates_file=candidate_file,
            candidate_id=candidate_id,
            mirna_db="toy_db",
            mirna_species=["human"],
            output_dir=tmp_path / "bwa",
            backend=MiRNASeedBackend.BWA,
        )
        optimized_output_dir = run_mirna_seed_analysis(
            candidates_file=candidate_file,
            candidate_id=candidate_id,
            mirna_db="toy_db",
            mirna_species=["human"],
            output_dir=tmp_path / "optimized",
            backend=MiRNASeedBackend.PYAHOCORASICK,
        )
    finally:
        monkeypatch.undo()

    bwa_identity = _load_semantic_mirna_hit_identities(bwa_output_dir, candidate_id)
    optimized_identity = _load_semantic_mirna_hit_identities(optimized_output_dir, candidate_id)

    assert bwa_identity, "toy fixture should produce BWA miRNA seed hits"
    assert optimized_identity, "toy fixture should produce optimized miRNA seed hits"

    _assert_semantic_hit_parity(
        baseline_hits=bwa_identity,
        candidate_hits=optimized_identity,
        baseline_label="bwa",
        candidate_label="pyahocorasick",
    )
