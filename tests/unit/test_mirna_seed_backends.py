"""Tests for in-process miRNA seed scanning backends."""

import json
from pathlib import Path

import pandas as pd

from sirnaforge.core.off_target import (
    MiRNASeedBackend,
    mirna_seed_hit_identity,
    parse_fasta_file,
    run_mirna_seed_analysis,
    scan_mirna_seed_matches,
)
from sirnaforge.models.schemas import MiRNAAlignmentSchema


def test_scan_mirna_seed_matches_pyahocorasick_matches_exhaustive_python() -> None:
    """Optimized seed scanning should match the exhaustive semantic hit set."""
    test_data_dir = Path(__file__).parent / "data"
    candidate_sequences = parse_fasta_file(test_data_dir / "toy_candidates.fasta")
    mirna_sequences = parse_fasta_file(test_data_dir / "toy_mirna_db.fasta")

    exhaustive_hits = scan_mirna_seed_matches(
        candidate_sequences,
        mirna_sequences,
        backend=MiRNASeedBackend.EXHAUSTIVE_PYTHON,
    )
    pyahocorasick_hits = scan_mirna_seed_matches(
        candidate_sequences,
        mirna_sequences,
        backend=MiRNASeedBackend.PYAHOCORASICK,
    )

    assert exhaustive_hits, "toy data should produce exhaustive seed hits"
    assert pyahocorasick_hits, "toy data should produce pyahocorasick seed hits"
    assert {mirna_seed_hit_identity(hit) for hit in exhaustive_hits} == {
        mirna_seed_hit_identity(hit) for hit in pyahocorasick_hits
    }


def test_run_mirna_seed_analysis_pyahocorasick_writes_schema_compatible_outputs(tmp_path: Path, monkeypatch) -> None:
    """In-process seed scanning should keep the existing miRNA output contract stable."""
    test_data_dir = Path(__file__).parent / "data"
    candidate_file = test_data_dir / "toy_candidates.fasta"
    mirna_db = test_data_dir / "toy_mirna_db.fasta"

    monkeypatch.setattr(
        "sirnaforge.data.mirna_manager.MiRNADatabaseManager.get_database",
        lambda *_args, **_kwargs: mirna_db,
    )

    output_dir = run_mirna_seed_analysis(
        candidates_file=candidate_file,
        candidate_id="toy",
        mirna_db="toy_db",
        mirna_species=["human"],
        output_dir=tmp_path,
        backend=MiRNASeedBackend.PYAHOCORASICK,
    )

    analysis_path = output_dir / "toy_mirna_analysis.tsv"
    summary_path = output_dir / "toy_mirna_summary.json"

    assert analysis_path.exists()
    assert summary_path.exists()

    analysis_df = pd.read_csv(analysis_path, sep="\t")
    validated_df = MiRNAAlignmentSchema.validate(analysis_df, lazy=True)
    assert not validated_df.empty

    summary = json.loads(summary_path.read_text())
    assert summary["total_hits"] == len(validated_df)
    assert summary["total_raw_alignments"] == len(validated_df)
    assert summary["hits_per_species"]["human"] == len(validated_df)
