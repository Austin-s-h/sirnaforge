"""Tests for in-process miRNA seed scanning backends."""

import json
from pathlib import Path

import pandas as pd
import pytest

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


def test_scan_mirna_seed_matches_preserves_seed_mismatch_and_score_semantics() -> None:
    """Seed mismatch counting and score calculation should remain stable across in-process scanning."""
    candidate_sequences = {"candidate_1": "AUGCUAACG"}
    mirna_sequences = {
        "mir_exact": "TGCTAAC",
        "mir_one_seed_mismatch": "TGCTAGC",
    }

    hits = scan_mirna_seed_matches(
        candidate_sequences,
        mirna_sequences,
        backend=MiRNASeedBackend.PYAHOCORASICK,
        seed_start=2,
        seed_end=8,
        max_mismatches=2,
        max_hits=100,
    )

    by_mirna = {hit["rname"]: hit for hit in hits}
    assert set(by_mirna) == {"mir_exact", "mir_one_seed_mismatch"}
    assert by_mirna["mir_exact"]["qseq"] == "ATGCTAACG"
    assert by_mirna["mir_exact"]["nm"] == 0
    assert by_mirna["mir_exact"]["seed_mismatches"] == 0
    assert by_mirna["mir_exact"]["offtarget_score"] == 0.0

    assert by_mirna["mir_one_seed_mismatch"]["nm"] == 1
    assert by_mirna["mir_one_seed_mismatch"]["seed_mismatches"] == 1
    assert by_mirna["mir_one_seed_mismatch"]["offtarget_score"] == 15.0


def test_scan_mirna_seed_matches_uses_full_sequence_when_shorter_than_seed_window() -> None:
    """Short queries should still produce deterministic hits using the full normalized sequence."""
    hits = scan_mirna_seed_matches(
        {"short_candidate": "AUGC"},
        {"short_mirna": "ATGC"},
        backend=MiRNASeedBackend.PYAHOCORASICK,
        seed_start=2,
        seed_end=8,
        max_mismatches=0,
        max_hits=10,
    )

    assert len(hits) == 1
    hit = hits[0]
    assert hit["qseq"] == "ATGC"
    assert hit["cigar"] == "4M"
    assert hit["coord"] == 0
    assert hit["nm"] == 0
    assert hit["seed_mismatches"] == 0
    assert hit["offtarget_score"] == 0.0


def test_scan_mirna_seed_matches_pyahocorasick_dependency_error_is_clear(monkeypatch) -> None:
    """Missing optional pyahocorasick dependency should fail with actionable guidance."""

    def _missing_dependency(name: str):
        if name == "ahocorasick":
            raise ImportError("simulated missing dependency")
        return __import__(name)

    monkeypatch.setattr("importlib.import_module", _missing_dependency)

    with pytest.raises(RuntimeError, match="required 'pyahocorasick' package"):
        scan_mirna_seed_matches(
            {"candidate_1": "AUGCUAACG"},
            {"mir_exact": "TGCTAAC"},
            backend=MiRNASeedBackend.PYAHOCORASICK,
        )
