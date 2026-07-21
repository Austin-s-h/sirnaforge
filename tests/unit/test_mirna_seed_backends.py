"""Tests for in-process miRNA seed scanning backends."""

import importlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from sirnaforge.core.off_target import (
    MiRNASeedBackend,
    mirna_seed_hit_identity,
    normalize_mirna_seed_hit,
    parse_fasta_file,
    run_mirna_seed_analysis,
    scan_mirna_seed_matches,
)
from sirnaforge.models.schemas import MiRNAAlignmentSchema

TEST_DATA_DIR = Path(__file__).parent / "data"
TOY_PARITY_CANDIDATES = TEST_DATA_DIR / "toy_mirna_parity_candidates.fasta"
TOY_PARITY_DB = TEST_DATA_DIR / "toy_mirna_parity_db.fasta"

# The parity fixtures exercise the seed-region contract explicitly. Guide seeds
# (guide positions 2-8) are scanned across each miRNA; ``coord`` is the 0-based start
# of that match in the miRNA. A match is only a real seed off-target when the guide seed
# lands on the miRNA's own seed region (0-based coord == seed_start - 1 == 1).
#
# Fixture design:
# - mir_seed_exact (UACGUACG...ACGUACG): carries cand_exact's seed at coord 1 (seed hit)
#   plus a perfect but NON-seed copy at coord 9 -> raw-only, must NOT count as a hit.
# - mir_seed_one: carries cand_one's seed at coord 1.
# - mir_nonseed: only contains cand_exact's seed away from the seed region (coord 5) ->
#   perfect match but raw-only, never a hit.
# - cand_pos2 differs from cand_exact only at guide position 2 (window position 1); it is
#   a regression guard for the seed-window -> guide coordinate mapping. A position-2
#   mismatch must be counted in seed_mismatches (=1), NOT silently dropped to a perfect seed.

# RAW alignments (every guide-seed placement within max_mismatches), keyed by semantic identity
# (qname, mirna_id, coord, nm, seed_mismatches, offtarget_score).
TOY_PARITY_RAW_SEMANTIC_IDENTITIES = {
    ("cand_exact", "mir_nonseed", 5, 0, 0, 0.0),
    ("cand_exact", "mir_seed_exact", 1, 0, 0, 0.0),
    ("cand_exact", "mir_seed_exact", 5, 1, 1, 15.0),
    ("cand_exact", "mir_seed_exact", 9, 0, 0, 0.0),
    ("cand_exact", "mir_seed_one", 1, 1, 1, 15.0),
    ("cand_one", "mir_nonseed", 5, 1, 1, 15.0),
    ("cand_one", "mir_seed_exact", 1, 1, 1, 15.0),
    ("cand_one", "mir_seed_exact", 5, 2, 2, 32.0),
    ("cand_one", "mir_seed_exact", 9, 1, 1, 15.0),
    ("cand_one", "mir_seed_one", 1, 0, 0, 0.0),
    ("cand_pos2", "mir_nonseed", 5, 1, 1, 15.0),
    ("cand_pos2", "mir_seed_exact", 1, 1, 1, 15.0),
    ("cand_pos2", "mir_seed_exact", 5, 2, 2, 30.0),
    ("cand_pos2", "mir_seed_exact", 9, 1, 1, 15.0),
    ("cand_pos2", "mir_seed_one", 1, 2, 2, 30.0),
}

# FILTERED hits: only the raw alignments whose guide seed lands on the miRNA seed region
# (0-based coord == 1). Perfect non-seed matches (coord 5 / coord 9) are excluded.
TOY_PARITY_SEMANTIC_IDENTITIES = {identity for identity in TOY_PARITY_RAW_SEMANTIC_IDENTITIES if identity[2] == 1}


def _run_toy_mirna_seed_analysis(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    backend: MiRNASeedBackend,
    candidate_id: str,
) -> Path:
    monkeypatch.setattr(
        "sirnaforge.data.mirna_manager.MiRNADatabaseManager.get_database",
        lambda *_args, **_kwargs: TOY_PARITY_DB,
    )
    return run_mirna_seed_analysis(
        candidates_file=TOY_PARITY_CANDIDATES,
        candidate_id=candidate_id,
        mirna_db="toy_db",
        mirna_species=["human"],
        output_dir=tmp_path / candidate_id,
        backend=backend,
    )


def _semantic_identities(rows: list[dict[str, Any]]) -> set[tuple[str, str, int, int, int, float]]:
    return {mirna_seed_hit_identity(row) for row in rows}


def _read_tsv_rows(path: Path) -> list[dict[str, Any]]:
    frame = pd.read_csv(path, sep="	")
    return frame.to_dict(orient="records")


def _read_json_rows(path: Path) -> list[dict[str, Any]]:
    return json.loads(path.read_text())


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


def test_deterministic_mirna_parity_fixture_locks_exhaustive_semantic_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The deterministic parity fixture should lock the exhaustive semantic miRNA contract."""
    output_dir = _run_toy_mirna_seed_analysis(
        tmp_path,
        monkeypatch,
        backend=MiRNASeedBackend.EXHAUSTIVE_PYTHON,
        candidate_id="toy_parity_exhaustive",
    )

    analysis_rows = _read_tsv_rows(output_dir / "toy_parity_exhaustive_mirna_analysis.tsv")
    raw_analysis_rows = _read_tsv_rows(output_dir / "toy_parity_exhaustive_mirna_analysis_raw.tsv")
    hits_rows = _read_json_rows(output_dir / "toy_parity_exhaustive_mirna_hits.json")
    raw_hits_rows = _read_json_rows(output_dir / "toy_parity_exhaustive_mirna_hits_raw.json")
    summary = json.loads((output_dir / "toy_parity_exhaustive_mirna_summary.json").read_text())

    # Filtered outputs keep only seed-region hits (coord == 1); raw outputs keep every alignment.
    assert _semantic_identities(analysis_rows) == TOY_PARITY_SEMANTIC_IDENTITIES
    assert _semantic_identities(hits_rows) == TOY_PARITY_SEMANTIC_IDENTITIES
    assert _semantic_identities(raw_analysis_rows) == TOY_PARITY_RAW_SEMANTIC_IDENTITIES
    assert _semantic_identities(raw_hits_rows) == TOY_PARITY_RAW_SEMANTIC_IDENTITIES

    # Perfect matches outside the seed region are retained as raw alignments but must NOT be hits.
    assert TOY_PARITY_RAW_SEMANTIC_IDENTITIES - TOY_PARITY_SEMANTIC_IDENTITIES, (
        "fixture must include non-seed alignments to exercise the seed-region filter"
    )
    assert all(identity[2] == 1 for identity in _semantic_identities(analysis_rows))

    assert summary["total_hits"] == len(TOY_PARITY_SEMANTIC_IDENTITIES)
    assert summary["total_raw_alignments"] == len(TOY_PARITY_RAW_SEMANTIC_IDENTITIES)
    # hits_per_species reports RAW alignment counts (documented on MiRNASummary).
    assert summary["hits_per_species"] == {"human": len(TOY_PARITY_RAW_SEMANTIC_IDENTITIES)}


@pytest.mark.parametrize("backend", [MiRNASeedBackend.EXHAUSTIVE_PYTHON, MiRNASeedBackend.PYAHOCORASICK])
def test_perfect_non_seed_match_is_raw_only_not_a_hit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    backend: MiRNASeedBackend,
) -> None:
    """A perfect guide-seed match outside the miRNA seed region must be raw-only, never a hit.

    ``mir_nonseed`` and the coord-9 copy inside ``mir_seed_exact`` are perfect (nm==0) matches
    that fall outside the seed region; they must appear in the raw outputs but must NOT be
    counted in the filtered analysis or ``total_hits``.
    """
    output_dir = _run_toy_mirna_seed_analysis(
        tmp_path, monkeypatch, backend=backend, candidate_id=f"nonseed_{backend.value}"
    )

    filtered = _read_tsv_rows(output_dir / f"nonseed_{backend.value}_mirna_analysis.tsv")
    raw = _read_tsv_rows(output_dir / f"nonseed_{backend.value}_mirna_analysis_raw.tsv")

    # The non-seed perfect match exists in raw output...
    assert any(row["mirna_id"] == "mir_nonseed" and int(row["nm"]) == 0 for row in raw)
    assert any(row["mirna_id"] == "mir_seed_exact" and int(row["coord"]) == 9 for row in raw)
    # ...but never survives into the filtered (seed-region) hits.
    assert all(int(row["coord"]) == 1 for row in filtered)
    assert not any(row["mirna_id"] == "mir_nonseed" for row in filtered)


@pytest.mark.parametrize("backend", [MiRNASeedBackend.EXHAUSTIVE_PYTHON, MiRNASeedBackend.PYAHOCORASICK])
def test_seed_position_two_mismatch_is_counted_not_treated_as_perfect(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    backend: MiRNASeedBackend,
) -> None:
    """A mismatch at guide position 2 (seed-window position 1) must count as a seed mismatch.

    Regression guard for the seed-window -> guide coordinate mapping: ``cand_pos2`` differs
    from the miRNA seed only at guide position 2. Its seed-region alignment must report
    seed_mismatches == nm == 1 (not 0), so it is never mislabeled as a perfect seed match.
    """
    output_dir = _run_toy_mirna_seed_analysis(
        tmp_path, monkeypatch, backend=backend, candidate_id=f"pos2_{backend.value}"
    )
    filtered = _read_tsv_rows(output_dir / f"pos2_{backend.value}_mirna_analysis.tsv")

    seed_hit = next(row for row in filtered if row["qname"] == "cand_pos2" and row["mirna_id"] == "mir_seed_exact")
    assert int(seed_hit["coord"]) == 1
    assert int(seed_hit["nm"]) == 1
    assert int(seed_hit["seed_mismatches"]) == 1, "position-2 mismatch must be counted in the seed region"
    assert float(seed_hit["offtarget_score"]) > 0.0, "a seed mismatch must not score as a perfect match"


def test_run_mirna_seed_analysis_pyahocorasick_matches_exhaustive_semantic_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The optimized analysis path should match the exhaustive oracle on semantic hit identity."""
    exhaustive_output_dir = _run_toy_mirna_seed_analysis(
        tmp_path,
        monkeypatch,
        backend=MiRNASeedBackend.EXHAUSTIVE_PYTHON,
        candidate_id="toy_parity_exhaustive",
    )
    pyahocorasick_output_dir = _run_toy_mirna_seed_analysis(
        tmp_path,
        monkeypatch,
        backend=MiRNASeedBackend.PYAHOCORASICK,
        candidate_id="toy_parity_pyaho",
    )

    output_pairs = [
        ("mirna_analysis.tsv", _read_tsv_rows, TOY_PARITY_SEMANTIC_IDENTITIES),
        ("mirna_analysis_raw.tsv", _read_tsv_rows, TOY_PARITY_RAW_SEMANTIC_IDENTITIES),
        ("mirna_hits.json", _read_json_rows, TOY_PARITY_SEMANTIC_IDENTITIES),
        ("mirna_hits_raw.json", _read_json_rows, TOY_PARITY_RAW_SEMANTIC_IDENTITIES),
    ]
    for suffix, reader, expected_identities in output_pairs:
        exhaustive_rows = reader(exhaustive_output_dir / f"toy_parity_exhaustive_{suffix}")
        pyahocorasick_rows = reader(pyahocorasick_output_dir / f"toy_parity_pyaho_{suffix}")
        assert _semantic_identities(pyahocorasick_rows) == expected_identities
        assert _semantic_identities(pyahocorasick_rows) == _semantic_identities(exhaustive_rows)

    exhaustive_summary = json.loads((exhaustive_output_dir / "toy_parity_exhaustive_mirna_summary.json").read_text())
    pyahocorasick_summary = json.loads((pyahocorasick_output_dir / "toy_parity_pyaho_mirna_summary.json").read_text())
    assert pyahocorasick_summary["total_hits"] == exhaustive_summary["total_hits"]
    assert pyahocorasick_summary["total_raw_alignments"] == exhaustive_summary["total_raw_alignments"]
    assert pyahocorasick_summary["hits_per_species"] == exhaustive_summary["hits_per_species"]


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
    raw_analysis_path = output_dir / "toy_mirna_analysis_raw.tsv"
    hits_path = output_dir / "toy_mirna_hits.json"
    raw_hits_path = output_dir / "toy_mirna_hits_raw.json"
    summary_path = output_dir / "toy_mirna_summary.json"
    schema_columns = list(MiRNAAlignmentSchema.to_schema().columns.keys())

    assert analysis_path.exists()
    assert raw_analysis_path.exists()
    assert hits_path.exists()
    assert raw_hits_path.exists()
    assert summary_path.exists()

    analysis_df = pd.read_csv(analysis_path, sep="\t")
    validated_df = MiRNAAlignmentSchema.validate(analysis_df, lazy=True)
    assert not validated_df.empty
    assert analysis_df.columns.tolist() == schema_columns

    raw_analysis_df = pd.read_csv(raw_analysis_path, sep="\t")
    raw_validated_df = MiRNAAlignmentSchema.validate(raw_analysis_df, lazy=True)
    assert raw_analysis_df.columns.tolist() == schema_columns

    # Filtered = seed-region hits only (0-based coord == 1); raw keeps every alignment.
    # The toy database intentionally contains non-seed matches, so filtered ⊂ raw.
    assert (validated_df["coord"] == 1).all(), "filtered analysis must contain only seed-region hits"
    assert len(validated_df) < len(raw_validated_df), "toy fixture should exercise the seed-region filter"
    seed_region_raw = raw_validated_df[raw_validated_df["coord"] == 1].reset_index(drop=True)
    pd.testing.assert_frame_equal(validated_df.reset_index(drop=True), seed_region_raw)

    hits = json.loads(hits_path.read_text())
    raw_hits = json.loads(raw_hits_path.read_text())
    assert hits
    assert raw_hits
    assert set(hits[0]) == set(schema_columns)
    assert set(raw_hits[0]) == set(schema_columns)

    summary = json.loads(summary_path.read_text())
    assert summary["total_hits"] == len(validated_df)
    assert summary["total_raw_alignments"] == len(raw_validated_df)
    # hits_per_species reports RAW alignment counts (documented on MiRNASummary).
    assert summary["hits_per_species"]["human"] == len(raw_validated_df)


def test_normalize_mirna_seed_hit_normalizes_bwa_and_in_process_hits() -> None:
    """Semantic normalization should align backend-specific hit metadata onto one contract."""
    backend_hit = {
        "qname": "cand-1",
        "qseq": "AUGCAUGCA",
        "mirna_id": "hsa-miR-test",
        "coord": 3,
        "strand": "+",
        "cigar": "7M",
        "mapq": 255,
        "as_score": 11,
        "nm": 1,
        "mismatch_positions": [3],
        "seed_mismatches": 1,
        "offtarget_score": 15.0,
    }
    bwa_hit = {
        "qname": "cand-1",
        "qseq": "AUGCAUGCA",
        "rname": "hsa-miR-test",
        "coord": "hsa-miR-test:4",
        "strand": "+",
        "cigar": "7M",
        "mapq": 255,
        "as_score": "11",
        "nm": 1,
        "mismatch_positions": [3],
        "seed_mismatches": 1,
        "offtarget_score": 15.0,
    }

    normalized_backend = normalize_mirna_seed_hit(backend_hit)
    normalized_bwa = normalize_mirna_seed_hit(bwa_hit, coord_is_one_based=True)

    assert normalized_backend == normalized_bwa
    assert normalized_backend.qseq == "ATGCATGCA"
    assert normalized_backend.semantic_identity() == mirna_seed_hit_identity(backend_hit)
    assert normalized_bwa.semantic_identity() == mirna_seed_hit_identity(bwa_hit, coord_is_one_based=True)


def test_run_mirna_seed_analysis_missing_pyahocorasick_dependency_fails_clearly(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The entrypoint should fail clearly when the optional fast backend dependency is unavailable."""
    test_data_dir = Path(__file__).parent / "data"
    candidate_file = test_data_dir / "toy_candidates.fasta"
    mirna_db = test_data_dir / "toy_mirna_db.fasta"

    monkeypatch.setattr(
        "sirnaforge.data.mirna_manager.MiRNADatabaseManager.get_database",
        lambda *_args, **_kwargs: mirna_db,
    )

    real_import_module = importlib.import_module

    def _fake_import_module(name: str, package: str | None = None):
        if name == "ahocorasick":
            raise ImportError("missing test dependency")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", _fake_import_module)

    with pytest.raises(RuntimeError, match="backend 'pyahocorasick' is unavailable"):
        run_mirna_seed_analysis(
            candidates_file=candidate_file,
            candidate_id="toy",
            mirna_db="toy_db",
            mirna_species=["human"],
            output_dir=tmp_path,
            backend=MiRNASeedBackend.PYAHOCORASICK,
        )
