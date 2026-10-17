"""Contracts for measured review datasets, independent of efficacy model fitting."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.build_review_benchmarks import build_dataset, normalize_record


@pytest.mark.unit
def test_passenger_core_is_oriented_without_inventing_full_strands() -> None:
    """A published sense core yields a guide core, not an invented full duplex."""
    source = {
        "collection_id": "katoh_2007",
        "source_record_id": "position1",
        "sequence": "AUGGUGAGCAAGGGCGAGG",
        "orientation": "sense",
        "representation": "paired_core_only",
        "paired_length": 19,
        "label_kind": "measured",
        "measurement_raw": "54.75",
        "endpoint": "normalized_egfp_knockdown",
        "units": "percent",
        "direction": "higher_is_better",
        "chemistry": "unmodified_core",
    }
    row = normalize_record(source)
    assert row["paired_guide_sequence"] == "CCTCGCCCTTGCTCACCAT"
    assert row["full_guide_sequence"] is None
    assert row["full_passenger_sequence"] is None
    assert row["measurement_raw"] == "54.75"
    assert row["measurement_value"] == "54.75"
    assert row["compatibility_status"] == "compatible_core_only"
    assert row["geometry_evidence"] == "declared_paired_core_only"
    assert row["target_identity_status"] == "not_established"


@pytest.mark.unit
def test_full_duplex_preserves_overhangs_censoring_and_chemistry() -> None:
    """Censored potency and modified strands are not unmodified inhibition labels."""
    guide = "UGGUCCAGAGGCCGCUGGUTT"
    passenger = "ACCAGCGGCCUCUGGACCATT"
    source = {
        "collection_id": "kliuchnikov_2024",
        "source_record_id": "siSER-1m",
        "representation": "full_duplex",
        "guide_sequence": guide,
        "passenger_sequence": passenger,
        "guide_sequence_original": "uGfgucCfagaggccGfcUfggudTdT",
        "passenger_sequence_original": "accaGfcGfGfCfcucuggaccadTdT",
        "paired_length": 19,
        "label_kind": "measured",
        "measurement_raw": ">100",
        "endpoint": "IC50",
        "units": "nM",
        "direction": "lower_is_better",
        "chemistry": "2OMe_2F_modified_core_dTdT_overhangs",
    }
    row = normalize_record(source)
    assert row["full_guide_sequence"] == guide.replace("U", "T")
    assert row["guide_overhang"] == row["passenger_overhang"] == "TT"
    assert row["paired_guide_sequence"] == "TGGTCCAGAGGCCGCTGGT"
    assert row["measurement_raw"] == ">100"
    assert (row["measurement_value"], row["measurement_relation"]) == ("100", ">")
    assert row["endpoint"] == "IC50"
    assert row["chemistry"] == source["chemistry"]
    assert row["guide_sequence_original"] == source["guide_sequence_original"]
    assert row["compatibility_status"] == "compatible_full_duplex"


@pytest.mark.unit
@pytest.mark.parametrize(
    ("updates", "status"),
    [
        ({"paired_length": 15, "sequence": "AUGGUGAGCAAGGGC"}, "incompatible_length"),
        ({"label_kind": "predicted"}, "excluded_not_measured"),
        ({"measurement_raw": "NaN"}, "excluded_invalid_measurement"),
    ],
)
def test_incompatible_or_unevidenced_examples_are_not_accepted(updates: dict, status: str) -> None:
    """A canonical sequence is not enough to qualify a measured training example."""
    source = {
        "collection_id": "control",
        "source_record_id": "1",
        "sequence": "AUGGUGAGCAAGGGCGAGG",
        "orientation": "guide",
        "representation": "paired_core_only",
        "paired_length": 19,
        "label_kind": "measured",
        "measurement_raw": "12.5",
        "endpoint": "residual_expression",
        "units": "percent",
        "direction": "lower_is_better",
        "chemistry": "not_stated",
        **updates,
    }
    assert normalize_record(source)["compatibility_status"] == status


@pytest.mark.unit
def test_builder_preserves_dose_rows_and_reports_exclusions(tmp_path: Path) -> None:
    """Same guide at two doses remains two assays; an unpaired record stays out."""
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source = {
        "collection_id": "test",
        "sequence": "AUGGUGAGCAAGGGCGAGG",
        "orientation": "guide",
        "representation": "paired_core_only",
        "paired_length": 19,
        "label_kind": "measured",
        "measurement_raw": "12.5",
        "endpoint": "residual_expression",
        "units": "percent",
        "direction": "lower_is_better",
        "chemistry": "not_stated",
    }
    rows = [
        {**source, "source_record_id": "1:3nM", "concentration": "3nM"},
        {**source, "source_record_id": "1:30nM", "concentration": "30nM"},
        {**source, "source_record_id": "bad", "label_kind": "predicted"},
    ]
    (source_dir / "extract.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
    (source_dir / "collections.json").write_text(
        json.dumps({"collections": [{"id": "test", "extract": "extract.jsonl", "expected_rows": len(rows)}]})
    )
    manifest = build_dataset(source_dir, tmp_path / "built")
    assert manifest["counts"] == {"source_rows": 3, "compatible_rows": 2, "excluded_rows": 1, "sequence_clusters": 1}
    written = [json.loads(line) for line in (tmp_path / "built" / "observations.jsonl").read_text().splitlines()]
    assert [row["concentration"] for row in written if row["compatibility_status"].startswith("compatible_")] == [
        "30nM",
        "3nM",
    ]
    assert len({row["record_id"] for row in written}) == len(written)
    assert len({row["sequence_cluster"] for row in written}) == 1
    build_dataset(source_dir, tmp_path / "again")
    assert (tmp_path / "built" / "records.csv").read_bytes() == (tmp_path / "again" / "records.csv").read_bytes()


@pytest.mark.unit
def test_vendored_review_sources_rebuild_with_exact_label_and_scope_accounting(tmp_path: Path) -> None:
    """Every source observation survives audit, and only compatible records enter the CSV."""
    root = Path(__file__).resolve().parents[2]
    data = root / "tests/data/benchmarks/ml_reviews"
    catalog = json.loads((data / "collections.json").read_text())
    manifest = build_dataset(data, tmp_path, root / "tests/data/benchmarks/oligogym/records.csv")
    assert manifest["counts"]["source_rows"] == sum(entry["expected_rows"] for entry in catalog["collections"])
    assert (
        manifest["counts"]["source_rows"] == manifest["counts"]["compatible_rows"] + manifest["counts"]["excluded_rows"]
    )
    originals = {
        f"{entry['id']}:{row['source_record_id']}": row
        for entry in catalog["collections"]
        for row in (json.loads(line) for line in (data / entry["extract"]).read_text().splitlines())
    }
    for line in (tmp_path / "observations.jsonl").read_text().splitlines():
        row = json.loads(line)
        original = originals[row["record_id"]]
        assert row["measurement_raw"] == original["measurement_raw"]
        assert row["endpoint"] == original["endpoint"]
        assert row["chemistry"] == original["chemistry"]
        assert row["target_identity_status"] == "not_established"
        if row["compatibility_status"].startswith("compatible_"):
            assert 19 <= row["paired_length"] <= 23
    assert (tmp_path / "records.csv").read_bytes() == (data / "built/records.csv").read_bytes()


@pytest.mark.unit
def test_database_katoh_cores_do_not_invent_overhang_chemistry() -> None:
    """A secondary database's core-only records cannot establish overhang chemistry."""
    root = Path(__file__).resolve().parents[2]
    data = root / "tests/data/benchmarks/ml_reviews"
    sources = [json.loads(line) for line in (data / "source_extracts/sirnaefficacydb.jsonl").read_text().splitlines()]
    katoh = [row for row in sources if row.get("origin_study") == "Takayuki Katoh"]
    assert katoh
    for source in katoh:
        row = normalize_record(source)
        assert row["chemistry"] == "not_stated_in_secondary_source"
        assert row["representation"] == "paired_core_only"
        assert row["full_guide_sequence"] is None
        assert row["full_passenger_sequence"] is None
        assert row["guide_overhang"] is None
        assert row["passenger_overhang"] is None
