"""Contracts for versioned, projected report artifacts."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pyarrow as pa
import pytest

from sirnaforge.reporting.artifacts import (
    ReportArtifactStore,
    ReportModeFacts,
    ReportModeThresholds,
    resolve_report_mode,
)
from tests.unit.test_reporting_payload import GUIDE, OTHER, _candidate_row, _hit_row, _write_run


@pytest.mark.unit
def test_auto_mode_uses_configurable_budgets_and_explicit_modes_win() -> None:
    """Mode selection is deterministic from source measurements and does not depend on biology."""
    thresholds = ReportModeThresholds(
        compact_max_input_bytes=100,
        compact_max_candidate_rows=10,
        compact_max_guides=5,
        paged_max_input_bytes=1_000,
        paged_max_candidate_rows=100,
        paged_max_guides=50,
    )
    assert resolve_report_mode("auto", ReportModeFacts(100, 10, 5), thresholds) == "compact_single_file"
    assert resolve_report_mode("auto", ReportModeFacts(101, 10, 5), thresholds) == "paged_single_file"
    assert resolve_report_mode("auto", ReportModeFacts(1_001, 10, 5), thresholds) == "data_backed"
    assert resolve_report_mode("data_backed", ReportModeFacts(1, 1, 1), thresholds) == "data_backed"


@pytest.mark.unit
def test_parquet_store_projects_and_streams_source_tables_and_traces_exports(tmp_path: Path) -> None:
    """The store contains queryable report columns, hashes, and export lineage without source-table tails."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST1", 10), _candidate_row("c2", OTHER, "ENST2", 20)],
        [_hit_row(GUIDE, "human", "ENST9", 1, "off_target", "GENE9")],
    )
    candidate_path = run / "sirnaforge/candidates_all.csv"
    candidate_path.write_text(
        candidate_path.read_text().replace("\n", ",unused_column\n", 1).replace(",PASS\n", ",PASS,not-for-report\n")
    )
    store = ReportArtifactStore.from_run(run, tmp_path / "report.artifacts")

    manifest = json.loads((store.root / "report_manifest.json").read_text())
    assert manifest["artifact_schema_version"] == 1
    assert manifest["report_mode"] == "compact_single_file"
    assert manifest["mode_facts"]["distinct_guides"] == 2
    assert "source_run_dir" not in manifest
    assert not list(store.root.glob("*.guides.sqlite"))
    assert manifest["tables"]["candidates"]["source_sha256"]
    assert "unused_column" not in manifest["tables"]["candidates"]["columns"]
    batches = list(store.iter_batches("candidates", columns=["guide_sequence"], batch_size=1))
    assert len(batches) == 2
    assert all(batch.schema.field(0).type == pa.string() for batch in batches)
    assert sum(batch.num_rows for batch in batches) == 2

    exported = store.export_table("candidates")
    with exported.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    updated = json.loads((store.root / "report_manifest.json").read_text())
    assert updated["exports"]["candidates"]["artifact_table_sha256"] == manifest["tables"]["candidates"]["sha256"]
    assert updated["exports"]["candidates"]["sha256"]
