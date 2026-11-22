"""Stable report evidence states preserve producer non-decisions across every reader."""

import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

from sirnaforge.config.run_policy import FILTER_SPECS
from sirnaforge.reporting import build_payload, render_html, write_report_from_run
from sirnaforge.reporting.evidence import (
    EVIDENCE_STATES,
    FILTER_UNITS,
    WARNING_REGISTRY,
    guide_mirna_evidence,
    mirna_evidence,
    normalize_evidence_state,
)
from sirnaforge.reporting.payload import reevaluate_gates
from tests.unit.test_reporting_payload import GUIDE, _execute_mirna_card, _write_fully_evidenced_run

pytestmark = pytest.mark.unit


def _run(tmp_path: Path, recorded: str | None) -> Path:
    run = _write_fully_evidenced_run(tmp_path)
    path = run / "sirnaforge/manifest.json"
    manifest = json.loads(path.read_text())
    manifest["reference_evidence"] = {
        "schema_version": 1,
        "mirna": {
            "source": "mirbase",
            "acquisition": {"status": "complete"},
            "screening": {"status": "complete", "guide_sequences": [GUIDE]},
        },
    }
    path.write_text(json.dumps(manifest))
    if recorded is not None:
        path = run / "sirnaforge/candidates_all.csv"
        frame = pd.read_csv(path)
        # The fallback value is zero and would pass: the explicit producer non-decision must win.
        frame["max_off_target_count_verdict"] = recorded
        frame.to_csv(path, index=False)
    return run


@pytest.mark.parametrize(
    ("recorded", "state", "warning_id"),
    [
        ("unavailable", "unavailable", "GATE_RUN_UNAVAILABLE"),
        ("not_established", "not_established", "GATE_RUN_NOT_ESTABLISHED"),
        ("unknown", "not_established", "GATE_RUN_UNKNOWN"),
        ("NOT-ESTABLISHED", "not_established", "GATE_RUN_NOT_ESTABLISHED"),
        ("not_evaluated", "not_established", "GATE_RUN_NOT_EVALUATED"),
        ("unrecognized_pending", "not_established", "GATE_RUN_NOT_ESTABLISHED"),
    ],
)
def test_explicit_nondecision_stays_unresolved_even_with_a_passing_value(
    tmp_path: Path, recorded: str, state: str, warning_id: str
) -> None:
    """No spelling of missing evidence can become zero/clean or be re-thresholded into a pass."""
    payload = build_payload(_run(tmp_path, recorded))
    guide = payload.guides[0]
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_off_target_count")
    assert payload.filters[index]["unevaluable_reason_id"] in WARNING_REGISTRY
    assert guide.status == "unknown"
    assert guide.gates[index][1] == 2
    assert guide.report_evidence["gates"][index]["state"] == state
    assert warning_id in guide.report_evidence["warning_ids"]
    assert payload.run["unresolved_guides"] == 1
    assert payload.run["passes_without_unresolved_evidence"] == 0
    assert reevaluate_gates(payload.filters, guide.gates, {"max_off_target_count": 1000})[index] == guide.gates[index]
    html = render_html(payload)
    card = _execute_mirna_card(html, "evidenceCard(G[0])")
    assert f'data-warning-id="{warning_id}"' in card
    assert WARNING_REGISTRY[warning_id]["remediation"] in card
    result = _execute_mirna_card(html, "JSON.stringify(PRESETS.unresolved.test(G[0],G[0]._live))")
    assert result == "true"


@pytest.mark.parametrize("mode", ["compact_single_file", "paged_single_file", "data_backed"])
def test_resolved_run_has_no_unresolved_card_and_keeps_policy_auditable(tmp_path: Path, mode: str) -> None:
    """A fully observed run presents cleanly while retaining the recorded source and scope."""
    run = _run(tmp_path, None)
    output = tmp_path / "report.html"
    result = write_report_from_run(run, output, mode=mode)
    if mode == "compact_single_file":
        summary = result.manifest["summary"]
        payload = build_payload(run)
        assert payload.run["unresolved_guides"] == 0
        assert payload.run["passes_without_unresolved_evidence"] == 1
        assert 'id="evidence-unresolved"' not in output.read_text()
        assert "Setting key" in _execute_mirna_card(output.read_text(), "thresholdMetadata(FILTERS[0])")
    else:
        summary = result.manifest["summary"]
        assert summary["unresolved_guides"] == 0
        assert summary["passes_without_unresolved_evidence"] == 1
        assert '<details id="evidence-notices">' not in output.read_text()
    assert summary["status_counts"]["pass"] == 1


@pytest.mark.parametrize("mode", ["paged_single_file", "data_backed"])
def test_staged_readers_keep_generic_gate_warning_identity(tmp_path: Path, mode: str) -> None:
    """The staged artifacts and reader inspector retain canonical states and registry IDs."""
    result = write_report_from_run(_run(tmp_path, "unavailable"), tmp_path / "report.html", mode=mode)
    assert result.artifact_dir is not None
    guides = pq.read_table(result.artifact_dir / "guides.parquet").to_pylist()
    diagnostic = json.loads(guides[0]["report_evidence"])
    assert "GATE_RUN_UNAVAILABLE" in diagnostic["warning_ids"]
    assert guides[0]["unresolved_evidence"] is True
    assert result.manifest["summary"]["evidence_warning_counts"]["GATE_RUN_UNAVAILABLE"] == 1
    assert result.manifest["summary"]["passes_without_unresolved_evidence"] == 0
    assert 'data-warning-id="GATE_RUN_UNAVAILABLE"' in result.output.read_text()


@pytest.mark.parametrize("value", [*EVIDENCE_STATES, "unknown", "NOT-EVALUATED", "pending", None])
def test_normalized_vocabulary_never_promotes_an_unknown(value: object) -> None:
    """Normalization accepts canonical values and conservative aliases, with no fallback pass."""
    normalized = normalize_evidence_state(value)
    assert normalized in EVIDENCE_STATES
    if value not in EVIDENCE_STATES:
        assert normalized != "pass"


def test_every_registered_gate_has_explicit_report_units() -> None:
    """Canonical policy settings must be inspectable without guessing their measurement unit."""
    assert {spec.setting_key for spec in FILTER_SPECS} <= set(FILTER_UNITS)
    assert FILTER_UNITS["gc_min"] == "% GC"
    assert FILTER_UNITS["max_repeat_transcript_fraction"] == "fraction (0–1)"
    assert FILTER_UNITS["max_transcript_seed_sites"] == "sites"


@pytest.mark.parametrize(
    ("acquisition", "expected"),
    [
        (" UNAVAILABLE ", "acquisition_failed"),
        ("FAIL", "acquisition_failed"),
        ("not-established", "unknown"),
        (None, "unknown"),
    ],
)
def test_incomplete_reference_phase_cannot_hide_behind_a_completed_screen(
    acquisition: str | None, expected: str
) -> None:
    """Source availability aliases and absent provenance cannot produce a completed-zero message."""
    contract = mirna_evidence(
        {
            "reference_evidence": {
                "schema_version": 1,
                "mirna": {
                    "acquisition": {"status": acquisition},
                    "screening": {"status": "complete", "guide_sequences": [GUIDE]},
                },
            }
        }
    )
    assert contract["state"] == expected
    diagnostic = guide_mirna_evidence(contract, GUIDE)
    assert diagnostic["warning_id"] is not None


@pytest.mark.parametrize("recorded", ["pass", "warn", "fail"])
def test_recorded_outcome_conflict_is_quarantined_instead_of_reinterpreted(tmp_path: Path, recorded: str) -> None:
    """Contradictory observations cannot silently erase a producer warning/failure or invent a pass."""
    run = _run(tmp_path, recorded)
    if recorded == "pass":
        path = run / "sirnaforge/candidates_all.csv"
        frame = pd.read_csv(path)
        frame["max_off_target_count_observed"] = 1000
        frame.to_csv(path, index=False)
    payload = build_payload(run)
    guide = payload.guides[0]
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_off_target_count")
    assert guide.gates[index][1] == 2
    assert "GATE_RECORDED_CONFLICT" in guide.report_evidence["warning_ids"]
    assert guide.report_evidence["gates"][index]["recorded_state_raw"] == recorded
    assert payload.run["passes_without_unresolved_evidence"] == 0
    assert reevaluate_gates(payload.filters, guide.gates, {"max_off_target_count": 0})[index] == guide.gates[index]
    card = _execute_mirna_card(render_html(payload), "evidenceCard(G[0])")
    assert 'data-warning-id="GATE_RECORDED_CONFLICT"' in card
    assert f"Recorded state: {recorded}" in card
