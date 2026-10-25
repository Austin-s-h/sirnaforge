"""Actual workflow exports use the selected report mode and preserve sidecar registration."""

import json
from pathlib import Path

import pytest

from sirnaforge.reporting import ReportMode
from sirnaforge.workflow import WorkflowConfig
from tests.unit.test_manifest_provenance import _qualified_workflow, _run_step6

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("mode", ["compact_single_file", "paged_single_file", "data_backed"])
def test_workflow_generates_selected_report_and_quilt_registration(tmp_path: Path, mode: ReportMode) -> None:
    """Workflow generates selected report and quilt registration."""
    workflow = _qualified_workflow(tmp_path, mode)
    workflow.config.report_mode = mode
    workflow.config.quilt_summarize = True

    base = _run_step6(workflow)

    assert (base / "report.html").is_file()
    diagnostic = json.loads((base / "report.diagnostics.json").read_text())
    assert diagnostic["mode_selection"]["mode"] == mode
    assert (workflow.config.output_dir / "quilt_summarize.json").is_file()
    attestation = json.loads((base / "report_manifest.json").read_text())
    assert attestation["report_html"]["exists"] is True
    assert attestation["report_html"]["render_error"] is None
    if mode == "data_backed":
        assert list((base / "report.details").glob("guide-*.js"))
    else:
        assert not (base / "report.details").exists()


def test_invalid_report_mode_fails_before_output_creation(tmp_path: Path) -> None:
    """Invalid report mode fails before output creation."""
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="Unsupported report mode"):
        WorkflowConfig(output_dir=output, gene_query="test", report_mode="unrecognized")  # type: ignore[arg-type]
    assert not output.exists()
