"""Assertions over the artifacts a real workflow run writes inside the image.

Not a test module (no ``test_`` prefix, so pytest does not collect it): shared helpers for the
container tier.

Why these exist. 0.7.1's two flagship subsystems -- #103's HTML report and #100's screening-evidence
contract -- execute end to end **only** here. Every host test of them either calls ``build_payload``
in process or stubs the screen (``tests/integration/test_offtarget_only_contract.py`` says so in its
docstring), and the ``.nf`` side is otherwise asserted by reading ``.nf`` files as text
(``tests/unit/test_evidence_plan_threading.py``). Yet no container test looked at ``report.html``,
``report_manifest.json`` or ``evidence.json`` -- and ``workflow.py`` wraps both the render and the
sidecar write in ``except Exception -> logger.warning``, so the CLI exits 0 either way. A total
failure of the report inside the image left the release stage green.

The artifacts are already written by runs the tier already performs, so pinning them costs no
runtime. Each helper asserts the machine-readable verdict the subsystem publishes about itself
rather than re-deriving it, because that verdict is the thing a consumer reads.
"""

import json
from pathlib import Path
from typing import Any

from sirnaforge.models.evidence import EvidenceStatus

#: Every status a screening unit is allowed to publish, read off the installed package rather than
#: retyped here so the container tier cannot drift from the vocabulary the image actually ships. A
#: value outside this set means producer and reader disagree, which is the #100 defect class (`hsa`
#: vs `human` made every miRNA unit reconcile FAILED on an ordinary run).
EVIDENCE_STATUSES = frozenset(status.value for status in EvidenceStatus)


def _load(path: Path, what: str) -> Any:
    assert path.is_file(), f"missing {what}: {path}"
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:  # pragma: no cover - a parse failure is the finding
        raise AssertionError(f"{what} at {path} is not valid JSON: {exc}") from None


def assert_report_rendered(output_dir: Path) -> dict[str, Any]:
    """The #103 report rendered, and its sidecar says so.

    ``report_manifest.json`` carries ``report_html.render_error``, which is exactly the field
    ``workflow.py`` fills in instead of failing the run. Asserting it is what turns a swallowed
    render exception from a green run into a red one.
    """
    manifest = _load(output_dir / "sirnaforge" / "report_manifest.json", "report sidecar")

    report_html = manifest.get("report_html", {})
    assert report_html.get("render_error") is None, f"report render failed: {report_html.get('render_error')}"
    assert report_html.get("exists") is True, f"sidecar reports no report.html: {report_html}"
    assert report_html.get("size_bytes", 0) > 0, f"report.html is empty: {report_html}"

    rendered = output_dir / "sirnaforge" / str(report_html.get("path", "report.html"))
    assert rendered.is_file(), f"sidecar names {rendered} but it is not on disk"

    return manifest


def assert_evidence_reconciled(output_dir: Path, *, expect_channels: set[str]) -> dict[str, Any]:
    """Every planned screening unit came back with a status, and the channels asked for are there.

    This is #100's contract in one assertion: the plan and the evidence are reconciled unit for
    unit. It deliberately does **not** require success. A unit may legitimately come back
    ``failed`` (an unreachable miRNA database for one species does exactly that, and #106 is about
    reporting such a run honestly rather than hiding it) -- what must never happen is a planned
    unit vanishing, or a status the reader's vocabulary does not contain.
    """
    envelope = _load(output_dir / "off_target" / "results" / "aggregated" / "evidence.json", "aggregated evidence")

    for key in ("plan", "evidence"):
        assert key in envelope, f"evidence envelope missing {key!r}: {sorted(envelope)}"

    def units(section: str) -> set[tuple[str, str]]:
        return {(e["channel"], e["species"]) for e in envelope[section]["entries"]}

    planned, observed = units("plan"), units("evidence")
    assert planned, "the run planned no screening units at all"
    assert planned == observed, (
        f"plan and evidence do not reconcile: planned-not-observed={sorted(planned - observed)}, "
        f"observed-not-planned={sorted(observed - planned)}"
    )

    channels = {channel for channel, _ in observed}
    assert expect_channels <= channels, f"expected channels {sorted(expect_channels)}, evidence has {sorted(channels)}"

    for entry in envelope["evidence"]["entries"]:
        status = entry.get("status")
        assert status in EVIDENCE_STATUSES, (
            f"unknown evidence status {status!r} for {entry.get('channel')}/{entry.get('species')}"
        )

    return envelope


def offtarget_status(output_dir: Path) -> str | None:
    """``offtarget_summary.status`` from the run's own summary.

    A test that claims to screen transcriptomes should assert this is ``completed``; one that only
    runs miRNA seeds should assert ``partial``, because that is what the run reports about itself.
    Reading it is how the tier tells 'screened clean' from 'never screened'.
    """
    summary = _load(output_dir / "logs" / "workflow_summary.json", "workflow summary")
    return summary.get("offtarget_summary", {}).get("status")
