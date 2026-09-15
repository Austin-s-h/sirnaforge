"""``quilt_summarize.json`` registration for the single-file HTML report (issue #103).

Before this module existed, ``sirnaforge report`` wrote exactly one file and nothing in the repo
registered it with Quilt: a package built from a run directory rendered no report in its package view
at all. The schema asserted against here (``QUILT_SUMMARIZE_SCHEMA`` below) is quiltdata/quilt's own
``shared/schemas/quilt_summarize.json``, pinned inline rather than fetched at test time so this suite
runs offline -- the document it describes is a JSON *array* of ``file | [file, ...]`` rows, not the
dict shape ``mcp__quilt__generate_quilt_summarize_json`` hands back, and a test that only checked "is
this valid JSON" would pass on that wrong shape too.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.reporting import ReportPayload, quilt_summarize_entries, write_quilt_summarize

#: quiltdata/quilt's shared/schemas/quilt_summarize.json, fetched 2026-09-14 and pinned so this test
#: does not depend on the network or on that repo staying reachable.
QUILT_SUMMARIZE_SCHEMA: dict[str, Any] = {
    "definitions": {
        "fileShortcut": {
            "type": "string",
            "description": "Path relative to quilt_summarize.json",
            "examples": ["file1.json"],
        },
        "fileExtended": {
            "type": "object",
            "properties": {
                "path": {"$ref": "#/definitions/fileShortcut"},
                "description": {"type": "string"},
                "expand": {"description": "Whether preview is expanded by default or not", "type": "boolean"},
                "title": {"type": "string"},
                "types": {
                    "type": "array",
                    "items": {
                        "oneOf": [{"$ref": "#/definitions/typeShorthand"}, {"$ref": "#/definitions/typeExtended"}]
                    },
                },
                "width": {"$ref": "#/definitions/width"},
            },
            "required": ["path"],
        },
        "typeShorthand": {
            "type": "string",
            "enum": ["echarts", "html", "igv", "json", "jupyter", "perspective", "text", "vega", "voila"],
        },
        "typeExtended": {
            "type": "object",
            "properties": {
                "name": {"$ref": "#/definitions/typeShorthand"},
                "style": {
                    "type": "object",
                    "properties": {
                        "height": {
                            "description": "Height as an absolute value (in `px`, `vh`, `em` etc.)",
                            "examples": ["1000px", "500vh"],
                            "type": "string",
                        }
                    },
                },
                "config": {
                    "description": "Perspective specific option. Restores renderer to a state previously "
                    "returned by saving config",
                    "type": "object",
                },
                "settings": {"description": "Perspective specific option. Sets config opened", "type": "boolean"},
            },
            "required": ["name"],
        },
        "file": {"oneOf": [{"$ref": "#/definitions/fileShortcut"}, {"$ref": "#/definitions/fileExtended"}]},
        "width": {
            "anyOf": [
                {"description": "Ratio number for flex-based width", "examples": [1.5], "type": "number"},
                {"description": "Width in pixels or percents", "examples": ["100px"], "type": "string"},
            ]
        },
        "row": {"anyOf": [{"$ref": "#/definitions/file"}, {"type": "array", "items": {"$ref": "#/definitions/file"}}]},
    },
    "type": "array",
    "items": {"$ref": "#/definitions/row"},
}


def _payload(run: dict[str, Any] | None = None) -> ReportPayload:
    """A minimal payload -- only ``run`` matters to this module; ``guides``/``filters`` do not."""
    return ReportPayload(
        schema_version="1.0.0",
        run={
            "gene_query": "TP53",
            "guides": 4159,
            "candidate_rows": 29605,
            "status_counts": {"fail": 1, "unknown": 2, "warn": 3, "pass": 4153},
            "embed_scope": "human, nm<=2",
            "agreement": {"comparable": True, "contradicted_run_pass": 0, "overruled_run_fail": 0},
            **(run or {}),
        },
        filters=[],
        guides=[],
        provenance={},
        caveats=[],
    )


def _write(run: Path, *rel_paths: str, content: str = "x") -> None:
    for rel in rel_paths:
        p = run / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)


@pytest.mark.unit
def test_the_document_is_an_array_not_the_mcp_tools_dict(tmp_path: Path) -> None:
    """The real file is a JSON array of rows. A dict is the wrong shape entirely, not a near miss."""
    run = tmp_path / "run"
    _write(run, "sirnaforge/manifest.json")
    entries = quilt_summarize_entries(_payload(), run / "report.html", run)
    assert isinstance(entries, list)
    jsonschema.validate(entries, QUILT_SUMMARIZE_SCHEMA)


@pytest.mark.unit
def test_report_html_is_first_full_width_and_expanded(tmp_path: Path) -> None:
    """The report is the obvious primary: alone in its row, expanded, and typed ``html``."""
    run = tmp_path / "run"
    entries = quilt_summarize_entries(_payload(), run / "report.html", run)

    first = entries[0]
    assert isinstance(first, dict), "a lone entry is a bare file, not wrapped in a row array"
    assert first["path"] == "report.html"
    assert first["types"] == ["html"]
    assert first["expand"] is True


@pytest.mark.unit
def test_the_description_names_gene_guides_status_and_agreement(tmp_path: Path) -> None:
    """Title and description are built from the run the report already describes, not restated."""
    run = tmp_path / "run"
    entries = quilt_summarize_entries(_payload(), run / "report.html", run)
    first = entries[0]

    assert "TP53" in first["title"]
    description = first["description"]
    assert "4,159" in description or "4159" in description
    assert "1 fail" in description
    assert "human, nm<=2" in description


@pytest.mark.unit
def test_an_absent_artifact_is_omitted_not_registered_broken(tmp_path: Path) -> None:
    """Absence is not an error the module raises on -- the entry for it simply does not exist.

    A run built with ``--skip-off-targets`` leaves ``off_target/`` empty. Registering
    ``combined_offtargets.tsv`` anyway would give the catalog a path it can never render, which the
    catalog shows as a broken preview -- worse than no entry at all.
    """
    run = tmp_path / "run"
    _write(run, "sirnaforge/candidates_all.csv", "sirnaforge/manifest.json", "report.html")
    entries = quilt_summarize_entries(_payload(), run / "report.html", run)

    blob = json.dumps(entries)
    assert "combined_offtargets" not in blob
    assert "combined_mirna_hits" not in blob
    assert "canonical" not in blob
    assert "orf_validation" not in blob
    # And everything that *is* registered must actually resolve to a real file.
    out_dir = (run / "report.html").parent
    for path in _all_paths(entries):
        assert out_dir.joinpath(path).resolve().is_file(), f"{path} does not resolve to a real file"
    jsonschema.validate(entries, QUILT_SUMMARIZE_SCHEMA)


@pytest.mark.unit
def test_every_registered_path_resolves_to_a_real_file(tmp_path: Path) -> None:
    """Every path in the document is relative to where the document is written, and it exists there."""
    run = tmp_path / "run"
    _write(
        run,
        "sirnaforge/candidates_all.csv",
        "sirnaforge/candidates_pass.csv",
        "sirnaforge/candidates_pass.fasta",
        "sirnaforge/manifest.json",
        "logs/workflow_summary.json",
        "off_target/results/aggregated/combined_offtargets.tsv",
        "off_target/results/aggregated/combined_mirna_hits.tsv",
        "orf_reports/orf_validation.txt",
        "transcripts/TP53_transcripts.fasta",
        "transcripts/TP53_canonical.fasta",
        "sub/report.html",
    )
    report_path = run / "sub" / "report.html"
    entries = quilt_summarize_entries(_payload(), report_path, run)

    paths = _all_paths(entries)
    assert len(paths) >= 10
    for rel in paths:
        assert report_path.parent.joinpath(rel).resolve().is_file(), f"{rel} does not resolve to a real file"
    jsonschema.validate(entries, QUILT_SUMMARIZE_SCHEMA)


@pytest.mark.unit
def test_multi_file_groups_are_rendered_as_row_arrays(tmp_path: Path) -> None:
    """Two candidate CSVs share one row (a multi-column row is an array of file, per the schema)."""
    run = tmp_path / "run"
    _write(run, "sirnaforge/candidates_all.csv", "sirnaforge/candidates_pass.csv")
    entries = quilt_summarize_entries(_payload(), run / "report.html", run)

    csv_row = next(row for row in entries if isinstance(row, list) and any("candidates_" in f["path"] for f in row))
    assert len(csv_row) == 2
    assert {f["types"][0] for f in csv_row} == {"perspective"}


@pytest.mark.unit
def test_write_quilt_summarize_lands_beside_the_report(tmp_path: Path) -> None:
    """Quilt only reads a summarize file at the package root beside the object it names."""
    run = tmp_path / "run"
    report_path = run / "report.html"
    report_path.parent.mkdir(parents=True)
    report_path.write_text("<html></html>")

    written = write_quilt_summarize(_payload(), report_path, run)

    assert written == report_path.parent / "quilt_summarize.json"
    assert written.exists()
    on_disk = json.loads(written.read_text())
    assert on_disk == quilt_summarize_entries(_payload(), report_path, run)
    jsonschema.validate(on_disk, QUILT_SUMMARIZE_SCHEMA)


@pytest.mark.unit
def test_paths_are_relative_to_the_summarize_file_not_the_report(tmp_path: Path) -> None:
    """When the two directories differ, Quilt resolves against the summarize file. So must we.

    The workflow writes ``sirnaforge/report.html`` but publishes the run directory, so its summarize
    file belongs at the run root. Relativising against the report's directory instead put every row one
    level too deep -- ``candidates_all.csv`` for a file that is really ``sirnaforge/candidates_all.csv``
    -- which the catalog renders as a broken preview on every single row.
    """
    run = tmp_path / "run"
    report_path = run / "sirnaforge" / "report.html"
    _write(run, "sirnaforge/report.html", "sirnaforge/candidates_all.csv", "sirnaforge/manifest.json")

    written = write_quilt_summarize(_payload(), report_path, run, out_path=run / "quilt_summarize.json")

    assert written == run / "quilt_summarize.json"
    paths = _all_paths(json.loads(written.read_text()))
    assert "sirnaforge/report.html" in paths
    assert "sirnaforge/candidates_all.csv" in paths
    for rel in paths:
        assert (run / rel).exists(), f"{rel} does not resolve from the summarize file's own directory"
    jsonschema.validate(json.loads(written.read_text()), QUILT_SUMMARIZE_SCHEMA)


def _all_paths(entries: list[Any]) -> list[str]:
    paths: list[str] = []
    for row in entries:
        for entry in row if isinstance(row, list) else [row]:
            paths.append(entry["path"] if isinstance(entry, dict) else entry)
    return paths


def _fake_write_report(payload: ReportPayload, out: Path | str) -> Path:
    """Stand-in for :func:`sirnaforge.reporting.write_report`.

    Writes the one file it promises, nothing more, so these CLI tests exercise only the
    ``--quilt-summarize`` wiring around it.
    """
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("<html></html>")
    return out_path


@pytest.mark.unit
def test_cli_writes_quilt_summarize_by_default(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``sirnaforge report`` writes ``quilt_summarize.json`` beside the HTML unless told not to."""
    monkeypatch.setattr("sirnaforge.cli.build_payload", lambda _run_dir: _payload())
    monkeypatch.setattr("sirnaforge.cli.write_report", _fake_write_report)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    output = tmp_path / "out" / "report.html"

    result = CliRunner().invoke(app, ["report", str(run_dir), "-o", str(output)])

    assert result.exit_code == 0, result.output
    assert (output.parent / "quilt_summarize.json").exists()


@pytest.mark.unit
def test_cli_no_quilt_summarize_skips_it(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The opt-out flag actually opts out, rather than the option being decorative."""
    monkeypatch.setattr("sirnaforge.cli.build_payload", lambda _run_dir: _payload())
    monkeypatch.setattr("sirnaforge.cli.write_report", _fake_write_report)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    output = tmp_path / "out" / "report.html"

    result = CliRunner().invoke(app, ["report", str(run_dir), "-o", str(output), "--no-quilt-summarize"])

    assert result.exit_code == 0, result.output
    assert not (output.parent / "quilt_summarize.json").exists()


@pytest.mark.unit
def test_cli_recognises_the_quilt_summarize_option(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``--quilt-summarize`` is a real option, not the ``exit(2) no such option`` of an unpatched CLI."""
    monkeypatch.setattr("sirnaforge.cli.build_payload", lambda _run_dir: _payload())
    monkeypatch.setattr("sirnaforge.cli.write_report", _fake_write_report)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    output = tmp_path / "report.html"

    result = CliRunner().invoke(app, ["report", str(run_dir), "-o", str(output), "--quilt-summarize"])

    assert result.exit_code != 2, result.output
    assert result.exit_code == 0, result.output
