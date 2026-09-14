"""Register the report and its evidence in ``quilt_summarize.json`` (issue #103).

Quilt's package view renders exactly what a package's own ``quilt_summarize.json`` names, and nothing
else -- it is a JSON *array* of ``file | [file, ...]`` rows (schema:
``quiltdata/quilt``'s ``shared/schemas/quilt_summarize.json``), not the dict some tools hand back. This
module builds that array from the same :class:`ReportPayload` the report itself renders, so a package
viewer never sees a claim the run didn't make.

Every path is written relative to where ``quilt_summarize.json`` itself lands, and an artifact that
does not exist under ``run_dir`` is omitted rather than registered and left dangling: a run built with
``--skip-off-targets`` leaves ``off_target/`` empty, and ``<GENE>_canonical.fasta`` exists only on the
gene-search path. Registering either unconditionally would give the catalog a path it can never render,
which reads as broken evidence rather than absent evidence.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from sirnaforge.reporting.payload import ReportPayload

#: File extension the schema calls out by name for the two big candidate tables and the two aggregated
#: hit tables -- a Perspective datagrid with filter, plot and pivot, loaded client-side by the catalog.
_GRID = ["perspective"]
_JSON = ["json"]
_TEXT = ["text"]


def _rel(path: Path | None, out_dir: Path) -> str | None:
    """``path`` relative to ``out_dir``, or ``None`` when ``path`` is unset or does not exist.

    Absence is the normal case for several of these -- not an error the caller should raise on.
    """
    if path is None or not path.exists():
        return None
    return os.path.relpath(path, out_dir)


def _file(
    path: str, *, title: str, types: list[Any], description: str | None = None, expand: bool | None = None
) -> dict[str, Any]:
    """One ``fileExtended`` entry.

    ``expand`` is omitted rather than set ``False``, matching the schema's own default, so a diff
    against a hand-written summarize file shows only real changes.
    """
    entry: dict[str, Any] = {"path": path, "title": title, "types": types}
    if description is not None:
        entry["description"] = description
    if expand is not None:
        entry["expand"] = expand
    return entry


def _row(*entries: dict[str, Any] | None) -> dict[str, Any] | list[dict[str, Any]] | None:
    """A summarize-file row: the lone survivor, several side by side, or nothing to show at all.

    A ``row`` is a ``file`` or an *array* of ``file`` -- never an empty array, which the schema does
    not forbid but no catalog build should ever emit, so a group that lost every member is dropped.
    """
    present = [e for e in entries if e is not None]
    if not present:
        return None
    return present[0] if len(present) == 1 else present


def _description(payload: ReportPayload) -> str:
    """The one line a reader sees before opening the report.

    Gene, guide count, verdict tally, embedded off-target scope and report/pipeline agreement, all
    read from the run the report already describes rather than restated.
    """
    run = payload.run
    counts = run.get("status_counts") or {}
    agreement = run.get("agreement") or {}
    tally = (
        f"{counts.get('pass', 0)} pass, {counts.get('warn', 0)} warn, "
        f"{counts.get('unknown', 0)} unknown, {counts.get('fail', 0)} fail"
    )
    if agreement.get("comparable"):
        agree = (
            f"; agrees with the run's own verdicts "
            f"({agreement.get('contradicted_run_pass', 0)} contradicted, "
            f"{agreement.get('overruled_run_fail', 0)} overruled -- both should read 0)"
        )
    else:
        agree = f"; {agreement.get('reason', 'agreement not comparable')}"
    return (
        f"{run.get('guides', 0):,} guides for {run.get('gene_query', 'unknown')}: {tally}. "
        f"Off-target evidence embedded at {run.get('embed_scope', 'unknown')}{agree}."
    )


def quilt_summarize_entries(payload: ReportPayload, report_path: Path | str, run_dir: Path | str) -> list[Any]:
    """Build the ``quilt_summarize.json`` array for one run, without writing it.

    Args:
        payload: The payload the report itself was rendered from.
        report_path: Where ``report.html`` is (or will be) written. ``quilt_summarize.json`` is always
            written beside it, so every path here is relative to ``report_path``'s directory.
        run_dir: The run directory the payload was built from -- may differ from ``report_path``'s
            directory when a reader points ``-o`` elsewhere.

    Returns:
        A JSON-serialisable list matching quiltdata/quilt's ``quilt_summarize.json`` schema.
    """
    out_dir = Path(report_path).parent
    run_dir = Path(run_dir)
    gene = str(payload.run.get("gene_query") or "run")
    sirnaforge_dir = run_dir / "sirnaforge"
    aggregated = run_dir / "off_target" / "results" / "aggregated"

    rows: list[Any] = [
        _file(
            os.path.relpath(Path(report_path), out_dir),
            title=f"siRNAforge report — {gene}",
            types=["html"],
            description=_description(payload),
            expand=True,
        )
    ]

    pass_csv = _rel(sirnaforge_dir / "candidates_pass.csv", out_dir)
    all_csv = _rel(sirnaforge_dir / "candidates_all.csv", out_dir)
    rows.append(
        _row(
            _file(pass_csv, title="Passing candidates", types=_GRID) if pass_csv else None,
            _file(all_csv, title="All candidates", types=_GRID) if all_csv else None,
        )
    )

    manifest = _rel(sirnaforge_dir / "manifest.json", out_dir)
    workflow_summary = _rel(run_dir / "logs" / "workflow_summary.json", out_dir)
    rows.append(
        _row(
            _file(manifest, title="Run manifest", types=_JSON) if manifest else None,
            _file(workflow_summary, title="Workflow summary", types=_JSON) if workflow_summary else None,
        )
    )

    offtargets = _rel(aggregated / "combined_offtargets.tsv", out_dir)
    mirna_hits = _rel(aggregated / "combined_mirna_hits.tsv", out_dir)
    rows.append(
        _row(
            _file(offtargets, title="Off-target hits (aggregated)", types=_GRID) if offtargets else None,
            _file(mirna_hits, title="miRNA seed hits (aggregated)", types=_GRID) if mirna_hits else None,
        )
    )

    orf_report = _rel(run_dir / "orf_reports" / "orf_validation.txt", out_dir)
    pass_fasta = _rel(sirnaforge_dir / "candidates_pass.fasta", out_dir)
    # Named by the gene query at write time, so glob rather than reconstruct the filename -- resolving
    # it the same way the workflow named it would duplicate that logic and could drift from it.
    transcripts_fasta = _rel(next(iter(sorted(run_dir.glob("transcripts/*_transcripts.fasta"))), None), out_dir)
    canonical_fasta = _rel(next(iter(sorted(run_dir.glob("transcripts/*_canonical.fasta"))), None), out_dir)
    rows.append(
        _row(
            _file(orf_report, title="ORF validation", types=_TEXT) if orf_report else None,
            _file(pass_fasta, title="Passing candidates (FASTA)", types=_TEXT) if pass_fasta else None,
            _file(transcripts_fasta, title=f"{gene} transcripts (FASTA)", types=_TEXT) if transcripts_fasta else None,
            _file(canonical_fasta, title=f"{gene} canonical transcript (FASTA)", types=_TEXT)
            if canonical_fasta
            else None,
        )
    )

    return [r for r in rows if r is not None]


def write_quilt_summarize(
    payload: ReportPayload,
    report_path: Path | str,
    run_dir: Path | str,
    out_path: Path | str | None = None,
) -> Path:
    """Write ``quilt_summarize.json``, returning the path written.

    Defaults to ``report_path``'s own directory -- Quilt only reads a summarize file it finds at the
    package root, so co-locating it with ``report.html`` is not a convenience, it is the contract.
    """
    entries = quilt_summarize_entries(payload, report_path, run_dir)
    out = Path(out_path) if out_path is not None else Path(report_path).parent / "quilt_summarize.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(entries, indent=2) + "\n", encoding="utf-8")
    return out
