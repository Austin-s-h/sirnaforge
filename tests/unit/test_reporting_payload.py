"""Contract tests for the single-file HTML report (issue #103).

Two of these pin the failure modes the issue names as most likely to ship a wrong report: the
evidence join keyed on the candidate id rather than the guide sequence, and non-liability alignments
counted as off-targets. The rest pin that an unevaluable gate stays unknown rather than becoming a
pass, which is the mistake the report exists to make visible.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS, unclassified_cells
from sirnaforge.reporting import ReportInputError, build_payload, render_html
from sirnaforge.reporting.payload import EMBED_MAX_NM

GUIDE = "ACGUACGUACGUACGUACGUA"
OTHER = "UUUUCCCCAAAAGGGGUUUUC"

_CANDIDATE_COLUMNS = (
    "id,guide_sequence,passenger_sequence,transcript_id,position,gc_content,asymmetry_score,"
    "empirical_score,paired_fraction,design_score,composite_score,weight_vector,passes_filters"
)


def _candidate_row(cid: str, guide: str, transcript: str, position: int, *, gc: float = 45.0) -> str:
    dna = guide.replace("U", "T")
    return f"{cid},{dna},{dna},{transcript},{position},{gc},0.8,0.5,0.3,55.0,60.0,postscreen_sirna_v4,PASS"


def _hit_row(guide: str, species: str, rname: str, nm: int, hit_class: str, symbol: str) -> dict[str, object]:
    cells = unclassified_cells()
    cells.update(hit_class=hit_class, hit_symbol=symbol, matched_symbol=symbol)
    return {
        "qname": "screen_0",
        "qseq": guide.replace("U", "T"),
        "species": species,
        "rname": rname,
        "coord": 100,
        "strand": "+",
        "cigar": "21M",
        "mapq": 60,
        "as_score": 42,
        "nm": nm,
        "seed_mismatches": 0,
        "offtarget_score": 1.0,
        **cells,
    }


def _write_run(tmp_path: Path, candidates: list[str], hits: list[dict[str, object]]) -> Path:
    run = tmp_path / "run"
    (run / "sirnaforge").mkdir(parents=True)
    (run / "off_target" / "results" / "aggregated").mkdir(parents=True)
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([_CANDIDATE_COLUMNS, *candidates]) + "\n")
    (run / "sirnaforge" / "manifest.json").write_text(json.dumps({"gene_query": "TP53", "tool_version": "test"}))

    table = run / "off_target" / "results" / "aggregated" / "combined_offtargets.tsv"
    if hits:
        header = list(hits[0].keys())
        lines = ["\t".join(header)] + ["\t".join(str(h[c]) for c in header) for h in hits]
        table.write_text("\n".join(lines) + "\n")
    else:
        base = [
            "qname",
            "qseq",
            "species",
            "rname",
            "coord",
            "strand",
            "cigar",
            "mapq",
            "as_score",
            "nm",
            "seed_mismatches",
            "offtarget_score",
        ]
        table.write_text("\t".join([*base, *CLASSIFICATION_COLUMNS]) + "\n")
    return run


@pytest.mark.unit
def test_evidence_is_keyed_by_guide_sequence_not_candidate_id(tmp_path: Path) -> None:
    """One guide on three transcripts shows identical evidence on all three.

    Joining on ``id`` would leave two of the three reading as zero-off-target. Screening runs once per
    distinct guide, so the guide sequence is the only key that can be right.
    """
    run = _write_run(
        tmp_path,
        [_candidate_row(f"c{i}", GUIDE, f"ENST0000000000{i}", 10 * i) for i in (1, 2, 3)],
        [_hit_row(GUIDE, "human", "ENST00000099999", 1, "off_target", "SOMEGENE")],
    )
    payload = build_payload(run)

    assert len(payload.guides) == 1, "three candidate rows are one guide"
    entry = payload.guides[0]
    assert entry.n_rows == 3
    assert len(entry.isoforms) == 3
    assert entry.liability_count == 1, "the single alignment reaches the guide, not just one of its rows"


@pytest.mark.unit
def test_on_target_and_ortholog_alignments_are_never_counted_as_off_targets(tmp_path: Path) -> None:
    """28-50% of alignments on real runs are not liabilities; counting them would be the core defect."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST00000000001", 10)],
        [
            _hit_row(GUIDE, "human", "ENST00000000001", 0, "on_target", "TP53"),
            _hit_row(GUIDE, "mouse", "ENSMUST00000000002", 1, "ortholog", "Trp53"),
            _hit_row(GUIDE, "human", "ENST00000077777", 2, "repeat", "REPEATY"),
            _hit_row(GUIDE, "human", "ENST00000099999", 2, "off_target", "REALLIABILITY"),
        ],
    )
    entry = build_payload(run).guides[0]

    assert entry.liability_count == 1, "only the off_target row is a liability"
    liabilities = {r["s"] for r in entry.offtarget_rows}
    assert liabilities == {"REALLIABILITY"}
    displayed = {s["symbol"] for s in entry.offtarget_by_symbol}
    assert {"TP53", "Trp53", "REPEATY"} <= displayed, "non-liabilities are still displayed, just not counted"
    assert sum(1 for s in entry.offtarget_by_symbol if s["is_liability"]) == 1


@pytest.mark.unit
def test_a_gate_whose_column_the_run_does_not_export_is_unknown_not_pass(tmp_path: Path) -> None:
    """An unknown must never read as clean.

    Six of the twelve active 0.7.1 gates read human-stratified counters the run does not export yet
    (#101). On the public baseline that leaves 450 guides the run failed but this cannot re-derive --
    they are reported ``not established``, never flipped to pass.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    entry = build_payload(run).guides[0]

    assert entry.n_gates_unknown > 0, "the synthetic run exports none of the screening counters"
    assert entry.status == "unknown", "a guide with unevaluable gates is not a pass"
    assert entry.n_gates_failed == 0, "and it is not a failure either"


@pytest.mark.unit
def test_a_hit_table_with_no_classification_is_refused(tmp_path: Path) -> None:
    """The report will not classify hits itself; it would then be able to disagree with the run."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    table = run / "off_target" / "results" / "aggregated" / "combined_offtargets.tsv"
    table.write_text("qname\tqseq\tspecies\trname\tnm\n")

    with pytest.raises(ReportInputError, match="hit_class"):
        build_payload(run)


@pytest.mark.unit
def test_out_of_scope_evidence_is_reported_rather_than_rendered_as_clean(tmp_path: Path) -> None:
    """A guide whose only hits lie outside the embedded scope must not look like a clean one."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST00000000001", 10)],
        [_hit_row(GUIDE, "human", "ENST00000099999", EMBED_MAX_NM + 3, "off_target", "FARGENE")],
    )
    entry = build_payload(run).guides[0]

    assert entry.offtarget_rows == [], "the hit is outside the embedded scope"
    assert entry.offtarget_embedded_scope_empty_but_counts_exist is True
    assert entry.offtarget_matrix, "but its counts are carried, completely"


@pytest.mark.unit
def test_a_register_neighbour_is_flagged_rather_than_looking_like_a_duplicate(tmp_path: Path) -> None:
    """Two designs one nucleotide apart share a window and can score identically."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST00000000001", 100), _candidate_row("c2", OTHER, "ENST00000000001", 101)],
        [],
    )
    payload = build_payload(run)
    flagged = [i for g in payload.guides for i in g.isoforms if i["register_neighbours"]]

    assert flagged, "a design 1 nt away on the same transcript is a register neighbour"


@pytest.mark.unit
def test_the_report_contradicts_no_run_pass_verdict(tmp_path: Path) -> None:
    """The agreement check is published so a reader can see the gate panel matches the pipeline."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    agreement = build_payload(run).run["agreement"]

    assert agreement["comparable"] is True
    assert agreement["contradicted_run_pass"] == 0


@pytest.mark.unit
def test_the_rendered_report_is_one_file_with_no_sidecars(tmp_path: Path) -> None:
    """Quilt's default sandbox withholds same-origin access, so a sidecar could never be fetched."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST00000000001", 10)],
        [_hit_row(GUIDE, "human", "ENST00000099999", 1, "off_target", "SOMEGENE")],
    )
    html = render_html(build_payload(run))

    assert html.startswith("<!DOCTYPE html>")
    for placeholder in ("PLOTLY_JS_PLACEHOLDER", "GUIDES_JSON_PLACEHOLDER", "FILTERS_JSON_PLACEHOLDER"):
        assert placeholder not in html, f"{placeholder} was not substituted"
    assert "<svg" in html, "the chart is hand-drawn inline SVG"
    assert "GUIDE_SEQUENCE_TODO" not in html


@pytest.mark.unit
def test_the_report_reaches_nothing_outside_itself(tmp_path: Path) -> None:
    """Quilt's default iframe sandbox withholds ``allow-same-origin``, so any reach outward fails.

    Enforced as a static check because that is what #103 asks for, and it is only enforceable now that
    plotly is gone: its bundle carried 45 external URLs and browser-storage references in map traces
    the report never invokes, and a static check cannot tell an inert string from a live call (D16).
    """
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST00000000001", 10)],
        [_hit_row(GUIDE, "human", "ENST00000099999", 1, "off_target", "SOMEGENE")],
    )
    html = render_html(build_payload(run))

    assert not re.findall(r"https?://", html), "the report contains an external URL"
    external = [m for m in re.findall(r"""(?:src|href)\s*=\s*["']([^"']+)""", html) if not m.startswith("#")]
    assert external == [], f"the report links a non-inline resource: {external}"
    for forbidden in ("fetch(", "XMLHttpRequest", "localStorage", "sessionStorage", "document.cookie"):
        assert forbidden not in html, f"{forbidden} cannot work in Quilt's default sandbox"
