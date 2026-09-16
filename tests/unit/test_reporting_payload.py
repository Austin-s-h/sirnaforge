"""Contract tests for the single-file HTML report (issue #103).

Two of these pin the failure modes the issue names as most likely to ship a wrong report: the
evidence join keyed on the candidate id rather than the guide sequence, and non-liability alignments
counted as off-targets. The rest pin that an unevaluable gate stays unknown rather than becoming a
pass, which is the mistake the report exists to make visible.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest

from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS, unclassified_cells
from sirnaforge.reporting import ReportInputError, build_payload, render_html
from sirnaforge.reporting.payload import (
    EMBED_MAX_NM,
    MIN_UNCOVERED_NT,
    REASON_FILTER_OFF,
    REASON_OK,
    REASON_RUN_NOT_EVALUATED,
    STATUSES,
    reevaluate_gates,
)

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


#: An observed value that satisfies each declared gate under the default profile. Twelve gates are
#: evaluated and four are off; a run that records all of them is the only shape in which the report
#: can legitimately call a guide clean.
_PASSING_OBSERVED = {
    "gc_content_min": 45.0,
    "gc_content_max": 45.0,
    "max_poly_runs": 2,
    "max_repeat_transcript_fraction": 0.0,
    "max_paired_fraction": 0.3,
    "min_asymmetry_score": 0.8,
    "min_empirical_score": 0.5,
    "min_isoform_coverage": "",
    "max_off_target_count": 0,
    "max_transcriptome_hits_0mm": 0,
    "max_transcriptome_hits_1mm": 0,
    "max_transcriptome_hits_2mm": 0,
    "max_transcriptome_seed_perfect": "",
    "max_mirna_perfect_seed": 0,
    "max_mirna_1mm_seed": "",
    "fail_on_high_risk_mirna": 0,
    "max_total_offtarget_hits": "",
}


def _write_fully_evidenced_run(
    tmp_path: Path,
    *,
    run_label: str = "PASS",
    recorded: Mapping[str, str] | None = None,
    guides: Sequence[str] = (GUIDE,),
) -> Path:
    """A run that records what every declared gate observed, as the fixed design path now does.

    ``recorded`` adds ``<filter_id>_verdict`` columns, which is the only way a fixture can express the
    verdict the *run* reached for one gate -- "in force, evidence unavailable" is a state no threshold
    comparison over the observed columns reconstructs. Every guide in ``guides`` carries the same
    observed values and the same recorded verdicts, so a gate recorded ``not_evaluated`` here is
    not evaluated for the whole run rather than for one row.
    """
    recorded = recorded or {}
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    columns = _CANDIDATE_COLUMNS + "".join(f",{k}_observed" for k in _PASSING_OBSERVED)
    columns += "".join(f",{k}_verdict" for k in recorded)
    rows = []
    for n, guide in enumerate(guides, start=1):
        row = _candidate_row(f"c{n}", guide, "ENST00000000001", 10 * n).rsplit(",", 1)[0] + f",{run_label}"
        row += "".join(f",{v}" for v in _PASSING_OBSERVED.values())
        row += "".join(f",{v}" for v in recorded.values())
        rows.append(row)
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([columns, *rows]) + "\n")
    return run


#: How each report status is produced from the observed columns, given the default gate panel. ``fail``
#: moves the GC floor's own observation below it, ``warn`` trips the warn-action asymmetry floor, and
#: ``unknown`` is the one state no observation can express -- it needs the run's own recorded non-verdict.
_GRADE_OBSERVED: Mapping[str, Mapping[str, object]] = {
    "pass": {},
    "warn": {"min_asymmetry_score": 0.1},
    "fail": {"gc_content_min": 10.0},
    "unknown": {},
}
_GRADE_RECORDED = {"unknown": ("max_repeat_transcript_fraction", "unknown")}
_GRADE_RUN_LABEL = {"fail": "GC_OUT_OF_RANGE"}


def _distinct_guide(index: int) -> str:
    """A distinct 21-mer per index. Gates read the ``_observed`` columns, so only distinctness matters."""
    tail = "".join("ACGU"[(index >> shift) & 3] for shift in (0, 2, 4, 6))
    return ("ACGUACGUACGUACGUA"[:17] + tail)[:21]


def _write_graded_run(tmp_path: Path, grades: Sequence[str]) -> Path:
    """A fully evidenced run with one guide per requested status, best-scoring first.

    The cap has to be tested against a run whose verdicts differ, because what it must never do is
    drop a guide by score alone. Every guide records what all seventeen gates observed, as the fixed
    design path does, and the run's own ``passes_filters`` label agrees with the status asked for.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    verdict_columns = sorted({column for column, _ in (_GRADE_RECORDED.get(g, ("", "")) for g in grades) if column})
    columns = (
        _CANDIDATE_COLUMNS
        + "".join(f",{k}_observed" for k in _PASSING_OBSERVED)
        + "".join(f",{c}_verdict" for c in verdict_columns)
    )
    rows = []
    for n, grade in enumerate(grades, start=1):
        observed = dict(_PASSING_OBSERVED, **_GRADE_OBSERVED[grade])
        recorded = _GRADE_RECORDED.get(grade)
        row = _candidate_row(f"c{n}", _distinct_guide(n), "ENST00000000001", 10 * n).split(",")
        row[10] = str(1000 - n)  # composite_score: the order the payload ranks them in
        row[12] = _GRADE_RUN_LABEL.get(grade, "PASS")
        cells = [*row, *(str(v) for v in observed.values())]
        cells += [recorded[1] if recorded and recorded[0] == c else "" for c in verdict_columns]
        rows.append(",".join(cells))
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([columns, *rows]) + "\n")
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

    A gate the run exports no readable column for -- neither ``<filter_id>_observed`` nor the
    descriptor's own column -- is reported ``not established``, never flipped to pass.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    entry = build_payload(run).guides[0]

    assert entry.n_gates_unknown > 0, "the synthetic run exports none of the screening counters"
    assert entry.status == "unknown", "a guide with unevaluable gates is not a pass"
    assert entry.n_gates_failed == 0, "and it is not a failure either"


@pytest.mark.unit
def test_an_integer_counter_is_read_rather_than_nulled(tmp_path: Path) -> None:
    """``numpy.int64`` is not a Python ``int``, and ``DataFrame.iloc`` hands back numpy scalars.

    An isinstance check against ``(int, float)`` therefore passed every float column and nulled every
    integer one, so ``max_off_target_count`` read ``unknown`` on all 5,706 guides of an internal run
    whose threshold was rejecting 65% of its candidates.
    """
    columns = _CANDIDATE_COLUMNS + ",off_target_count,max_off_target_count_observed"
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10) + ",7,7"], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text(
        f"{columns}\n{_candidate_row('c1', GUIDE, 'ENST00000000001', 10)},7,7\n"
    )
    payload = build_payload(run)
    entry = payload.guides[0]
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_off_target_count")

    assert entry.metrics["off_target_count"] == 7, "an int64 counter is a number, not a blank"
    assert entry.gates[index][0] == 7, "and the gate reads it"
    assert entry.gates[index][2] == REASON_OK


@pytest.mark.unit
def test_a_gate_reads_the_value_the_run_compared_not_a_wider_counter(tmp_path: Path) -> None:
    """``<filter_id>_observed`` wins over the descriptor's column, because the run compared it.

    The human-stratified gates read counters 0.7.1 does not export under the descriptor's name, but it
    does export ``<filter_id>_observed``. Preferring the same-named all-species column instead would
    let the report fail a guide on a scope wider than the gate's: on a four-species internal run the
    two disagree 17,600 hits against 63,801.
    """
    columns = _CANDIDATE_COLUMNS + ",transcriptome_hits_1mm,max_transcriptome_hits_1mm_observed"
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    # 40 all-species hits would fail the ceiling of 10; the 2 the gate actually counted pass it.
    (run / "sirnaforge" / "candidates_all.csv").write_text(
        f"{columns}\n{_candidate_row('c1', GUIDE, 'ENST00000000001', 10)},40,2\n"
    )
    payload = build_payload(run)
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_transcriptome_hits_1mm")

    assert payload.filters[index]["read_column"] == "max_transcriptome_hits_1mm_observed"
    assert payload.guides[0].gates[index][0] == 2
    assert payload.guides[0].gates[index][1] == 0, "2 <= 10 passes; reading the wider counter would fail it"


@pytest.mark.unit
def test_an_empty_observed_column_falls_back_to_the_descriptor_column(tmp_path: Path) -> None:
    """A present-but-empty ``_observed`` column is a gate that recorded no verdict, not evidence.

    Preferring it unconditionally turned ``gc_content_min``/``gc_content_max`` -- answerable from the
    exported ``gc_content`` -- into unknowns.
    """
    columns = _CANDIDATE_COLUMNS + ",gc_content_max_observed"
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text(
        f"{columns}\n{_candidate_row('c1', GUIDE, 'ENST00000000001', 10, gc=45.0)},\n"
    )
    payload = build_payload(run)
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "gc_content_max")

    assert payload.filters[index]["read_column"] == "gc_content"
    assert payload.guides[0].gates[index][0] == 45.0


@pytest.mark.unit
def test_the_gate_panel_comes_from_the_run_not_from_library_defaults(tmp_path: Path) -> None:
    """A report must publish the thresholds its run applied.

    Resolving a fresh default policy published a GC ceiling of 60 against a run that set 65, which
    failed 227 guides on a threshold the run never applied and contradicted 41 of its own PASSes.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10, gc=62.0)], [])
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        stated={"gc_max": 65.0},
        query_species="human",
        screen_species=["human"],
    )
    manifest = json.loads((run / "sirnaforge" / "manifest.json").read_text())
    manifest["run_policy"] = policy.as_manifest()
    (run / "sirnaforge" / "manifest.json").write_text(json.dumps(manifest))

    payload = build_payload(run)
    ceiling = next(f for f in payload.filters if f["filter_id"] == "gc_content_max")

    assert ceiling["threshold"] == 65.0, "the run's ceiling, not the profile default"
    assert payload.provenance["policy_source"] == "the run's own manifest"
    index = payload.filters.index(ceiling)
    assert payload.guides[0].gates[index][1] == 0, "62% GC passes the ceiling this run actually set"


@pytest.mark.unit
def test_a_run_with_no_published_policy_says_the_gates_are_defaults(tmp_path: Path) -> None:
    """Falling back is allowed. Falling back silently is what shipped the wrong ceiling."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    payload = build_payload(run)

    assert payload.provenance["policy_source"] == "library defaults"
    assert any("library defaults" in c for c in payload.caveats)


@pytest.mark.unit
def test_a_rejection_no_declared_gate_expresses_is_not_overruled(tmp_path: Path) -> None:
    """The report may report less than the run. It may not report more.

    ``REPEAT_ELEMENT`` is stamped by the pipeline and declared by no filter, so no descriptor can
    re-derive it. Calling such a guide clean would have published 185 guides of one internal run as
    passing that the run threw out -- the fabricated-evidence direction the reverse metric now guards.
    """
    payload = build_payload(_write_fully_evidenced_run(tmp_path, run_label="REPEAT_ELEMENT"))
    entry = payload.guides[0]

    assert entry.n_gates_unknown == 0, "every declared gate is evidenced in this fixture"
    assert entry.run_verdict == "REPEAT_ELEMENT"
    assert entry.undeclared_run_rejection is True
    assert entry.status == "unknown", "a rejection the registry cannot express is not a pass"
    assert payload.run["agreement"]["overruled_run_fail"] == 0
    assert any("REPEAT_ELEMENT" in c for c in payload.caveats)


@pytest.mark.unit
def test_a_fully_evidenced_run_can_reach_a_pass(tmp_path: Path) -> None:
    """The headline defect: no guide of any run could be called clean.

    Three gates decided during enumeration and recorded no verdict, so every guide carried an
    unevaluable gate and one internal report published ``0 pass`` across 5,706 guides. With the
    verdicts recorded, a guide that satisfies every gate reads as one.
    """
    payload = build_payload(_write_fully_evidenced_run(tmp_path))
    entry = payload.guides[0]

    assert entry.n_gates_unknown == 0
    assert entry.status == "pass"
    assert payload.run["status_counts"]["pass"] == 1
    assert sum(payload.run["status_counts"].values()) == len(payload.guides), "every guide is counted once"


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
def test_the_payload_carries_a_design_map_and_the_report_draws_it(tmp_path: Path) -> None:
    """Where designs sit on the transcript is the axis the guide index cannot show.

    The map is drawn server-side by ``reporting.tracks`` so a notebook figure and this card come out of
    one function; only the selection marker is client-side, because only it depends on the selection.
    """
    run = _write_run(
        tmp_path,
        [_candidate_row(f"c{i}", GUIDE, "ENST00000000001", 100 * i) for i in (1, 2, 3)],
        [],
    )
    (run / "orf_reports").mkdir()
    (run / "orf_reports" / "orf_validation.txt").write_text(
        "transcript_id\tsequence_length\tlongest_orf_start\tlongest_orf_end\nENST00000000001\t1000\t100\t900\n"
    )
    payload = build_payload(run)
    entry = next(t for t in payload.run["transcripts"] if t["transcript_id"] == "ENST00000000001")

    assert (entry["cds_start"], entry["cds_end"], entry["length"]) == (100, 900, 1000)
    assert entry["windows"] == 3
    assert sum(len(v) for v in entry["series"].values()) == 3

    html = render_html(payload)
    assert "DESIGN MAP" in html.upper() or "Design map" in html
    assert "CDS 100-900" in html, "the region bar is drawn from the run's own ORF call"


@pytest.mark.unit
@pytest.mark.parametrize(
    ("drop", "maps", "why"),
    [
        ((), True, "everything present"),
        (("composite_score",), True, "design_score is the documented fallback"),
        (("composite_score", "design_score"), False, "no value column left to plot"),
        (("position",), False, "no x axis"),
        (("passes_filters",), False, "no run verdict to class a point by"),
    ],
)
def test_the_map_needs_a_value_column_and_says_so_when_it_has_none(
    tmp_path: Path, drop: tuple[str, ...], maps: bool, why: str
) -> None:
    """``composite_score`` does not exist before post-screen scoring, so the map cannot require it.

    A design-only run, a run whose screening failed, and the "keeping design-time score" path all
    produce candidates with no composite. Requiring it dropped the whole map for those runs with no
    caveat; requiring a *value* column and naming the one used keeps the panel and explains any loss.
    """
    names = _CANDIDATE_COLUMNS.split(",")
    values = dict(zip(names, _candidate_row("c1", GUIDE, "ENST1", 10).split(","), strict=True))
    keep = [c for c in names if c not in drop]
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text(
        ",".join(keep) + "\n" + ",".join(values[c] for c in keep) + "\n"
    )

    payload = build_payload(run)
    assert bool(payload.run["transcripts"]) is maps, why
    if maps:
        assert payload.run["transcripts"][0]["value_column"] == ("design_score" if drop else "composite_score")
    else:
        assert any("no design map" in c for c in payload.caveats), "a lost panel is explained, not silent"
    assert render_html(payload).startswith("<!DOCTYPE html>")


@pytest.mark.unit
def test_the_map_does_not_colour_an_unestablished_window_as_passing(tmp_path: Path) -> None:
    """The map is subject to the same rule as the status column: report less, never more.

    Classifying a run-PASS row as "passes every gate" regardless of whether the report could establish
    it drew all 11,520 passing windows of one internal run under that legend while the report itself
    called every one of their guides *not established*.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    payload = build_payload(run)
    classes = {k for t in payload.run["transcripts"] for k in t["series"]}

    assert payload.guides[0].status == "unknown", "the fixture exports none of the screening counters"
    assert classes == {"unknown"}, "a window the report cannot establish is not drawn as passing"


@pytest.mark.unit
def test_the_structure_is_laid_out_from_the_published_dot_bracket(tmp_path: Path) -> None:
    """The picture must be of the fold the gates used, so the dot-bracket travels with the guide.

    Layouts are keyed by structure and computed once per distinct one: 40,079 candidates on an internal
    run carry 1,333 distinct structures, which is the difference between embedding them and not.
    """
    fold = ".....((((....))))....."
    columns = _CANDIDATE_COLUMNS + ",structure,mfe"
    rows = [
        f"{_candidate_row('c1', GUIDE, 'ENST00000000001', 10)},{fold},-1.9",
        f"{_candidate_row('c2', OTHER, 'ENST00000000001', 40)},{fold},-1.9",
    ]
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([columns, *rows]) + "\n")

    payload = build_payload(run)

    assert {g.structure for g in payload.guides} == {fold}, "each guide carries its own published fold"
    assert len(payload.run["structure_layouts"]) == 1, "two guides, one distinct structure, one layout"


@pytest.mark.unit
def test_a_malformed_structure_costs_its_own_panel_and_nothing_else(tmp_path: Path) -> None:
    """A dot-bracket that cannot describe the guide must not take the rest of the detail pane down.

    The JS drawer indexed a coordinate array by sequence position; a structure shorter than its guide
    made that undefined, threw inside the template literal, and blanked gates, isoforms and off-targets
    along with the structure.
    """
    columns = _CANDIDATE_COLUMNS + ",structure"
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text(
        f"{columns}\n{_candidate_row('c1', GUIDE, 'ENST00000000001', 10)},..((..\n"
    )
    html = render_html(build_payload(run))

    assert "function balanced(" in html, "the client checks the same things structure.py raises on"
    assert "db.length!==g.guide.length" in html
    assert "function card(render, g)" in html, "and one panel raising cannot cost the pane"


@pytest.mark.unit
def test_the_map_note_states_the_threshold_the_code_uses(tmp_path: Path) -> None:
    """A caption naming a hardcoded 40 nt outlives the constant it was copied from."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST1", 1), _candidate_row("c2", OTHER, "ENST1", 900)],
        [],
    )
    (run / "orf_reports").mkdir()
    (run / "orf_reports" / "orf_validation.txt").write_text(
        "transcript_id\tsequence_length\tlongest_orf_start\tlongest_orf_end\nENST1\t1000\t100\t900\n"
    )
    html = render_html(build_payload(run))

    assert f"{MIN_UNCOVERED_NT} nt or more carry no candidate in this table" in html
    assert "no enumerated window" not in html, "the table's coverage is not the transcript's"


@pytest.mark.unit
def test_the_mirna_panel_names_the_mirna_it_matched(tmp_path: Path) -> None:
    """A seed match with no name is not a finding.

    The aggregate publishes the name as ``mirna_id``; the payload read ``rname`` or ``mirna``, neither
    of which the table has, so every row rendered with an empty name -- 18,078 anonymous seed matches
    on one internal run. Which miRNA is mimicked is the whole question, and ``coord`` travels with it
    because a motif matching away from position 1 is what a real defect in this scanner once counted
    as a perfect seed hit.
    """
    hit = _hit_row(GUIDE, "hsa", "unused", 0, "off_target", "X")
    hit = {k: v for k, v in hit.items() if k not in CLASSIFICATION_COLUMNS}
    hit.update(mirna_id="Hsa-Mir-24-P2_3p", database="mirgenedb", coord=1, seed_mismatches=0)
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1", 10)], [])
    table = run / "off_target" / "results" / "aggregated" / "combined_mirna_hits.tsv"
    header = list(hit)
    table.write_text("\t".join(header) + "\n" + "\t".join(str(hit[c]) for c in header) + "\n")

    entry = build_payload(run).guides[0]

    assert len(entry.mirna) == 1
    assert entry.mirna[0]["mirna"] == "Hsa-Mir-24-P2_3p", "the name the table published"
    assert entry.mirna[0]["database"] == "mirgenedb"
    assert entry.mirna[0]["coord"] == 1, "the seed offset travels with the hit"


@pytest.mark.unit
def test_cross_species_conservation_is_published_per_species_with_its_mismatches(tmp_path: Path) -> None:
    """The ortholog class is conservation evidence, and the summary of it cannot answer the question.

    ``conservation_score`` is (species hit)/3 and counts a species conserved at up to 8 mismatches; on
    one internal run mouse ortholog hits ran 1,413 at nm=0 against 1,493 at nm>=3. Publishing the best
    ``(nm, seed_mismatches)`` per species is what lets a reader require "perfect in macaque, seed-intact
    in mouse" -- which took that run's cross-reactive pool from 49 guides to 113.

    The symbol is a synthetic placeholder: this path groups by species, so nothing here reads it. The
    ``GENEX``/``Genex`` casing only mirrors the primate/rodent symbol convention the hits would carry.
    """
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST1", 10)],
        [
            _hit_row(GUIDE, "macaque", "ENSMMUT1", 0, "ortholog", "GENEX"),
            _hit_row(GUIDE, "mouse", "ENSMUST2", 3, "ortholog", "Genex"),
            _hit_row(GUIDE, "mouse", "ENSMUST1", 1, "ortholog", "Genex"),
        ],
    )
    entry = build_payload(run).guides[0]

    assert entry.ortholog["macaque"] == {"nm": 0, "seed_mismatches": 0}
    assert entry.ortholog["mouse"]["nm"] == 1, "the best alignment, not the first or the worst"
    assert "rat" not in entry.ortholog, "a species with no ortholog hit is absent, not zero"

    html = render_html(build_payload(run))
    assert "function passesConservation(" in html, "and it is a reader's threshold, not a gate"
    assert "allow 1 mm outside the seed" in html


@pytest.mark.unit
def test_the_renderer_knows_every_verdict_the_payload_emits(tmp_path: Path) -> None:
    """A verdict code the renderer cannot name is a TypeError in the browser, not a blank cell.

    ``VERDICT`` listed four codes while the payload emits five: a warn-action gate the guide exceeds
    is code 4, so ``VERDICT[4]`` was undefined and ``v.replace`` threw for every guide carrying one.
    On one internal run that was 3,098 of 5,000 embedded guides -- the whole gate panel, gone.
    """
    observed = dict(_PASSING_OBSERVED, min_asymmetry_score=0.1)  # below the 0.65 warn-action floor
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    columns = _CANDIDATE_COLUMNS + "".join(f",{k}_observed" for k in observed)
    row = _candidate_row("c1", GUIDE, "ENST00000000001", 10) + "".join(f",{v}" for v in observed.values())
    (run / "sirnaforge" / "candidates_all.csv").write_text(f"{columns}\n{row}\n")

    payload = build_payload(run)
    entry = payload.guides[0]
    assert entry.n_gates_warned == 1, "the fixture must produce a warn verdict"
    assert entry.status == "warn", "a warn is not a fail and not a plain pass"

    html = render_html(payload)
    codes = re.search(r"const VERDICT=\[(.*?)\];", html)
    assert codes is not None
    names = [c.strip().strip("'") for c in codes.group(1).split(",")]
    assert len(names) == 5 and names[4] == "warn", "every emitted verdict code needs a name"
    assert ".v-warn{" in html, "and a style, or the pill renders unstyled"
    assert "warn:1" in html.replace(" ", ""), "and a sort rank, or sorting by status yields NaN"


@pytest.mark.unit
def test_the_cart_exports_a_tsv_and_degrades_where_downloads_are_blocked(tmp_path: Path) -> None:
    """Export is a real .tsv download, with a stated fallback rather than a silent failure.

    Quilt's default iframe sandbox can withhold ``allow-downloads``, so the Blob path is wrapped and
    the textarea is the declared way out -- a button that quietly does nothing is worse than one that
    says why. ``navigator.clipboard`` stays out entirely: it needs a permission the sandbox also
    withholds, and there is no way to tell a refusal from a success.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1", 10)], [])
    html = render_html(build_payload(run))

    for control in ('id="filters"', 'id="addtop"', 'id="cartcard"', 'id="carttsv"', 'id="cartdl"'):
        assert control in html, f"{control} is missing"
    assert "function passesFilters(" in html, "thresholds are applied client-side"
    assert "function togglePick(" in html

    assert "URL.createObjectURL" in html and "a.download=name" in html, "export writes a real file"
    assert "revokeObjectURL" in html, "and does not leak the blob url"
    assert "download refused by this viewer" in html, "a blocked download is reported, not swallowed"
    assert "navigator.clipboard" not in html, "needs a permission the sandbox withholds"
    assert "GENE_JSON_PLACEHOLDER" not in html, "the filename's gene was substituted"


@pytest.mark.unit
def test_a_gate_that_is_off_still_reports_what_it_measured(tmp_path: Path) -> None:
    """An off gate is not evaluated; that does not mean it measured nothing.

    Three of the four off gates on the reference run carry a real per-guide number, and
    ``max_mirna_1mm_seed``'s own policy definition says "the number is reported and nothing acts on
    it" -- which returning None for the value made untrue. The verdict stays ``not_evaluated``.
    """
    columns = _CANDIDATE_COLUMNS + ",mirna_hits_1mm_seed"
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text(f"{columns}\n{_candidate_row('c1', GUIDE, 'ENST1', 10)},4\n")
    payload = build_payload(run)
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_mirna_1mm_seed")
    value, verdict, reason = payload.guides[0].gates[index]

    assert payload.filters[index]["action"] == "off"
    assert value == 4, "the number the gate would have compared"
    assert verdict == 3, "and it is still not_evaluated, not a pass"
    assert reason == REASON_FILTER_OFF


@pytest.mark.unit
def test_an_in_force_gate_the_run_did_not_evaluate_is_unknown_rather_than_a_pass(tmp_path: Path) -> None:
    """A gate that is in force and undecided is not evidence of cleanliness (#103).

    ``n_gates_unknown`` counted only the UNKNOWN code, so a gate the run recorded ``not_evaluated``
    landed on verdict code 3 and was neither failed nor unknown: ``status`` fell through to ``pass``,
    the not-evaluated banner was suppressed, and the guide entered the Passing preset -- against
    #103's own rule that a guide failing no gate is clean only when every gate could be evaluated.
    ``max_repeat_transcript_fraction`` reaches this state on any run whose repeat scan did not happen:
    the gate stays in force, and ``repeat_transcript_fraction`` is exported unconditionally from a
    0.0-default field, so a confident PASS was re-derived for a gate nobody evaluated.

    Verdict code 3 is therefore reserved for a gate that was never in force at all -- off, or with no
    declared threshold -- and the run's own ``not_evaluated`` on an in-force gate is published as
    UNKNOWN with :data:`REASON_RUN_NOT_EVALUATED` naming which non-decision it was.
    """
    payload = build_payload(
        _write_fully_evidenced_run(tmp_path, recorded={"max_repeat_transcript_fraction": "not_evaluated"})
    )
    entry = payload.guides[0]
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_repeat_transcript_fraction")
    descriptor = payload.filters[index]

    assert (descriptor["action"], descriptor["threshold"]) == ("fail", 0.001), "the gate is in force on this run"
    assert tuple(entry.gates[index]) == (0, 2, REASON_RUN_NOT_EVALUATED), "in force and undecided is unknown"
    assert entry.n_gates_unknown == 1, "the tally must see it, or the status below cannot"
    assert entry.n_gates_failed == 0
    assert entry.run_verdict == "PASS"
    assert entry.undeclared_run_rejection is False, "the run passed the guide; only the gate is undecided"
    assert entry.status == "unknown", "an unevaluated in-force gate cannot make a guide clean"
    assert payload.run["status_counts"]["pass"] == 0

    # The browser agrees by construction rather than by a second tally rule it would have to be taught:
    # the shipped counter counts the UNKNOWN code, which is the code this state is now published under.
    assert "t[1]===V_UNKNOWN" in render_html(payload), "the client's unknown tally reads the UNKNOWN code"


@pytest.mark.unit
def test_a_gate_no_guide_was_evaluated_for_is_frozen_rather_than_given_a_dead_slider(tmp_path: Path) -> None:
    """``evaluable`` means one thing: moving this control can change a verdict (#103).

    A gate every one of whose rows the run recorded ``not_evaluated`` still exports the number it
    measured, so "this run produced at least one value" called it evaluable and the report shipped a
    live slider that cannot decide anything: a reader moves it, every row stays frozen on its reason
    code, nothing changes, and nothing says why. It is published frozen with its own reason instead --
    the reason the panel already prints beside the control it withholds.
    """
    payload = build_payload(
        _write_fully_evidenced_run(
            tmp_path, guides=(GUIDE, OTHER), recorded={"max_repeat_transcript_fraction": "not_evaluated"}
        )
    )
    frozen = next(f for f in payload.filters if f["filter_id"] == "max_repeat_transcript_fraction")
    index = payload.filters.index(frozen)

    assert len(payload.guides) == 2
    assert {tuple(g.gates[index]) for g in payload.guides} == {(0, 2, REASON_RUN_NOT_EVALUATED)}
    assert frozen["evaluable"] is False, "no slider position can change a verdict this run never reached"
    assert frozen["control"] is None
    assert frozen["n_values"] == 2, "the values were measured; the run simply applied no verdict to them"
    assert frozen["unevaluable_reason"] in render_html(payload), "and the panel says why it is frozen"

    # A gate the run decided for at least one guide stays movable, or `evaluable` would mean
    # "sometimes decidable" here and "decidable" everywhere else.
    live = next(f for f in payload.filters if f["filter_id"] == "gc_content_max")
    assert live["evaluable"] is True
    assert live["control"] is not None


@pytest.mark.unit
def test_the_guide_cap_is_a_published_count_and_every_number_counts_the_embedded_guides(tmp_path: Path) -> None:
    """A cap the counts do not know about is a report advertising guides it does not contain.

    One internal report headlined 5,706 guides, pill-counted 448/976/0/4,282 and embedded 5,000 of them:
    448 pass, 935 warn, 3,617 fail. The 706 it dropped included 41 pass-warned guides -- shippable
    candidates that cannot be searched, carted or exported -- and its agreement line claimed coverage of
    1,424 run PASSes while holding 1,383 of them. The cap now lives above every count.
    """
    # The fails score HIGHEST here, so keeping the top 5 by score and keeping the 5 a reader can use
    # are different answers -- which is the whole distinction the old cap could not make.
    grades = ["fail", "fail", "fail", "fail", "pass", "pass", "warn", "unknown"]
    payload = build_payload(_write_graded_run(tmp_path, grades), max_guides=5)
    run = payload.run

    assert len(payload.guides) == 5
    assert (run["guides"], run["guides_total"], run["guides_dropped"]) == (5, 8, 3)
    assert run["guide_embed_limit"] == 5
    assert sum(run["status_counts"].values()) == len(payload.guides), "the pills count what is embedded"
    assert run["status_counts"] == {"fail": 1, "unknown": 1, "warn": 1, "pass": 2}
    assert run["guides_dropped_by_status"] == {"fail": 3, "unknown": 0, "warn": 0, "pass": 0}
    assert run["agreement"]["run_pass_guides"] == 4, "the denominator is embedded run PASSes, not the run's"
    assert run["candidate_rows"] == sum(g.n_rows for g in payload.guides)
    assert run["candidate_rows_total"] == 8, "the run's own row count is still published, under its own name"
    assert any("embeds 5 of the run's 8 guides" in c for c in payload.caveats), "a dropped guide is stated"
    assert render_html(payload).startswith("<!DOCTYPE html>")


@pytest.mark.unit
def test_the_cap_drops_the_verdicts_the_run_rejected_before_anything_a_reader_could_ship(tmp_path: Path) -> None:
    """Dropping by score alone is what lost 41 shippable guides; dropping by verdict cannot.

    A ``fail`` is a guide the run and the report agree to reject, so it is the only kind whose absence
    costs a reader nothing they could act on. When even that is not enough room, the shortfall is named.
    """
    # Again the fail is the best-scoring guide in the run, so a score cap would have kept it and thrown
    # away one of the three a reader came for.
    run = _write_graded_run(tmp_path, ["fail", "pass", "warn", "unknown"])
    kept = build_payload(run, max_guides=3)

    assert [g.status for g in kept.guides] == ["pass", "warn", "unknown"], "the fail went first"
    assert kept.run["guides_dropped_by_status"]["fail"] == 1

    squeezed = build_payload(run, max_guides=1)
    assert [g.status for g in squeezed.guides] == ["pass"], "then unknown, then warn, and a pass last"
    assert squeezed.run["guides_dropped_by_status"] == {"fail": 1, "unknown": 1, "warn": 1, "pass": 0}
    assert any("a reader might have acted on" in c for c in squeezed.caveats), "and it says so plainly"
    # Three of this run's four guides carry a run PASS; one of them is embedded. The agreement line has
    # to say one, or it claims coverage of two guides a reader cannot open.
    assert squeezed.run["agreement"]["run_pass_guides"] == 1


@pytest.mark.unit
def test_an_uncapped_payload_embeds_the_run_whole_and_says_nothing_about_a_cap(tmp_path: Path) -> None:
    """``max_guides=None`` is the escape hatch the cap's own caveat points a reader at."""
    grades = ["pass", "warn", "fail", "fail"]
    payload = build_payload(_write_graded_run(tmp_path, grades), max_guides=None)

    assert len(payload.guides) == payload.run["guides_total"] == 4
    assert payload.run["guides_dropped"] == 0
    assert payload.run["guides_dropped_by_status"] == dict.fromkeys(STATUSES, 0)
    assert payload.run["guide_embed_limit"] is None
    assert not [c for c in payload.caveats if "not embedded" in c]


@pytest.mark.unit
def test_a_slider_domain_is_snapped_to_its_step_and_clamped_to_the_settings_own_bounds(tmp_path: Path) -> None:
    """An unsnapped, unclamped domain published thresholds the pipeline would refuse.

    ``gc_content_min`` opened at 27.804348 with a step of 0.1, so every drag landed on ...704348, 40 was
    unreachable, and 39.704348 is what went into the URL fragment and the cart TSV. Padding an observed
    minimum by 5% also published "at most -59 off-targets" and a ``max_paired_fraction`` floor of
    -0.03913 -- thresholds that fail every guide, on quantities their own fields declare non-negative.
    """
    observed = [
        dict(_PASSING_OBSERVED, gc_content_min=27.804348, max_off_target_count=0, max_paired_fraction=0.0),
        dict(_PASSING_OBSERVED, gc_content_min=64.7, max_off_target_count=59, max_paired_fraction=0.783),
    ]
    columns = _CANDIDATE_COLUMNS + "".join(f",{k}_observed" for k in _PASSING_OBSERVED)
    rows = [
        _candidate_row(f"c{n}", _distinct_guide(n), "ENST00000000001", 10 * n)
        + "".join(f",{v}" for v in values.values())
        for n, values in enumerate(observed, start=1)
    ]
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([columns, *rows]) + "\n")
    payload = build_payload(run)
    by_id = {f["filter_id"]: f for f in payload.filters}

    gc = by_id["gc_content_min"]["control"]
    assert gc["step"] == 0.1
    assert round((40.0 - gc["min"]) / gc["step"], 6).is_integer(), "a reader must be able to land on 40"
    for f in payload.filters:
        control = f["control"]
        if control is None:
            continue
        values = [g.gates[payload.filters.index(f)][0] for g in payload.guides]
        values = [v for v in values if v is not None]
        assert control["min"] <= min(values) and control["max"] >= max(values), f"{f['filter_id']} hides a value"
        assert control["min"] <= f["threshold"] <= control["max"], f"{f['filter_id']} cannot be reset"
        if f["bound_min"] is not None:
            assert control["min"] >= min(f["bound_min"], *values), f"{f['filter_id']} goes below its field"
        if f["bound_max"] is not None:
            assert control["max"] <= max(f["bound_max"], *values), f"{f['filter_id']} goes above its field"
    assert by_id["max_off_target_count"]["control"]["min"] == 0, "there is no such thing as -59 off-targets"
    assert by_id["max_paired_fraction"]["control"]["min"] == 0


@pytest.mark.unit
def test_an_inert_gates_domain_can_still_reach_a_threshold_that_changes_a_verdict(tmp_path: Path) -> None:
    """The one question worth asking of an inert gate is what raising it would do, so it must be askable.

    ``min_empirical_score`` sat at 0.4 against an attainable range of 0.4-0.6, and the domain padded the
    observed values to 0.39-0.61 -- either side of the two thresholds that answer the question, and both
    of them values the field itself rejects. Clamping to the field's own bounds makes 0.6 reachable.
    """
    observed = [dict(_PASSING_OBSERVED, min_empirical_score=score) for score in (0.4, 0.5, 0.6)]
    columns = _CANDIDATE_COLUMNS + "".join(f",{k}_observed" for k in _PASSING_OBSERVED)
    rows = [
        _candidate_row(f"c{n}", _distinct_guide(n), "ENST00000000001", 10 * n)
        + "".join(f",{v}" for v in values.values())
        for n, values in enumerate(observed, start=1)
    ]
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([columns, *rows]) + "\n")
    payload = build_payload(run)
    empirical = next(f for f in payload.filters if f["filter_id"] == "min_empirical_score")
    index = payload.filters.index(empirical)

    assert (empirical["bound_min"], empirical["bound_max"]) == (0.4, 0.6)
    assert (empirical["control"]["min"], empirical["control"]["max"]) == (0.4, 0.6)

    at_run = [reevaluate_gates(payload.filters, g.gates)[index][1] for g in payload.guides]
    at_top = [
        reevaluate_gates(payload.filters, g.gates, {"min_empirical_score": empirical["control"]["max"]})[index][1]
        for g in payload.guides
    ]
    assert at_run != at_top, "a domain that cannot change a verdict is a control over nothing"


@pytest.mark.unit
def test_liability_is_published_per_species_as_well_as_summed_over_them(tmp_path: Path) -> None:
    """The gate that rejects the most guides counts every species at once, and nothing decomposed it.

    Run-wide liability alignments on one internal run: 92,658 human, 91,536 mouse, 25,243 macaque, 23,427
    rat. ``max_off_target_count`` is scoped to all species and failed 3,317 guides, 782 of them on nothing
    else -- and 458 of those 782 pass the same ceiling on human liabilities alone, against 448 passes in
    the whole report. Both numbers are now published, because which one a reader takes decides whether
    they agree with the run.
    """
    hits = [_hit_row(GUIDE, "human", f"ENST0000001{i:04d}", 1, "off_target", f"HUMANGENE{i}") for i in range(3)]
    hits += [_hit_row(GUIDE, "mouse", f"ENSMUST0000{i:04d}", 1, "off_target", f"Mousegene{i}") for i in range(17)]
    hits += [_hit_row(GUIDE, "mouse", "ENSMUST9999999", 0, "ortholog", "Genex")]
    columns = _CANDIDATE_COLUMNS + ",off_target_count,max_off_target_count_observed,total_offtarget_hits_query"
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], hits)
    (run / "sirnaforge" / "candidates_all.csv").write_text(
        f"{columns}\n{_candidate_row('c1', GUIDE, 'ENST00000000001', 10)},20,20,3\n"
    )
    payload = build_payload(run)
    entry = payload.guides[0]
    index = next(i for i, f in enumerate(payload.filters) if f["filter_id"] == "max_off_target_count")
    ceiling = payload.filters[index]["threshold"]

    assert entry.liability_count == 20, "the all-species number the gate compared"
    assert entry.liability_by_species == {"human": 3, "mouse": 17}, "the ortholog hit is not a liability"
    assert sum(entry.liability_by_species.values()) == entry.liability_count
    assert entry.gates[index][1] == 1, "20 liabilities fail a ceiling of 15"
    assert entry.liability_by_species["human"] <= ceiling, "and on human alone the same ceiling passes it"
    assert entry.metrics["total_offtarget_hits_query"] == 3, "the query-species gate input #101 exports"
    assert payload.run["liability_rows_by_species"] == {"human": 3, "mouse": 17}
    assert payload.run["screened_species"] == ["human", "mouse"]
    # The run-wide totals stay in one scope: the species split sums to the liability rows, and those
    # plus the non-liabilities are the alignments the header prints.
    assert sum(payload.run["liability_rows_by_species"].values()) == payload.run["liability_rows"] == 20
    assert payload.run["liability_rows"] + payload.run["non_liability_rows"] == payload.run["hit_rows"] == 21


@pytest.mark.unit
def test_the_panel_says_which_gates_did_the_work_and_which_decided_nothing(tmp_path: Path) -> None:
    """Seventeen equal-looking sliders, and nothing saying that two of them decide everything.

    Sole cause of rejection over one internal payload: ``gc_content_min`` 983, ``max_off_target_count``
    782, ``max_repeat_transcript_fraction`` 171, ``max_paired_fraction`` 31,
    ``max_transcriptome_hits_0mm`` 2 -- and ``min_empirical_score`` inert by construction, its floor
    being the lowest value its own field permits.
    """
    payload = build_payload(_write_graded_run(tmp_path, ["pass", "warn", "fail", "fail"]))
    by_id = {f["filter_id"]: f for f in payload.filters}

    assert (by_id["gc_content_min"]["rejects"], by_id["gc_content_min"]["sole_rejects"]) == (2, 2)
    assert by_id["gc_content_min"]["inert"] is False
    assert by_id["min_asymmetry_score"]["warns"] == 1, "a warn is a finding, and not a rejection"
    assert by_id["min_asymmetry_score"]["inert"] is False

    empirical = by_id["min_empirical_score"]
    assert (empirical["rejects"], empirical["warns"]) == (0, 0)
    assert empirical["inert"] is True
    assert "declared minimum of 0.4" in empirical["inert_reason"], "why it cannot reject, not just that it did not"
    assert by_id["max_off_target_count"]["inert"] is True, "nothing in this fixture exceeds the ceiling"
    assert by_id["max_off_target_count"]["inert_reason"] is not None
    assert "declared" not in by_id["max_off_target_count"]["inert_reason"], "a ceiling of 15 could still reject"

    frozen = next(f for f in payload.filters if not f["evaluable"])
    assert frozen["inert"] is False and frozen["inert_reason"] is None, "frozen is not inert; it is undecided"

    # The ranking a renderer builds from this: sole_rejects can never exceed the guides that failed.
    assert sum(f["sole_rejects"] for f in payload.filters) <= payload.run["status_counts"]["fail"]


@pytest.mark.unit
def test_canonical_isoform_status_comes_from_the_run_or_is_unknown(tmp_path: Path) -> None:
    """The isoform picker is ordered by window count, which is length: the canonical one sat tenth.

    Ordering it properly needs the length and the canonical flag, and canonical status is a fact only the
    run can supply -- so it is read from the transcript FASTA the workflow writes, and is ``None``
    everywhere when the run wrote none. Nothing here infers it from being the longest.
    """
    rows = [_candidate_row("c1", GUIDE, "ENST00000000001", 10), _candidate_row("c2", GUIDE, "ENST00000000002", 20)]
    run = _write_run(tmp_path, rows, [])
    (run / "transcripts").mkdir()
    (run / "transcripts" / "TP53_canonical.fasta").write_text(
        ">ENST00000000002 TP53 type:protein_coding length:1200 canonical:true\nACGU\n"
    )
    (run / "orf_reports").mkdir()
    (run / "orf_reports" / "orf_validation.txt").write_text(
        "transcript_id\tsequence_length\tlongest_orf_start\tlongest_orf_end\n"
        "ENST00000000001\t7988\t100\t900\nENST00000000002\t1200\t100\t900\n"
    )
    payload = build_payload(run)
    isoforms = {i["transcript"]: i for i in payload.guides[0].isoforms}

    assert payload.run["canonical_transcript_ids"] == ["ENST00000000002"]
    assert payload.run["canonical_source"] == "TP53_canonical.fasta"
    assert (isoforms["ENST00000000002"]["canonical"], isoforms["ENST00000000002"]["length"]) == (True, 1200)
    assert (isoforms["ENST00000000001"]["canonical"], isoforms["ENST00000000001"]["length"]) == (False, 7988)

    # What the renderer needs: canonical first, length only as a tie-break -- and the longest is not it.
    ordered = sorted(payload.run["transcripts"], key=lambda t: (not t["canonical"], -t["length"]))
    assert [t["transcript_id"] for t in ordered] == ["ENST00000000002", "ENST00000000001"]


@pytest.mark.unit
def test_a_run_that_records_no_canonical_status_reports_unknown_rather_than_false(tmp_path: Path) -> None:
    """Not recorded and not canonical are different claims, and only one of them is this run's."""
    payload = build_payload(_write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], []))

    assert payload.run["canonical_source"] is None
    assert payload.run["canonical_transcript_ids"] == []
    assert all(i["canonical"] is None for g in payload.guides for i in g.isoforms)
    assert all(t["canonical"] is None for t in payload.run["transcripts"])


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
