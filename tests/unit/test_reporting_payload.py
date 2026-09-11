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

from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS, unclassified_cells
from sirnaforge.reporting import ReportInputError, build_payload, render_html
from sirnaforge.reporting.payload import EMBED_MAX_NM, REASON_OK

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


def _write_fully_evidenced_run(tmp_path: Path, *, run_label: str = "PASS") -> Path:
    """A run that records what every declared gate observed, as the fixed design path now does."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    columns = _CANDIDATE_COLUMNS + "".join(f",{k}_observed" for k in _PASSING_OBSERVED)
    row = _candidate_row("c1", GUIDE, "ENST00000000001", 10).rsplit(",", 1)[0] + f",{run_label}"
    row += "".join(f",{v}" for v in _PASSING_OBSERVED.values())
    (run / "sirnaforge" / "candidates_all.csv").write_text(f"{columns}\n{row}\n")
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
    integer one, so ``max_off_target_count`` read ``unknown`` on all 5,706 guides of an MSH3 run whose
    threshold was rejecting 65% of its candidates.
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
    let the report fail a guide on a scope wider than the gate's: on a four-species MSH3 run the two
    disagree 17,600 hits against 63,801.
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
    re-derive it. Calling such a guide clean would have published 185 guides of one MSH3 run as
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
    unevaluable gate and one MSH3 report published ``0 pass`` across 5,706 guides. With the verdicts
    recorded, a guide that satisfies every gate reads as one.
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
def test_the_map_does_not_colour_an_unestablished_window_as_passing(tmp_path: Path) -> None:
    """The map is subject to the same rule as the status column: report less, never more.

    Classifying a run-PASS row as "passes every gate" regardless of whether the report could establish
    it drew all 11,520 passing windows of one MSH3 run under that legend while the report itself called
    every one of their guides *not established*.
    """
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST00000000001", 10)], [])
    payload = build_payload(run)
    classes = {k for t in payload.run["transcripts"] for k in t["series"]}

    assert payload.guides[0].status == "unknown", "the fixture exports none of the screening counters"
    assert classes == {"unknown"}, "a window the report cannot establish is not drawn as passing"


@pytest.mark.unit
def test_the_structure_is_laid_out_from_the_published_dot_bracket(tmp_path: Path) -> None:
    """The picture must be of the fold the gates used, so the dot-bracket travels with the guide.

    Layouts are keyed by structure and computed once per distinct one: 40,079 candidates on an MSH3 run
    carry 1,333 distinct structures, which is the difference between embedding them and not.
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
def test_the_renderer_knows_every_verdict_the_payload_emits(tmp_path: Path) -> None:
    """A verdict code the renderer cannot name is a TypeError in the browser, not a blank cell.

    ``VERDICT`` listed four codes while the payload emits five: a warn-action gate the guide exceeds
    is code 4, so ``VERDICT[4]`` was undefined and ``v.replace`` threw for every guide carrying one.
    On one MSH3 run that was 3,098 of 5,000 embedded guides -- the whole gate panel, gone.
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
