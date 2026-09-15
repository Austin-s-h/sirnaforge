"""Client-side re-thresholding, and the fields a preset view needs (issue #103).

The report lets a reader move a gate's threshold and watch it re-decided in the browser.
:func:`reevaluate_gates` is the Python reference that shipped evaluator has to agree with, and the
property no threshold may ever defeat: a gate the run itself did not decide -- off, undeclared, a
missing column, or its own recorded ``unknown``/``not_evaluated`` -- freezes untouched, whatever a
reader asks for. The draft this was salvaged from got that property right in spirit but wrong in
code: it froze per-*filter*, not per-*(guide, gate)*, so a row the run marked ``not_evaluated`` while
still recording a real number fell through to a fresh pass/fail. The tests below pin the fix.

The middle of the file pins the per-filter and per-guide fields a preset view (all / passing /
near-miss / off-target-clean / register-deduplicated) cannot be built without, and that the draft
either lacked or got from the wrong column.

The last section is the other half of the same contract: that the RENDERED page says what the payload
computed. Those numbers were audited in a browser, and three of them were wrong on screen while every
test here passed -- a guide total counting a set the file did not hold, a design-map legend frozen at
the run's thresholds while every other number went live, and a scope warning on 92% of guides. They are
asserted against the rendered document rather than through a JavaScript engine, because a test that
needs node cannot be the thing that stops a renderer publishing a number it did not compute.
"""

from __future__ import annotations

import json
import re
from html import unescape
from pathlib import Path

import pytest

from sirnaforge.core.hit_annotation import unclassified_cells
from sirnaforge.reporting import build_payload, render_html
from sirnaforge.reporting.payload import (
    EMBED_MAX_NM,
    EMBED_SPECIES,
    REASON_EMPTY_VALUE,
    REASON_OK,
    REASON_RUN_NOT_EVALUATED,
    VERDICT_WARN,
    reevaluate_gates,
)
from sirnaforge.reporting.render import _GAPS_KEY, _NOT_EMBEDDED

GUIDE_A = "ACGUACGUACGUACGUACGUA"
#: One nucleotide off GUIDE_A's position on the same transcript, so the two form a register cluster.
GUIDE_B = "UUUUCCCCAAAAGGGGUUUUC"

_COLUMNS = (
    "id,guide_sequence,passenger_sequence,transcript_id,position,gc_content,asymmetry_score,"
    "design_score,composite_score,weight_vector,passes_filters,off_target_screened,screen_query_id"
)


def _row(cid: str, guide: str, tx: str, pos: int, *, composite: float, screened: bool = True) -> str:
    dna = guide.replace("U", "T")
    return f"{cid},{dna},{dna},{tx},{pos},45.0,0.8,55.0,{composite},postscreen_sirna_v4,PASS,{screened},screen_{cid}"


def _write_run(tmp_path: Path) -> Path:
    """Two guides one nucleotide apart on one transcript -- the smallest register cluster there is."""
    run = tmp_path / "run"
    (run / "sirnaforge").mkdir(parents=True)
    rows = [
        _row("c1", GUIDE_A, "ENST1", 100, composite=90.0),
        _row("c2", GUIDE_B, "ENST1", 101, composite=80.0),
    ]
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([_COLUMNS, *rows]) + "\n")
    (run / "sirnaforge" / "manifest.json").write_text(json.dumps({"gene_query": "TP53", "tool_version": "test"}))
    return run


@pytest.fixture
def payload(tmp_path: Path):
    """The payload built from :func:`_write_run`, for tests that need a real run."""
    return build_payload(_write_run(tmp_path))


# ---------------------------------------------------------------------------
# reevaluate_gates: the freeze rule, per-(guide, gate) not per-filter
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_not_evaluated_gate_with_a_measured_value_cannot_be_decided_by_any_threshold() -> None:
    """The draft's own hole: NOT_EVALUATED carrying a real number must stay code 3, always.

    A filter can be globally re-thresholdable while one guide's own row still recorded
    ``not_evaluated`` with a measured value -- the run applied no verdict for *this row*, and no
    threshold a reader chooses can conjure the decision it declined to make.
    """
    filters = [{"filter_id": "f", "comparator": "le", "threshold": 5, "action": "fail"}]
    gate = [3, 3, REASON_RUN_NOT_EVALUATED]  # value=3 measured, verdict=NOT_EVALUATED
    for threshold in (-100, 0, 3, 5, 100):
        out = reevaluate_gates(filters, [gate], {"f": threshold})
        assert out == [gate], f"threshold {threshold} decided a gate the run never applied"


@pytest.mark.unit
def test_a_warn_gate_moved_past_its_floor_stays_warn_not_fail() -> None:
    """A warn gate a reader thresholds past must stay warn, not fail.

    Folding it into fail would flip ``contradicted_run_pass`` for every guide the run merely warned
    about -- the one metric that proves the report agrees with the pipeline.
    """
    filters = [{"filter_id": "w", "comparator": "ge", "threshold": 0.5, "action": "warn"}]
    gate = [0.6, 0, REASON_OK]  # passes at the run's own floor
    out = reevaluate_gates(filters, [gate], {"w": 0.9})  # moved past this guide's own value

    assert out == [[0.6, VERDICT_WARN, REASON_OK]], "a warn floor exceeded must stay warn, never fail"


@pytest.mark.unit
def test_an_empty_value_stays_unknown_however_the_threshold_moves() -> None:
    """REASON_EMPTY_VALUE's value is already None; re-comparing it must invent nothing."""
    filters = [{"filter_id": "e", "comparator": "le", "threshold": 5, "action": "fail"}]
    gate = [None, 2, REASON_EMPTY_VALUE]
    out = reevaluate_gates(filters, [gate], {"e": 999})

    assert out == [[None, 2, REASON_EMPTY_VALUE]]


@pytest.mark.unit
def test_a_run_unknown_and_an_off_gate_are_both_frozen_untouched() -> None:
    """Every non-OK, non-empty-value reason freezes -- not only the two named above."""
    filters = [
        {"filter_id": "u", "comparator": "le", "threshold": 5, "action": "fail"},
        {"filter_id": "o", "comparator": "le", "threshold": 5, "action": "off"},
    ]
    gates = [[None, 2, 5], [4, 3, 3]]  # REASON_RUN_UNKNOWN=5, REASON_FILTER_OFF=3
    out = reevaluate_gates(filters, gates, {"u": -1, "o": -1})

    assert out == gates


@pytest.mark.unit
def test_a_passing_gate_can_still_be_rethresholded_into_a_fail() -> None:
    """The one path that must move: a REASON_OK gate answers a new threshold, not the old one."""
    filters = [{"filter_id": "p", "comparator": "le", "threshold": 5, "action": "fail"}]
    gate = [3, 0, REASON_OK]
    out = reevaluate_gates(filters, [gate], {"p": 1})

    assert out == [[3, 1, REASON_OK]]


# ---------------------------------------------------------------------------
# Per-filter preset fields: evaluable / unevaluable_reason / control / n_values
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_every_filter_carries_the_rethreshold_fields(payload) -> None:  # noqa: ANN001
    """Every filter carries the four fields a preset view's slider UI reads."""
    for f in payload.filters:
        assert "evaluable" in f
        assert "unevaluable_reason" in f
        assert "control" in f
        assert "n_values" in f
        if f["evaluable"]:
            assert f["unevaluable_reason"] is None
            assert set(f["control"]) == {"min", "max", "step"}
            assert f["n_values"] > 0
        else:
            assert f["unevaluable_reason"]
            assert f["control"] is None


@pytest.mark.unit
def test_an_off_filter_is_frozen_and_carries_no_control(payload) -> None:  # noqa: ANN001
    """An off filter is never evaluable, so it never gets a slider."""
    off = [f for f in payload.filters if f["action"] == "off"]
    assert off, "the resolved policy declares at least one off filter"
    for f in off:
        assert f["evaluable"] is False
        assert f["control"] is None


# ---------------------------------------------------------------------------
# GuideEntry: off_target_screened / screen_query_id
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_off_target_screened_and_the_query_id_ride_along_on_the_guide(payload) -> None:  # noqa: ANN001
    """Both fields are on ``GuideEntry``, read straight from the fixture's screened candidates."""
    for guide in payload.guides:
        assert guide.off_target_screened is True
        assert guide.screen_query_id


@pytest.mark.unit
def test_an_unscreened_guide_does_not_read_as_screened(tmp_path: Path) -> None:
    """A never-screened guide must not silently pass an off-target-clean preset."""
    run = tmp_path / "run"
    (run / "sirnaforge").mkdir(parents=True)
    rows = [_row("c1", GUIDE_A, "ENST1", 100, composite=90.0, screened=False)]
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([_COLUMNS, *rows]) + "\n")
    (run / "sirnaforge" / "manifest.json").write_text(json.dumps({"gene_query": "TP53", "tool_version": "test"}))

    payload = build_payload(run)

    assert payload.guides[0].off_target_screened is False


# ---------------------------------------------------------------------------
# Isoform register clusters: a candidate/guide back-reference for dedup
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_register_neighbours_carry_a_cluster_id_a_reader_can_deduplicate_on(payload) -> None:  # noqa: ANN001
    """``register_neighbours`` alone is bare ints with no guide behind them.

    A register-deduplicated preset needs to know which of several near-duplicate designs to keep;
    the cluster id plus a representative flag answers that without resolving a position back to a
    guide by hand.
    """
    isoforms = [i for g in payload.guides for i in g.isoforms]
    for i in isoforms:
        assert "register_cluster" in i
        assert "register_representative" in i

    clustered = [i for i in isoforms if i["register_neighbours"]]
    assert clustered, "the fixture's two guides sit 1 nt apart on the same transcript"
    cluster_ids = {i["register_cluster"] for i in clustered}
    assert None not in cluster_ids

    representatives = [i for i in clustered if i["register_representative"]]
    assert len(representatives) == 1, "the higher-composite-score guide is the cluster's sole keeper"


# ---------------------------------------------------------------------------
# What the reader sees: the renderer's own numbers (#103, phase 2)
# ---------------------------------------------------------------------------
# The payload became the authority on how many guides a report holds and on what each gate did; these
# pin that the rendered page says so. Every assertion below reads the real rendered document or the
# real payload -- no JavaScript engine -- because a test that needs node cannot be the thing that stops
# a renderer shipping a number it did not compute.

_LIAB_COLUMNS = _COLUMNS + ",off_target_count,max_off_target_count_observed"

#: A guide per verdict the cap has to choose between, plus the liabilities the species split needs.
#: ``gc_content`` is what carries the verdict: the default panel's GC floor rejects 20.0 and passes 45.0.
_GRADED = (
    ("g1", "ACGUACGUACGUACGUACGUA", 45.0, "PASS", 3),
    ("g2", "UUUUCCCCAAAAGGGGUUUUC", 45.0, "PASS", 1),
    ("g3", "GGGGAAAACCCCUUUUGGGGA", 20.0, "GC_OUT_OF_RANGE", 2),
    ("g4", "CCCCGGGGUUUUAAAACCCCG", 20.0, "GC_OUT_OF_RANGE", 0),
    ("g5", "AAGGCCUUAAGGCCUUAAGGC", 45.0, "PASS", 4),
    ("g6", "UUGGAACCUUGGAACCUUGGA", 45.0, "PASS", 0),
)

#: Species and mismatch count per liability alignment, cycled. The first is inside the embedded scope
#: (query species, nm <= EMBED_MAX_NM); the rest are outside it, one for each reason it can be outside.
_LIAB_SHAPE = ((EMBED_SPECIES, 1), ("mouse", 0), (EMBED_SPECIES, EMBED_MAX_NM + 2), ("mouse", 4))


def _graded_run(tmp_path: Path) -> Path:
    """Six guides over two transcripts, with liabilities inside and outside the embedded scope."""
    run = tmp_path / "graded"
    (run / "sirnaforge").mkdir(parents=True)
    (run / "off_target" / "results" / "aggregated").mkdir(parents=True)
    rows, hits = [], []
    for n, (cid, guide, gc, verdict, n_liab) in enumerate(_GRADED, start=1):
        dna = guide.replace("U", "T")
        for t, tx in enumerate(("ENST00000000001", "ENST00000000002")):
            rows.append(
                f"{cid}_{t},{dna},{dna},{tx},{60 * n + t},{gc},0.8,55.0,{1000 - n},"
                f"postscreen_sirna_v4,{verdict},True,screen_{cid},{n_liab},{n_liab}"
            )
        for k in range(n_liab):
            species, nm = _LIAB_SHAPE[k % len(_LIAB_SHAPE)]
            cells = unclassified_cells()
            cells.update(hit_class="off_target", hit_symbol=f"GENE{k}", matched_symbol=f"GENE{k}")
            hits.append(
                {
                    "qname": f"screen_{cid}",
                    "qseq": dna,
                    "species": species,
                    "rname": f"ENSTOFF{n}{k}",
                    "coord": 100 + k,
                    "strand": "+",
                    "cigar": "21M",
                    "mapq": 60,
                    "as_score": 42,
                    "nm": nm,
                    "seed_mismatches": 0,
                    "offtarget_score": 1.0,
                    **cells,
                }
            )
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join([_LIAB_COLUMNS, *rows]) + "\n")
    (run / "sirnaforge" / "manifest.json").write_text(json.dumps({"gene_query": "TP53", "tool_version": "test"}))
    header = list(hits[0].keys())
    (run / "off_target" / "results" / "aggregated" / "combined_offtargets.tsv").write_text(
        "\n".join(["\t".join(header)] + ["\t".join(str(h[c]) for c in header) for h in hits]) + "\n"
    )
    return run


def _flat(markup: str) -> str:
    """Markup as the reader's eye takes it: entities resolved, line breaks and indentation collapsed.

    A sentence in a Jinja template is wrapped for the source file, so asserting on the raw markup pins
    where the template happens to break its lines rather than what the header says.
    """
    return re.sub(r"\s+", " ", unescape(re.sub(r"<[^>]+>", " ", markup))).strip()


def _embedded(html: str, name: str) -> object:
    """One ``const NAME = <json>;`` line back out of the rendered document, parsed.

    The renderer substitutes single-line JSON into a named const, so the document itself is readable as
    data. Asserting against this rather than against ``payload`` is the point: it is the only way to
    catch the renderer changing a number on its way to the page.
    """
    match = re.search(rf"^const {name} = (.*);$", html, re.MULTILINE)
    assert match, f"the rendered document declares no {name}"
    return json.loads(match.group(1))


@pytest.mark.unit
def test_the_renderer_embeds_every_guide_the_payload_holds_and_no_more(tmp_path: Path) -> None:
    """No second cap in the renderer. The one in ``build_payload`` sits above every published count.

    A slice here is what made the header lie: the payload counted 5,706 guides, the renderer embedded
    5,000, and the pills, the agreement line and the index denominator all counted the wrong set. It is
    not enough to have made the two numbers equal by default -- ``max_guides=None`` asks for the run
    whole, and a slice downstream of that would silently truncate it again.
    """
    payload = build_payload(_graded_run(tmp_path), max_guides=None)
    guides = _embedded(render_html(payload), "G")

    assert len(guides) == len(payload.guides) == len(_GRADED)
    assert [g["guide"] for g in guides] == [g.guide for g in payload.guides]


@pytest.mark.unit
def test_the_header_counts_the_guides_in_the_file_and_states_what_the_cap_cost(tmp_path: Path) -> None:
    """The dropped guides are on screen, by verdict, with the limit that dropped them.

    ``guides_dropped_by_status`` carries every status including the zeros, so the reader can see that
    the cap took rejections and not shortlist candidates -- or that it did take some, which is the case
    worth acting on and the one the old silent slice hid.
    """
    payload = build_payload(_graded_run(tmp_path), max_guides=3)
    html = render_html(payload)
    header = _flat(html[: html.index("<main>")])

    assert len(_embedded(html, "G")) == 3
    assert "3 guides of the run's 6" in header
    assert "3 of the run's 6 guides are not in this file" in header
    assert "The embed limit is 3 guides" in header
    for status, n in payload.run["guides_dropped_by_status"].items():
        assert f"{n} {status}" in header, f"the breakdown omits {status}"
    # And every count the header prints is over the embedded three, never the run's six.
    assert f"over the {payload.run['guides']} guides in this file" in header
    assert sum(payload.run["status_counts"].values()) == 3


@pytest.mark.unit
def test_an_uncapped_report_makes_no_claim_about_a_cap(tmp_path: Path) -> None:
    """The dropped-guide statement is conditional: a report holding the run whole must not imply loss."""
    html = render_html(build_payload(_graded_run(tmp_path), max_guides=None))

    assert "are not in this file" not in html
    assert "The embed limit is" not in html


@pytest.mark.unit
def test_the_candidate_map_ships_dots_to_colour_rather_than_a_frozen_legend(tmp_path: Path) -> None:
    """The largest visual on the page was the last thing on it still frozen at the run's thresholds.

    A legend pre-rendered in Python cannot follow a reader's threshold, so it read
    ``PASS, all gates evaluated (394)`` over 394 green dots after a gate move that left 39 guides
    passing. The fix is structural rather than a wording change: the map hands over points with the
    guide behind each one and no legend at all, so there is nothing left that *could* be stale.
    """
    payload = build_payload(_graded_run(tmp_path))
    html = render_html(payload)
    maps = _embedded(html, "MAPS")
    assert maps, "the fixture must produce at least one transcript map"

    by_guide = {g.guide: i for i, g in enumerate(payload.guides)}
    for tid, entry in maps.items():
        assert "legend" not in entry, "a pre-rendered legend is the frozen legend, whatever it says"
        assert entry["points"], f"{tid} carries no points for the client to colour"
        for position, _value, index in entry["points"]:
            assert index == -1 or 0 <= index < len(payload.guides)
            if index >= 0:
                guide = payload.guides[index]
                assert by_guide[guide.guide] == index
                assert any(row["transcript"] == tid and row["position"] == position for row in guide.isoforms), (
                    f"point {position} on {tid} was resolved to a guide not enumerated there"
                )

    # The legend is data the client paints from the classes it just drew, and it covers every class the
    # points can carry -- including the one no live verdict can answer for.
    keys = {entry["key"] for entry in _embedded(html, "MAP_LEGEND")}
    assert keys >= {"pass", "warn", "unknown", "fail"} | {_NOT_EMBEDDED, _GAPS_KEY}
    assert _embedded(html, "NOT_EMBEDDED") == _NOT_EMBEDDED


@pytest.mark.unit
def test_every_map_property_the_client_reads_is_one_the_renderer_emits(tmp_path: Path) -> None:
    """The defect class that shipped a page printing the word "undefined" where the legend belongs.

    ``MAPS[t].legend`` and ``MAPS[t].note`` outlived the Python that built them: the renderer stopped
    emitting both, the client kept reading both, and every test passed because nothing executed the
    document. A map entry is bound to ``mp`` and only to ``mp``, so every property read off one can be
    checked against the keys :func:`_design_maps` actually produces.
    """
    payload = build_payload(_graded_run(tmp_path))
    html = render_html(payload)
    maps = _embedded(html, "MAPS")
    keys = set(next(iter(maps.values())))

    read = set(re.findall(r"\bmp\.([A-Za-z_]+)", html))
    assert read, "the map functions no longer bind a map entry to `mp`; this test cannot see them"
    assert read <= keys, f"the client reads map properties the renderer does not emit: {sorted(read - keys)}"


@pytest.mark.unit
def test_the_maps_geometry_belongs_to_the_svg_it_is_shipped_beside(tmp_path: Path) -> None:
    """Client-drawn dots land in the plot area of the backdrop they are drawn over, or they are lies.

    The backdrop is measured once with the points on it -- the only way to learn the value bounds the
    primitive rounds its gridlines to -- and drawn again with those bounds and no points. If the two
    calls differ in anything that moves the plot area, a title being the live case, every dot lands off
    the axes by that offset while looking perfectly plausible.
    """
    payload = build_payload(_graded_run(tmp_path))
    maps = _embedded(render_html(payload), "MAPS")

    for tid, entry in maps.items():
        geometry, svg = entry["geometry"], entry["svg"]
        # The primitive draws a gridline at each of five fractions of the plot area, so the top and the
        # bottom of the geometry it handed over must both be gridlines of the SVG it handed over with it.
        for edge in ("y0", "y1"):
            assert f'y1="{format(geometry[edge], "g")}"' in svg, (
                f"{tid}: the geometry's {edge} is not an edge of this SVG's plot area, so every dot the "
                "client places in it is offset from the axes it is drawn against"
            )
        assert geometry["value_max"] > geometry["value_min"]


@pytest.mark.unit
def test_a_window_whose_guide_the_cap_dropped_keeps_the_runs_own_colour(tmp_path: Path) -> None:
    """Every candidate row is plotted whether or not its guide was embedded (payload's own choice).

    A dropped guide has no live gate table, so its windows are the one class on the map a moved
    threshold cannot answer for. They must not borrow a verdict colour -- and must not silently vanish
    either, which would make a capped map look like a smaller run.
    """
    payload = build_payload(_graded_run(tmp_path), max_guides=2)
    maps = _embedded(render_html(payload), "MAPS")

    orphans = [pt for entry in maps.values() for pt in entry["points"] if pt[2] == -1]
    resolved = [pt for entry in maps.values() for pt in entry["points"] if pt[2] >= 0]
    assert orphans, "the cap dropped four of six guides; their windows are still plotted"
    assert resolved, "and the two embedded guides' windows still resolve to a live verdict"

    legend = {entry["key"]: entry for entry in _embedded(render_html(payload), "MAP_LEGEND")}
    assert "run" in legend[_NOT_EMBEDDED]["label"], "the key has to say whose verdict that colour is"


@pytest.mark.unit
def test_the_header_states_the_off_target_scope_and_how_much_sits_outside_it(tmp_path: Path) -> None:
    """The scope, once, at run level -- not as a warning on nearly every guide.

    On one internal run 197,727 of 232,864 liability alignments sit at nm >= 3, outside the embedded
    scope, so 4,616 of 5,000 guides raised the per-guide warning. A warning on 92% of rows is not a
    warning: the reader learns to skip it and the informative cases go with it. The size of the scope is
    a run-level fact and is stated as one, with the arithmetic that closes it.
    """
    payload = build_payload(_graded_run(tmp_path))
    html = render_html(payload)
    header = _flat(html[: html.index("<main>")])

    liabilities = sum(g.liability_count for g in payload.guides)
    in_scope = sum(len(g.offtarget_rows) for g in payload.guides)
    assert 0 < in_scope < liabilities, "the fixture must put liability both inside and outside the scope"

    assert payload.run["embed_scope"] in header
    assert f"{in_scope:,} of the {liabilities:,} liability" in header
    assert f"The other {liabilities - in_scope:,}" in header
    # And the split that says whether what sits outside is a human liability or another index at all.
    for species, n in payload.run["liability_rows_by_species"].items():
        assert f"{species} {n:,}" in header


@pytest.mark.unit
def test_the_per_guide_scope_note_is_reserved_for_a_guide_the_report_is_not_rejecting(
    tmp_path: Path,
) -> None:
    """Structural, because which branch a browser takes is not a property of the markup.

    The banner is the informative case -- a guide nothing here rejects, carrying liability this file
    cannot itemise, which is exactly when the missing rows could change a decision. A guide already
    being rejected gets the same fact stated plainly instead of a warning it does not need.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    assert "const shortlistable = g._live.status==='pass' || g._live.status==='warn';" in html
    banner = "Nothing here rejects this guide"
    assert html.count(banner) == 1, "one banner, in one branch"
    guard = html.index("const shortlistable")
    assert guard < html.index(banner) < html.index("counted, charted and gated, but not itemised here"), (
        "the banner must sit in the shortlistable arm, ahead of the plain statement for a rejected guide"
    )
    # And the payload flag that used to trigger it -- true for 92% of guides -- no longer decides
    # anything on screen. It still rides on the guide; nothing reads it.
    assert "g.offtarget_embedded_scope_empty_but_counts_exist" not in html


@pytest.mark.unit
def test_the_liability_column_and_its_bound_can_be_narrowed_to_one_species(tmp_path: Path) -> None:
    """One all-species count cannot tell a human liability from a non-human-index artefact.

    ``max_off_target_count`` compares the total over every screened species, and on one internal run
    about half of that total is non-human -- so a guide it rejected may carry no human liability at all.
    The decomposition is on the guide (``liability_by_species``, which sums to ``liability_count``), the
    index cell shows it, and one selector re-scopes the column, its sort and the liability bound
    together rather than adding a second control over the same number.
    """
    payload = build_payload(_graded_run(tmp_path))
    html = render_html(payload)

    assert _embedded(html, "SCREENED_SPECIES") == payload.run["screened_species"]
    assert _embedded(html, "EMBED_SPECIES") == EMBED_SPECIES
    assert len(payload.run["screened_species"]) > 1, "the fixture must screen more than one species"
    for guide in _embedded(html, "G"):
        assert sum(guide["liability_by_species"].values()) == guide["liability_count"]

    # The selector is built by mapping the run's own species list; no species is spelled in this path.
    assert "SCREENED_SPECIES.map(sp=>opt(sp,sp))" in html
    # The bound reads the scoped count, so narrowing the scope narrows the bound with it.
    assert "{k:'liab',      label:'liabilities',     dir:'max', get:g=>liabIn(g)}" in html
    assert "function liabIn(g){ return liabScope===null ? g.liability_count : liabOf(g, liabScope); }" in html
    # And the fragment carries the scope, or a copied URL restores a bound without its question.
    assert "parts.push('sp='+encodeURIComponent(liabScope))" in html


@pytest.mark.unit
def test_the_gate_panel_ranks_the_gates_and_says_which_decided_nothing(tmp_path: Path) -> None:
    """Seventeen equal-looking sliders invited a reader to move the ones that decide nothing.

    On one internal run five of the seventeen reject anything at all, two account for 1,765 of the
    rejections, and two cannot reject anybody because their threshold sits on the limit their own
    setting declares. The counts are the payload's; what is pinned here is that they reach the reader,
    ranked, with the inert ones labelled and their reason given rather than left looking consequential.
    """
    payload = build_payload(_graded_run(tmp_path))
    html = render_html(payload)
    filters = _embedded(html, "FILTERS")

    for entry in filters:
        for key in ("rejects", "sole_rejects", "warns", "unknowns", "inert", "inert_reason"):
            assert key in entry, f"{entry['filter_id']} reaches the panel without {key}"
    assert any(f["rejects"] for f in filters), "the fixture must have a gate that rejected something"
    inert = [f for f in filters if f["inert"]]
    assert inert, "and one that decided nothing, or the label below has nothing to say"
    assert all(f["inert_reason"] for f in inert), "an inert gate must say why it is inert"

    assert "inert in this run" in html
    assert "sole_rejects" in html and "by this gate alone" in html
    assert "gates decided nothing for the guides in this file" in re.sub(r"\s+", " ", html)
    # Ranked on the gate each rejection is solely attributable to, with filter order as the tiebreak so
    # the panel reads the same on every reload.
    assert "((b[0].sole_rejects||0)-(a[0].sole_rejects||0))" in html


@pytest.mark.unit
def test_the_header_counts_are_re_rendered_from_the_live_verdicts(tmp_path: Path) -> None:
    """Four frozen pills beside a re-coloured map is the frozen legend again with the numbers swapped.

    The template seeds them from the payload so they are right before any script runs; ``renderStatusPills``
    then owns them. The agreement line is deliberately *not* re-derived -- an agreement with the run has
    to be measured at the run's own thresholds -- so it says which thresholds it was computed at, and
    says when those are no longer the ones in force.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    for status in ("pass", "warn", "unknown", "fail"):
        assert f'id="n_{status}"' in html, f"the {status} pill has no handle to re-render"
    assert "function renderStatusPills()" in html
    assert html.count("renderStatusPills();") == 2, "on every threshold move, and on the first paint"
    assert "at the run's own thresholds" in re.sub(r"\s+", " ", html), (
        "the agreement line has to name the thresholds it was computed at"
    )
    assert 'id="agreestale"' in html and "Not re-derived at your thresholds" in html


@pytest.mark.unit
def test_the_phase_two_markup_still_reaches_nothing_outside_itself(tmp_path: Path) -> None:
    """The contract every panel added here has to keep: Quilt's default sandbox withholds same-origin.

    Re-asserted over a run with liabilities, several species and a design map, because that is the
    markup this phase changed -- the payload module's own copy of this check runs on one guide and no
    map, and would not have seen a species selector or a client-drawn overlay.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    assert not re.findall(r"https?://", html), "the report contains an external URL"
    external = [m for m in re.findall(r"""(?:src|href)\s*=\s*["']([^"']+)""", html) if not m.startswith("#")]
    assert external == [], f"the report links a non-inline resource: {external}"
    for forbidden in ("fetch(", "XMLHttpRequest", "localStorage", "sessionStorage", "document.cookie"):
        assert forbidden not in html, f"{forbidden} cannot work in Quilt's default sandbox"


# ---------------------------------------------------------------------------
# The controls, and the state a reader can send (#103, phase 3)
# ---------------------------------------------------------------------------
# Four defects the maintainer reproduced by driving the rendered page in a browser: both refusal
# banners shipped inside a `<details>` that is closed, so `checkVisibility()` was false and
# `elementFromPoint` at the banner's own box returned the table header behind it; a threshold box and
# its slider could show different numbers, after which one touch of the slider silently dropped the
# real threshold to its own maximum; the fragment was read once at load with no `hashchange` listener,
# while `syncHash` pushed a history entry per keystroke; and half the visible state -- the search box,
# the status set, the conservation controls and the cart -- was never in the URL at all.
#
# Asserted against the rendered document, and against the payload that feeds it, rather than through a
# JavaScript engine. What that cannot see is stated in the report for this phase: the behaviours below
# are pinned as wiring -- which function calls which -- and the browser run that confirmed them is not
# reproducible here.


def _js_body(html: str, signature: str) -> str:
    """One JavaScript function body out of the rendered document, matched on its braces.

    Where a call sits is the whole claim in this section: `openFilters()` in the refusal path and the
    same call somewhere harmless are indistinguishable to a substring search, and "the search box calls
    syncHash" is a statement about one handler rather than about the file.
    """
    start = html.index(signature)
    open_at = html.index("{", start)
    depth = 0
    for i in range(open_at, len(html)):
        if html[i] == "{":
            depth += 1
        elif html[i] == "}":
            depth -= 1
            if depth == 0:
                return html[open_at : i + 1]
    raise AssertionError(f"{signature} has no balanced body in the rendered document")


@pytest.mark.unit
def test_a_refusal_opens_the_container_that_was_swallowing_it(tmp_path: Path) -> None:
    """Both refusal banners were unreachable, and both exist only to be seen.

    Loading `#t=<gate>:<value>,<unknown gate>:<value>` set `refusedFilters`, cleared the banner's
    `display` -- and `checkVisibility()` still returned false, because the banner lives inside the
    thresholds `<details>` and that element ships closed. The two banners were written so that a refused
    filter or a refused keystroke "is visible rather than swallowed"; the container was doing the
    swallowing. So anything that has to be seen opens it: a refused name, a refused number, or a
    fragment that set any control at all.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    assert 'id="refused"' in html and 'id="badnum"' in html
    assert 'class="filters" id="filters"' in html and "<details open" not in html, (
        "the panel still ships closed, which is the point: it is opened when there is something to see"
    )
    assert "d.open=true" in _js_body(html, "function openFilters()")
    assert "openFilters();" in _js_body(html, "function renderRefused()"), "the URL-refusal banner"
    assert "openFilters();" in _js_body(html, "function markControl("), "the refused-keystroke banner"

    # And a fragment that moved anything opens the panel even when nothing was refused: the reader has
    # to be able to see the control it moved.
    fragment = _js_body(html, "function applyFragmentState(")
    assert "if(applied.controls || applied.refused.length) openFilters();" in fragment
    # `controls` is a property of the fragment, set by every key that reaches a control and false for the
    # bare-guide form an earlier report wrote.
    parse = _js_body(html, "function applyHashFragment(")
    assert "refused, controls:false}" in parse, "a bare guide sequence sets no control"
    assert parse.count("controls=true;") == 8, "every control-bearing key must say so"


@pytest.mark.unit
def test_a_threshold_the_slider_cannot_show_is_kept_and_labelled_not_pinned(tmp_path: Path) -> None:
    """The box and the slider are one threshold, and the code asserted they could never disagree.

    They could. The slider's domain is this run's own values, snapped and clamped to the setting's
    declared range, so a threshold the pipeline would accept can sit outside it -- and a range input
    cannot represent that: it pins at its own maximum. Typing 90 into a gate whose domain ends at
    76.104348 left the box at 90 and the slider at the pin, and one touch of the slider then dropped the
    real threshold to 76.104348 without the reader asking.

    The domain is deliberately NOT widened to make the disagreement go away -- the clamping is what
    stopped this report publishing "at most -59 off-targets" on a count-valued gate -- so the threshold
    stays where it was typed, the slider says it cannot show it, and it is disabled so it cannot answer
    for a number it does not hold.
    """
    payload = build_payload(_graded_run(tmp_path))
    html = render_html(payload)

    # The state is reachable, not hypothetical: at least one gate's slider stops short of a value its own
    # declared range allows, so a legal threshold exists that the slider cannot show.
    narrow = [
        f
        for f in payload.filters
        if f["evaluable"] and f["bound_max"] is not None and f["control"]["max"] < f["bound_max"]
    ]
    assert narrow, "the fixture must contain a gate whose slider domain is narrower than its own limits"

    assert "v<c.min || v>c.max" in _js_body(html, "function offScale("), "off-scale is a domain question"
    assert "f.control" in _js_body(html, "function offScale("), "and the domain is the payload's, not re-derived"
    row = _js_body(html, "function gateRow(")
    assert 'min="${c.min}" max="${c.max}" step="${c.step}"' in row, "the published domain, unchanged"
    assert "${offScale(f)?' disabled':''}" in row, "an off-scale slider cannot be dragged"
    assert 'id="os_${esc(f.filter_id)}">${offScaleNote(f)}' in row, "and says so where the reader is looking"
    note = _js_body(html, "function offScaleNote(")
    assert "Off this slider's scale." in note
    assert "the box above is\n    the threshold in force" in note, "which of the two controls is in force"

    # Nothing writes the pin back: the box always follows, the slider only while it can show T, and an
    # input arriving from a disabled slider is ignored rather than allowed to drop the threshold.
    paint = _js_body(html, "function paintControl(")
    assert "if(box && box!==moved) box.value = T[id];" in paint
    assert "rng.disabled = off; if(!off && rng!==moved) rng.value = T[id];" in paint
    assert "if(el.id.startsWith('rng_') && offScale(f)){ paintControl(f, null); return; }" in html


@pytest.mark.unit
def test_the_fragment_is_re_read_when_it_changes_and_costs_no_history(tmp_path: Path) -> None:
    """The URL claimed one state while the page showed another, and Back could not escape the report.

    `location.hash = encodeHash()` is a navigation: it pushed one entry per keystroke on every
    threshold, so `history.length` climbed during light use and Back then changed the URL without
    changing the view -- while nothing anywhere listened for `hashchange`, so the fragment was read
    exactly once, at load. `replaceState` keeps the URL current without owning the reader's Back button,
    and the listener makes Back, Forward, an edited address bar and a pasted link all reach the same
    code a fresh load does.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    assert "location.hash = encodeHash();" not in html, "the per-keystroke history push"
    sync = _js_body(html, "function syncHash()")
    assert "history.replaceState(null, '', '#'+h);" in sync
    # pushState and replaceState throw on an opaque origin, which is the sandbox this report is read in,
    # so the fallback is the old behaviour rather than losing the URL.
    assert "catch(e){ if(location.hash.slice(1)!==h) location.hash = h; }" in sync

    wire = _js_body(html, "function wireUI()")
    assert "window.addEventListener('hashchange'" in wire
    assert "if(location.hash.slice(1) === encodeHash()) return;" in wire, (
        "the report must ignore its own writes, or the fallback path would fight the reader's typing"
    )
    assert "applyFragmentState(location.hash.slice(1));" in wire
    # One path for the first load and for every later change, so the two can never drift.
    assert wire.count("applyFragmentState(") == 2
    fragment = _js_body(html, "function applyFragmentState(")
    assert fragment.index("resetControls();") < fragment.index("applyHashFragment(raw)"), (
        "a fragment describes the whole state, so what it omits is a default and not what was on screen"
    )


@pytest.mark.unit
def test_every_control_a_reader_can_see_travels_in_the_url(tmp_path: Path) -> None:
    """The promise that a reader who sends the URL sends the rows they saw was false as shipped.

    `t=`, `r=`, `sp=`, `preset=` and `g=` were encoded; the search box, the four status checkboxes, the
    conservation species and the seed-tolerance toggle were not, and `applyFilters()` from the search
    box never called `syncHash()` at all. The status set alone hides every fail and unknown row by
    default -- 3,617 of them on one internal run -- so the URL was carrying the verdicts without the
    rows. Each key is omitted when its control holds the default, which is what keeps `Reset` on an
    empty fragment.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))
    encode = _js_body(html, "function encodeHash()")
    parse = _js_body(html, "function applyHashFragment(")

    for key in ("g=", "t=", "r=", "sp=", "preset=", "q=", "s=", "c=", "k="):
        assert f"'{key}'" in encode, f"{key} is not written into the fragment"
        assert f"k==='{key[:-1]}'" in parse, f"{key} is written but never read back"

    # Omitted at the default, so a reset URL is empty and a shared one carries only reader choices.
    assert "if(F.q) parts.push('q='" in encode
    assert "if(!statusIsDefault()) parts.push('s='" in encode
    assert "if(cons.length || !F.consSeedIntact)" in encode
    assert "F.status = new Set(DEFAULT_STATUS)" in _js_body(html, "function resetControls()")
    assert "F.q = ''" in _js_body(html, "function resetControls()")
    # A set key REPLACES its set rather than adding to it, or a re-applied fragment would union with
    # whatever was on screen.
    assert "F.status=new Set();" in parse and "F.cons=new Set(); F.consSeedIntact=true;" in parse

    # The search text lives in F, not in the input, so the box is a view of it: a fragment can restore
    # it and Reset can clear it.
    assert "const raw=(F.q||'').trim();" in _js_body(html, "function applyFilters()")
    assert "if(qbox) qbox.value = F.q;" in _js_body(html, "function buildFilterUI()")

    # And every control that changes any of them writes the URL.
    assert "F.q = e.target.value; applyFilters(); syncHash();" in _js_body(html, "function wireUI()")
    build = _js_body(html, "function buildFilterUI()")
    assert build.count("syncHash();") == 3, "the status set, the conservation species and the seed toggle"
    assert "syncHash();" in _js_body(html, "function togglePick(")


@pytest.mark.unit
def test_the_cart_persists_in_the_fragment_because_it_has_nowhere_else_to_live(tmp_path: Path) -> None:
    """The audit's own suggestion for the cart was browser storage, which this report may not use.

    It reaches nothing outside itself, and the sandbox it is read in withholds same-origin, so no
    storage API is reachable there in the first place. The URL is the only in-page place a pick list can
    persist -- and it has the property storage does not: a reader can send it. The cost is a cap, and
    the card states on screen both that the address is the cart and when it is carrying fewer guides
    than the cart holds.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    for forbidden in ("localStorage", "sessionStorage", "document.cookie", "indexedDB"):
        assert forbidden not in html, f"{forbidden} is not available to this report at all"

    # Sequences, never indices: an index would resolve to a different guide in a differently built
    # report, where a sequence this file does not hold can be refused like any other invented name.
    assert "function cartInUrl(){ return G.filter(g=>CART.has(g.guide)).map(g=>g.guide)" in html
    assert ".slice(0, CART_IN_URL_MAX); }" in html
    assert "if(carted.length) parts.push('k='+carted.join(','));" in _js_body(html, "function encodeHash()")
    parse = _js_body(html, "function applyHashFragment(")
    assert "if(G.some(g=>g.guide===key)) CART.add(key); else refused.push(key);" in parse

    # Said on screen, in both directions.
    assert 'id="carturl"' in html
    cart = _js_body(html, "function renderCart()")
    assert "The cart is in this page's URL" in cart and "Nothing is stored in\n       this browser" in cart
    assert "This page's URL carries ${carried} of these ${rows.length} guides" in cart

    # A fragment is the whole state, so it owns the cart; the Reset button beside the thresholds does
    # not, because it must never destroy a pick list nothing else on this page can rebuild.
    assert "CART.clear();" in _js_body(html, "function applyFragmentState(")
    assert "CART" not in _js_body(html, "function resetControls()")


def _details_span(html: str, marker: str) -> tuple[int, int]:
    """Offsets of one `<details>` element's own content, nesting counted.

    Containment is the whole mechanism of the banner defect, so it is measured rather than inferred
    from an id appearing somewhere in the same document.
    """
    start = html.index(marker)
    depth = 0
    i = start
    while i < len(html):
        if html.startswith("<details", i):
            depth += 1
            i += 8
        elif html.startswith("</details>", i):
            depth -= 1
            if depth == 0:
                return start, i
            i += 10
        else:
            i += 1
    raise AssertionError(f"{marker} is never closed")


@pytest.mark.unit
def test_both_banners_sit_inside_the_container_that_has_to_be_opened(tmp_path: Path) -> None:
    """The defect was containment, so containment is what this asserts.

    `#refused` and `#badnum` were reachable in the DOM with `display` cleared and still failed
    `checkVisibility()`, because both live inside the thresholds `<details>` and it ships closed --
    `elementFromPoint` at the banner's own box returned the index's table header behind it. A test that
    only checks the two ids exist and that `openFilters` is called would stay green if a later edit
    moved either banner out of the panel, which would make the call pointless, or moved a third thing
    that must be seen in without wiring one. Both offsets are inside the element's own span, and the
    element still ships closed, which is why opening it has to be someone's job.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    open_at, close_at = _details_span(html, '<details class="filters" id="filters">')
    for banner in ('id="refused"', 'id="badnum"'):
        at = html.index(banner)
        assert open_at < at < close_at, f"{banner} is no longer inside #filters; openFilters() cannot reveal it"
    assert " open>" not in html[open_at : open_at + 60], "the panel ships closed"
    # And the lever names that element, so it is the same one measured above.
    assert "getElementById('filters')" in _js_body(html, "function openFilters()")


@pytest.mark.unit
def test_no_fragment_key_is_written_without_being_read_back(tmp_path: Path) -> None:
    """Every key the URL carries has to survive the round trip, derived rather than listed.

    The shareable-state defect was asymmetry: the encoder wrote five keys and the page held nine
    controls. A hardcoded list of the nine cannot catch the next one -- a tenth control encoded and
    never parsed would ship a URL that silently restores less than it carries, which is the same defect
    with a different name. So the two halves are read out of the rendered source and compared as sets.
    """
    html = render_html(build_payload(_graded_run(tmp_path)))

    written = set(re.findall(r"parts\.push\('([a-z]+)='", _js_body(html, "function encodeHash()")))
    read = set(re.findall(r"k==='([a-z]+)'", _js_body(html, "function applyHashFragment(")))

    assert written, "encodeHash writes no keys at all, so the regex no longer matches the source"
    assert written == read, (
        f"written but never read back: {sorted(written - read)}; read but never written: {sorted(read - written)}"
    )
