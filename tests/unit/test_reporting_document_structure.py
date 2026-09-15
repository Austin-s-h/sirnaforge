"""Static assertions over the rendered report document and over the report template's text (#103).

What these tests assert, and all they assert, are properties of the *text* the renderer produces and of
the template it produces it from: that the payload fixture still exercises every reason code
:mod:`sirnaforge.reporting.payload` can freeze on; that no declared ``filter_id`` appears as a literal
in the template, so the shipped evaluator has to read ``FILTERS`` rather than branch per gate; that the
control panel keeps its two labelled groups and that every threshold box is a real ``type="number"``;
that the ``<script>`` slice ends where a browser ends a script; and that a payload string spelling
``</script>`` is escaped, escaped *reversibly*, and never reaches the document unescaped.

**The JavaScript this report ships is no longer executed by any test.** These tests read it as text;
nothing runs it. A divergence between :func:`sirnaforge.reporting.payload.reevaluate_gates` and the
evaluator in ``render._TEMPLATE`` would therefore **not** fail this suite: the browser could publish a
verdict the pipeline would not reproduce and every assertion in this file would still pass. That is
coverage this repository used to have and now does not.

It was given up deliberately. The harness these tests replace rendered the report, sliced the
substituted ``<script>`` out of it and ran the shipped evaluator under ``node`` -- every guide x every
filter at five threshold sets, plus the presets, the reader filters, the gate panel's Why column, the
cart export and the numeric boundary. The maintainer chose to delete it rather than bundle a JavaScript
engine into this project's test dependencies: less machinery, honestly labelled. From here the shipped
evaluator is covered by **review**, not by tests -- a change to the evaluator, the presets, the reader
filters, ``gatesCard``, ``cartTsv``, ``encodeHash``/``applyHashFragment`` or ``parseControlValue`` has
to be read, because nothing will run it.

One half of the contract does survive. ``tests/unit/test_reporting_rethreshold.py`` still tests the
**Python** side of the re-threshold rule in 10 tests: a ``not_evaluated`` gate carrying a measured
value cannot be decided by any threshold, a warn gate moved past its floor stays ``warn``, an empty
value stays ``unknown``, a run-``unknown`` gate and an ``off`` gate are both frozen untouched, a
passing gate can still be re-thresholded into a fail. What is gone is specifically "**and the browser
does the same**".

Marked ``unit`` rather than ``integration``: ``tests/conftest.py`` promotes ``integration`` to the
release-only tier, and these must run in ``make test-dev`` on every iteration, not once at release.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from sirnaforge.core.hit_annotation import unclassified_cells
from sirnaforge.models.policy import DECLARED_FILTER_IDS
from sirnaforge.reporting import build_payload, render_html
from sirnaforge.reporting.render import _TEMPLATE

pytestmark = pytest.mark.unit

_META_COLUMNS = (
    "id",
    "guide_sequence",
    "passenger_sequence",
    "transcript_id",
    "position",
    "design_score",
    "composite_score",
    "weight_vector",
    "passes_filters",
)

#: Nine synthetic gates, none sharing a name with a real filter (#103's own DECLARED_FILTER_IDS), so a
#: coincidental match with production code is never mistaken for coverage of it. Comparators le/lt/
#: ge/gt and actions fail/warn/off are all represented, over columns each guide below sets up its own
#: reason code through -- REASON_OK (a real comparison), REASON_EMPTY_VALUE, REASON_FILTER_OFF,
#: REASON_NO_THRESHOLD, REASON_RUN_UNKNOWN and REASON_RUN_NOT_EVALUATED. REASON_MISSING_COLUMN is
#: deliberately not exercised here: the column a gate reads is missing for the whole run, not for one
#: guide, so it would put every guide's ``n_gates_unknown`` above zero.
_FILTER_SPECS = (
    {"filter_id": "ok_fail_le", "column": "metric_a", "comparator": "le", "threshold": 10, "action": "fail"},
    {"filter_id": "ok_warn_ge", "column": "metric_b", "comparator": "ge", "threshold": 5, "action": "warn"},
    {"filter_id": "ok_fail_lt", "column": "metric_f", "comparator": "lt", "threshold": 10, "action": "fail"},
    {"filter_id": "ok_warn_gt", "column": "metric_g", "comparator": "gt", "threshold": 5, "action": "warn"},
    {"filter_id": "empty_value_le", "column": "metric_c", "comparator": "le", "threshold": 3, "action": "fail"},
    {"filter_id": "filter_off_le", "column": "metric_off", "comparator": "le", "threshold": 2, "action": "off"},
    {
        "filter_id": "no_threshold_le",
        "column": "metric_nothresh",
        "comparator": "le",
        "threshold": None,
        "action": "fail",
    },
    {"filter_id": "run_unknown_gate", "column": "metric_d", "comparator": "le", "threshold": 100, "action": "fail"},
    {"filter_id": "run_not_evaluated_gate", "column": "metric_e", "comparator": "ge", "threshold": 1, "action": "fail"},
)

assert not {f["filter_id"] for f in _FILTER_SPECS} & set(DECLARED_FILTER_IDS), (
    "the fixture's synthetic filter ids must stay disjoint from the real registry"
)

_SCOPE = {"species": [], "max_mismatches": None, "hit_classes": []}


def _manifest_filters() -> list[dict[str, object]]:
    return [
        {
            **spec,
            "scope": _SCOPE,
            "stage": "design",
            "setting_key": spec["filter_id"],
            "definition": f"synthetic gate for report-structure testing ({spec['filter_id']})",
            "evidence_exported": True,
        }
        for spec in _FILTER_SPECS
    ]


#: guide label -> {column: value}. Blank/absent means the cell is empty in the CSV. Chosen so every
#: reason code in payload.py's _evaluate (187-303) fires on at least one guide, and so the guides with
#: no composite score, with real liabilities and with no screen at all are each present.
_GUIDES: dict[str, dict[str, object]] = {
    # ALL_PASS and WARN_ONLY are the two guides _HITS gives real liabilities to (4 and 2), so they are
    # marked screened: a guide with alignments in the hit table was plainly submitted to the aligner.
    "ALL_PASS": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
        "off_target_screened": True,
    },
    "WARN_ONLY": {
        "metric_a": 1,
        "metric_b": 2,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
        "off_target_screened": True,
    },
    "FAIL_ONE": {
        "metric_a": 15,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "FAIL",
    },
    "EMPTY_UNKNOWN": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": "",
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "FAIL",
    },
    "RUN_UNKNOWN_GATE": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "FAIL",
        "run_unknown_gate_verdict": "unknown",
    },
    "RUN_NOTEVAL_GATE": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
        "run_not_evaluated_gate_verdict": "not_evaluated",
    },
    "UNDECLARED_REJECT": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "REPEAT_ELEMENT",
    },
    "REGISTER_PAIR_A": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
    },
    "REGISTER_PAIR_B": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
    },
    "OFFTARGET_CLEAN": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
        "off_target_screened": True,
    },
    "NEVER_SCREENED": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
        "off_target_screened": False,
    },
    # Every gate decides normally; only the composite score is absent. design_score is kept, so the
    # guide still sorts and still renders.
    "NO_COMPOSITE": {
        "metric_a": 1,
        "metric_b": 8,
        "metric_f": 1,
        "metric_g": 10,
        "metric_c": 1,
        "metric_off": 5,
        "metric_nothresh": 5,
        "metric_d": 50,
        "metric_e": 5,
        "passes_filters": "PASS",
        "composite_score": "",
    },
}

#: id -> (transcript, position, composite_score). REGISTER_PAIR_A/B share a transcript 1 nt apart, so
#: payload.py's cross-guide register index (payload.py:733-745) links them; A scores higher and is
#: the cluster's representative.
_PLACEMENT = {
    "ALL_PASS": ("ENST00000000001", 100, 90.0),
    "WARN_ONLY": ("ENST00000000002", 100, 80.0),
    "FAIL_ONE": ("ENST00000000003", 100, 70.0),
    "EMPTY_UNKNOWN": ("ENST00000000004", 100, 60.0),
    "RUN_UNKNOWN_GATE": ("ENST00000000005", 100, 50.0),
    "RUN_NOTEVAL_GATE": ("ENST00000000006", 100, 40.0),
    "UNDECLARED_REJECT": ("ENST00000000007", 100, 30.0),
    "REGISTER_PAIR_A": ("ENST00000000008", 200, 95.0),
    "REGISTER_PAIR_B": ("ENST00000000008", 201, 20.0),
    "OFFTARGET_CLEAN": ("ENST00000000009", 100, 65.0),
    "NEVER_SCREENED": ("ENST00000000010", 100, 10.0),
    "NO_COMPOSITE": ("ENST00000000011", 100, 5.0),
}

#: (guide label, hit_class, symbol), one aggregated alignment each. ``off_target`` is a liability class
#: and ``on_target`` is not, so these decide each guide's ``liability_count``: ALL_PASS carries 4 and
#: WARN_ONLY 2. Every other guide stays at 0.
_HITS = (
    ("OFFTARGET_CLEAN", "on_target", "TP53"),
    *(("ALL_PASS", "off_target", "SOMEGENE") for _ in range(4)),
    *(("WARN_ONLY", "off_target", "OTHERGENE") for _ in range(2)),
)

_EXTRA_COLUMNS = (
    "metric_a",
    "metric_b",
    "metric_c",
    "metric_off",
    "metric_nothresh",
    "metric_d",
    "metric_e",
    "metric_f",
    "metric_g",
    "run_unknown_gate_verdict",
    "run_not_evaluated_gate_verdict",
    "off_target_screened",
)


def _seq(label: str) -> str:
    """A distinct, valid-looking 21-mer per label -- content is irrelevant, only distinctness is."""
    digest = abs(hash(label)) % (4**21)
    letters = "ACGT"
    return "".join(letters[(digest // 4**k) % 4] for k in range(21))


def _write_fixture_run(tmp_path: Path, gene_query: str = "TP53") -> Path:
    """A run exercising every reason code and every comparator.

    ``empty_value_le`` is blank on exactly one guide's row (REASON_EMPTY_VALUE); every other filter's
    column is real for every row, so the difference between the two is a fact about one guide, not
    about the run. The same discipline gives the rendered panel its edge cases: one guide with no
    composite score at all, and two guides carrying 4 and 2 real liabilities (``_HITS``).

    ``gene_query`` is a parameter only so one test can hand the renderer a string that spells
    ``</script>``; every other caller takes the default.
    """
    run = tmp_path / "run"
    (run / "sirnaforge").mkdir(parents=True)
    (run / "off_target" / "results" / "aggregated").mkdir(parents=True)

    columns = [*_META_COLUMNS, *_EXTRA_COLUMNS]
    lines = [",".join(columns)]
    for label, (transcript, position, score) in _PLACEMENT.items():
        seq = _seq(label)
        overrides = _GUIDES[label]
        cells: dict[str, object] = {
            "id": f"CAND_{label}",
            "guide_sequence": seq,
            "passenger_sequence": seq,
            "transcript_id": transcript,
            "position": position,
            "design_score": score,
            "composite_score": score,
            "weight_vector": "postscreen_sirna_v4",
        }
        # composite_score is overridable so one guide can have none at all (NO_COMPOSITE).
        for column in _EXTRA_COLUMNS + ("passes_filters", "composite_score"):
            if column in overrides:
                cells[column] = overrides[column]
        lines.append(",".join(str(cells.get(c, "")) for c in columns))
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join(lines) + "\n")

    manifest = {
        "gene_query": gene_query,
        "tool_version": "test",
        "run_policy": {
            "run_mode": "qualified",
            "profile": {"name": "document-structure-fixture"},
            "filters": _manifest_filters(),
        },
    }
    (run / "sirnaforge" / "manifest.json").write_text(json.dumps(manifest))

    base_columns = [
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
    hit_columns = [*base_columns, *unclassified_cells()]
    rows = []
    for index, (label, hit_class, symbol) in enumerate(_HITS):
        cells = unclassified_cells()
        cells.update(hit_class=hit_class, hit_symbol=symbol, matched_symbol=symbol)
        rows.append(
            [
                f"screen_{index}",
                _seq(label),
                "human",
                f"ENST0000009{index:04d}",
                str(100 + index),
                "+",
                "21M",
                "60",
                "42",
                "0",
                "0",
                "1.0",
                *(str(v) for v in cells.values()),
            ]
        )
    table = run / "off_target" / "results" / "aggregated" / "combined_offtargets.tsv"
    table.write_text("\n".join("\t".join(row) for row in [hit_columns, *rows]) + "\n")
    return run


@pytest.fixture(scope="module")
def fixture_run(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The portable run directory every assertion in this module is made over."""
    run = tmp_path_factory.mktemp("reporting_document_structure")
    return _write_fixture_run(run)


@pytest.fixture(scope="module")
def payload(fixture_run: Path):  # noqa: ANN201 - ReportPayload, kept unannotated to avoid the import
    """The payload built once for the module; every test reads it, none of them mutate it."""
    return build_payload(fixture_run)


_SCRIPT_OPEN, _SCRIPT_CLOSE = "<script>", "</script>"


def _slice_script(html: str) -> str:
    """The one ``<script>`` body, cut where a browser cuts it: at the FIRST ``</script>``.

    Stated as a rule rather than a regex, because a greedy ``<script>(.*)</script>`` runs to the LAST
    ``</script>`` in the document: a payload string spelling ``</script>`` would then be read as part
    of the script body, while the browser ends the script at the first occurrence and drops the rest of
    the document. The renderer escapes ``<`` in every embedded JSON blob (``render._embed``) so this can
    never happen; the count below is what keeps that true (#103, finding 7).
    """
    assert html.count(_SCRIPT_CLOSE) == 1, (
        "the rendered document closes its own <script> more than once: a payload string spells "
        "</script> and every browser would truncate the document there"
    )
    open_at = html.index(_SCRIPT_OPEN) + len(_SCRIPT_OPEN)
    body = html[open_at : html.index(_SCRIPT_CLOSE, open_at)]
    assert "function reevaluateGates(" in body, "the slice missed the script body the report ships"
    return body


@pytest.mark.unit
def test_the_fixture_covers_every_reason_code_the_evaluator_must_freeze_on(payload) -> None:  # noqa: ANN001
    """Guard the comparison against becoming vacuous: every reason in payload.py must fire at least once."""
    seen = {reason for g in payload.guides for _, _, reason in g.gates}
    assert seen == {0, 2, 3, 4, 5, 6}, f"the fixture stopped exercising a reason code: {sorted(seen)}"


@pytest.mark.unit
def test_no_declared_filter_id_or_column_is_a_branch_in_the_template() -> None:
    """The evaluator is generic (#103): it reads FILTERS, it never branches on a specific filter.

    A static grep over the six hand-written v1 metric boxes is what the issue's own review caught --
    ``gc_content``, ``asymmetry_score`` and ``off_target_count`` were literal strings in the old
    template, and those three are gates a reader now moves generically. None of the real registry's 17
    ids may appear as a literal here. The three restored reader filters are not a counterexample: they
    read guide fields no gate covers, under their own keys, and no `filter_id` among them.
    """
    for filter_id in DECLARED_FILTER_IDS:
        assert f"'{filter_id}'" not in _TEMPLATE and f'"{filter_id}"' not in _TEMPLATE, (
            f"{filter_id} is hardcoded in the template; the evaluator must read it from FILTERS instead"
        )


@pytest.mark.unit
def test_the_panel_separates_gate_controls_from_reader_filters(payload) -> None:  # noqa: ANN001
    """A gate control changes a verdict; a reader filter only selects rows. The UI must say which.

    Structural, because the distinction is what the removal of these boxes cost and what restoring them
    has to preserve: two labelled groups, three reader inputs in the second one, and the search box
    still applied on top of ``passesFilters`` rather than beside it.
    """
    html = render_html(payload)

    assert 'id="frows"' in html and 'id="rrows"' in html, "the two groups are separate containers"
    assert "Gate thresholds" in html and "re-derives the verdict this run computed" in html
    assert "Reader filters" in html and "only select among rows; no verdict changes" in html
    for key in ("composite", "isoforms", "liab"):
        assert f"{{k:'{key}'," in html, f"the {key} reader filter is missing"
    assert "if(!passesReaderFilters(g)) return false;" in html, "they must compose inside passesFilters"
    assert "view = matching().filter(" in html, "and the search box applies on top of that, not beside it"


@pytest.mark.unit
def test_every_threshold_box_refuses_non_numeric_input_at_the_control(payload) -> None:  # noqa: ANN001
    """The DOM half of the numeric rule: the boxes are typed, and both handlers refuse what fnum refuses.

    Structural because the handlers read ``validity.badInput``, which only a browser sets, so no test in
    this repository can drive them. ``type="number"`` is what makes that flag exist; ``inputmode`` never
    constrained anything.
    """
    html = render_html(payload)

    assert 'id="ctl_${esc(f.filter_id)}" type="number" step="any"' in html, "the gate box is typed"
    assert 'id="rf_${esc(f.k)}" type="number" step="any"' in html, "the reader box is typed"
    assert "inputmode" not in html, "inputmode is a keyboard hint, not a constraint, and it read as one"
    assert html.count("if(v===undefined){ markControl(el,false); return; }") == 2, (
        "both handlers -- the gate thresholds and the reader bounds -- must refuse it, not just one"
    )
    assert 'id="badnum"' in html and "not a finite number, so nothing moved" in html, "and it says so"
    assert "el.validity.badInput" in html, "the only trace a type=number box leaves of unparseable text"
    # No min/max on the box, deliberately: payload._control_domain bounds the slider, and the number box
    # is what keeps every threshold reachable. What is refused is not an unusual number but a non-number.
    assert 'type="number" step="any" value="${T[f.filter_id]}"' in html


@pytest.mark.unit
def test_the_script_slice_ends_where_a_browser_would() -> None:
    """A document that closes its own script twice must be refused (#103, finding 7).

    A greedy ``<script>(.*)</script>`` slice runs to the LAST ``</script>``, so the forged document
    below reads as one long script that parses cleanly, while a browser stops at the first close tag and
    drops everything after it. ``_slice_script`` refuses it instead, which is what lets the test below
    read the embedded literal out of a document a browser would actually have executed whole.
    """
    forged = '<!DOCTYPE html><body><script>const G=["</script>"];\nrun()</script></body></html>'
    assert re.search(r"<script>(.*)</script>", forged, re.DOTALL).group(1).count("</script>") == 1, (
        "the greedy regex swallowed the injected close tag and would have read it as live script"
    )
    with pytest.raises(AssertionError, match="closes its own"):
        _slice_script(forged)


@pytest.mark.unit
def test_a_payload_spelling_a_script_close_cannot_truncate_the_document(tmp_path: Path) -> None:
    """A gene query, transcript id or symbol may contain anything; the document must survive it.

    ``json.dumps`` does not escape ``<``, and no layer between a run's own metadata and the rendered
    file sanitises it, so a ``</script>`` anywhere in the payload ended the script block early and took
    every panel with it. The renderer escapes it (``render._embed``); this asserts the document still has
    exactly one script, and that the escape is **reversible** -- what the renderer writes is an ordinary
    JSON unicode escape, so parsing the embedded literal back yields the original string, which is what
    the browser's own parser would do with it. That is the step this test used to take through node.
    """
    hostile = "TP53</script><script>run()</script>"
    run = _write_fixture_run(tmp_path, gene_query=hostile)
    html = render_html(build_payload(run))

    assert html.count("</script>") == 1, "a payload string closed the script block"
    assert "\\u003c/script\\u003e" in html, "the close tag must be escaped, not stripped or dropped"
    assert hostile not in html, "and it must not appear unescaped anywhere, markup included"

    literal = re.search(r"^const GENE = (.*);$", _slice_script(html), re.MULTILINE)
    assert literal is not None, "the embedded GENE literal moved; this test can no longer read it"
    assert json.loads(literal.group(1)) == hostile, "escaping must be reversible: the export filename reads GENE"
