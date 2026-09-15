"""Parity between the report's shipped JS gate evaluator and the Python filter (issue #103).

The report now lets a reader move a gate threshold and re-reads every gate in the browser. If that
client-side evaluator and the pipeline's own gate application ever disagree, the report can publish a
verdict the pipeline would not reproduce -- exactly the fabricated-evidence failure the descriptor-
driven design exists to prevent (payload.py's module docstring).

No Playwright: it is not a dependency of this repo and its browser binary is a separate install, which
would make this a skip rather than a gate. Node is not a new tool here either -- ``.nvmrc`` pins 20 and
CI's ``ubuntu-latest`` ships it -- so the harness below does the one thing Playwright would otherwise
be needed for (real placeholder substitution) in Python, and hands the resulting real JS text to node
directly: :func:`render_html` runs first, the ``<script>`` body is sliced out of the rendered
document, and a small driver appended to it drives the evaluator the file actually ships.

Marked ``unit`` rather than ``integration``: ``tests/conftest.py`` promotes ``integration`` to the
release-only tier, and this must run in ``make test-dev`` on every iteration, not once at release.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.core.hit_annotation import unclassified_cells
from sirnaforge.models.policy import DECLARED_FILTER_IDS
from sirnaforge.reporting import build_payload, render_html
from sirnaforge.reporting.payload import REASON_OK, reevaluate_gates
from sirnaforge.reporting.render import _TEMPLATE

pytestmark = pytest.mark.unit

NODE_INSTALL_HINT = "install node >=20 (.nvmrc pins 20); ubuntu-latest and this dev box both ship it"

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
#: guide, so it would put every guide's ``n_gates_unknown`` above zero and make near-miss -- which
#: requires ``n_gates_unknown===0`` -- unreachable by construction. The evaluator's freeze check does
#: not distinguish reason codes (it is one ``reason !== REASON_OK`` branch for all six), so exercising
#: the other five already covers the code path this one would take too.
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
            "definition": f"synthetic gate for parity testing ({spec['filter_id']})",
            "evidence_exported": True,
        }
        for spec in _FILTER_SPECS
    ]


#: guide label -> {column: value}. Blank/absent means the cell is empty in the CSV. Chosen so every
#: reason code in payload.py's _evaluate (187-303) fires on at least one guide, and so near-miss,
#: off-target-clean and register-deduplicated each have a genuine positive and a genuine trap case.
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
    # Every gate decides normally; only the composite score is absent. This is the guide a composite
    # floor must EXCLUDE rather than admit -- the reader-filter half of "an absent value cannot
    # satisfy a threshold". design_score is kept, so the guide still sorts and still renders.
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
#: payload.py's cross-guide register index (payload.py:733-745) links them; A scores higher and must be
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
#: WARN_ONLY 2, which is what makes "at most 3 liabilities" a request with a real answer -- and one the
#: off-target-clean preset, which covers exactly zero, cannot express. Every other guide stays at 0.
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


def _guide_of(label: str) -> str:
    """The payload's key for a fixture guide.

    ``payload._normalise_guide`` upper-cases and spells the guide in RNA (T -> U); the fixture's CSV
    rows store the DNA spelling, like every other ``candidates_all.csv``, so a lookup key must too.
    """
    return _seq(label).replace("T", "U")


def _write_fixture_run(tmp_path: Path, gene_query: str = "TP53") -> Path:
    """A run exercising every reason code, every comparator, and the presets' and filters' edge cases.

    ``empty_value_le`` is blank on exactly one guide's row (REASON_EMPTY_VALUE); every other filter's
    column is real for every row, so the difference between the two is a fact about one guide, not
    about the run. The same discipline gives the reader filters their edge cases: one guide with no
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
        "run_policy": {"run_mode": "qualified", "profile": {"name": "parity-fixture"}, "filters": _manifest_filters()},
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
    """The portable run directory every comparison in this module is made over."""
    run = tmp_path_factory.mktemp("client_evaluator_parity")
    return _write_fixture_run(run)


@pytest.fixture(scope="module")
def payload(fixture_run: Path):  # noqa: ANN201 - ReportPayload, kept unannotated to avoid the import
    """The payload built once for the module; every test reads it, none of them mutate it."""
    return build_payload(fixture_run)


_SCRIPT_OPEN, _SCRIPT_CLOSE = "<script>", "</script>"


def _slice_script(html: str) -> str:
    """The one ``<script>`` body, cut where a browser cuts it: at the FIRST ``</script>``.

    This slice is the harness's own trustworthiness, so it is stated as a rule rather than a regex. A
    greedy ``<script>(.*)</script>`` ran to the LAST ``</script>`` in the document, so a payload string
    spelling ``</script>`` was handed to node whole -- inside a JavaScript string literal it parses and
    runs perfectly -- while the browser ended the script at the first occurrence and dropped the rest of
    the document. Parity then passed on a page no browser could run, and a gate that can pass on a
    broken document is not a gate (#103, finding 7). The renderer escapes ``<`` in every embedded JSON
    blob (``render._embed``) so this can never happen; the count below is what keeps that true.
    """
    assert html.count(_SCRIPT_CLOSE) == 1, (
        "the rendered document closes its own <script> more than once: a payload string spells "
        "</script> and every browser would truncate the document there"
    )
    open_at = html.index(_SCRIPT_OPEN) + len(_SCRIPT_OPEN)
    body = html[open_at : html.index(_SCRIPT_CLOSE, open_at)]
    assert "function reevaluateGates(" in body, "the slice missed the script body the report ships"
    return body


@pytest.fixture(scope="module")
def rendered_script(payload) -> str:  # noqa: ANN001
    """The rendered report's ``<script>`` body -- placeholder substitution has already happened."""
    return _slice_script(render_html(payload))


def _node_or_fail() -> str:
    """The ``node`` binary, or a loud failure -- a skip here is not a gate (#103's own instruction)."""
    node = shutil.which("node")
    if node is None:
        pytest.fail(f"node is required for this parity test and was not found on PATH; {NODE_INSTALL_HINT}")
    return node


def _python_reevaluate_gates(guide, filters: list[dict[str, object]], overrides: dict[str, float]) -> dict[str, object]:
    """The Python side of the parity contract: re-threshold one guide exactly as the JS evaluator does.

    Delegates the per-gate freeze/re-decide rule to :func:`sirnaforge.reporting.payload.reevaluate_gates`
    itself rather than a second copy of it -- that function *is* the contract the shipped JS has to
    agree with, so this test compares against it directly. Only the guide-level aggregation
    (status/n_gates_*) is this file's own, and it restates :class:`GuideEntry`'s own logic
    (payload.py:87-134), not a new rule.
    """
    gates = reevaluate_gates(filters, guide.gates, overrides)
    n_failed = sum(1 for _, v, _ in gates if v == 1)
    n_unknown = sum(1 for _, v, _ in gates if v == 2)
    n_warned = sum(1 for _, v, _ in gates if v == 4)
    undeclared = guide.run_verdict not in (None, "PASS") and not (n_failed or n_unknown)
    if n_failed:
        status = "fail"
    elif n_unknown or undeclared:
        status = "unknown"
    elif n_warned:
        status = "warn"
    else:
        status = "pass"
    return {
        "gates": gates,
        "n_gates_failed": n_failed,
        "n_gates_unknown": n_unknown,
        "n_gates_warned": n_warned,
        "status": status,
    }


def _override_sets(payload) -> dict[str, dict[str, float]]:  # noqa: ANN001
    """The run's own thresholds, and four sets that move every movable control that exists."""
    evaluable = [f for f in payload.filters if f["action"] != "off" and f["threshold"] is not None]
    frozen = next(f for f in payload.filters if not (f["action"] != "off" and f["threshold"] is not None))
    return {
        "the run's own": {},
        "every evaluable control at zero": {f["filter_id"]: 0 for f in evaluable},
        "every evaluable control one below its threshold": {f["filter_id"]: f["threshold"] - 1 for f in evaluable},
        "every evaluable control far above its threshold": {f["filter_id"]: f["threshold"] + 1000 for f in evaluable},
        "a frozen filter and an unknown id asked for too": {
            **{f["filter_id"]: f["threshold"] for f in evaluable},
            frozen["filter_id"]: 999999,
            "no_such_filter_id": 1,
        },
    }


def _run_node_driver(
    script: str, override_sets: dict[str, dict[str, float]], frozen_filter_id: str
) -> dict[str, object]:
    """Append a driver to the shipped script and run it under node, returning its parsed stdout.

    The driver only ever calls functions the report ships (``reevaluateGates``, ``PRESETS``,
    ``applyHashFragment``); it adds no logic of its own that the parity check could be fooled by.
    """
    driver = f"""
const __OVERRIDE_SETS = {json.dumps(override_sets)};
const __RESULTS = [];
for (const [label, overrides] of Object.entries(__OVERRIDE_SETS)) {{
  for (const g of G) {{
    const live = reevaluateGates(g, overrides);
    __RESULTS.push({{label, guide: g.guide, overrides, live, passing: PRESETS.passing.test(g, live)}});
  }}
}}
const __PRESETS = G.map(g => ({{
  guide: g.guide,
  status: g._live.status,
  near_miss: PRESETS.near_miss.test(g, g._live),
  off_target_clean: PRESETS.off_target_clean.test(g, g._live),
  register_dedup: PRESETS.register_dedup.test(g, g._live),
}}));
const __REFUSAL = applyHashFragment('t=' + {json.dumps(frozen_filter_id)} + ':999999,no_such_filter_id:1');
process.stdout.write(JSON.stringify({{results: __RESULTS, presets: __PRESETS, refusal: __REFUSAL}}));
"""
    return _node_run(script, driver)


def _node_run(script: str, driver: str) -> dict[str, object]:
    """Run the shipped script plus a driver under node and return its parsed stdout.

    A temp ``.mjs``, not ``node -e``: the substituted script is ~60 KB and this is the shape the issue's
    own instructions ask for -- the placeholder substitution happens for real, in Python, before node
    ever sees a byte of it.
    """
    node = _node_or_fail()
    with tempfile.TemporaryDirectory() as tmp:
        driver_path = Path(tmp) / "evaluator.mjs"
        driver_path.write_text(script + driver, encoding="utf-8")
        completed = subprocess.run(
            [node, str(driver_path)],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
    assert completed.returncode == 0, f"the shipped evaluator raised under node:\n{completed.stderr}"
    return json.loads(completed.stdout)


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
def test_the_shipped_evaluator_agrees_with_the_python_filter_on_every_gate(
    payload,
    rendered_script: str,  # noqa: ANN001
) -> None:
    """Every guide x every filter, at five threshold sets. Not a sample.

    Also the property test the issue asks for: no threshold this drives can move a gate whose reason
    is not OK, and no preset can promote an ``unknown`` guide into the ``passing`` view.
    """
    overrides_sets = _override_sets(payload)
    evaluable = [f for f in payload.filters if f["action"] != "off" and f["threshold"] is not None]
    frozen = next(f for f in payload.filters if not (f["action"] != "off" and f["threshold"] is not None))
    output = _run_node_driver(rendered_script, overrides_sets, frozen["filter_id"])

    by_guide = {g.guide: g for g in payload.guides}
    compared = 0
    for entry in output["results"]:
        guide = by_guide[entry["guide"]]
        expected = _python_reevaluate_gates(guide, payload.filters, entry["overrides"])
        assert entry["live"]["gates"] == expected["gates"], (entry["label"], entry["guide"])
        assert entry["live"]["n_gates_failed"] == expected["n_gates_failed"], (entry["label"], entry["guide"])
        assert entry["live"]["n_gates_unknown"] == expected["n_gates_unknown"], (entry["label"], entry["guide"])
        assert entry["live"]["n_gates_warned"] == expected["n_gates_warned"], (entry["label"], entry["guide"])
        assert entry["live"]["status"] == expected["status"], (entry["label"], entry["guide"])

        # Property: a frozen reason never moves, at any threshold this evaluator was handed.
        for i, (f, triple) in enumerate(zip(payload.filters, guide.gates, strict=True)):
            _value, _verdict, reason = triple
            if reason != REASON_OK:
                assert entry["live"]["gates"][i] == triple, (entry["label"], entry["guide"], f["filter_id"])

        # Property: the 'passing' preset never includes a guide this evaluator itself calls unknown.
        assert not (entry["live"]["status"] == "unknown" and entry["passing"]), (entry["label"], entry["guide"])
        compared += 1

    assert compared == len(overrides_sets) * len(payload.guides), compared
    assert len(evaluable) >= 6, "the fixture must exercise more than a couple of movable gates"


@pytest.mark.unit
def test_near_miss_off_target_clean_and_register_dedup_presets(payload, rendered_script: str) -> None:  # noqa: ANN001
    """The three presets that are not a plain status filter, each pinned to its trap case.

    ``FAIL_ONE`` fails exactly one gate and is unknown on none -- a genuine near miss.
    ``OFFTARGET_CLEAN`` carries ``off_target_screened=True`` and a real liability_count of 0;
    ``NEVER_SCREENED`` has the same liability_count==0 but ``off_target_screened=False``, which is the
    trap ``liability_count==0`` alone would fall into. ``REGISTER_PAIR_A`` outscores its 1-nt neighbour
    ``REGISTER_PAIR_B`` and must be the cluster's sole representative.
    """
    output = _run_node_driver(rendered_script, {"the run's own": {}}, DECLARED_FILTER_IDS[0])
    by_guide = {row["guide"]: row for row in output["presets"]}
    guide_of = {label: _guide_of(label) for label in _GUIDES}

    assert by_guide[guide_of["FAIL_ONE"]]["near_miss"] is True
    assert by_guide[guide_of["ALL_PASS"]]["near_miss"] is False
    assert by_guide[guide_of["EMPTY_UNKNOWN"]]["near_miss"] is False, "unknown can never read as near-miss"

    assert by_guide[guide_of["OFFTARGET_CLEAN"]]["off_target_clean"] is True
    assert by_guide[guide_of["NEVER_SCREENED"]]["off_target_clean"] is False, (
        "liability_count==0 alone must not be read as clean -- this guide was simply never screened"
    )

    assert by_guide[guide_of["REGISTER_PAIR_A"]]["register_dedup"] is True
    assert by_guide[guide_of["REGISTER_PAIR_B"]]["register_dedup"] is False


@pytest.mark.unit
def test_a_frozen_filter_named_in_the_url_fragment_is_refused_visibly(payload, rendered_script: str) -> None:  # noqa: ANN001
    """The URL fragment can name a frozen gate or a nonexistent one; both must be listed, not dropped."""
    frozen = next(f for f in payload.filters if not (f["action"] != "off" and f["threshold"] is not None))
    output = _run_node_driver(rendered_script, {"the run's own": {}}, frozen["filter_id"])
    refused = output["refusal"]["refused"]
    assert frozen["filter_id"] in refused
    assert "no_such_filter_id" in refused


# --- reader filters: composite score, isoforms hit, liabilities (#103 remainder) --------------------
# Three of the v1 report's six metric boxes are restored beside the gate controls, because no gate
# covers them: composite_score and transcript_hits are gated by nothing, and liability_count is
# reachable only as "exactly zero" through the off-target-clean preset. They select rows; they never
# move a verdict. The bounds below are chosen against _PLACEMENT's composite scores (5..95).
_READER_BOUND_SETS = {
    "nothing set": {},
    "composite at least 60": {"composite": 60},
    "isoforms at least 1": {"isoforms": 1},
    "at most 0 liabilities": {"liab": 0},
    "at most 3 liabilities": {"liab": 3},
    "composite floor and liability ceiling together": {"composite": 30, "liab": 3},
}

#: reader filter key -> (direction, how to read the value off a payload guide). The Python side of the
#: reader-filter contract, stated once; the JS reads its own READER_FILTERS table.
_READER_ACCESS: dict[str, tuple[str, Callable[[Any], float | None]]] = {
    "composite": ("min", lambda g: g.composite_score),
    "isoforms": ("min", lambda g: g.transcript_hits),
    "liab": ("max", lambda g: g.liability_count),
}


def _python_reader_selection(payload, bounds: dict[str, float]) -> list[str]:  # noqa: ANN001
    """The guides a reader filter set admits, by the same rule: an absent value satisfies nothing."""
    kept = []
    for guide in payload.guides:
        ok = True
        for key, bound in bounds.items():
            direction, read = _READER_ACCESS[key]
            value = read(guide)
            if value is None or (value < bound if direction == "min" else value > bound):
                ok = False
                break
        if ok:
            kept.append(guide.guide)
    return kept


def _run_reader_filter_driver(script: str, bound_sets: dict[str, dict[str, float]]) -> dict[str, object]:
    """Drive the shipped reader filters, their composition, and a fragment round-trip, under node.

    Every call is to a function or a binding the report itself ships -- ``passesReaderFilters``,
    ``passesFilters``, ``PRESETS``, ``encodeHash``, ``applyHashFragment``, ``resetControls``, ``R``,
    ``T``, ``activePreset`` and ``F``. The only logic here is the bookkeeping that records what they
    returned; the reset below is the same ``resetControls`` the ``#freset`` button calls, not a
    restatement of it.
    """
    driver = f"""
const __BOUND_SETS = {json.dumps(bound_sets)};
const __SELECTION = Object.entries(__BOUND_SETS).map(([label, bounds]) => ({{
  label, bounds, guides: G.filter(g => passesReaderFilters(g, bounds)).map(g => g.guide),
}}));

// Composition: a reader bound, a preset and the status checkboxes at once must be exactly the
// intersection of the three taken separately -- driven through passesFilters itself, not re-derived.
const __COMPOSE = (() => {{
  const bounds = {{composite: 30}}, preset = 'register_dedup', statuses = ['pass'];
  const readerOnly = G.filter(g => passesReaderFilters(g, bounds)).map(g => g.guide);
  const presetOnly = G.filter(g => PRESETS[preset].test(g, g._live)).map(g => g.guide);
  const statusOnly = G.filter(g => statuses.includes(g._live.status)).map(g => g.guide);
  Object.assign(R, bounds); activePreset = preset; F.status = new Set(statuses);
  const combined = G.filter(passesFilters).map(g => g.guide);
  resetControls();
  return {{readerOnly, presetOnly, statusOnly, combined}};
}})();

// Fragment round-trip and reset. The reset is the #freset button's own resetControls; the reload is
// the fragment alone, as a fresh load of the copied URL would be.
const __TRIP = (() => {{
  const gate = FILTERS.find(f => evaluable(f));
  const moved = gate.threshold + 1;
  T[gate.filter_id] = moved; R.composite = 30; R.liab = 3; activePreset = 'near_miss';
  recomputeLive();
  const before = {{
    hash: encodeHash(), T: {{...T}}, R: {{...R}}, preset: activePreset,
    rows: G.filter(passesFilters).map(g => g.guide),
  }};
  resetControls();
  recomputeLive();
  const reset = {{
    hash: encodeHash(), T: {{...T}}, R: {{...R}}, preset: activePreset,
    rows: G.filter(passesFilters).map(g => g.guide),
  }};
  const applied = applyHashFragment(before.hash);
  recomputeLive();
  const reloaded = {{
    hash: encodeHash(), T: {{...T}}, R: {{...R}}, preset: activePreset, refused: applied.refused,
    rows: G.filter(passesFilters).map(g => g.guide),
  }};
  resetControls(); recomputeLive();
  return {{gate: gate.filter_id, moved, before, reset, reloaded}};
}})();

// An unknown reader-filter name in the fragment is refused, exactly as an unknown gate id is.
const __UNKNOWN = (() => {{
  const applied = applyHashFragment('r=no_such_reader_filter:7');
  return {{refused: applied.refused, R: {{...R}}}};
}})();

process.stdout.write(JSON.stringify(
  {{selection: __SELECTION, compose: __COMPOSE, trip: __TRIP, unknown: __UNKNOWN}}));
"""
    return _node_run(script, driver)


@pytest.fixture(scope="module")
def reader_filter_output(rendered_script: str) -> dict[str, object]:
    """One node invocation shared by the reader-filter assertions below."""
    return _run_reader_filter_driver(rendered_script, _READER_BOUND_SETS)


@pytest.mark.unit
def test_each_reader_filter_selects_the_rows_the_python_side_would(payload, reader_filter_output) -> None:  # noqa: ANN001
    """Composite floor, isoform floor and liability ceiling, each against the payload's own numbers.

    ``at most 3 liabilities`` is the request the off-target-clean preset cannot express: that preset
    covers exactly zero, so without this box a reader has no way to ask for a small, tolerable number.
    """
    by_label = {entry["label"]: entry for entry in reader_filter_output["selection"]}
    assert set(by_label) == set(_READER_BOUND_SETS), "the driver dropped a bound set"

    for label, bounds in _READER_BOUND_SETS.items():
        expected = _python_reader_selection(payload, bounds)
        assert by_label[label]["guides"] == expected, label

    # Guard against a vacuous comparison: each bound set must actually cut something, or nothing above
    # would notice a filter that had stopped being applied at all.
    everything = by_label["nothing set"]["guides"]
    assert everything == [g.guide for g in payload.guides]
    for label in ("composite at least 60", "at most 0 liabilities", "at most 3 liabilities"):
        assert 0 < len(by_label[label]["guides"]) < len(everything), label

    # And the ceiling is a real ceiling, not just "zero or everything": ALL_PASS carries 4 liabilities
    # and WARN_ONLY 2, so 3 admits one of them and refuses the other.
    at_most_3 = by_label["at most 3 liabilities"]["guides"]
    assert _guide_of("WARN_ONLY") in at_most_3 and _guide_of("ALL_PASS") not in at_most_3
    at_most_0 = by_label["at most 0 liabilities"]["guides"]
    assert _guide_of("WARN_ONLY") not in at_most_0 and _guide_of("ALL_PASS") not in at_most_0


@pytest.mark.unit
def test_a_guide_with_no_composite_score_is_excluded_by_a_composite_floor(payload, reader_filter_output) -> None:  # noqa: ANN001
    """An absent value cannot satisfy a threshold -- the rule the gate evaluator already holds to.

    Admitting ``NO_COMPOSITE`` under a floor it was never measured against would publish it as having
    cleared a bar nobody applied to it, which is the fabricated-evidence direction.
    """
    no_composite = _guide_of("NO_COMPOSITE")
    absent = next(g for g in payload.guides if g.guide == no_composite)
    assert absent.composite_score is None, "the fixture stopped exercising an absent composite score"

    by_label = {entry["label"]: entry for entry in reader_filter_output["selection"]}
    assert no_composite in by_label["nothing set"]["guides"], "it is a real row when nothing is asked"
    assert no_composite not in by_label["composite at least 60"]["guides"]
    assert no_composite not in by_label["composite floor and liability ceiling together"]["guides"]
    # Its liability count is a real 0, so a ceiling that says nothing about the composite keeps it.
    assert no_composite in by_label["at most 0 liabilities"]["guides"]


@pytest.mark.unit
def test_a_reader_filter_composes_with_the_presets_and_the_status_checkboxes(reader_filter_output) -> None:  # noqa: ANN001
    """All three narrow together through ``passesFilters``; none of them overrides another.

    A composite floor, the register-deduplicated preset and a ``pass``-only status set, chosen so each
    of the three cuts rows the other two keep -- so ``combined`` being strictly smaller than every one
    of them is evidence that all three were actually applied, not just the last one wired in.
    """
    compose = reader_filter_output["compose"]
    intersection = [
        guide
        for guide in compose["readerOnly"]
        if guide in set(compose["presetOnly"]) and guide in set(compose["statusOnly"])
    ]
    assert compose["combined"] == intersection
    assert compose["combined"], "the fixture must leave at least one row, or this proves nothing"
    for key in ("readerOnly", "presetOnly", "statusOnly"):
        assert len(compose["combined"]) < len(compose[key]), f"{key} alone already decided the view"


@pytest.mark.unit
def test_the_reader_filters_survive_the_url_fragment_and_the_reset_clears_them(reader_filter_output) -> None:  # noqa: ANN001
    """A copied URL restores the same rows; ``Reset`` puts every control back, gates and filters alike."""
    trip = reader_filter_output["trip"]
    assert "r=composite:30,liab:3" in trip["before"]["hash"], trip["before"]["hash"]
    assert f"t={trip['gate']}:{trip['moved']}" in trip["before"]["hash"]
    assert "preset=near_miss" in trip["before"]["hash"]

    # Reset: no r=, no t=, no preset= left, and the rows genuinely change -- a reset that cleared
    # nothing would satisfy every assertion below it.
    assert trip["reset"]["hash"] == ""
    assert all(bound is None for bound in trip["reset"]["R"].values())
    assert trip["reset"]["preset"] == "all"
    assert trip["reset"]["rows"] != trip["before"]["rows"]

    # Reload from the fragment alone: same bounds, same thresholds, same preset, same rows.
    assert trip["reloaded"]["refused"] == []
    assert trip["reloaded"]["R"] == {"composite": 30, "isoforms": None, "liab": 3}
    assert trip["reloaded"]["T"][trip["gate"]] == trip["moved"]
    assert trip["reloaded"]["preset"] == "near_miss"
    assert trip["reloaded"]["rows"] == trip["before"]["rows"]
    assert trip["reloaded"]["hash"] == trip["before"]["hash"], "the fragment must round-trip byte for byte"


@pytest.mark.unit
def test_an_unknown_reader_filter_name_in_the_fragment_is_refused_not_invented(reader_filter_output) -> None:  # noqa: ANN001
    """The same treatment a frozen gate id gets: listed as refused, and nothing set behind the reader."""
    unknown = reader_filter_output["unknown"]
    assert unknown["refused"] == ["no_such_reader_filter"]
    assert all(bound is None for bound in unknown["R"].values())


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
def test_node_absent_fails_loudly_rather_than_skipping(monkeypatch: pytest.MonkeyPatch) -> None:
    """A skip is not a gate: without node this must fail with the install hint, never pass quietly."""
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    with pytest.raises(pytest.fail.Exception, match="node"):
        _node_or_fail()


# --- the gate panel's Why column, the cart's provenance, and the numeric boundary (#103 remainder) ---
#: The fixture gate the panel and cart assertions move, and a guide the move flips. ``ok_fail_le`` reads
#: ``metric_a`` with ``le 10``; ``FAIL_ONE`` measures 15, so it fails at the run's threshold and passes
#: that gate at 20. This is the reviewer's own reproduction, kept as the fixture case.
_MOVED_GATE = "ok_fail_le"
_MOVED_TO = 20
_FLIPPED_GUIDE = "FAIL_ONE"

#: The cart export's columns and their order, frozen. They are the file's contract with every ticket and
#: every script that has ever read one, so adding threshold provenance above them must not disturb them;
#: a reorder or a rename has to fail here rather than in someone's spreadsheet.
_CART_COLUMNS = (
    "guide",
    "passenger",
    "status",
    "composite_score",
    "design_score",
    "isoforms_hit",
    "isoforms_in_run",
    "gc_content",
    "asymmetry_score",
    "off_target_count",
    "liabilities",
    "structure",
    "mouse_nm",
    "mouse_seed_mm",
    "macaque_nm",
    "macaque_seed_mm",
    "rat_nm",
    "rat_seed_mm",
)

#: One rendered gate row: the filter id, the verdict on its pill, the value, the whole threshold cell and
#: the Why sentence. Parsed out of the panel the report actually builds rather than re-derived, because
#: the defect was a disagreement *between* two cells of this row.
_GATE_ROW = re.compile(
    r'<tr><td class="mono">(?P<filter_id>[^<]+)</td>\s*'
    r'<td><span class="pill v-(?P<verdict>[a-z_]+)">[^<]*</span></td>\s*'
    r"<td>(?P<value>[^<]*)</td>"
    r'<td class="mono">(?P<threshold_cell>.*?)</td>\s*'
    r"<td>[^<]*</td><td>[^<]*</td><td>(?P<why>[^<]*)</td></tr>",
    re.DOTALL,
)

#: The comparison a Why sentence opens with, when the gate was actually compared: value, polarity,
#: comparator, threshold. ``fmt`` groups thousands, so the digits are read with separators stripped.
_WHY_COMPARISON = re.compile(
    r"^(?P<value>-?[\d,.]+) (?P<negated>not )?(?P<comparator>le|lt|ge|gt) (?P<threshold>-?[\d,.]+)"
)

_COMPARE = {
    "le": lambda v, t: v <= t,
    "lt": lambda v, t: v < t,
    "ge": lambda v, t: v >= t,
    "gt": lambda v, t: v > t,
}


def _num(text: str) -> float:
    return float(text.replace(",", ""))


def _run_panel_driver(script: str, guide: str) -> dict[str, Any]:
    """Build the gate panel for one guide, at the run's thresholds and with one gate loosened.

    ``gatesCard`` touches no DOM -- it returns the panel's HTML as a string -- so the harness can drive
    the panel the report actually renders rather than a restatement of it, which is the only way to
    catch two cells of one row disagreeing.
    """
    driver = f"""
const __GUIDE = {json.dumps(guide)}, __GATE = {json.dumps(_MOVED_GATE)}, __TO = {json.dumps(_MOVED_TO)};
const __PANEL = (() => {{
  const g = G.find(x => x.guide === __GUIDE);
  const atRun = {{html: gatesCard(g), thresholds: {{...T}}}};
  T[__GATE] = __TO;
  recomputeLive();
  const atMoved = {{html: gatesCard(g), thresholds: {{...T}}}};
  resetControls(); recomputeLive();
  return {{atRun, atMoved}};
}})();
process.stdout.write(JSON.stringify(__PANEL));
"""
    return _node_run(script, driver)


def _panel_rows(html: str) -> dict[str, dict[str, str]]:
    rows = {m.group("filter_id"): m.groupdict() for m in _GATE_ROW.finditer(html)}
    assert rows, "no gate row parsed out of the panel; the row markup changed shape"
    return rows


@pytest.mark.unit
def test_the_gate_panel_states_the_comparison_that_produced_the_verdict_beside_it(
    payload,  # noqa: ANN001
    rendered_script: str,
) -> None:
    """The Why sentence and the pill next to it must come from ONE threshold (#103, finding 1).

    The Why column built its sentence from ``f.threshold`` -- the run's -- while the verdict beside it
    came from the live re-thresholded table, so any reader move made the sentence arithmetically false
    in both directions: "10 not le 10 so FAIL" beside a PASS pill, or a satisfied comparison beside a
    FAIL. A report whose stated arithmetic does not produce its own verdict is publishing a verdict the
    pipeline would not reproduce, which is the failure this whole report exists to prevent.

    Checked as a property over every compared row of both panels -- the threshold the sentence cites is
    the threshold in force, and the sentence's polarity agrees with the pill -- not just on the row that
    was moved.
    """
    panel = _run_panel_driver(rendered_script, _guide_of(_FLIPPED_GUIDE))
    comparators = {f["filter_id"]: f["comparator"] for f in payload.filters}
    run_thresholds = {f["filter_id"]: f["threshold"] for f in payload.filters}

    compared = 0
    for label in ("atRun", "atMoved"):
        thresholds = panel[label]["thresholds"]
        for filter_id, row in _panel_rows(panel[label]["html"]).items():
            stated = _WHY_COMPARISON.match(row["why"])
            if stated is None:  # a frozen or undecided gate states its missing input, not a comparison
                continue
            in_force = thresholds.get(filter_id, run_thresholds[filter_id])
            assert _num(stated.group("threshold")) == in_force, (label, filter_id, row["why"])
            assert stated.group("comparator") == comparators[filter_id], (label, filter_id)
            satisfied = _COMPARE[comparators[filter_id]](_num(stated.group("value")), in_force)
            assert satisfied is (stated.group("negated") is None), (label, filter_id, row["why"])
            assert satisfied is (row["verdict"] == "pass"), (
                f"{label}: {filter_id}'s Why sentence and its pill disagree: {row['why']} / {row['verdict']}"
            )
            compared += 1
    assert compared >= 8, f"only {compared} rows stated a comparison; the property is nearly vacuous"

    # The move really flips this gate -- without that, everything above holds trivially.
    at_run = _panel_rows(panel["atRun"]["html"])[_MOVED_GATE]
    at_moved = _panel_rows(panel["atMoved"]["html"])[_MOVED_GATE]
    assert (at_run["verdict"], at_moved["verdict"]) == ("fail", "pass"), (at_run, at_moved)

    # And the reader's threshold is never silently substituted for the run's: the moved row names both,
    # the unmoved panel names neither, and the panel says up front whose thresholds these verdicts are.
    assert f"le {_MOVED_TO}" in at_moved["threshold_cell"]
    assert f"run {run_thresholds[_MOVED_GATE]}" in at_moved["threshold_cell"], at_moved["threshold_cell"]
    assert f"the run used le {run_thresholds[_MOVED_GATE]}" in at_moved["why"], at_moved["why"]
    assert "1 threshold moved from the run's" in panel["atMoved"]["html"]
    assert "yours" not in at_run["threshold_cell"] and "your threshold" not in at_run["why"]
    assert "moved from the run's" not in panel["atRun"]["html"]


def _run_panels_driver(script: str, guides: list[str]) -> dict[str, Any]:
    """Build the gate panel for several guides, all at the run's own thresholds."""
    driver = f"""
const __GUIDES = {json.dumps(guides)};
process.stdout.write(JSON.stringify(Object.fromEntries(
  __GUIDES.map(k => [k, gatesCard(G.find(x => x.guide === k))]))));
"""
    return _node_run(script, driver)


#: The two fixture guides whose rows the RUN itself left undecided on a gate it held in force: one
#: recorded ``unknown`` (reason 5, no value), one recorded ``not_evaluated`` (reason 6, value 5 present).
_UNDECIDED = {"RUN_UNKNOWN_GATE": ("run_unknown_gate", 5), "RUN_NOTEVAL_GATE": ("run_not_evaluated_gate", 6)}


@pytest.mark.unit
def test_a_gate_the_run_left_undecided_states_its_non_decision_not_a_comparison(
    payload,  # noqa: ANN001
    rendered_script: str,
) -> None:
    """Reasons 5 and 6 must say which non-decision they were, never invent a comparison (#103).

    ``gateReason`` special-cased reasons 1-4 and let 5 and 6 fall through to the comparison line, so a
    gate the run held in force and then declined to decide got a sentence describing a comparison
    nobody made. Reason 6 is the sharp case: payload.py keeps the value the run recorded and applies no
    verdict to it, so the fallthrough printed ``5 not ge 1`` -- arithmetically false, and beside an
    ``unknown`` pill -- which is the fabricated-evidence direction stated in words.

    The two states must also stay distinguishable in the panel. Both now ship under the UNKNOWN verdict
    code (payload.py ``_evaluate``, #103), because an in-force gate the run could not decide is
    "applied, evidence unavailable" and not "not applied", so the pill alone can no longer tell a
    reader which of the two happened; only the reason, and therefore only this column, can.
    """
    reasons = {
        f["filter_id"]: reason
        for g in payload.guides
        if g.guide == _guide_of("RUN_NOTEVAL_GATE")
        for f, (_v, _verdict, reason) in zip(payload.filters, g.gates, strict=True)
    }
    assert reasons["run_not_evaluated_gate"] == 6, "the fixture stopped exercising reason 6 on this guide"

    panels = _run_panels_driver(rendered_script, [_guide_of(label) for label in _UNDECIDED])
    for label, (filter_id, reason) in _UNDECIDED.items():
        html = panels[_guide_of(label)]
        row = _panel_rows(html)[filter_id]
        assert row["verdict"] == "unknown", (label, row)
        assert _WHY_COMPARISON.match(row["why"]) is None, (
            f"{label}: {filter_id} states a comparison the run never made: {row['why']}"
        )
        if reason == 5:
            assert row["why"] == "the run recorded unknown for this gate; no value to re-compare", row
            assert row["value"] == "—", f"reason 5 nulls the value; the panel showed {row['value']!r}"
        else:
            assert row["why"] == ("the run did not evaluate this gate; value measured, no verdict applied"), row
            assert _num(row["value"]) == 5, "reason 6 keeps the number the run recorded"
        # The banner over these rows names the set it counts. It counts UNKNOWN verdicts, which after #103
        # include this row, while the pills that literally read "not evaluated" are the ones excluded --
        # so calling the count "gates not evaluated" named the wrong set.
        assert re.search(r"<b>\d+ of \d+ gates undecided\.</b>", html), label
        assert "gates not evaluated" not in html, label


def _run_cart_driver(script: str, guides: list[str]) -> dict[str, Any]:
    """Export the cart as TSV twice: at the run's thresholds, and with one gate and one bound moved."""
    driver = f"""
const __PICKS = {json.dumps(guides)}, __GATE = {json.dumps(_MOVED_GATE)}, __TO = {json.dumps(_MOVED_TO)};
const __CART = (() => {{
  __PICKS.forEach(k => CART.add(k));
  const atRun = cartTsv(cartRows());
  T[__GATE] = __TO; R.composite = 30; activePreset = 'near_miss';
  recomputeLive();
  const atMoved = cartTsv(cartRows());
  CART.clear(); resetControls(); recomputeLive();
  return {{atRun, atMoved}};
}})();
process.stdout.write(JSON.stringify(__CART));
"""
    return _node_run(script, driver)


def _tsv_parts(text: str) -> tuple[list[str], list[str], list[list[str]]]:
    """A cart export split into its comment lines, its header cells and its body rows."""
    lines = text.split("\n")
    comments = [line for line in lines if line.startswith("#")]
    rest = [line for line in lines if not line.startswith("#")]
    return comments, rest[0].split("\t"), [line.split("\t") for line in rest[1:]]


@pytest.mark.unit
def test_the_cart_export_records_the_thresholds_its_status_column_came_from(
    payload,  # noqa: ANN001
    rendered_script: str,
) -> None:
    """A ``status`` a reader re-thresholded cannot travel without the thresholds (#103, finding 3).

    The export's whole header was one row of column names, and its ``status`` column carries the live
    re-thresholded status -- so a TSV pasted into a ticket claimed a verdict that could not be
    reproduced from the run, and nothing in the file said so. The provenance rides as ``#`` comment
    lines above an unchanged header, which is why the column assertions below are frozen literals.
    """
    picks = [_guide_of(_FLIPPED_GUIDE), _guide_of("ALL_PASS")]
    cart = _run_cart_driver(rendered_script, picks)
    by_guide = {g.guide: g for g in payload.guides}

    at_run_comments, header, at_run_rows = _tsv_parts(cart["atRun"])
    at_moved_comments, moved_header, at_moved_rows = _tsv_parts(cart["atMoved"])
    assert tuple(header) == _CART_COLUMNS, "the export's columns or their order changed"
    assert moved_header == header, "the provenance lines must not disturb the columns"
    assert len(at_run_rows) == len(picks) and len(at_moved_rows) == len(picks)

    # At the run's own thresholds the file says so, and its statuses are the run's own.
    assert any(line == "#status_basis=run_thresholds\tmoved_gates=0" for line in at_run_comments), at_run_comments
    assert not [line for line in at_run_comments if line.startswith("#moved_gate=")]
    status = header.index("status")
    for row in at_run_rows:
        assert row[status] == by_guide[row[0]].status, row

    # Moved: the basis flips, every moved gate is named with both values, the reader bound and the
    # preset that decided which guides `Add top n` would have picked ride along, and a status really did
    # change -- so a reader holding the old file could not have told the two apart.
    assert any(line == "#status_basis=reader_rethresholded\tmoved_gates=1" for line in at_moved_comments)
    assert (
        f"#moved_gate={_MOVED_GATE}\tcomparator=le\trun_threshold="
        f"{next(f['threshold'] for f in payload.filters if f['filter_id'] == _MOVED_GATE)}"
        f"\treader_threshold={_MOVED_TO}"
    ) in at_moved_comments, at_moved_comments
    assert "#reader_filter=composite\tdirection=min\tbound=30" in at_moved_comments
    assert any(line.startswith("#sirnaforge_cart\t") and "preset=near_miss" in line for line in at_moved_comments)
    flipped = _guide_of(_FLIPPED_GUIDE)
    before = next(row[status] for row in at_run_rows if row[0] == flipped)
    after = next(row[status] for row in at_moved_rows if row[0] == flipped)
    assert before != after, "the fixture must actually re-threshold a status, or this proves nothing"


def _run_numeric_boundary_driver(script: str) -> dict[str, Any]:
    """Drive the one numeric boundary every control and the fragment parser go through."""
    driver = """
const __LABEL = v => v === undefined ? 'undefined' : (v === null ? 'null' : v);
const __GATE = FILTERS.find(f => evaluable(f)).filter_id;
const __RUN = FILTERS.find(f => f.filter_id === __GATE).threshold;

const __PARSED = Object.fromEntries(
  ['', '   ', '0', '12', ' 12.5 ', '-3', '1e3', 'abc', '5px', 'Infinity', '-Infinity', 'NaN', '1e999']
    .map(s => [s, __LABEL(parseControlValue(s))]));

// The defect, reproduced: NaN in T fails every re-decidable gate on every guide, because every
// comparison against NaN is false -- and the fragment it writes cannot be read back.
T[__GATE] = Number('abc');
recomputeLive();
const __POISONED = {statuses: G.map(g => g._live.status), hash: encodeHash()};
resetControls(); recomputeLive();
R.composite = Number('abc');
const __POISONED_BOUND = {hash: encodeHash()};
resetControls(); recomputeLive();

// And the fragment parser refuses the same values the controls do, including an empty one -- which
// Number() reads as 0, and which would therefore have set a real threshold of 0.
const __FRAGMENTS = ['t=' + __GATE + ':abc', 't=' + __GATE + ':Infinity', 't=' + __GATE + ':',
                     'r=composite:abc', 'r=composite:'].map(frag => {
  const applied = applyHashFragment(frag);
  const out = {frag, refused: applied.refused, gate: __LABEL(T[__GATE]), composite: __LABEL(R.composite)};
  resetControls(); recomputeLive();
  return out;
});

process.stdout.write(JSON.stringify({parsed: __PARSED, poisoned: __POISONED,
  poisonedBound: __POISONED_BOUND, fragments: __FRAGMENTS, gate: __GATE, run: __RUN}));
"""
    return _node_run(script, driver)


@pytest.mark.unit
def test_a_value_that_is_not_a_number_is_refused_rather_than_poisoning_every_verdict(
    rendered_script: str,
) -> None:
    """``abc`` and ``Infinity`` must never reach a comparison (#103, findings 4 and 6).

    The threshold boxes carried only ``inputmode="decimal"``, a soft-keyboard hint rather than a
    constraint, so both reached the handler. ``Number('abc')`` is NaN and every comparison against NaN
    is false, so one junk keystroke failed every re-decidable gate on every guide -- and ``encodeHash``
    then wrote a fragment that ``applyHashFragment`` itself refuses, so the reader's own reload could
    not reproduce the page they were looking at. One boundary answers both halves: nothing non-finite
    enters ``T`` or ``R``, and nothing non-finite leaves in the fragment.
    """
    out = _run_numeric_boundary_driver(rendered_script)

    assert out["parsed"] == {
        "": "null",
        "   ": "null",  # blank is "no bound", never 0
        "0": 0,
        "12": 12,
        " 12.5 ": 12.5,
        "-3": -3,
        "1e3": 1000,
        "abc": "undefined",
        "5px": "undefined",
        "Infinity": "undefined",
        "-Infinity": "undefined",
        "NaN": "undefined",
        "1e999": "undefined",  # overflows to Infinity, which is not a threshold either
    }

    # The fragment never carries what the parser would refuse -- not even when T or R is poisoned by
    # some other route, which is what makes this a property of the writer rather than of the controls.
    assert "t=" not in out["poisoned"]["hash"], out["poisoned"]["hash"]
    assert "r=" not in out["poisonedBound"]["hash"], out["poisonedBound"]["hash"]
    # And this is what a NaN threshold does, which is why it is refused at the control.
    assert set(out["poisoned"]["statuses"]) == {"fail"}, "the fixture must show NaN failing everything"

    for entry in out["fragments"]:
        expected = "composite" if entry["frag"].startswith("r=") else out["gate"]
        assert expected in entry["refused"], entry
        assert entry["gate"] == out["run"], f"{entry['frag']} moved a threshold it should have refused"
        assert entry["composite"] == "null", f"{entry['frag']} set a bound it should have refused"


@pytest.mark.unit
def test_every_threshold_box_refuses_non_numeric_input_at_the_control(payload) -> None:  # noqa: ANN001
    """The DOM half of the same rule: the boxes are typed, and both handlers refuse what fnum refuses.

    Structural because the handlers are the one part of the panel node cannot drive -- they read
    ``validity.badInput``, which only a browser sets. ``type="number"`` is what makes that flag exist;
    ``inputmode`` never constrained anything.
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
    """The harness's own trustworthiness (#103, finding 7): a truncated document must not pass parity.

    A greedy ``<script>(.*)</script>`` slice runs to the LAST ``</script>``, so the forged document
    below hands node a script that parses and runs while a browser stops at the first one and drops
    everything after it. The slice refuses it instead.
    """
    forged = '<!DOCTYPE html><body><script>const G=["</script>"];\nrun()</script></body></html>'
    assert re.search(r"<script>(.*)</script>", forged, re.DOTALL).group(1).count("</script>") == 1, (
        "the greedy regex swallowed the injected close tag and would have handed node a live script"
    )
    with pytest.raises(AssertionError, match="closes its own"):
        _slice_script(forged)


@pytest.mark.unit
def test_a_payload_spelling_a_script_close_cannot_truncate_the_document(tmp_path: Path) -> None:
    """A gene query, transcript id or symbol may contain anything; the document must survive it.

    ``json.dumps`` does not escape ``<``, and no layer between a run's own metadata and the rendered
    file sanitises it, so a ``</script>`` anywhere in the payload ended the script block early and took
    every panel with it. The renderer escapes it (``render._embed``); this asserts the document still has
    exactly one script, and that the string arrives in the browser's JavaScript unchanged.
    """
    hostile = "TP53</script><script>run()</script>"
    run = _write_fixture_run(tmp_path, gene_query=hostile)
    html = render_html(build_payload(run))

    assert html.count("</script>") == 1, "a payload string closed the script block"
    assert "\\u003c/script\\u003e" in html, "the close tag must be escaped, not stripped or dropped"
    assert hostile not in html, "and it must not appear unescaped anywhere, markup included"

    echoed = _node_run(_slice_script(html), "process.stdout.write(JSON.stringify({gene: GENE}));")
    assert echoed["gene"] == hostile, "escaping must be reversible: the export filename reads GENE"
