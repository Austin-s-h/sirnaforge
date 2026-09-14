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
from pathlib import Path

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
}

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


def _write_fixture_run(tmp_path: Path) -> Path:
    """A run directory exercising every reason code, every comparator and three presets' edge cases.

    ``empty_value_le`` is blank on exactly one guide's row (REASON_EMPTY_VALUE); every other filter's
    column is real for every row, so the difference between the two is a fact about one guide, not
    about the run.
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
        for column in _EXTRA_COLUMNS + ("passes_filters",):
            if column in overrides:
                cells[column] = overrides[column]
        lines.append(",".join(str(cells.get(c, "")) for c in columns))
    (run / "sirnaforge" / "candidates_all.csv").write_text("\n".join(lines) + "\n")

    manifest = {
        "gene_query": "TP53",
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
    cells = unclassified_cells()
    cells.update(hit_class="on_target", hit_symbol="TP53", matched_symbol="TP53")
    hit_columns = [*base_columns, *cells]
    hit_row = [
        "screen_0",
        _seq("OFFTARGET_CLEAN"),
        "human",
        "ENST00000000009",
        "100",
        "+",
        "21M",
        "60",
        "42",
        "0",
        "0",
        "1.0",
        *cells.values(),
    ]
    table = run / "off_target" / "results" / "aggregated" / "combined_offtargets.tsv"
    table.write_text("\t".join(hit_columns) + "\n" + "\t".join(hit_row) + "\n")
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


@pytest.fixture(scope="module")
def rendered_script(payload) -> str:  # noqa: ANN001
    """The rendered report's ``<script>`` body -- placeholder substitution has already happened."""
    html = render_html(payload)
    match = re.search(r"<script>(.*)</script>", html, re.DOTALL)
    assert match is not None, "the rendered report has no <script> block"
    return match.group(1)


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
    node = _node_or_fail()
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
    # A temp .mjs, not `node -e`: the substituted script is ~60 KB and this is the shape the issue's
    # own instructions ask for -- the placeholder substitution happens for real, in Python, before
    # node ever sees a byte of it.
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

    A static grep over the six real, hand-written FILTERS_UI boxes this replaced is what the issue's
    own review caught -- ``gc_content``, ``asymmetry_score`` and ``off_target_count`` were literal
    strings in the old template. None of the real registry's 17 ids may appear as a literal here.
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
    # payload._normalise_guide upper-cases and spells the guide in RNA (T -> U); the fixture's CSV
    # rows store the DNA spelling, like every other candidates_all.csv, so the lookup key must too.
    guide_of = {label: _seq(label).replace("T", "U") for label in _GUIDES}

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


@pytest.mark.unit
def test_node_absent_fails_loudly_rather_than_skipping(monkeypatch: pytest.MonkeyPatch) -> None:
    """A skip is not a gate: without node this must fail with the install hint, never pass quietly."""
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    with pytest.raises(pytest.fail.Exception, match="node"):
        _node_or_fail()
