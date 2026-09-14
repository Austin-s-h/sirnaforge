"""Client-side re-thresholding, and the fields a preset view needs (issue #103).

The report lets a reader move a gate's threshold and watch it re-decided in the browser.
:func:`reevaluate_gates` is the Python reference that shipped evaluator has to agree with, and the
property no threshold may ever defeat: a gate the run itself did not decide -- off, undeclared, a
missing column, or its own recorded ``unknown``/``not_evaluated`` -- freezes untouched, whatever a
reader asks for. The draft this was salvaged from got that property right in spirit but wrong in
code: it froze per-*filter*, not per-*(guide, gate)*, so a row the run marked ``not_evaluated`` while
still recording a real number fell through to a fresh pass/fail. The tests below pin the fix.

The rest of the file pins the per-filter and per-guide fields a preset view (all / passing /
near-miss / off-target-clean / register-deduplicated) cannot be built without, and that the draft
either lacked or got from the wrong column.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sirnaforge.reporting import build_payload
from sirnaforge.reporting.payload import (
    REASON_EMPTY_VALUE,
    REASON_OK,
    REASON_RUN_NOT_EVALUATED,
    VERDICT_WARN,
    reevaluate_gates,
)

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
