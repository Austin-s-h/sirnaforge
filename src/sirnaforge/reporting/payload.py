"""Build the report payload from a finished run directory.

Reads ``candidates_all.csv``, the aggregated hit tables and ``manifest.json`` into a typed payload.
Performs **no** alignment, folding, network access or classification: ``hit_class`` is read from the
table the run published, so the report can never disagree with the run.

The report unit is the **guide sequence**, not the candidate id. A guide is enumerated once per
transcript it was found on (median 6, up to 34 rows on the public 0.7.1 baseline), while off-target
screening runs once per distinct guide; joining on ``id`` instead leaves most candidates falsely
reading as having zero off-targets. The guide sequence is the join key, and per-transcript rows
collapse into an isoform sub-table.

This module also owns **how many guides the report contains**, and every count it publishes is a
count of those guides: see :data:`DEFAULT_MAX_EMBEDDED_GUIDES` and :func:`_select_embedded`.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from collections.abc import Container, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from sirnaforge.config.run_policy import (
    SETTING_BY_KEY,
    EntryPoint,
    ResolvedRunPolicy,
    RunPolicyError,
    filters_from_manifest,
    resolve_run_policy,
)
from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS, hit_class_of, is_annotated
from sirnaforge.core.hit_classification import HitClass
from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation
from sirnaforge.models.sirna import FilterCriteria, OffTargetFilterCriteria
from sirnaforge.reporting.structure import layouts_for
from sirnaforge.reporting.tracks import TranscriptRegions, transcript_regions, uncovered_stretches

#: Bump when the payload's shape changes, so a report and the run it describes can never be
#: silently mismatched. 1.2.0 moved the guide cap here, published as counts (``guides_total``/
#: ``guides_dropped``/``guides_dropped_by_status``), and added per-filter ``rejects``/
#: ``sole_rejects``/``inert`` and per-guide ``liability_by_species``.
PAYLOAD_SCHEMA_VERSION = "1.2.0"

#: Hit rows embedded per guide. Everything outside this scope is carried as counts only -- a
#: deliberate scope decision, not a limitation, so a guide with hits only outside it must render its
#: counts rather than an empty panel.
EMBED_SPECIES = "human"
EMBED_MAX_NM = 2

#: Classes that are a liability. On-target isoforms, orthologues and repeat-mediated hits are
#: displayed, but never counted as off-targets.
LIABILITY = frozenset({HitClass.OFF_TARGET, HitClass.UNDETERMINED})

#: Guides one report embeds. Kept because the file is read in a browser and one internal run of
#: 5,706 guides already renders to 31 MB, applied **here**, above every count the report prints, and
#: never below them -- the only cap on how many guides a report holds. Set ``max_guides=None`` to
#: embed the run whole.
DEFAULT_MAX_EMBEDDED_GUIDES = 5000

#: Statuses the cap drops last, in this order: a ``fail`` guide is one the run and report agree to
#: reject, so it is the only kind whose absence costs a reader nothing actionable. Anything a reader
#: might ship, or that the report cannot establish either way, outranks it -- which is why the
#: dropped set is chosen by verdict, not by taking the first N of a score sort.
CAP_DROP_ORDER = ("fail", "unknown", "warn", "pass")

_UNKNOWN_SYMBOL = "unknown"


class ReportInputError(ValueError):
    """A run directory that cannot produce a truthful report."""


@dataclass
class GuideEntry:
    """Everything the report shows for one guide sequence."""

    guide: str
    passenger: str | None
    overhang: str | None
    modifications: str | None
    n_rows: int
    design_score: float | None
    composite_score: float | None
    weight_vector: str | None
    metrics: dict[str, Any]
    gates: list[list[Any]]
    isoforms: list[dict[str, Any]]
    offtarget_by_symbol: list[dict[str, Any]]
    offtarget_matrix: list[dict[str, Any]]
    offtarget_rows: list[dict[str, Any]]
    offtarget_embedded_scope_empty_but_counts_exist: bool
    liability_count: int
    mirna: list[dict[str, Any]]
    #: species -> liability alignments in that species, summing to :attr:`liability_count`. The
    #: all-species number is what ``max_off_target_count`` gates; this is what a reader needs to ask
    #: whether it would still reject the guide on one species (:func:`_liability_by_species`).
    liability_by_species: dict[str, int] = field(default_factory=dict)
    #: species -> {nm, seed_mismatches} for the best ortholog alignment. Empty when none was found.
    ortholog: dict[str, dict[str, Any]] = field(default_factory=dict)
    run_verdict: str | None = None
    structure: str | None = None
    transcript_hits: int | None = None
    #: Whether screening reached the aligner for this guide (``off_target_screened``, models/sirna.py).
    #: A clean preset built on ``liability_count == 0`` alone would publish a never-screened guide as
    #: clean -- 0 hits and no screen look identical unless this rides along.
    off_target_screened: bool = False
    #: The query id this guide's rows were screened under (``screen_query_id``, #103's join key).
    #: Carried for provenance; None when the candidate was never submitted to the aligner.
    screen_query_id: str | None = None

    @property
    def undeclared_run_rejection(self) -> bool:
        """The run rejected this guide everywhere and no declared gate accounts for it.

        ``REPEAT_ELEMENT`` is the live case: the pipeline stamps it, but the 17-filter registry
        declares no repeat gate, so nothing here can re-derive the rejection. Must not read as clean --
        on one internal run this would have published 185 guides as passing that the run rejected.
        """
        return self.run_verdict not in (None, "PASS") and not (self.n_gates_failed or self.n_gates_unknown)

    @property
    def n_gates_failed(self) -> int:
        """Gates this guide independently fails -- all of them, not only the first to fire."""
        return sum(1 for g in self.gates if g[1] == _VERDICT_CODE[FilterEvaluation.FAIL.value])

    @property
    def n_gates_unknown(self) -> int:
        """Gates in force that reached no verdict; never folded into a pass.

        Every in-force non-decision -- a missing column, an empty value, the run's own ``unknown``, or
        its own ``not_evaluated`` -- is emitted as the UNKNOWN verdict by :func:`_evaluate`, with the
        reason code saying which, so this tally and the browser's (which counts UNKNOWN alone) cannot
        disagree about the same guide.
        """
        return sum(1 for g in self.gates if g[1] == _VERDICT_CODE[FilterEvaluation.UNKNOWN.value])

    @property
    def n_gates_warned(self) -> int:
        """Gates this guide exceeds whose action is ``warn``: reported, and not a rejection."""
        return sum(1 for g in self.gates if g[1] == VERDICT_WARN)

    @property
    def status(self) -> str:
        """``fail`` / ``unknown`` / ``warn`` / ``pass``.

        A guide failing no gate is clean only when every gate could actually be evaluated: on the
        public 0.7.1 baseline, 450 guides fail nothing here while the run failed them, because those
        gates read counters the run does not export. Collapsing that into ``pass`` is the
        fabricated-evidence mistake this report exists to make visible.

        ``warn`` stays separate from ``fail`` for the mirror-image reason (:data:`VERDICT_WARN`). A
        rejection no declared gate can express is ``unknown`` too (:attr:`undeclared_run_rejection`) --
        the report may report less than the run, never more.
        """
        if self.n_gates_failed:
            return "fail"
        if self.n_gates_unknown or self.undeclared_run_rejection:
            return "unknown"
        return "warn" if self.n_gates_warned else "pass"


@dataclass
class ReportPayload:
    """The whole report, ready to render."""

    schema_version: str
    run: dict[str, Any]
    filters: list[dict[str, Any]]
    guides: list[GuideEntry]
    provenance: dict[str, Any]
    caveats: list[str] = field(default_factory=list)


def _normalise_guide(seq: object) -> str:
    """Uppercase RNA spelling, so a U/T difference cannot split one guide into two."""
    return str(seq or "").strip().upper().replace("T", "U")


def _num(value: object) -> float | int | None:
    """A number, or None for NaN/blank -- never a fabricated zero.

    Coerces through ``float`` rather than ``isinstance(value, (int, float))``: ``numpy.int64`` is
    **not** a Python ``int`` while ``numpy.float64`` *is* a Python ``float``, and ``DataFrame.iloc``
    hands back numpy scalars. The isinstance test nulled every integer-dtype column, which made
    ``max_off_target_count`` read ``unknown`` for every guide on a real run while its threshold was
    rejecting 65% of them.
    """
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, str) and not value.strip():
        return None
    try:
        number = float(value)  # type: ignore[arg-type]  # non-numeric raises, which is the answer
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):  # NaN or infinite is not a value
        return None
    return int(number) if number.is_integer() else round(number, 6)


def _flag(value: object) -> bool:
    """A boolean read from a CSV cell, which may already be a bool or the string it was printed as.

    ``pandas.read_csv`` infers a bool dtype only when a column is *exclusively* ``True``/``False``
    tokens; one blank cell downgrades the whole column to ``object`` and hands this the literal
    string ``"False"``, which is truthy under a bare ``bool()``.
    """
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    return str(value).strip().lower() == "true"


#: Why a gate reached its verdict, as a code. The descriptor is emitted once for the whole report, so
#: a per-guide outcome carries only what varies -- value, verdict and this. Repeating filter_id,
#: column, comparator, threshold and scope on every guide cost 9 MB of the first render's 38 MB.
REASON_OK = 0
REASON_MISSING_COLUMN = 1
REASON_EMPTY_VALUE = 2
REASON_FILTER_OFF = 3
REASON_NO_THRESHOLD = 4
#: The run itself recorded UNKNOWN for this gate. A distinct reason from an empty or absent column,
#: because the run measured nothing on purpose and said so, rather than the report failing to find it.
REASON_RUN_UNKNOWN = 5
#: The run recorded NOT_EVALUATED for this gate, so it made no claim and the report makes none either.
#: A *reason*, not a verdict. A gate that is in force and undecided is published under the UNKNOWN
#: verdict code and this reason (#103); see :func:`_evaluate`.
REASON_RUN_NOT_EVALUATED = 6

#: Verdict codes, keyed by the run's own vocabulary. ``not_evaluated`` (3) means the gate was never in
#: force for this run -- off, or with no declared threshold -- and nothing else: an in-force gate the
#: run declined to decide is UNKNOWN (2), because "not applied" and "applied, evidence unavailable" are
#: different claims and only the second one bears on whether a guide is clean (#103).
_VERDICT_CODE = {
    FilterEvaluation.PASS.value: 0,
    FilterEvaluation.FAIL.value: 1,
    FilterEvaluation.UNKNOWN.value: 2,
    FilterEvaluation.NOT_EVALUATED.value: 3,
}

#: A gate whose threshold was exceeded but whose action is ``warn``: a real finding the reader should
#: see, and not a rejection. It needs its own code because folding it into ``fail`` would make the
#: report contradict the run on every warn-flagged guide -- inverting the one metric
#: (``contradicted_run_pass``) that exists to prove the report and the pipeline agree.
VERDICT_WARN = 4

#: Guide statuses, worst first. ``warn`` belongs in the header tally like the rest: counting only
#: pass/unknown/fail dropped every warn-status guide out of a total that claims to be all of them.
STATUSES = ("fail", "unknown", "warn", "pass")


def _select_embedded(guides: Sequence[GuideEntry], max_guides: int | None) -> tuple[list[GuideEntry], dict[str, int]]:
    """The guides this report will contain, and what the cap cost, by verdict.

    ``guides`` arrives best-scoring first. What is dropped is chosen by verdict in
    :data:`CAP_DROP_ORDER` and, within a verdict, worst-scoring first -- so a run larger than the cap
    loses the guides its own gates rejected before it loses one a reader could ship. Taking the first
    N of the score sort instead dropped 41 pass/warn guides from one internal report while its header
    still counted them, uncartable and unexported.

    Returns:
        ``(embedded, dropped_by_status)`` with ``embedded`` still in score order and every status
        present in the second, zeros included, so a consumer never has to guard a missing key.
    """
    empty = dict.fromkeys(STATUSES, 0)
    if max_guides is None or len(guides) <= max_guides:
        return list(guides), empty
    worst_first = sorted(range(len(guides)), key=lambda i: (CAP_DROP_ORDER.index(guides[i].status), -i))
    dropped = set(worst_first[: len(guides) - max_guides])
    counts = Counter(guides[i].status for i in dropped)
    return [g for i, g in enumerate(guides) if i not in dropped], {**empty, **counts}


def observed_column(descriptor: Any, populated: Container[str]) -> str | None:
    """The column carrying the value this gate compared, or None if the run exports neither.

    ``<filter_id>_observed`` is preferred over the descriptor's own ``column``, because it is the
    number the pipeline itself compared. It is the answer for gates that read human-stratified
    counters: 0.7.1 does not export ``transcriptome_hits_1mm_human`` under that name (#101), but does
    export ``max_transcriptome_hits_1mm_observed`` -- on a four-species internal run the observed
    column reproduces each gate's own verdict on 100% of 40,081 rows while the same-named all-species
    column disagrees (17,600 hits against 63,801). The same-named column holds a wider-scoped quantity
    than the gate compared, so reading it instead would let the report contradict the run while
    looking like it agreed.

    ``populated`` must hold only columns that carry at least one value. A column present but empty for
    every row is the shape a gate takes when the run never recorded its verdict, and preferring it
    would turn a gate the descriptor column *can* answer into an unknown.
    """
    for candidate in (f"{descriptor.filter_id}_observed", descriptor.column):
        if candidate in populated:
            return candidate
    return None


def _recorded_verdict(row: pd.Series, filter_id: str) -> str | None:
    """The verdict the RUN recorded for one gate on one row, from its ``<filter_id>_verdict`` column.

    Distinct from :func:`_run_verdict`, which reads the single ``passes_filters`` label over a guide's
    rows. Only this per-gate record can say "in force, evidence unavailable" -- a state no threshold
    comparison reconstructs.
    """
    value = row.get(f"{filter_id}_verdict")
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    return str(value)


def _verdict_code(passed: bool, action: str) -> int:
    """The verdict code one comparison earns, given the gate's action.

    Both :func:`_evaluate` and :func:`reevaluate_gates` decide this, so they decide it here: a
    ``warn`` gate the guide exceeds must come back as :data:`VERDICT_WARN` and never the fail code
    (see its definition for what folding the two together would cost).
    """
    if passed:
        return _VERDICT_CODE[FilterEvaluation.PASS.value]
    if action == FilterAction.WARN.value:
        return VERDICT_WARN
    return _VERDICT_CODE[FilterEvaluation.FAIL.value]


def _evaluate(descriptor: Any, row: pd.Series, column: str | None) -> tuple[float | int | None, int, int]:
    """Evaluate one descriptor against one candidate row, independently of every other gate.

    Returns ``unknown`` rather than a pass when the run exports no column the threshold can read, and
    equally when the run kept the gate in force but recorded no verdict for it. Reporting either as a
    pass would be the fabricated-evidence failure this report exists to make visible.

    Args:
        descriptor: The gate as configured for this run.
        row: One candidate row.
        column: Column to read, from :func:`observed_column`; None when the run exports neither.

    Returns:
        ``(value, verdict_code, reason_code)``. The descriptor itself is emitted once per report.
    """
    # Read the value regardless of verdict: an off or threshold-less gate still *measured* something --
    # three of the four off gates on the reference run carry a real per-guide number, and
    # `max_mirna_1mm_seed`'s own policy definition says "the number is reported and nothing acts on it".
    value = _num(row.get(column)) if column is not None else None
    not_evaluated = _VERDICT_CODE[FilterEvaluation.NOT_EVALUATED.value]
    if descriptor.action.value == FilterAction.OFF.value:
        return value, not_evaluated, REASON_FILTER_OFF
    if descriptor.threshold is None:
        return value, not_evaluated, REASON_NO_THRESHOLD
    # The run's own non-decision wins before any comparison, never toward a pass: `observed_column`
    # falls back to the descriptor's own column when the observed one is empty for every row, and
    # three of those fallbacks default to 0, so without this an unscreened guide's
    # `max_off_target_count` would read as a pass at 0.
    #
    # Both non-decisions land on UNKNOWN, not NOT_EVALUATED: a gate reaching this line is in force with
    # a threshold -- off and threshold-less gates returned above -- and NOT_EVALUATED sits on a verdict
    # code the unknown tally does not count, so `status` would fall through to `pass` with the banner
    # suppressed. `max_repeat_transcript_fraction` is the live case: a repeat scan that never ran
    # leaves `repeat_transcript_fraction` at its 0.0 default with verdict `not_evaluated`. The reason
    # code keeps the two apart along with the value: UNKNOWN has no measurement to show, while
    # `not_evaluated` keeps the number the run recorded.
    recorded = _recorded_verdict(row, descriptor.filter_id)
    if recorded in (FilterEvaluation.UNKNOWN.value, FilterEvaluation.NOT_EVALUATED.value):
        run_unknown = recorded == FilterEvaluation.UNKNOWN.value
        return (
            None if run_unknown else value,
            _VERDICT_CODE[FilterEvaluation.UNKNOWN.value],
            REASON_RUN_UNKNOWN if run_unknown else REASON_RUN_NOT_EVALUATED,
        )
    if column is None:
        return None, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_MISSING_COLUMN
    if value is None:
        return None, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_EMPTY_VALUE

    # The one implementation of the comparator table in Python -- collapsed onto
    # ``FilterComparator.passes`` (#103) so :func:`reevaluate_gates` and the report's JS restate
    # exactly this table, not a second copy that could drift from it.
    passed = descriptor.comparator.passes(value, descriptor.threshold)
    return value, _verdict_code(passed, descriptor.action.value), REASON_OK


def reevaluate_gates(
    filters: Sequence[Mapping[str, Any]],
    gates: Sequence[Sequence[Any]],
    thresholds: Mapping[str, float] | None = None,
) -> list[list[Any]]:
    """Re-apply one guide's gates at reader-chosen thresholds.

    The contract the report's client-side evaluator implements, kept here so it can be asserted
    without a browser and the browser checked against it. Two rules:

    * A gate this guide's own row did not decide is returned **untouched**. Freezing is keyed on the
      row's ``reason`` code, not on whether the *filter* is globally evaluable: a filter can export a
      column for most guides and still leave one row at :data:`REASON_RUN_NOT_EVALUATED` with a real
      measured value, and moving a threshold cannot conjure the decision the run declined to make.
      Every reason except :data:`REASON_OK` and :data:`REASON_EMPTY_VALUE` freezes the triple exactly
      as the run produced it.
    * :data:`REASON_EMPTY_VALUE` stays ``unknown``. Its value is already ``None``, so nothing is
      re-compared; this is what stops a reader re-thresholding their way out of a non-decision by
      supplying a value the run never observed.

    A ``warn``-action gate the guide exceeds comes back as :data:`VERDICT_WARN`, never ``fail`` (see
    its definition for why). Comparisons go through :meth:`FilterComparator.passes`, the single Python
    implementation of the comparator table, so a client-side evaluator restating it in JavaScript has
    exactly one table to restate.

    Args:
        filters: Filter entries as emitted in the payload, in payload order.
        gates: That guide's ``[value, verdict, reason]`` triples, positionally matching ``filters``.
        thresholds: ``filter_id`` -> threshold. Absent ids keep the run's own threshold.

    Returns:
        Fresh triples in the same order.
    """
    chosen = thresholds or {}
    out: list[list[Any]] = []
    for f, gate in zip(filters, gates, strict=True):
        value, _verdict, reason = gate
        if reason == REASON_EMPTY_VALUE:
            out.append([value, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_EMPTY_VALUE])
            continue
        if reason != REASON_OK:
            out.append(list(gate))
            continue
        threshold = chosen.get(f["filter_id"], f["threshold"])
        passed = FilterComparator(f["comparator"]).passes(value, threshold)
        out.append([value, _verdict_code(passed, f["action"]), REASON_OK])
    return out


#: Slider granularity: a continuous domain is cut into roughly this many steps, then snapped to a
#: power of ten so the readout is legible. Every control also carries a number box, so the step
#: bounds the slider's resolution, never the reachable thresholds.
CONTROL_STEPS = 50

#: The models a filter threshold is validated on, keyed as ``SettingSpec.model`` names them. Only the
#: two filter models can carry a gate threshold; the other setting targets hold no gate.
_SETTING_MODELS: Mapping[str, Any] = {"filters": FilterCriteria, "offtarget_filters": OffTargetFilterCriteria}


@dataclass(frozen=True)
class _DeclaredBounds:
    """What values one setting is allowed to take, read off its own field metadata.

    Derived from the model rather than a hand-written table here, so a clamp cannot drift from the
    declaration it restates. ``exclusive`` records a ``gt``/``lt`` bound, so a slider edge lands one
    step inside a limit the model would reject rather than on it.
    """

    low: float | None = None
    high: float | None = None
    low_exclusive: bool = False
    high_exclusive: bool = False

    def clamp(self, value: float, *, upper: bool, step: float) -> float:
        """``value`` pulled inside the declared bound on one side, if the model declares one."""
        limit = self.high if upper else self.low
        if limit is None:
            return value
        if self.high_exclusive if upper else self.low_exclusive:
            limit = limit - step if upper else limit + step
        return min(value, limit) if upper else max(value, limit)


def _declared_bounds(setting_key: str) -> _DeclaredBounds:
    """The declared range of the setting a gate's threshold lives on, or an empty one if it has none."""
    spec = SETTING_BY_KEY.get(setting_key)
    model = _SETTING_MODELS.get(spec.model) if spec is not None else None
    field_info = model.model_fields.get(spec.field) if model is not None and spec is not None else None
    if field_info is None:
        return _DeclaredBounds()
    bounds = _DeclaredBounds()
    for meta in field_info.metadata:
        for name, exclusive in (("ge", False), ("gt", True)):
            limit = getattr(meta, name, None)
            if limit is not None and (bounds.low is None or float(limit) > bounds.low):
                bounds = _DeclaredBounds(float(limit), bounds.high, exclusive, bounds.high_exclusive)
        for name, exclusive in (("le", False), ("lt", True)):
            limit = getattr(meta, name, None)
            if limit is not None and (bounds.high is None or float(limit) < bounds.high):
                bounds = _DeclaredBounds(bounds.low, float(limit), bounds.low_exclusive, exclusive)
    return bounds


def _snap(value: float, step: float, *, up: bool) -> float:
    """``value`` moved outward onto the step grid, so every slider position is a legible number.

    Unsnapped bounds were the defect: ``gc_content_min`` opened at 27.804348 with a step of 0.1, so
    every drag landed on ...704348 and a reader could not reach 40 -- 39.704348 is what ended up in the
    URL fragment and the cart TSV's ``reader_threshold``.
    """
    scaled = value / step
    grid = math.ceil(round(scaled, 9)) if up else math.floor(round(scaled, 9))
    return round(grid * step, 9)


def _control_domain(threshold: float, values: Sequence[float], bounds: _DeclaredBounds) -> dict[str, Any]:
    """Slider bounds for one filter: this run's own values, snapped to the grid and inside the model's.

    The span comes from the data because the descriptors carry no display range, and a domain wider
    than the run's values spends its travel where no verdict changes. Two rules keep it honest:

    * Bounds are snapped to the step, so the grid contains round numbers a reader means to type.
    * Bounds are clamped to the setting's **own** declared range. Padding an observed minimum by 5%
      once published "at most -59 off-targets" on a ``ge 0``, count-valued gate -- a threshold that
      fails every guide and the pipeline would refuse -- and locked ``min_empirical_score`` to
      0.39-0.61 either side of a declared 0.4-0.6 range, exactly its two verdict-changing extremes.
      Clamping never crosses an observed value or the run's own threshold, so every measurement this
      run produced stays reachable, and the run's threshold sits on an edge only when the run set it
      at the setting's own limit.

    The declared range itself is published beside the control, not inside it: it is a property of the
    setting, and the slider is one of several things that read it.
    """
    lo, hi = min([*values, threshold]), max([*values, threshold])
    integral = all(float(v).is_integer() for v in [*values, threshold])
    span = hi - lo
    pad = max(span * 0.05, 1e-6) if span else (1.0 if integral else max(abs(hi) * 0.1, 0.1))
    step = 1.0 if integral else 10 ** math.floor(math.log10((span + 2 * pad) / CONTROL_STEPS))
    low = min(bounds.clamp(_snap(lo - pad, step, up=False), upper=False, step=step), lo)
    high = max(bounds.clamp(_snap(hi + pad, step, up=True), upper=True, step=step), hi)
    return {
        "min": int(low) if integral else round(low, 6),
        "max": int(math.ceil(high)) if integral else round(high, 6),
        "step": int(step) if integral else step,
    }


@dataclass(frozen=True)
class _GateEffect:
    """What one gate actually did to the guides this report embedded.

    Seventeen sliders with no indication of which decided anything is the panel the audit found: over
    one internal run five gates reject at all, two of the seventeen account for 1,765 of the
    rejections, and two more are inert by construction. ``sole`` ranks them -- a guide rejected by one
    gate alone is a guide that gate is solely responsible for losing.
    """

    rejects: int = 0
    sole_rejects: int = 0
    warns: int = 0
    unknowns: int = 0

    @property
    def decides_nothing(self) -> bool:
        """No embedded guide is failed or flagged by this gate at the run's own threshold."""
        return not self.rejects and not self.warns


def _gate_effects(guides: Sequence[GuideEntry], n_filters: int) -> list[_GateEffect]:
    """Per-filter reject/sole-reject/warn/unknown tallies over the embedded guides, in filter order."""
    fail_code = _VERDICT_CODE[FilterEvaluation.FAIL.value]
    unknown_code = _VERDICT_CODE[FilterEvaluation.UNKNOWN.value]
    rejects, sole, warns, unknowns = ([0] * n_filters for _ in range(4))
    for guide in guides:
        failed = [i for i, gate in enumerate(guide.gates) if gate[1] == fail_code]
        for i in failed:
            rejects[i] += 1
        if len(failed) == 1:
            sole[failed[0]] += 1
        for i, gate in enumerate(guide.gates):
            if gate[1] == VERDICT_WARN:
                warns[i] += 1
            elif gate[1] == unknown_code:
                unknowns[i] += 1
    return [_GateEffect(*counts) for counts in zip(rejects, sole, warns, unknowns, strict=True)]


def _inert_reason(descriptor: Any, effect: _GateEffect, bounds: _DeclaredBounds, setting_key: str) -> str | None:
    """Why this gate decided nothing in this run, or None when it decided something.

    Two inertnesses, and the second is the one worth acting on: a gate that happens to reject nobody
    here, and one that *cannot* reject anybody because its threshold sits on the limit its own setting
    declares. ``min_empirical_score`` is the live case -- a ``ge`` floor at 0.4 against a declared
    0.4-0.6 range -- and the report gave it a slider as prominent as the two gates that decided 1,765
    rejections.
    """
    if not effect.decides_nothing:
        return None
    floor = descriptor.comparator in (FilterComparator.GE, FilterComparator.GT)
    limit = bounds.low if floor else bounds.high
    # Only a non-strict comparator is inert *at* its limit: ``gt``/``lt`` still reject the limit value.
    at_limit = (
        limit is not None
        and descriptor.threshold is not None
        and descriptor.comparator in (FilterComparator.GE, FilterComparator.LE)
        and float(limit) == float(descriptor.threshold)
    )
    reason = "no guide embedded in this report is failed or flagged by this gate at the run's own threshold"
    if at_limit:
        reason += (
            f"; its threshold is {setting_key}'s own declared "
            f"{'minimum' if floor else 'maximum'} of {limit}, so no permitted value can fail it"
        )
    return reason


def _filter_view(
    descriptor: Any,
    *,
    column: str | None,
    gates: Sequence[Sequence[Any]],
    setting_key: str,
    definition: str,
    effect: _GateEffect,
) -> dict[str, Any]:
    """One filter as the report carries it: the descriptor, plus whether a reader may re-threshold it.

    ``evaluable`` is decided from this run's own gate triples and from ``column`` -- the
    ``<filter_id>_observed``-or-descriptor column :func:`observed_column` resolved, never the
    descriptor's bare column, which would freeze every human-stratified gate ``observed_column`` can
    in fact answer. A gate is re-thresholdable only when it was applied, has a threshold, and this run
    left at least one guide's row on :data:`REASON_OK`, the one reason :func:`reevaluate_gates` and the
    report's JS will re-compare; anything else stays frozen, since no slider position can answer a
    question the run has no evidence for. ``evaluable`` therefore means exactly one thing: moving this
    control can change a verdict. ``n_values`` still counts the measurements, because they were
    measured -- it is the deciding that never happened.

    ``rejects``/``sole_rejects``/``warns``/``unknowns`` say what the gate *did* at the run's own
    threshold, and ``inert`` says it decided nothing, so a panel of seventeen equal-looking sliders can
    be ranked by which of them the run's verdicts actually turned on. All four count the guides this
    report embedded, not the run's -- see :func:`_select_embedded`.
    """
    observed = [value for value, _verdict, _reason in gates if value is not None]
    decided = sum(1 for _value, _verdict, reason in gates if reason == REASON_OK)
    reason: str | None = None
    if descriptor.action.value == FilterAction.OFF.value:
        reason = "this run has the filter off, so it reached no verdict to re-threshold"
    elif descriptor.threshold is None:
        reason = "no threshold is declared for this filter"
    elif column is None:
        reason = f"the run exports neither {descriptor.filter_id}_observed nor {descriptor.column}"
    elif not observed:
        reason = f"{column} is exported but empty for every guide in this run"
    elif not decided:
        reason = (
            f"{column} is exported, but the run evaluated this gate for none of the {len(gates)} guides "
            "in this run, so no threshold can decide it"
        )
    bounds = _declared_bounds(setting_key)
    return {
        **descriptor.model_dump(mode="json"),
        "setting_key": setting_key,
        "definition": definition,
        "scope_label": _scope_label(descriptor),
        "read_column": column,
        "evaluable": reason is None,
        "unevaluable_reason": reason,
        "control": None if reason is not None else _control_domain(descriptor.threshold, observed, bounds),
        "n_values": len(observed),
        # The setting's own declared limits, or None where the model states none. Beside the control
        # rather than in it: a consumer refusing a typed threshold the pipeline would reject needs these
        # even for a gate that has no slider.
        "bound_min": bounds.low,
        "bound_max": bounds.high,
        "rejects": effect.rejects,
        "sole_rejects": effect.sole_rejects,
        "warns": effect.warns,
        "unknowns": effect.unknowns,
        "inert": reason is None and effect.decides_nothing,
        "inert_reason": None if reason is not None else _inert_reason(descriptor, effect, bounds, setting_key),
    }


def _transcript_hits(rows: pd.DataFrame) -> int | None:
    """How many distinct transcripts carry this guide -- the numerator of isoform coverage.

    Not the row count: a guide's site can occur twice in one transcript, so enumerations exceed
    isoforms (14 guides on one internal run have more rows than transcripts, one with 19 rows over
    10). Prefers the run's own ``transcript_hit_count`` and falls back to counting, so the report
    agrees with the column when it exists.
    """
    stated = _num(rows["transcript_hit_count"].iloc[0]) if "transcript_hit_count" in rows.columns else None
    if stated is not None:
        return int(stated)
    return int(rows["transcript_id"].nunique()) if "transcript_id" in rows.columns else None


def _run_verdict(rows: pd.DataFrame) -> str | None:
    """The run's own verdict for a guide: PASS if any of its rows passed, else its first label.

    Guide-level because screening is: one passing enumeration is enough for the run to keep the
    guide, so the report must read the group rather than the best row.
    """
    if "passes_filters" not in rows.columns:
        return None
    labels = [str(v).split(" (")[0] for v in rows["passes_filters"].dropna()]
    if not labels:
        return None
    return "PASS" if "PASS" in labels else labels[0]


def _liability_by_species(hits: pd.DataFrame) -> dict[str, int]:
    """Liability alignments per species -- the decomposition the headline number hides.

    Half of every liability total on a multi-species run is non-human: 92,658 human against 91,536
    mouse, 25,243 macaque and 23,427 rat on one internal run of 232,864 alignments.
    ``max_off_target_count`` is scoped to *all* species and rejected 3,317 of that run's guides, 782 of
    them on nothing else -- and 458 of those 782 pass the same ceiling on human liabilities alone,
    against 448 total passes in the whole report. Which number a consumer reads decides whether it
    agrees with the run, so both are published and neither is called *the* liability count.

    Derived from the published hit table, not the candidate row: the row's ``off_target_count`` is
    all-species, and the ``*_query`` counters #101 added (``total_offtarget_hits_query``,
    ``transcriptome_hits_{0,1,2}mm_query``, ``mirna_hits_0mm_seed_query``, ``mirna_hits_high_risk_query``
    -- all carried in ``metrics`` now) decompose the *stratified* gate inputs and the query species
    only. Neither gives a per-species split of the genuine-liability count.

    A hit whose species cell is blank is counted as ``query``, matching :func:`_offtarget_views`.
    """
    if hits.empty or "_class" not in hits.columns:
        return {}
    liabilities = hits[hits["_class"].isin({c.value for c in LIABILITY})]
    if liabilities.empty:
        return {}
    if "species" not in liabilities.columns:
        return {"query": int(len(liabilities))}
    labels = liabilities["species"].astype(str).str.strip()
    labels = labels.mask(labels.eq("") | labels.eq("nan"), "query")
    return {str(species): int(n) for species, n in sorted(labels.value_counts().items())}


@dataclass(frozen=True)
class _TranscriptFacts:
    """What the run itself recorded about its transcripts, for ordering and labelling the isoform picker.

    The picker orders by candidate-window count, which is length in disguise: on one internal run it
    opened on a 7,988 nt transcript with the canonical one tenth down the list. Ordering it properly
    needs length and which transcript is canonical -- the latter a fact only the run can supply.
    Canonical status is read from the FASTA the workflow writes (``transcripts/<gene>_canonical.fasta``,
    or a ``canonical:true`` header) and is ``None`` for every transcript when the run wrote no such
    file: unknown, never inferred from sequence presence or from being the longest.
    """

    lengths: Mapping[str, int] = field(default_factory=dict)
    canonical_ids: frozenset[str] = frozenset()
    #: Which file said so, or None when this run recorded canonical status nowhere.
    source: str | None = None

    def length(self, transcript: str) -> int | None:
        """Transcript length as the run recorded it, or None."""
        return self.lengths.get(transcript)

    def canonical(self, transcript: str) -> bool | None:
        """True/False when the run recorded canonical status, None when it recorded none at all."""
        return None if self.source is None else transcript in self.canonical_ids


def _transcript_facts(run_dir: Path, regions: Mapping[str, TranscriptRegions]) -> _TranscriptFacts:
    """Transcript lengths and canonical ids, from the run's own ORF report and transcript FASTAs."""
    lengths = {tid: r.length for tid, r in regions.items()}
    canonical: set[str] = set()
    sources: list[str] = []
    for path in sorted(run_dir.glob("**/transcripts/*.fasta")):
        by_filename = path.name.endswith("_canonical.fasta")
        named_here = False
        with path.open() as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                fields = header.split()
                if not fields:
                    continue
                transcript = fields[0]
                stated = next((f.split(":", 1)[1] for f in fields if f.startswith("length:")), None)
                declared = _num(stated)
                if declared is not None and transcript not in lengths:
                    lengths[transcript] = int(declared)
                if by_filename or "canonical:true" in header.lower():
                    canonical.add(transcript)
                    named_here = True
        if named_here or by_filename:
            sources.append(path.name)
    return _TranscriptFacts(lengths, frozenset(canonical), ", ".join(sources) or None)


def _scope_label(descriptor: Any) -> str:
    """Human-readable scope, emitted once per filter rather than once per guide."""
    label = ",".join(sorted(descriptor.scope.species)) or "all species"
    if descriptor.scope.max_mismatches is not None:
        label += f", nm<={descriptor.scope.max_mismatches}"
    return label


def _gene_query(manifest: Mapping[str, Any], candidates: pd.DataFrame) -> str:
    """The run's target: the manifest's own ``gene_query``, else a gene column, else unknown.

    The manifest is authoritative: ``candidates_all.csv`` carries no gene column at all.
    """
    stated = str(manifest.get("gene_query") or "").strip()
    if stated:
        return stated
    for column in ("gene_name", "gene_symbol", "gene_id"):
        if column in candidates.columns:
            values = candidates[column].dropna().astype(str)
            values = values[values.str.strip().ne("")]
            if not values.empty:
                names = sorted(values.unique())
                return names[0] if len(names) == 1 else f"{names[0]} +{len(names) - 1} more"
    return "unknown"


def _verdict_agreement(candidates: pd.DataFrame, guides: list[GuideEntry]) -> dict[str, Any]:
    """Compare the report's independent verdicts against the run's own ``passes_filters`` column.

    Published in the report header, because a reader has no other way to tell whether the gate panel
    agrees with the pipeline. The asymmetry is the point: the report must never contradict a run PASS,
    while a run failure it cannot re-derive is reported ``unknown`` rather than quietly flipped.

    Every count here is over the guides the report **embedded**: comparing against the run's whole
    candidate table published "0/1424 run PASS contradicted" from a report holding 1,383 of those
    1,424 guides -- a coverage claim over 41 guides the reader cannot open.
    """
    if "passes_filters" not in candidates.columns:
        return {"comparable": False, "reason": "the run exports no passes_filters column"}
    embedded = {g.guide for g in guides}
    run_pass = set(candidates.loc[candidates["passes_filters"] == "PASS", "_guide"]) & embedded
    by_status = {g.guide: g.status for g in guides}
    contradicted = sorted(g for g in run_pass if by_status.get(g) == "fail")
    not_rederivable = sorted(g for g, s in by_status.items() if s == "unknown" and g not in run_pass)
    overruled = sorted(g.guide for g in guides if g.status in ("pass", "warn") and g.guide not in run_pass)
    return {
        "comparable": True,
        "run_pass_guides": len(run_pass),
        "contradicted_run_pass": len(contradicted),
        "run_failed_not_rederivable": len(not_rederivable),
        # Must stay 0. It counts guides the report calls clean while the run rejected them -- the
        # fabricated-evidence direction, which no earlier metric could see.
        "overruled_run_fail": len(overruled),
        "undeclared_run_rejections": sum(1 for g in guides if g.undeclared_run_rejection),
    }


def _read_hits(path: Path, caveats: list[str]) -> pd.DataFrame:
    """Read a published hit table, requiring the classification columns the report displays."""
    empty: pd.DataFrame = pd.DataFrame()
    if not path.exists():
        caveats.append(f"no hit table at {path.name}; every off-target panel is empty for want of input")
        return empty
    hits: pd.DataFrame = pd.read_csv(path, sep="\t", low_memory=False)
    missing = [c for c in CLASSIFICATION_COLUMNS if c not in hits.columns]
    if missing:
        raise ReportInputError(
            f"{path.name} is missing {missing}. The report will not classify hits itself -- it would then "
            "be able to disagree with the run. Re-publish the table through the workflow, which writes "
            "these columns, and report again."
        )
        # Deliberately fatal: a table without hit_class would render on-target isoform alignments as
        # off-targets -- 50.1% of rows on the public baseline.
    unannotated = int((~hits.apply(is_annotated, axis=1)).sum()) if len(hits) else 0
    if unannotated:
        caveats.append(f"{unannotated} hit rows carry no classification and are shown as undetermined")
    return hits


@dataclass(frozen=True)
class _GatePanel:
    """The gates the report will apply, and where they came from."""

    filters: tuple[Any, ...]
    profile_name: str
    run_mode: str
    source: str


def _gate_panel(policy: ResolvedRunPolicy | None, manifest: Mapping[str, Any], caveats: list[str]) -> _GatePanel:
    """The run's own gates: the live policy, else the manifest, else library defaults with a caveat.

    Library defaults are the last resort and must be declared, because they are not this run's gates:
    defaulting silently once published a GC ceiling of 60 against a run that set 65, failing 227
    guides on a threshold the run never applied and contradicting 41 of its PASSes.
    """
    if policy is not None:
        return _GatePanel(tuple(policy.filters), policy.profile.name, policy.run_mode.value, "the live run policy")

    block = manifest.get("run_policy")
    if isinstance(block, Mapping):
        try:
            filters = filters_from_manifest(block)
        except RunPolicyError as exc:
            caveats.append(f"the manifest's gate registry is unusable ({exc}); gates below are library defaults")
        else:
            if filters:
                profile = block.get("profile") or {}
                return _GatePanel(
                    filters,
                    str(profile.get("name") or "unknown") if isinstance(profile, Mapping) else "unknown",
                    str(block.get("run_mode") or "unknown"),
                    "the run's own manifest",
                )
            caveats.append("the manifest declares no gates; gates below are library defaults")
    else:
        caveats.append("this run published no policy in its manifest; gates below are library defaults")

    fallback = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human", screen_species=["human"]
    )
    return _GatePanel(tuple(fallback.filters), fallback.profile.name, fallback.run_mode.value, "library defaults")


def _has_transcript_context(manifest: Mapping[str, Any], caveats: list[str]) -> bool:
    """Whether this run's candidate rows sit on real transcripts, per the run's own declared entry point.

    ``sirnaforge offtarget`` screens pre-designed guides, so it enumerates nothing: every row carries a
    placeholder transcript id and position 1 (#100). Read as an enumeration, that fabricates an isoform
    table, a design map and a coverage denominator for target context this entry point does not have --
    so all three are withheld for such a run.

    Decided from the manifest's declared entry point, not by sniffing the placeholder value: the entry
    point is a fact the run published about itself, while a sentinel string is a coincidence a reader
    cannot verify. A run with no readable manifest declares nothing, and is read as a design run -- the
    same assumption every other provenance field already makes.
    """
    entry_point = str((manifest.get("run_policy") or {}).get("entry_point") or "")
    if entry_point != EntryPoint.OFFTARGET_ONLY.value:
        return True
    caveats.append(
        "this run screened pre-designed guides, which carry no transcript context: no isoform table, "
        "no design map and no coverage denominator are reported"
    )
    return False


def build_payload(
    run_dir: Path | str,
    *,
    policy: ResolvedRunPolicy | None = None,
    max_guides: int | None = DEFAULT_MAX_EMBEDDED_GUIDES,
) -> ReportPayload:
    """Build the payload for a finished run directory.

    Every count in ``run`` describes the guides this payload **contains**, and ``guides_total`` /
    ``guides_dropped`` / ``guides_dropped_by_status`` say what a cap cost. That is the whole point of
    capping here rather than downstream: the counts and the guides cannot disagree.

    Args:
        run_dir: A completed output directory. No live pipeline state is required.
        policy: Resolved policy supplying the gate descriptors. Read from the run's own
            ``manifest.json`` when omitted.
        max_guides: Guides to embed, worst verdicts dropped first
            (:data:`DEFAULT_MAX_EMBEDDED_GUIDES`). ``None`` embeds the run whole.

    Returns:
        A :class:`ReportPayload`.

    Raises:
        ReportInputError: The directory has no candidate table, or a hit table with no
            classification columns.
    """
    run_dir = Path(run_dir)
    caveats: list[str] = []

    candidates_csv = next(iter(sorted(run_dir.glob("**/candidates_all.csv"))), None)
    if candidates_csv is None:
        raise ReportInputError(f"no candidates_all.csv under {run_dir}")
    candidates = pd.read_csv(candidates_csv, low_memory=False)

    agg = run_dir / "off_target" / "results" / "aggregated"
    hits = _read_hits(agg / "combined_offtargets.tsv", caveats)
    mirna_path = agg / "combined_mirna_hits.tsv"
    mirna = pd.read_csv(mirna_path, sep="\t", low_memory=False) if mirna_path.exists() else pd.DataFrame()

    manifest_path = next(iter(sorted(run_dir.glob("**/manifest.json"))), None)
    manifest: dict[str, Any] = {}
    if manifest_path is not None:
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:  # a bad manifest costs provenance, not the report
            caveats.append(f"manifest.json could not be read ({exc}); the header shows less provenance")

    panel = _gate_panel(policy, manifest, caveats)
    transcript_context = _has_transcript_context(manifest, caveats)

    candidates["_guide"] = candidates["guide_sequence"].map(_normalise_guide)
    if not hits.empty:
        hits["_guide"] = hits["qseq"].map(_normalise_guide)
        hits["_nm"] = pd.to_numeric(hits["nm"], errors="coerce")
        hits["_class"] = hits.apply(lambda r: hit_class_of(r).value, axis=1)
    if not mirna.empty and "qseq" in mirna.columns:
        mirna["_guide"] = mirna["qseq"].map(_normalise_guide)

    hits_by_guide = dict(tuple(hits.groupby("_guide"))) if not hits.empty else {}
    mirna_by_guide = dict(tuple(mirna.groupby("_guide"))) if not mirna.empty and "_guide" in mirna else {}

    descriptors = [f.descriptor for f in panel.filters]
    populated = {c for c in candidates.columns if candidates[c].notna().any()}
    gate_columns = [observed_column(d, populated) for d in descriptors]
    register = _register_index(candidates)
    clusters = _register_clusters(candidates)
    regions = transcript_regions(run_dir) if transcript_context else {}
    facts = _transcript_facts(run_dir, regions) if transcript_context else _TranscriptFacts()
    guides: list[GuideEntry] = []

    # Sorted so `best` is the guide's best-scoring enumeration rather than whichever row the CSV
    # happens to list first; the gates read constant-per-guide columns either way.
    by = [c for c in ("composite_score", "design_score") if c in candidates.columns]
    ranked = candidates.sort_values(by, ascending=False, na_position="last") if by else candidates
    for guide, rows in ranked.groupby("_guide", sort=False):
        guides.append(
            _build_guide(
                guide=str(guide),
                rows=rows,
                best=rows.iloc[0],
                descriptors=descriptors,
                gate_columns=gate_columns,
                hits=hits_by_guide.get(guide, pd.DataFrame()),
                mirna=mirna_by_guide.get(guide, pd.DataFrame()),
                register=register,
                clusters=clusters,
                facts=facts,
                transcript_context=transcript_context,
            )
        )

    guides.sort(key=lambda g: (-(g.composite_score or g.design_score or -1), g.guide))

    # From here on `embedded` is the report. `guides` stays available for the design map, which plots
    # every candidate row whether or not its guide fits, and for the cap's own arithmetic.
    embedded, dropped_by_status = _select_embedded(guides, max_guides)
    n_dropped = len(guides) - len(embedded)
    if n_dropped:
        lost = ", ".join(f"{n} {s}" for s, n in dropped_by_status.items() if n)
        shippable = sum(n for s, n in dropped_by_status.items() if s != "fail")
        caveats.append(
            f"this report embeds {len(embedded)} of the run's {len(guides)} guides, and every count in it "
            f"-- header, pills, gate tallies and the agreement line -- counts those {len(embedded)}; the "
            f"{n_dropped} not embedded ({lost}) are the lowest-scoring of their verdict"
            + (
                f", including {shippable} a reader might have acted on: raise max_guides to see them"
                if shippable
                else ""
            )
        )

    undeclared = sorted({g.run_verdict for g in embedded if g.undeclared_run_rejection if g.run_verdict})
    if undeclared:
        caveats.append(
            f"{sum(1 for g in embedded if g.undeclared_run_rejection)} guides were rejected by the run as "
            f"{', '.join(undeclared)}, which no declared gate expresses; they are reported not "
            "established rather than clean"
        )

    # Read off the hit table rather than summed over guides, so it stays a run-wide alignment count in
    # the same scope as `hit_rows` even when the cap dropped guides that carried some of those rows.
    n_liab = int(hits["_class"].isin({c.value for c in LIABILITY}).sum()) if not hits.empty else 0
    effects = _gate_effects(embedded, len(panel.filters))
    agreement = _verdict_agreement(candidates, embedded)
    run = {
        "candidate_rows": sum(g.n_rows for g in embedded),
        "candidate_rows_total": int(len(candidates)),
        "guides": len(embedded),
        "guides_total": len(guides),
        "guides_dropped": n_dropped,
        "guides_dropped_by_status": dropped_by_status,
        "guide_embed_limit": max_guides,
        "hit_rows": int(len(hits)),
        "liability_rows": n_liab,
        "non_liability_rows": int(len(hits)) - n_liab,
        # Run-wide, and roughly half non-human on a four-species screen. Published because the gate that
        # rejects the most guides counts every species at once; see _liability_by_species.
        "liability_rows_by_species": _liability_by_species(hits),
        "screened_species": sorted(
            {s for s in hits["species"].astype(str).str.strip().unique() if s and s != "nan"}
            if not hits.empty and "species" in hits.columns
            else set()
        ),
        "mirna_rows": int(len(mirna)),
        "embed_scope": f"{EMBED_SPECIES}, nm<={EMBED_MAX_NM}",
        "gene_query": _gene_query(manifest, candidates),
        "status_counts": {s: sum(1 for g in embedded if g.status == s) for s in STATUSES},
        "agreement": agreement,
        # The isoform-coverage denominator. Every transcript the design step enumerated on, which is
        # not the same as the transcripts the map can plot -- a transcript can be present with no
        # scoreable row and still be an isoform this guide either does or does not reach.
        "transcript_ids": sorted(candidates["transcript_id"].dropna().astype(str).unique())
        if transcript_context and "transcript_id" in candidates.columns
        else [],
        "canonical_transcript_ids": sorted(facts.canonical_ids),
        # None means this run recorded canonical status nowhere, so every `canonical` field is unknown
        # rather than false. Nothing here guesses it from length or from sequence presence.
        "canonical_source": facts.source,
        "transcripts": _transcript_maps(candidates, guides, regions, facts, caveats) if transcript_context else [],
        # Keyed by dot-bracket and computed once per distinct structure, which is what makes the
        # layouts small enough to embed: 40,079 candidates carry 1,333 distinct structures.
        "structure_layouts": layouts_for(g.structure for g in embedded),
    }
    return ReportPayload(
        schema_version=PAYLOAD_SCHEMA_VERSION,
        run=run,
        filters=[
            _filter_view(
                f.descriptor,
                column=column,
                gates=[g.gates[i] for g in embedded],
                setting_key=f.setting_key,
                definition=f.definition,
                effect=effects[i],
            )
            for i, (f, column) in enumerate(zip(panel.filters, gate_columns, strict=True))
        ],
        guides=embedded,
        provenance={
            "run_dir": str(run_dir),
            "candidates_csv": str(candidates_csv.relative_to(run_dir)),
            "manifest_present": manifest_path is not None,
            "tool_version": str(manifest.get("tool_version") or "unknown"),
            "run_timestamp": str(manifest.get("run_timestamp") or "unknown"),
            "policy_profile": panel.profile_name,
            "run_mode": panel.run_mode,
            "policy_source": panel.source,
            "payload_schema_version": PAYLOAD_SCHEMA_VERSION,
        },
        caveats=caveats,
    )


#: Transcript stretches shorter than this are not reported as carrying no candidate: a 23-mer window
#: cannot start in the last 22 nt of a transcript, so short tails are arithmetic, not a design choice.
MIN_UNCOVERED_NT = 40

#: Value plotted on the design map, in order of preference. ``composite_score`` does not exist until
#: post-screen scoring, so a design-only run or one whose scoring failed has only ``design_score``.
#: Without the fallback those runs lose the whole map to an all-NaN column.
_MAP_VALUE_COLUMNS = ("composite_score", "design_score")


def _transcript_maps(
    candidates: pd.DataFrame,
    guides: list[GuideEntry],
    regions: Mapping[str, TranscriptRegions],
    facts: _TranscriptFacts,
    caveats: list[str],
) -> list[dict[str, Any]]:
    """Per-transcript position series for the design map, most-enumerated transcript first.

    A row the run rejected is plotted as a rejection; every other row carries its guide's report
    status, so a window the report cannot establish is never drawn as passing. Regions come from the
    run's ORF report; a transcript missing from it still gets a map, without a region bar.

    ``canonical`` rides along since window count (length) is the wrong ordering key for the isoform
    picker built from this list -- see :class:`_TranscriptFacts`. The order here is unchanged; a
    consumer now has what it needs to reorder it.

    Statuses are read from **all** the run's guides, not only the embedded ones: every candidate row is
    plotted either way, so dropping a guide from the payload must not recolour its windows.
    """
    required = {"position", "transcript_id", "passes_filters"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        caveats.append(f"no design map: this run exports no {', '.join(missing)}")
        return []
    value_column = next((c for c in _MAP_VALUE_COLUMNS if c in candidates.columns), None)
    if value_column is None:
        caveats.append(f"no design map: this run exports none of {', '.join(_MAP_VALUE_COLUMNS)}")
        return []

    if not regions:
        caveats.append("no ORF report in this run, so the design map cannot label CDS and UTR")

    status_by_guide = {g.guide: g.status for g in guides}
    ids = candidates["id"].astype(str) if "id" in candidates.columns else pd.Series("", index=candidates.index)
    rows = candidates[~ids.str.contains("DIRTY", na=False)].copy()
    rows["_pos"] = pd.to_numeric(rows["position"], errors="coerce")
    rows["_val"] = pd.to_numeric(rows[value_column], errors="coerce")
    rows = rows.dropna(subset=["_pos", "_val"])
    run_pass = rows["passes_filters"].astype(str).str.split(" (", regex=False).str[0].eq("PASS")
    rows["_class"] = [
        "fail" if not passed else status_by_guide.get(guide, "unknown")
        for guide, passed in zip(rows["_guide"], run_pass, strict=True)
    ]

    # The run exports no window-length column, so take it from the guides themselves.
    window = int(rows["guide_sequence"].astype(str).str.len().max() or 0) if len(rows) else 0
    out: list[dict[str, Any]] = []
    for tid, group in rows.groupby("transcript_id", sort=False):
        region = regions.get(str(tid))
        length = region.length if region else int(group["_pos"].max()) + max(window - 1, 0)
        series = {
            name: [[int(p), round(float(v), 3)] for p, v in zip(part["_pos"], part["_val"], strict=True)]
            for name, part in group.groupby("_class", sort=False)
        }
        out.append(
            {
                "transcript_id": str(tid),
                "length": int(length),
                "length_stated": facts.length(str(tid)) is not None,
                "canonical": facts.canonical(str(tid)),
                "cds_start": region.cds_start if region else None,
                "cds_end": region.cds_end if region else None,
                "windows": int(len(group)),
                "value_column": value_column,
                "series": series,
                "gaps": [
                    [start, end]
                    for start, end in uncovered_stretches(
                        (int(p) for p in group["_pos"]), int(length), window=max(window, 1), min_nt=MIN_UNCOVERED_NT
                    )
                ]
                if window
                else [],
            }
        )
    out.sort(key=lambda d: (-d["windows"], d["transcript_id"]))
    return out


def _build_guide(
    *,
    guide: str,
    rows: pd.DataFrame,
    best: pd.Series,
    descriptors: list[Any],
    gate_columns: list[str | None],
    hits: pd.DataFrame,
    mirna: pd.DataFrame,
    register: dict[str, list[int]],
    clusters: dict[str, dict[str, Any]],
    facts: _TranscriptFacts,
    transcript_context: bool = True,
) -> GuideEntry:
    gates = [list(_evaluate(d, best, c)) for d, c in zip(descriptors, gate_columns, strict=True)]

    # No transcript context means no enumeration to report; see _has_transcript_context.
    isoforms = _isoform_table(rows, register, clusters, facts) if transcript_context else []
    by_symbol, matrix, embedded, liability = _offtarget_views(hits)

    counts_exist = bool(len(hits)) and not embedded
    return GuideEntry(
        run_verdict=_run_verdict(rows),
        transcript_hits=_transcript_hits(rows),
        off_target_screened=_flag(best.get("off_target_screened")),
        screen_query_id=(str(best.get("screen_query_id")) if pd.notna(best.get("screen_query_id")) else None),
        structure=(str(best.get("structure")) if pd.notna(best.get("structure")) else None),
        guide=guide,
        passenger=(str(best.get("passenger_sequence")) if pd.notna(best.get("passenger_sequence")) else None),
        overhang=(str(best.get("passenger_overhang")) if pd.notna(best.get("passenger_overhang")) else None),
        modifications=(
            str(best.get("passenger_modifications")) if pd.notna(best.get("passenger_modifications")) else None
        ),
        n_rows=int(len(rows)),
        design_score=_num(best.get("design_score")),
        composite_score=_num(best.get("composite_score")),
        weight_vector=(str(best.get("weight_vector")) if pd.notna(best.get("weight_vector")) else None),
        metrics={
            k: _num(best.get(k))
            for k in (
                "gc_content",
                "asymmetry_score",
                "paired_fraction",
                "mfe",
                "empirical_score",
                "target_accessibility_p",
                "melting_temperature",
                "off_target_count",
                "undetermined_hits",
                "transcriptome_hits_total",
                "mirna_hits_0mm_seed",
                # The query-species gate inputs #101 added: the only species decomposition the row
                # carries, and what the human-stratified gates actually compared -- the all-species
                # columns above answer a different question. Absent from an older run's CSV, hence the
                # conditional (see _liability_by_species for the rest).
                "total_offtarget_hits_query",
                "transcriptome_hits_0mm_query",
                "transcriptome_hits_1mm_query",
                "transcriptome_hits_2mm_query",
                "transcriptome_hits_seed_0mm",
                "mirna_hits_0mm_seed_query",
                "mirna_hits_high_risk_query",
            )
            if k in best.index
        },
        gates=gates,
        isoforms=isoforms,
        offtarget_by_symbol=by_symbol,
        offtarget_matrix=matrix,
        offtarget_rows=embedded,
        offtarget_embedded_scope_empty_but_counts_exist=counts_exist,
        liability_count=liability,
        liability_by_species=_liability_by_species(hits),
        mirna=_mirna_table(mirna),
        ortholog=_ortholog_conservation(hits),
    )


#: Positions this close on one transcript share a window, so their designs can score identically.
REGISTER_NEIGHBOUR_NT = 2


def _register_index(candidates: pd.DataFrame) -> dict[str, list[int]]:
    """Transcript -> every enumerated position, across **all** guides.

    Register neighbours are cross-guide by nature: the point is that two *different* designs one
    nucleotide apart share a window. Built once from the whole candidate table; a per-guide map would
    see only its own rows and report no neighbours at all.
    """
    index: dict[str, list[int]] = defaultdict(list)
    for tx, pos in zip(candidates.get("transcript_id", []), candidates.get("position", []), strict=False):
        p = _num(pos)
        if p is not None:
            index[str(tx or "")].append(int(p))
    return index


def _register_clusters(candidates: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Cross-guide register clusters: candidate id -> its cluster and whether it is the representative.

    ``register_neighbours`` (below) names positions, not designs, so a reader cannot map a neighbour
    back to the guide that owns it or pick which member of a near-duplicate cluster to keep for a
    register-deduplicated view. This runs the same "positions within :data:`REGISTER_NEIGHBOUR_NT`"
    clustering over the whole candidate table, and names each run's best-scoring member as the
    representative a dedup view should keep. A candidate the run never enumerated a position for gets
    no entry.
    """
    score_column = next((c for c in _MAP_VALUE_COLUMNS if c in candidates.columns), None)
    n = len(candidates)
    ids = candidates["id"] if "id" in candidates.columns else pd.Series([""] * n, index=candidates.index)
    txs = candidates["transcript_id"] if "transcript_id" in candidates.columns else pd.Series([""] * n)
    positions = candidates["position"] if "position" in candidates.columns else pd.Series([None] * n)
    scores = candidates[score_column] if score_column else pd.Series([None] * n)

    by_tx: dict[str, list[tuple[int, str, float]]] = defaultdict(list)
    for cid, tx, pos, score in zip(ids, txs, positions, scores, strict=True):
        p = _num(pos)
        if p is None:
            continue
        s = _num(score)
        by_tx[str(tx or "")].append((int(p), str(cid or ""), s if s is not None else float("-inf")))

    out: dict[str, dict[str, Any]] = {}
    for tx, entries in by_tx.items():
        entries.sort(key=lambda e: e[0])
        window: list[tuple[int, str, float]] = []
        for entry in entries:
            if window and entry[0] - window[-1][0] > REGISTER_NEIGHBOUR_NT:
                _flush_cluster(tx, window, out)
                window = []
            window.append(entry)
        if window:
            _flush_cluster(tx, window, out)
    return out


def _flush_cluster(tx: str, members: list[tuple[int, str, float]], out: dict[str, dict[str, Any]]) -> None:
    """One connected run of within-window positions on one transcript; its best scorer is the keeper."""
    cluster_id = f"{tx}:{members[0][0]}"
    representative = max(members, key=lambda m: m[2])[1]
    for _pos, cid, _score in members:
        out[cid] = {"register_cluster": cluster_id, "register_representative": cid == representative}


def _isoform_table(
    rows: pd.DataFrame,
    register: dict[str, list[int]],
    clusters: dict[str, dict[str, Any]],
    facts: _TranscriptFacts,
) -> list[dict[str, Any]]:
    """One entry per transcript the guide was enumerated on, flagging register neighbours.

    A shorter guide starting at *s* sits inside the longer window at *s-1*, so two designs one
    nucleotide apart can share a window and receive an identical score -- a property of the
    enumeration, true for any target.

    Do not add a measured-knockdown ratio to the rendered card: a mislabelled one shipped here before
    ("1.4x", stated as knockdown but actually the *fraction-remaining* ratio between two designs of one
    internal reference set at positions 1982/1983 -- the knockdown ratio for that pair is 1.27x, and
    it wasn't even the set's strongest case, which differs 1.97x at 2733/2735). Neither number
    generalises to a target that reference set says nothing about.

    ``length`` and ``canonical`` travel with each row so a consumer can order these by canonical, then
    length, instead of by whichever transcript happens to carry the most windows. ``canonical`` is
    ``None`` when the run recorded no canonical status at all -- see :class:`_TranscriptFacts`.
    """
    out: list[dict[str, Any]] = []
    for _, r in rows.iterrows():
        tx = str(r.get("transcript_id") or "")
        pos = _num(r.get("position"))
        cid = str(r.get("id") or "")
        neighbours = [
            p
            for p in register.get(tx, [])
            if pos is not None and p != int(pos) and abs(p - int(pos)) <= REGISTER_NEIGHBOUR_NT
        ]
        # A candidate absent from `clusters` (no numeric position) is trivially its own cluster: never
        # dropped from a register-deduplicated view for want of a flag.
        cluster = clusters.get(cid, {"register_cluster": None, "register_representative": True})
        out.append(
            {
                "candidate_id": cid,
                "transcript": tx,
                "position": pos,
                "length": facts.length(tx),
                "canonical": facts.canonical(tx),
                "register_neighbours": sorted(neighbours),
                "register_cluster": cluster["register_cluster"],
                "register_representative": cluster["register_representative"],
            }
        )
    out.sort(key=lambda d: (d["transcript"], d["position"] if d["position"] is not None else -1))
    return out


def _offtarget_views(
    hits: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], int]:
    """Group hits by gene symbol, build the complete count matrix, and embed the in-scope rows."""
    if hits.empty:
        return [], [], [], 0

    liability = int(hits["_class"].isin({c.value for c in LIABILITY}).sum())

    grouped: Counter[tuple[str, str, str]] = Counter()
    min_nm: dict[tuple[str, str, str], float] = {}
    for _, r in hits.iterrows():
        symbol = str(r.get("hit_symbol") or r.get("matched_symbol") or "").strip() or _UNKNOWN_SYMBOL
        key = (symbol, str(r["_class"]), str(r.get("species") or "").strip() or "query")
        grouped[key] += 1
        nm = r["_nm"]
        if pd.notna(nm):
            min_nm[key] = min(min_nm.get(key, float(nm)), float(nm))
    by_symbol = [
        {
            "symbol": s,
            "hit_class": c,
            "species": sp,
            "n": n,
            "min_nm": min_nm.get((s, c, sp)),
            "is_liability": c in {x.value for x in LIABILITY},
        }
        for (s, c, sp), n in grouped.most_common()
    ]

    matrix_counter: Counter[tuple[str, str, str]] = Counter()
    for _, r in hits.iterrows():
        nm = r["_nm"]
        band = "unknown" if pd.isna(nm) else ("0" if nm == 0 else "1" if nm == 1 else "2" if nm == 2 else ">=3")
        matrix_counter[(str(r.get("species") or "query"), band, str(r["_class"]))] += 1
    matrix = [
        {"species": sp, "nm": band, "hit_class": cls, "n": n} for (sp, band, cls), n in sorted(matrix_counter.items())
    ]

    # Row-level drill-down exists to inspect liabilities, so only liability rows in the embedded scope
    # carry per-row detail. On-target isoform, ortholog and repeat alignments are still fully
    # displayed -- by gene in offtarget_by_symbol and completely in offtarget_matrix -- just not
    # per-row. Embedding every class cost 14.7 MB of the first render's 38 MB.
    in_scope = hits[
        (hits["species"].astype(str).str.strip().isin({EMBED_SPECIES, ""}))
        & (hits["_nm"].notna())
        & (hits["_nm"] <= EMBED_MAX_NM)
        & (hits["_class"].isin({c.value for c in LIABILITY}))
    ]
    embedded = [
        {
            "t": str(r.get("rname") or ""),
            "s": str(r.get("hit_symbol") or "").strip() or _UNKNOWN_SYMBOL,
            "c": str(r["_class"]),
            "nm": _num(r.get("nm")),
            "sm": _num(r.get("seed_mismatches")),
            "cig": str(r.get("cigar") or ""),
            "st": str(r.get("strand") or ""),
        }
        for _, r in in_scope.iterrows()
    ]
    embedded.sort(key=lambda d: (d["nm"] if d["nm"] is not None else 99, d["s"]))
    return by_symbol, matrix, embedded, liability


#: Columns the miRNA name may arrive under. ``mirna_id`` is what the aggregate actually publishes
#: (``Hsa-Mir-24-P2_3p``); the table has no ``rname``, so reading only that returned an empty string
#: for every row and listed 18,078 anonymous seed matches. Which miRNA is mimicked is the whole
#: question -- a match to a cardiac or neuronal family is a different finding from one to an
#: unexpressed paralogue, and a retraction in this programme turned on exactly that.
_MIRNA_NAME_COLUMNS = ("mirna_id", "rname", "mirna", "name")


def _ortholog_conservation(hits: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Best ortholog alignment per species: how conserved this guide's site is, per species.

    Cross-species conservation is already in the screen and is thrown away twice over. The screen
    aligns each guide against every requested species' transcriptome and classifies a hit on the
    orthologous gene as ``ortholog`` -- explicitly *not* a liability -- and the candidate table then
    summarises those hits as ``conservation_score``, i.e. ``(species hit) / 3``.

    That summary cannot answer what a cross-species programme asks, for two reasons. It is
    **mismatch-blind**: on one internal run it counted a species as conserved on alignments up to 8
    mismatches, and mouse ortholog hits ran 1,413 at nm=0 against 1,493 at nm>=3. And it is
    **seed-blind**: mouse nm=1 split 96 seed-intact against 79 seed-hit, and one mismatch outside
    positions 2-8 is a different molecule from one inside them -- allowing it took the
    mouse-and-macaque pool from 49 guides to 113.

    So the per-species best ``(nm, seed_mismatches)`` is published instead and thresholding is left to
    the reader: requiring "perfect in mouse and macaque" is a programme decision, not a property of
    the target, and hard-coding it as a gate would put one programme's requirement in every run.
    """
    if hits.empty or "_class" not in hits.columns:
        return {}
    ortholog = hits[hits["_class"].eq(HitClass.ORTHOLOG.value)]
    if ortholog.empty:
        return {}
    out: dict[str, dict[str, Any]] = {}
    for species, group in ortholog.groupby(ortholog["species"].astype(str).str.strip()):
        nm = pd.to_numeric(group["nm"], errors="coerce")
        seed = pd.to_numeric(group.get("seed_mismatches"), errors="coerce")
        order = nm.fillna(99) * 100 + seed.fillna(99)
        best = order.idxmin()
        out[species or "query"] = {"nm": _num(nm.get(best)), "seed_mismatches": _num(seed.get(best))}
    return out


def _mirna_table(mirna: pd.DataFrame) -> list[dict[str, Any]]:
    """Named miRNA seed hits at 0 and 1 seed mismatch, with the database and match coordinate."""
    if mirna.empty:
        return []
    name_column = next((c for c in _MIRNA_NAME_COLUMNS if c in mirna.columns), None)
    out: list[dict[str, Any]] = []
    for _, r in mirna.iterrows():
        sm = _num(r.get("seed_mismatches"))
        if sm is not None and sm > 1:
            continue
        out.append(
            {
                "mirna": str(r.get(name_column) or "") if name_column else "",
                "source": str(r.get("species") or r.get("source") or ""),
                "database": str(r.get("database") or ""),
                "seed_mismatches": sm,
                "nm": _num(r.get("nm")),
                # The seed offset. A v0.5.1 defect counted a guide-seed motif matching anywhere on a
                # miRNA as a perfect seed hit, and coord != 1 is what distinguished the artifact from
                # the real thing, so the report shows it rather than asking a reader to trust the fix.
                "coord": _num(r.get("coord")),
            }
        )
    out.sort(key=lambda d: (d["seed_mismatches"] if d["seed_mismatches"] is not None else 99, d["mirna"]))
    return out
