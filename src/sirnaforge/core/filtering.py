"""Pure gate evaluator: one place that turns (threshold, action, observed, evidence) into a verdict.

Five rules, applied once instead of scattered across call sites: an off filter never decides, a
filter with no threshold cannot decide either, a missing measurement is unknown rather than a pass, a
lower-bound count that would otherwise pass is undecidable on incomplete evidence, and only a
``FAIL``-action gate rejects.

Pure by construction: imports only :mod:`sirnaforge.models.policy`, takes plain values in, returns
plain records out, mutates nothing. ``SiRNACandidate.record_filter_verdict`` stays the single writer
to the candidate -- this module decides what to write, not where.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TypeVar

from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation, ScreeningChannel

#: Display/precedence order for the four verdicts. Not read within this module -- exists so a
#: consumer sorting/grouping verdicts has one canonical order. NOT_EVALUATED/UNKNOWN sit first since
#: neither is a claim.
EVALUATION_ORDER: tuple[FilterEvaluation, ...] = (
    FilterEvaluation.NOT_EVALUATED,
    FilterEvaluation.UNKNOWN,
    FilterEvaluation.PASS,
    FilterEvaluation.FAIL,
)


@dataclass(frozen=True, slots=True)
class GateSpec:
    """One gate as configured, independent of any candidate.

    ``evidence_pairs`` is the evidence a pass claim depends on: channel x species pairs whose
    completion a "did not exceed the threshold" reading must have before it can be trusted. A
    design-stage gate (no screening input) declares none, so evaluation never depends on
    ``complete_pairs``.

    Attributes:
        filter_id: Stable machine identity, also the key ``observed``/``status_for`` are read by.
        threshold: Value compared against. ``None`` means the gate has nothing to decide with.
        action: What the gate does with its own verdict.
        comparator: The passing direction.
        channels: Screening channels this gate's count is drawn from, for the caller's bookkeeping.
        evidence_pairs: Channel x species pairs whose completion a pass claim depends on.
    """

    filter_id: str
    threshold: float | None
    action: FilterAction
    comparator: FilterComparator
    channels: frozenset[ScreeningChannel] = frozenset()
    evidence_pairs: frozenset[tuple[str, str]] = frozenset()


@dataclass(frozen=True, slots=True)
class FilterOutcome:
    """What one gate decided about one observation.

    Attributes:
        filter_id: The gate this outcome belongs to.
        evaluation: The verdict.
        observed: The value the gate compared, or ``None`` when nothing was observed or the
            evaluation was suppressed to ``UNKNOWN`` on incomplete evidence.
        action: The action this gate was configured with.
        rejects: Whether this outcome, on its own, rejects the candidate.
        undecidable: Whether the evaluation is ``UNKNOWN`` -- the gate was in force but could not
            be decided, as distinct from ``NOT_EVALUATED``, which makes no claim at all.
    """

    filter_id: str
    evaluation: FilterEvaluation
    observed: float | None
    action: FilterAction
    rejects: bool
    undecidable: bool


def evaluate_gate(
    spec: GateSpec,
    observed: float | None,
    *,
    complete_pairs: frozenset[tuple[str, str]] | None = None,
) -> FilterOutcome:
    """Apply one gate's rules to one observation.

    Four ways to not reach a decided PASS/FAIL, checked in order: no threshold or an OFF action means
    no claim (``NOT_EVALUATED``); no measurement means nothing to compare (``UNKNOWN``); a
    naively-passing count whose required evidence is incomplete is a lower bound that cannot show a
    ceiling was respected (``UNKNOWN``, observed cleared); everything else is decided from the
    comparator, and only a decided ``FAIL`` on a ``FAIL``-action gate rejects. A naively-*failing*
    count is decided even on incomplete evidence -- a known failure is still a failure.
    """
    if spec.threshold is None or spec.action is FilterAction.OFF:
        return FilterOutcome(spec.filter_id, FilterEvaluation.NOT_EVALUATED, observed, spec.action, False, False)

    if observed is None:
        return FilterOutcome(spec.filter_id, FilterEvaluation.UNKNOWN, None, spec.action, False, True)

    exceeds = not spec.comparator.passes(observed, spec.threshold)
    if not exceeds and complete_pairs is not None and not spec.evidence_pairs <= complete_pairs:
        return FilterOutcome(spec.filter_id, FilterEvaluation.UNKNOWN, None, spec.action, False, True)

    evaluation = FilterEvaluation.FAIL if exceeds else FilterEvaluation.PASS
    rejects = exceeds and spec.action is FilterAction.FAIL
    return FilterOutcome(spec.filter_id, evaluation, observed, spec.action, rejects, False)


def evaluate_gates(
    specs: Sequence[GateSpec],
    observed: Mapping[str, float | None],
    *,
    complete_pairs: frozenset[tuple[str, str]] | None = None,
) -> tuple[FilterOutcome, ...]:
    """Evaluate every gate against its own observation, in declared order.

    Every gate is evaluated regardless of an earlier rejection: returning early would make a later
    gate's reported count depend on list order rather than what actually happened.
    """
    return tuple(evaluate_gate(spec, observed.get(spec.filter_id), complete_pairs=complete_pairs) for spec in specs)


def first_rejection(outcomes: Sequence[FilterOutcome]) -> FilterOutcome | None:
    """The first outcome that rejects on its own, in the sequence's own order, or ``None``."""
    for outcome in outcomes:
        if outcome.rejects:
            return outcome
    return None


_S = TypeVar("_S")


def derive_passes_filters(
    outcomes: Sequence[FilterOutcome],
    status_for: Mapping[str, _S],
    *,
    current: _S | bool = True,
) -> _S | bool:
    """Reproduce ``record_filter_verdict``'s first-label-wins rule over a batch of outcomes.

    ``passes_filters`` is a single label, not a count: the first ``FAIL``-action gate to reject owns
    it, and every later rejection is recorded only in ``filter_verdicts``. ``current`` is unchanged
    when nothing rejects, and is overwritten only while still passing -- a gate cannot demote a
    candidate an earlier gate already demoted.
    """
    label: _S | bool = current
    for outcome in outcomes:
        if outcome.evaluation is not FilterEvaluation.FAIL or outcome.action is not FilterAction.FAIL:
            continue
        if not is_passing(label):
            continue
        mapped = status_for.get(outcome.filter_id)
        if mapped is not None:
            label = mapped
    return label


def is_passing(passes_filters: object) -> bool:
    """Whether a ``passes_filters`` value is a passing representation.

    Three spellings must read as passing: bare ``True``, ``SiRNACandidate.FilterStatus.PASS``, and
    the CSV string ``"PASS"``. ``FilterStatus`` is a ``str`` subclass, so comparing against the
    literal covers both without importing ``models.sirna`` (this module's purity constraint).
    """
    return passes_filters is True or passes_filters == "PASS"


def unknown_filter_ids(outcomes: Sequence[FilterOutcome]) -> frozenset[str]:
    """Filter ids whose evaluation is ``UNKNOWN`` -- in force, but undecidable."""
    return frozenset(outcome.filter_id for outcome in outcomes if outcome.evaluation is FilterEvaluation.UNKNOWN)
