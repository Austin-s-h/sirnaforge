"""Pure gate evaluator: one place that turns (threshold, action, observed, evidence) into a verdict.

``workflow.py::_check_offtarget_filters`` / ``_apply_isoform_coverage_gate`` and
``core/design.py::_record_enumeration_verdicts`` / ``_flag_excess_pairing`` / ``_apply_score_filters``
each re-derive the same five rules independently: an off filter never decides, a filter with no
threshold cannot decide either, a missing measurement is unknown rather than a pass, a lower-bound
count that would otherwise pass is undecidable, and only a ``FAIL``-action gate rejects. This module
is that evaluator, so every call site applies the rules once instead of five times.

Pure by construction: it imports only :mod:`sirnaforge.models.policy`, takes plain values in and
returns plain records out, and mutates nothing. The candidate is still mutated by
``SiRNACandidate.record_filter_verdict``, which stays the single writer -- this module decides what
to write, not where.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TypeVar

from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation, ScreeningChannel

#: Display/precedence order for the four verdicts. Not read by anything in this module -- it exists
#: because a consumer sorting or grouping verdicts needs one canonical order instead of inventing its
#: own, and NOT_EVALUATED/UNKNOWN sit before the two decided outcomes because neither is a claim.
EVALUATION_ORDER: tuple[FilterEvaluation, ...] = (
    FilterEvaluation.NOT_EVALUATED,
    FilterEvaluation.UNKNOWN,
    FilterEvaluation.PASS,
    FilterEvaluation.FAIL,
)


@dataclass(frozen=True, slots=True)
class GateSpec:
    """One gate as configured, independent of any candidate.

    ``channels``/``evidence_pairs`` are the evidence a pass claim depends on: the channel x species
    pairs whose completion this gate's count requires before a "did not exceed the threshold" reading
    may be trusted. A design-stage gate (no screening input) declares neither, so its evaluation never
    depends on ``complete_pairs``.

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

    Four ways to not reach a decided PASS/FAIL, checked in order: no threshold or an OFF action
    means the gate makes no claim (``NOT_EVALUATED``); no measurement means nothing to compare
    (``UNKNOWN``); a naively-passing count whose required evidence is incomplete is a lower bound,
    which cannot show a ceiling was respected (``UNKNOWN``, observed cleared); everything else is
    decided from the comparator, and only a decided ``FAIL`` on a ``FAIL``-action gate rejects. A
    naively-*failing* count is decided even on incomplete evidence -- a known failure is still a
    failure regardless of what else was never measured.
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

    Every gate is evaluated regardless of an earlier rejection -- returning early on the first
    rejection is the defect this module exists to retire, since it left every later gate's reported
    count a function of list order rather than of what actually happened.
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

    ``passes_filters`` is a single label, not a count, so the first ``FAIL``-action gate to reject
    owns it and every later rejection is recorded only in ``filter_verdicts``. ``current`` is returned
    unchanged when nothing rejects, and is only ever overwritten while it is still a passing label --
    a gate cannot demote a candidate that an earlier gate already demoted.
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

    Three spellings must all read as passing: the bare bool ``True``, ``SiRNACandidate.FilterStatus
    .PASS``, and the CSV string ``"PASS"``. ``FilterStatus`` is a ``str`` subclass, so its ``PASS``
    member already compares equal to the literal ``"PASS"`` -- comparing against the string covers
    both without importing ``models.sirna``, which would break this module's purity.
    """
    return passes_filters is True or passes_filters == "PASS"


def unknown_filter_ids(outcomes: Sequence[FilterOutcome]) -> frozenset[str]:
    """Filter ids whose evaluation is ``UNKNOWN`` -- in force, but undecidable."""
    return frozenset(outcome.filter_id for outcome in outcomes if outcome.evaluation is FilterEvaluation.UNKNOWN)
