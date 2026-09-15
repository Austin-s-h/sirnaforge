"""core/filtering.py is the pure evaluator #100 factors the five gate rules into.

``workflow.py`` and ``core/design.py`` each re-derive, independently, whether a threshold-bearing
gate with an OFF action makes a claim, whether a missing measurement is a pass, and whether a
naively-passing lower-bound count on incomplete evidence may report clean. This module is the one
place those rules live; this file pins the table the issue's acceptance criteria describe so a
future caller (W2) can be routed through it with no behaviour change.
"""

from enum import Enum

import pytest

from sirnaforge.core.filtering import (
    Comparator,
    FilterOutcome,
    GateSpec,
    derive_passes_filters,
    evaluate_gate,
    evaluate_gates,
    first_rejection,
    is_passing,
    unknown_filter_ids,
)
from sirnaforge.models.policy import FilterAction, FilterEvaluation, ScreeningChannel


class _LegacyFilterStatus(str, Enum):
    """A minimal str-subclass enum standing in for ``SiRNACandidate.FilterStatus`` without the import."""

    PASS = "PASS"


def _spec(
    filter_id: str = "max_off_target_count",
    threshold: float | None = 5,
    action: FilterAction = FilterAction.FAIL,
    comparator: Comparator = Comparator.AT_MOST,
    channels: frozenset[ScreeningChannel] = frozenset({ScreeningChannel.TRANSCRIPTOME}),
    evidence_pairs: frozenset[tuple[str, str]] = frozenset(),
) -> GateSpec:
    """Build a GateSpec with sensible defaults, overridden per test."""
    return GateSpec(
        filter_id=filter_id,
        threshold=threshold,
        action=action,
        comparator=comparator,
        channels=channels,
        evidence_pairs=evidence_pairs,
    )


class TestEvaluateGate:
    """The five rules, table-driven against the issue's own examples."""

    def test_warn_action_records_fail_but_never_rejects(self) -> None:
        """A WARN-action gate that exceeds its threshold still records FAIL, but ``rejects`` stays False."""
        outcome = evaluate_gate(_spec(threshold=5, action=FilterAction.WARN), observed=9)
        assert outcome.evaluation is FilterEvaluation.FAIL
        assert outcome.rejects is False
        assert outcome.observed == 9

    def test_fail_action_exceeding_threshold_rejects(self) -> None:
        """The same breach under a FAIL action rejects."""
        outcome = evaluate_gate(_spec(threshold=5, action=FilterAction.FAIL), observed=9)
        assert outcome.evaluation is FilterEvaluation.FAIL
        assert outcome.rejects is True

    def test_no_threshold_is_not_evaluated_and_reports_observed(self) -> None:
        """An undeclared threshold makes no claim, but the measured value is still reported."""
        outcome = evaluate_gate(_spec(threshold=None), observed=42)
        assert outcome.evaluation is FilterEvaluation.NOT_EVALUATED
        assert outcome.observed == 42
        assert outcome.rejects is False
        assert outcome.undecidable is False

    def test_off_action_is_not_evaluated_even_with_a_threshold(self) -> None:
        """An OFF action overrides a declared threshold: still no claim, still no rejection."""
        outcome = evaluate_gate(_spec(threshold=5, action=FilterAction.OFF), observed=42)
        assert outcome.evaluation is FilterEvaluation.NOT_EVALUATED
        assert outcome.observed == 42
        assert outcome.rejects is False

    def test_missing_observed_is_unknown(self) -> None:
        """No measurement at all is UNKNOWN, never a pass."""
        outcome = evaluate_gate(_spec(threshold=5), observed=None)
        assert outcome.evaluation is FilterEvaluation.UNKNOWN
        assert outcome.observed is None
        assert outcome.rejects is False

    def test_incomplete_evidence_on_a_naive_pass_is_unknown_not_pass(self) -> None:
        """A naively-passing lower bound on incomplete evidence is undecidable, with the count cleared."""
        spec = _spec(
            threshold=5,
            action=FilterAction.FAIL,
            evidence_pairs=frozenset({("transcriptome", "human")}),
        )
        outcome = evaluate_gate(spec, observed=2, complete_pairs=frozenset())
        assert outcome.evaluation is FilterEvaluation.UNKNOWN
        assert outcome.observed is None
        assert outcome.undecidable is True
        assert outcome.rejects is False

    def test_incomplete_evidence_on_a_naive_fail_still_fails(self) -> None:
        """A known failure fails on incomplete evidence -- only the pass direction needs completeness."""
        spec = _spec(
            threshold=5,
            action=FilterAction.FAIL,
            evidence_pairs=frozenset({("transcriptome", "human")}),
        )
        outcome = evaluate_gate(spec, observed=9, complete_pairs=frozenset())
        assert outcome.evaluation is FilterEvaluation.FAIL
        assert outcome.observed == 9
        assert outcome.rejects is True
        assert outcome.undecidable is False

    def test_complete_evidence_lets_a_naive_pass_through(self) -> None:
        """Once the required pairs are complete, the same naive pass is decided as PASS."""
        spec = _spec(
            threshold=5,
            action=FilterAction.FAIL,
            evidence_pairs=frozenset({("transcriptome", "human")}),
        )
        outcome = evaluate_gate(spec, observed=2, complete_pairs=frozenset({("transcriptome", "human")}))
        assert outcome.evaluation is FilterEvaluation.PASS
        assert outcome.observed == 2
        assert outcome.rejects is False

    def test_no_complete_pairs_argument_means_every_pair_is_treated_as_evidenced(self) -> None:
        """``complete_pairs=None`` is the caller declaring it holds no per-species record at all."""
        spec = _spec(
            threshold=5,
            action=FilterAction.FAIL,
            evidence_pairs=frozenset({("transcriptome", "human")}),
        )
        outcome = evaluate_gate(spec, observed=2, complete_pairs=None)
        assert outcome.evaluation is FilterEvaluation.PASS
        assert outcome.observed == 2

    def test_at_least_comparator_fails_below_the_floor(self) -> None:
        """AT_LEAST is the min_* direction: below the floor exceeds."""
        outcome = evaluate_gate(_spec(threshold=0.5, comparator=Comparator.AT_LEAST), observed=0.2)
        assert outcome.evaluation is FilterEvaluation.FAIL
        assert outcome.rejects is True

    def test_at_least_comparator_passes_at_the_floor(self) -> None:
        """AT_LEAST passes at the floor itself, not only strictly above it."""
        outcome = evaluate_gate(_spec(threshold=0.5, comparator=Comparator.AT_LEAST), observed=0.5)
        assert outcome.evaluation is FilterEvaluation.PASS


class TestEvaluateGates:
    """Batch evaluation over a declared gate order."""

    def test_every_gate_is_evaluated_regardless_of_an_earlier_rejection(self) -> None:
        """No early return: the whole point is that later gates are not a function of list order."""
        specs = (
            _spec(filter_id="a", threshold=1, action=FilterAction.FAIL),
            _spec(filter_id="b", threshold=1, action=FilterAction.FAIL),
        )
        outcomes = evaluate_gates(specs, {"a": 9, "b": 9})
        assert [o.filter_id for o in outcomes] == ["a", "b"]
        assert all(o.evaluation is FilterEvaluation.FAIL for o in outcomes)
        assert all(o.rejects for o in outcomes)

    def test_a_filter_id_missing_from_observed_reads_as_unmeasured(self) -> None:
        """A gate whose id is absent from the observed mapping gets None, i.e. UNKNOWN."""
        outcomes = evaluate_gates((_spec(filter_id="a", threshold=1),), {})
        assert outcomes[0].evaluation is FilterEvaluation.UNKNOWN


class TestFirstRejection:
    """The single rejection callers use for the legacy ``(should_fail, status)`` shape."""

    def test_returns_the_first_rejecting_outcome_in_order(self) -> None:
        """The first rejecting outcome in sequence order wins, matching the old early-return call site."""
        outcomes = (
            FilterOutcome("a", FilterEvaluation.PASS, 1, FilterAction.FAIL, False, False),
            FilterOutcome("b", FilterEvaluation.FAIL, 9, FilterAction.FAIL, True, False),
            FilterOutcome("c", FilterEvaluation.FAIL, 9, FilterAction.FAIL, True, False),
        )
        rejection = first_rejection(outcomes)
        assert rejection is not None
        assert rejection.filter_id == "b"

    def test_none_when_nothing_rejects(self) -> None:
        """A batch with no rejecting outcome reports None, not a false rejection."""
        outcomes = (FilterOutcome("a", FilterEvaluation.PASS, 1, FilterAction.FAIL, False, False),)
        assert first_rejection(outcomes) is None


class TestDerivePassesFilters:
    """The first-label-wins rule ``record_filter_verdict`` already applies to ``passes_filters``."""

    def test_first_fail_action_rejection_wins_the_label(self) -> None:
        """The earliest rejecting outcome's mapped status owns the label; a later one does not overwrite it."""
        outcomes = (
            FilterOutcome("a", FilterEvaluation.FAIL, 9, FilterAction.FAIL, True, False),
            FilterOutcome("b", FilterEvaluation.FAIL, 9, FilterAction.FAIL, True, False),
        )
        label = derive_passes_filters(outcomes, {"a": "GC_OUT_OF_RANGE", "b": "POLY_RUNS"})
        assert label == "GC_OUT_OF_RANGE"

    def test_warn_action_never_overwrites_the_label(self) -> None:
        """A WARN-action FAIL is recorded per-gate but never becomes the candidate-level label."""
        outcomes = (FilterOutcome("a", FilterEvaluation.FAIL, 9, FilterAction.WARN, False, False),)
        label = derive_passes_filters(outcomes, {"a": "GC_OUT_OF_RANGE"})
        assert label is True

    def test_returns_current_unchanged_when_nothing_rejects(self) -> None:
        """Nothing rejects: the caller's starting label survives untouched."""
        outcomes = (FilterOutcome("a", FilterEvaluation.PASS, 1, FilterAction.FAIL, False, False),)
        assert derive_passes_filters(outcomes, {"a": "GC_OUT_OF_RANGE"}, current=True) is True

    def test_unknown_never_overwrites_the_label(self) -> None:
        """Whether UNKNOWN rejects is a selection-time policy decision, never this module's to make."""
        outcomes = (FilterOutcome("a", FilterEvaluation.UNKNOWN, None, FilterAction.FAIL, False, True),)
        assert derive_passes_filters(outcomes, {"a": "GC_OUT_OF_RANGE"}) is True

    def test_a_filter_id_with_no_mapped_status_leaves_the_label_alone(self) -> None:
        """A rejecting outcome with no entry in ``status_for`` cannot overwrite the label with nothing."""
        outcomes = (FilterOutcome("a", FilterEvaluation.FAIL, 9, FilterAction.FAIL, True, False),)
        assert derive_passes_filters(outcomes, {}) is True


class TestIsPassing:
    """The three legacy spellings of "passing" that must all keep working."""

    @pytest.mark.parametrize("value", [True, "PASS"])
    def test_accepts_every_legacy_passing_representation(self, value: object) -> None:
        """Bare bool True and the CSV string 'PASS' both read as passing."""
        assert is_passing(value) is True

    @pytest.mark.parametrize("value", [False, "FAIL", "GC_OUT_OF_RANGE", None])
    def test_rejects_everything_else(self, value: object) -> None:
        """Anything else, including a specific failure reason string, is not passing."""
        assert is_passing(value) is False

    def test_accepts_a_str_enum_member_that_compares_equal_to_pass(self) -> None:
        """FilterStatus.PASS is a str subclass, so it must read as passing without importing it here."""
        assert is_passing(_LegacyFilterStatus.PASS) is True


class TestUnknownFilterIds:
    """The set of filter ids a caller must treat as an evidence shortfall."""

    def test_collects_only_the_undecidable_ones(self) -> None:
        """Only UNKNOWN outcomes are collected; PASS/FAIL/NOT_EVALUATED are not shortfalls."""
        outcomes = (
            FilterOutcome("a", FilterEvaluation.UNKNOWN, None, FilterAction.FAIL, False, True),
            FilterOutcome("b", FilterEvaluation.PASS, 1, FilterAction.FAIL, False, False),
            FilterOutcome("c", FilterEvaluation.NOT_EVALUATED, None, FilterAction.OFF, False, False),
        )
        assert unknown_filter_ids(outcomes) == frozenset({"a"})

    def test_empty_for_no_outcomes(self) -> None:
        """An empty batch has no shortfall."""
        assert unknown_filter_ids(()) == frozenset()
