"""Pure eligibility and shortlist assembly, extracted from ``workflow.py``.

``workflow.py::_apply_post_screen_ranking`` / ``_evidence_shortfall`` / ``_required_evidence_shortfall``
/ ``_missing_required_evidence`` each held a piece of the same decision: which candidates a run can
stand behind, in what order, and why the rest were held out. This module is that decision, taking
plain views in and plain records out so it can be tested without a workflow instance and so the
Nextflow-only entry points (#100's ``run_offtarget_only_workflow``) can reach the same answer.

Adds :class:`~sirnaforge.models.sirna.SelectionState`'s ``PROVISIONAL`` member, the one behaviour
change this extraction makes: an ``EXPLORATORY`` run used to label an evidence-incomplete candidate
identically to a fully evidenced one (``eligible``). It now carries its own state, so a report can
tell the two apart and a provisional score is never silently indistinguishable from a qualified one.

Pure by construction: it imports :mod:`sirnaforge.models.policy` for the run-policy vocabulary and
:mod:`sirnaforge.models.sirna` for :class:`~sirnaforge.models.sirna.SelectionState` only, takes plain
values in, and mutates nothing. Nothing here reads a ``ScreeningPlan`` or a ``ScreeningEvidence``
directly -- it takes ``completed_pairs``, the frozenset the evidence producer derives.
"""

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

from sirnaforge.models.policy import EvidenceRequirements, RunMode, ScreeningChannel, UnknownEvidenceAction
from sirnaforge.models.sirna import SelectionState


@dataclass(frozen=True, slots=True)
class CandidateView:
    """The read-only facts about one candidate this module needs, independent of its live object.

    ``unknown_filter_ids`` is precomputed by the caller (typically
    :func:`sirnaforge.core.filtering.unknown_filter_ids` over that candidate's own gate outcomes)
    rather than re-derived here, so this module never has to know how a verdict was decided -- only
    that it was undecided.

    Attributes:
        candidate_id: Stable candidate identity. Not unique on its own: duplicate ids can reach
            selection through step5's ``finally`` block, which is why ordering and exclusion below
            key on ``ordinal``, not on this.
        ordinal: Position of this view in the batch handed to :func:`select`. The identity ordering
            and exclusion tracking actually use, because a dict keyed by ``candidate_id`` would
            collapse two candidates that happen to share one.
        passes_filters: Whether the candidate's gates, taken together, still call it passing.
        repeat_flagged: Whether repeat detection flagged this candidate.
        scored_after_screening: Whether a post-screen composite score was computed.
        off_target_screened: Whether the query-species transcriptome channel screened this candidate.
        ranking_score: The score :func:`select` ranks on.
        unknown_filter_ids: Filter ids whose evaluation on this candidate is ``UNKNOWN``.
    """

    candidate_id: str
    ordinal: int
    passes_filters: bool
    repeat_flagged: bool
    scored_after_screening: bool
    off_target_screened: bool
    ranking_score: float
    unknown_filter_ids: frozenset[str] = frozenset()


@dataclass(frozen=True, slots=True)
class SelectionInputs:
    """Everything about the run that is not a property of one candidate.

    Attributes:
        run_mode: How much evidence this run claims to have. ``None`` behaves like neither
            ``QUALIFIED`` nor ``EXPLORATORY``: nothing disqualifies and nothing is labelled
            provisional, matching a run with no resolved policy at all.
        requirements: Which screening evidence this run requires. ``None`` means no shortfall can be
            computed, matching the resolver's own guard.
        completed_pairs: Channel x species pairs with completed evidence. ``None`` means the run
            produced no evidence record at all, which reads as incomplete everywhere it is checked.
        query_species: The species whose transcriptome pair is answered from
            ``CandidateView.off_target_screened`` rather than from ``completed_pairs``, because that
            is the finest-grained record the evidence producers offer per candidate today.
        filter_channels: Which screening channels back each gate's count (``POST_SCREEN_FILTER_CHANNELS``
            in ``workflow.py``). A filter id absent from this mapping reads design-time or annotation
            evidence, so an incomplete channel never excuses it.
        repeat_rejects: Whether a repeat-flagged candidate is rejected outright, decided by the
            ``max_repeat_transcript_fraction`` gate's own action.
        top_n: How many eligible candidates the shortlist keeps. ``None`` means every one of them --
            the same convention ``DesignParameters.top_n`` carries, because that is where the value
            comes from and ``list[:None]`` needs no special-casing.
        design_input_shortfalls: Transcript ids a design-stage batch failure dropped, mapped to why.
            Additive: empty by default so a caller that has not wired #100's design-input tracking
            (W5) gets today's behaviour unchanged.
    """

    run_mode: RunMode | None
    requirements: EvidenceRequirements | None
    completed_pairs: frozenset[tuple[str, str]] | None
    query_species: str
    filter_channels: Mapping[str, frozenset[ScreeningChannel]]
    repeat_rejects: bool
    top_n: int | None
    design_input_shortfalls: Mapping[str, str] = MappingProxyType({})


@dataclass(frozen=True, slots=True)
class CandidateSelection:
    """What one candidate was decided to be, and the evidence reasons behind it.

    ``shortfall`` and ``ordering_shortfall`` are both populated regardless of ``state`` -- a
    candidate excluded by a gate before the evidence check still had an evidence position, and the
    run summary reads both independently of which reason actually applied.

    Attributes:
        candidate_id: The view's own candidate identity.
        ordinal: The view's own ordinal.
        state: The decided :class:`~sirnaforge.models.sirna.SelectionState`.
        shortfall: Required-evidence reasons this candidate is short of, whether or not that shortfall
            disqualified it (only a ``QUALIFIED`` run disqualifies on it).
        ordering_shortfall: Every declared-pair/undecided-gate reason this candidate is short of,
            whatever its requiredness -- what the sort key and the ``PROVISIONAL`` split read.
    """

    candidate_id: str
    ordinal: int
    state: SelectionState
    shortfall: tuple[str, ...]
    ordering_shortfall: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SelectionResult:
    """The resolved shortlist: ordering, membership and the summary the CLI reads.

    Attributes:
        order: Every ordinal, ranked (evidence-complete first, then scored, then by score) -- the
            order CSV export and FASTA writers use.
        eligible_ordinals: Ordinals that reached ranking (``ELIGIBLE`` or ``PROVISIONAL``), in rank
            order -- the shortlist ``top_ordinals`` is drawn from.
        top_ordinals: The first ``SelectionInputs.top_n`` of ``eligible_ordinals``.
        per_candidate: One :class:`CandidateSelection` per view, in ``order``.
        required_evidence_missing: Required channel/species pairs and design-input shortfalls this
            run produced no completed evidence for -- a property of the run, not of any candidate.
        summary: The exact key set ``cli.py::_fail_if_nothing_could_qualify`` and
            ``_SELECTION_COUNTERS`` read, plus ``provisional_candidates`` and ``design_input_excluded``.
    """

    order: tuple[int, ...]
    eligible_ordinals: tuple[int, ...]
    top_ordinals: tuple[int, ...]
    per_candidate: tuple[CandidateSelection, ...]
    required_evidence_missing: tuple[str, ...]
    summary: dict[str, Any]


def _tally(values: Iterable[str]) -> dict[str, int]:
    """Count each distinct value, ordered most frequent first, then by name for a stable output.

    Mirrors ``workflow.py::_tally`` exactly; kept as a private copy rather than an import so this
    module stays importable without ``workflow.py``, which is #100's largest file and not this
    slice's to touch.
    """
    counts = Counter(values)
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def evidence_shortfall(view: CandidateView, inputs: SelectionInputs, *, required_only: bool) -> tuple[str, ...]:
    """The evidence this candidate does not have, named so a summary can say which.

    ``required_only=True`` is what disqualifies in ``qualified`` mode. ``False`` is every declared
    pair and every undecided gate whatever its requiredness, and it is what **orders** candidates --
    a separate question, because nothing is required in ``exploratory`` mode, so the required set is
    always empty there and could not put an incomplete candidate below a complete one.

    The unknown rule is scoped to gates *all* of whose channels are required: a subset, and
    deliberately not an intersection. ``max_total_offtarget_hits`` counts the transcriptome and miRNA
    channels together, so an absent miRNA aggregate alone makes it ``UNKNOWN``, and an intersection put
    it in scope -- emptying the shortlist of a complete screen over a channel the policy calls
    exploratory, the mirror of the defect this scoping exists to prevent (#101). A gate reading no
    channel at all (``min_isoform_coverage``, the design-stage gates) is always in scope: its unknown
    is the run's.
    """
    if inputs.requirements is None:
        return ()
    requirements = inputs.requirements
    completed = inputs.completed_pairs
    in_scope_pairs = (
        requirements.required_pairs
        if required_only
        else frozenset(entry.key for entry in requirements.channel_requirements)
    )

    missing: list[str] = []
    for channel, species in sorted(in_scope_pairs):
        if channel == ScreeningChannel.TRANSCRIPTOME.value and species == inputs.query_species:
            if not view.off_target_screened:
                missing.append(f"no_evidence:{channel}:{species}")
            continue
        if completed is None or (channel, species) not in completed:
            missing.append(f"no_evidence:{channel}:{species}")

    # An UNKNOWN verdict is a shortfall for ordering purposes always; for *disqualifying* purposes
    # only when the run said an undecided gate should fail.
    unknown_counts = not required_only or requirements.unknown_evidence_action is UnknownEvidenceAction.FAIL
    if unknown_counts:
        in_scope_channels = {channel for channel, _ in in_scope_pairs}
        for filter_id in sorted(view.unknown_filter_ids):
            channels = inputs.filter_channels.get(filter_id)
            if channels is None or not required_only or {c.value for c in channels} <= in_scope_channels:
                missing.append(f"unknown:{filter_id}")

    return tuple(dict.fromkeys(missing))


def missing_required_evidence(inputs: SelectionInputs) -> tuple[str, ...]:
    """Required channel/species pairs and design-input shortfalls this run has no completed record for.

    A property of the run, not of a candidate, so an empty shortlist can be explained even when every
    candidate was excluded by something else first. Empty in design-only and exploratory modes, where
    nothing is required, because ``inputs.requirements.required_pairs`` is empty there by construction
    upstream. ``design_input:<transcript_id>`` entries are additive (#100/W5): a design-stage batch
    failure named here means a required target set cannot appear complete after it.
    """
    if inputs.requirements is None:
        channel_gaps: tuple[str, ...] = ()
    else:
        completed = inputs.completed_pairs
        channel_gaps = tuple(
            f"{channel}:{species}"
            for channel, species in sorted(inputs.requirements.required_pairs)
            if completed is None or (channel, species) not in completed
        )
    design_input_gaps = tuple(
        f"design_input:{transcript_id}" for transcript_id in sorted(inputs.design_input_shortfalls)
    )
    return tuple(dict.fromkeys((*channel_gaps, *design_input_gaps)))


def select(views: Sequence[CandidateView], inputs: SelectionInputs) -> SelectionResult:
    """Rank, exclude and label a batch of candidates -- the eligibility decision, made once.

    Ports ``workflow.py::_apply_post_screen_ranking``'s loop verbatim onto plain views, with one
    behaviour change: an ``EXPLORATORY`` candidate whose ordering shortfall is non-empty is labelled
    ``PROVISIONAL`` rather than ``ELIGIBLE``, splitting what used to be one undifferentiated state.
    Every other exclusion (repeat, gate, unscored-against-a-mixed-scale) keeps its own
    ``NOT_ELIGIBLE`` outcome, and a ``QUALIFIED`` run's required-evidence shortfall still withholds
    rather than merely reorders.
    """
    qualified = inputs.run_mode is RunMode.QUALIFIED
    exploratory = inputs.run_mode is RunMode.EXPLORATORY

    scored_count = sum(1 for view in views if view.scored_after_screening)
    # Mixed scales: some candidates carry a post-screen score and some do not, so the unscored ones
    # cannot be ranked against them. A wholly unscored batch is not mixed; the evidence rule (below)
    # answers that case instead.
    mixed_scales = 0 < scored_count < len(views)

    any_shortfall = {view.ordinal: evidence_shortfall(view, inputs, required_only=False) for view in views}
    required_shortfall = {view.ordinal: evidence_shortfall(view, inputs, required_only=True) for view in views}

    ordered = sorted(
        views,
        key=lambda view: (not any_shortfall[view.ordinal], view.scored_after_screening, view.ranking_score),
        reverse=True,
    )

    per_candidate: list[CandidateSelection] = []
    eligible_ordinals: list[int] = []
    shortfalls: list[tuple[int, tuple[str, ...]]] = []
    repeat_excluded = 0
    filter_excluded = 0
    unscored_excluded = 0

    for view in ordered:
        shortfall = required_shortfall[view.ordinal]
        ordering_shortfall = any_shortfall[view.ordinal]

        if view.repeat_flagged and inputs.repeat_rejects:
            repeat_excluded += 1
            state = SelectionState.NOT_ELIGIBLE
        elif not view.passes_filters:
            filter_excluded += 1
            state = SelectionState.NOT_ELIGIBLE
        elif shortfall and qualified:
            shortfalls.append((view.ordinal, shortfall))
            state = SelectionState.WITHHELD
        elif mixed_scales and not view.scored_after_screening:
            unscored_excluded += 1
            state = SelectionState.NOT_ELIGIBLE
        else:
            state = SelectionState.PROVISIONAL if (exploratory and ordering_shortfall) else SelectionState.ELIGIBLE
            eligible_ordinals.append(view.ordinal)

        per_candidate.append(CandidateSelection(view.candidate_id, view.ordinal, state, shortfall, ordering_shortfall))

    top_ordinals = tuple(eligible_ordinals[: inputs.top_n])
    required_evidence_missing = missing_required_evidence(inputs)
    eligible_count = sum(1 for selection in per_candidate if selection.state is SelectionState.ELIGIBLE)
    provisional_count = sum(1 for selection in per_candidate if selection.state is SelectionState.PROVISIONAL)

    summary: dict[str, Any] = {
        "run_mode": inputs.run_mode.value if inputs.run_mode is not None else None,
        "eligible_candidates": eligible_count,
        "top_candidates": len(top_ordinals),
        "scored_after_screening": scored_count,
        "repeat_excluded": repeat_excluded,
        "filter_excluded": filter_excluded,
        "evidence_excluded": len(shortfalls),
        "unscored_excluded": unscored_excluded,
        # Every distinct reason with its cost, not one example: a missing species alignment and a
        # missing coverage annotation are different problems with different fixes.
        "evidence_shortfall_reasons": _tally(reason for _, reasons in shortfalls for reason in reasons),
        # Run-level, and not derivable from the per-candidate reasons above: the loop moves to the next
        # exclusion test before it reaches the evidence check, so a run whose screen produced nothing
        # AND whose candidates also fail a gate reports no shortfall reasons at all.
        "required_evidence_missing": required_evidence_missing,
        "provisional_candidates": provisional_count,
        "design_input_excluded": len(inputs.design_input_shortfalls),
    }

    return SelectionResult(
        order=tuple(view.ordinal for view in ordered),
        eligible_ordinals=tuple(eligible_ordinals),
        top_ordinals=top_ordinals,
        per_candidate=tuple(per_candidate),
        required_evidence_missing=required_evidence_missing,
        summary=summary,
    )
