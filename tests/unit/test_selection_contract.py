"""core/selection.py is the pure eligibility/shortlist decision #100 factors out of ``workflow.py``.

``workflow.py::_apply_post_screen_ranking`` decided ranking, exclusion and the ``selection_state``
label in one place with no test surface of its own. This file pins that decision against plain
views so a future caller (W2, and the offtarget-only entry point) can be routed through it with no
behaviour change except the one #100 adds deliberately: an ``EXPLORATORY`` candidate with an
evidence shortfall is now ``PROVISIONAL``, not indistinguishable from a fully evidenced ``ELIGIBLE``
one.
"""

from sirnaforge.core.selection import (
    CandidateView,
    SelectionInputs,
    evidence_shortfall,
    missing_required_evidence,
    select,
)
from sirnaforge.models.policy import (
    ChannelRequirement,
    EvidenceRequirements,
    RunMode,
    ScreeningChannel,
    UnknownEvidenceAction,
)
from sirnaforge.models.sirna import SelectionState


def _make_view(
    candidate_id: str, ordinal: int, *, off_target_screened: bool, passes_filters: bool = True
) -> CandidateView:
    """One evidence-incomplete, otherwise-clean view: passes its gates, was never screened."""
    return CandidateView(
        candidate_id=candidate_id,
        ordinal=ordinal,
        passes_filters=passes_filters,
        repeat_flagged=False,
        scored_after_screening=False,
        off_target_screened=off_target_screened,
        ranking_score=1.0,
    )


def _human_transcriptome_requirements(unknown_evidence_action: UnknownEvidenceAction) -> EvidenceRequirements:
    return EvidenceRequirements(
        channel_requirements=(
            ChannelRequirement(channel=ScreeningChannel.TRANSCRIPTOME, species="human", requiredness="required"),
        ),
        unknown_evidence_action=unknown_evidence_action,
    )


class TestQualifiedWithholds:
    """(a) QUALIFIED + a query-species transcriptome shortfall withholds every candidate."""

    def test_withholds_every_candidate_for_incomplete_required_evidence(self) -> None:
        """No completed transcriptome evidence and QUALIFIED mode withholds the whole batch."""
        views = [_make_view(f"c{i}", i, off_target_screened=False) for i in range(3)]
        inputs = SelectionInputs(
            run_mode=RunMode.QUALIFIED,
            requirements=_human_transcriptome_requirements(UnknownEvidenceAction.WARN),
            completed_pairs=frozenset(),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        result = select(views, inputs)

        assert result.summary["eligible_candidates"] == 0
        assert result.summary["evidence_excluded"] == len(views)
        assert all(candidate.state is SelectionState.WITHHELD for candidate in result.per_candidate)


class TestExploratoryIsProvisional:
    """(b) EXPLORATORY with the same shortfall labels PROVISIONAL and ranks it below evidenced rows."""

    def test_provisional_state_and_summary_count(self) -> None:
        """Same shortfall, EXPLORATORY mode: every candidate is PROVISIONAL, none ELIGIBLE."""
        views = [_make_view(f"c{i}", i, off_target_screened=False) for i in range(3)]
        inputs = SelectionInputs(
            run_mode=RunMode.EXPLORATORY,
            requirements=_human_transcriptome_requirements(UnknownEvidenceAction.WARN),
            completed_pairs=frozenset(),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        result = select(views, inputs)

        assert all(candidate.state is SelectionState.PROVISIONAL for candidate in result.per_candidate)
        assert all(candidate.state is not SelectionState.ELIGIBLE for candidate in result.per_candidate)
        assert result.summary["provisional_candidates"] == len(views)

    def test_provisional_sorts_strictly_after_fully_evidenced_rows(self) -> None:
        """A mixed batch: one candidate is fully evidenced, the rest are short. Order must not mix them."""
        views = [
            _make_view("evidenced", 0, off_target_screened=True),
            _make_view("short-a", 1, off_target_screened=False),
            _make_view("short-b", 2, off_target_screened=False),
        ]
        inputs = SelectionInputs(
            run_mode=RunMode.EXPLORATORY,
            requirements=_human_transcriptome_requirements(UnknownEvidenceAction.WARN),
            completed_pairs=frozenset(),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        result = select(views, inputs)

        states = {candidate.ordinal: candidate.state for candidate in result.per_candidate}
        assert states[0] is SelectionState.ELIGIBLE
        assert states[1] is SelectionState.PROVISIONAL
        assert states[2] is SelectionState.PROVISIONAL

        eligible_position = result.order.index(0)
        provisional_positions = [result.order.index(1), result.order.index(2)]
        assert all(eligible_position < position for position in provisional_positions)


class TestDesignOnlyIsAlwaysEligible:
    """(c) DESIGN_ONLY: nothing was screened, so nothing is missing; every row is ELIGIBLE."""

    def test_every_state_is_eligible(self) -> None:
        """DESIGN_ONLY claims nothing, so an evidence shortfall never withholds or provisions."""
        views = [_make_view(f"c{i}", i, off_target_screened=False) for i in range(3)]
        inputs = SelectionInputs(
            run_mode=RunMode.DESIGN_ONLY,
            requirements=None,
            completed_pairs=None,
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        result = select(views, inputs)

        assert all(candidate.state is SelectionState.ELIGIBLE for candidate in result.per_candidate)
        assert result.summary["eligible_candidates"] == len(views)


class TestSummaryKeySuperset:
    """(d) summary.keys() is a superset of the ten keys cli.py reads, plus the two new ones."""

    #: cli.py::_fail_if_nothing_could_qualify and _SELECTION_COUNTERS's own ten keys, spelled out here
    #: rather than imported, so this test does not depend on cli.py's file ownership.
    _CLI_READ_KEYS = (
        "run_mode",
        "eligible_candidates",
        "top_candidates",
        "scored_after_screening",
        "repeat_excluded",
        "filter_excluded",
        "evidence_excluded",
        "unscored_excluded",
        "evidence_shortfall_reasons",
        "required_evidence_missing",
    )

    def test_summary_keys_superset(self) -> None:
        """The summary keeps every key the CLI already reads and adds the two new ones."""
        views = [_make_view("c0", 0, off_target_screened=True)]
        inputs = SelectionInputs(
            run_mode=RunMode.DESIGN_ONLY,
            requirements=None,
            completed_pairs=None,
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        result = select(views, inputs)

        assert set(self._CLI_READ_KEYS) <= result.summary.keys()
        assert "provisional_candidates" in result.summary
        assert "design_input_excluded" in result.summary


class TestEvidenceShortfallHelpers:
    """Unit coverage for the two pure helpers ``select`` builds on top of."""

    def test_evidence_shortfall_required_only_true_scopes_to_required_pairs(self) -> None:
        """required_only=True reports only the declared required pair, named by channel and species."""
        view = _make_view("c0", 0, off_target_screened=False)
        inputs = SelectionInputs(
            run_mode=RunMode.QUALIFIED,
            requirements=_human_transcriptome_requirements(UnknownEvidenceAction.WARN),
            completed_pairs=frozenset(),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        shortfall = evidence_shortfall(view, inputs, required_only=True)

        assert shortfall == ("no_evidence:transcriptome:human",)

    def test_evidence_shortfall_none_requirements_is_empty(self) -> None:
        """No policy at all cannot compute a shortfall, in either scoring mode."""
        view = _make_view("c0", 0, off_target_screened=False)
        inputs = SelectionInputs(
            run_mode=RunMode.DESIGN_ONLY,
            requirements=None,
            completed_pairs=None,
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        assert evidence_shortfall(view, inputs, required_only=True) == ()
        assert evidence_shortfall(view, inputs, required_only=False) == ()

    def test_unknown_verdict_disqualifies_only_when_action_is_fail(self) -> None:
        """The UNKNOWN-verdict rule is gated on UnknownEvidenceAction.FAIL, not on requiredness alone."""
        view = CandidateView(
            candidate_id="c0",
            ordinal=0,
            passes_filters=True,
            repeat_flagged=False,
            scored_after_screening=False,
            off_target_screened=True,
            ranking_score=1.0,
            unknown_filter_ids=frozenset({"min_isoform_coverage"}),
        )
        requirements_warn = _human_transcriptome_requirements(UnknownEvidenceAction.WARN)
        requirements_fail = _human_transcriptome_requirements(UnknownEvidenceAction.FAIL)
        inputs_warn = SelectionInputs(
            run_mode=RunMode.QUALIFIED,
            requirements=requirements_warn,
            completed_pairs=frozenset({("transcriptome", "human")}),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )
        inputs_fail = SelectionInputs(
            run_mode=RunMode.QUALIFIED,
            requirements=requirements_fail,
            completed_pairs=frozenset({("transcriptome", "human")}),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
        )

        assert evidence_shortfall(view, inputs_warn, required_only=True) == ()
        assert evidence_shortfall(view, inputs_fail, required_only=True) == ("unknown:min_isoform_coverage",)

    def test_missing_required_evidence_includes_design_input_shortfalls(self) -> None:
        """A dropped transcript batch is named in required_evidence_missing even with evidence complete."""
        inputs = SelectionInputs(
            run_mode=RunMode.QUALIFIED,
            requirements=_human_transcriptome_requirements(UnknownEvidenceAction.WARN),
            completed_pairs=frozenset({("transcriptome", "human")}),
            query_species="human",
            filter_channels={},
            repeat_rejects=False,
            top_n=10,
            design_input_shortfalls={"ENST00000000001": "worker crashed"},
        )

        missing = missing_required_evidence(inputs)

        assert missing == ("design_input:ENST00000000001",)
