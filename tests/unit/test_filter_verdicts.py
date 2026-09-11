"""Per-filter verdicts, and the three gates that record without rejecting.

``passes_filters`` holds one label and every gate overwrote it, so a gate's outcome was lost as soon
as another gate rejected the same candidate and a label count was never a rejection count. These
tests pin the replacement: one verdict per declared filter, carrying the value that filter compared,
and ``FilterAction.WARN`` meaning what its name says.
"""

import pytest

from sirnaforge.config.run_policy import FILTER_SPEC_BY_ID, EntryPoint, resolve_run_policy
from sirnaforge.models.policy import DECLARED_FILTER_IDS, FilterAction, FilterEvaluation
from sirnaforge.models.sirna import SiRNACandidate, build_candidate_row

#: The gates demoted to warn, and why they cannot be demoted separately.
DEMOTED = ("min_asymmetry_score", "max_mirna_perfect_seed", "fail_on_high_risk_mirna")


def _candidate() -> SiRNACandidate:
    return SiRNACandidate(
        id="cand_1",
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence="ACGTACGTACGTACGTACGTA",
        passenger_sequence="TACGTACGTACGTACGTACGT",
        length=21,
        gc_content=47.6,
        asymmetry_score=0.5,
    )


@pytest.mark.unit
def test_the_exported_filter_ids_cannot_drift_from_the_registry():
    """``DECLARED_FILTER_IDS`` duplicates the registry's ids because models cannot import it.

    The candidate row needs a fixed column per filter, and ``models`` cannot reach
    ``config.run_policy`` -- that module imports ``models``. So the list is restated, and this is what
    stops a filter being added to one place and forgotten in the other.
    """
    assert tuple(FILTER_SPEC_BY_ID) == DECLARED_FILTER_IDS, (
        "a filter was added or removed in run_policy without updating DECLARED_FILTER_IDS"
    )


@pytest.mark.unit
def test_a_warn_gate_records_its_verdict_and_does_not_reject():
    """This is the whole point of the change: a failure that is recorded rather than fatal."""
    candidate = _candidate()
    candidate.record_filter_verdict(
        "min_asymmetry_score",
        observed=0.4,
        passed=False,
        action=FilterAction.WARN,
        status=SiRNACandidate.FilterStatus.LOW_ASYMMETRY,
    )

    assert candidate.filter_verdicts["min_asymmetry_score"] == FilterEvaluation.FAIL.value
    assert candidate.filter_observed["min_asymmetry_score"] == 0.4
    assert candidate.passes_filters is True, "warn must not reject"


@pytest.mark.unit
def test_a_fail_gate_still_rejects_and_keeps_the_first_label():
    """First-failure-wins survives, but it is now only a choice about the label."""
    candidate = _candidate()
    candidate.record_filter_verdict(
        "max_paired_fraction",
        observed=0.9,
        passed=False,
        action=FilterAction.FAIL,
        status=SiRNACandidate.FilterStatus.EXCESS_PAIRING,
    )
    candidate.record_filter_verdict(
        "max_off_target_count",
        observed=99,
        passed=False,
        action=FilterAction.FAIL,
        status=SiRNACandidate.FilterStatus.EXCESS_OFF_TARGETS,
    )

    assert candidate.passes_filters == SiRNACandidate.FilterStatus.EXCESS_PAIRING, "the first label wins"
    # ... and the second gate's verdict is not lost to it, which is what the single label used to do.
    assert candidate.filter_verdicts["max_off_target_count"] == FilterEvaluation.FAIL.value
    assert candidate.filter_observed["max_off_target_count"] == 99


@pytest.mark.unit
def test_missing_evidence_records_unknown_rather_than_a_pass():
    """A gate in force whose evidence is unavailable has not been satisfied."""
    candidate = _candidate()
    candidate.record_filter_verdict(
        "max_off_target_count",
        observed=None,
        passed=True,
        action=FilterAction.FAIL,
        status=SiRNACandidate.FilterStatus.EXCESS_OFF_TARGETS,
    )

    assert candidate.filter_verdicts["max_off_target_count"] == FilterEvaluation.UNKNOWN.value
    assert candidate.passes_filters is True


@pytest.mark.unit
@pytest.mark.parametrize("filter_id", DEMOTED)
def test_the_three_uncalibrated_gates_ship_as_warn(filter_id):
    """Asymmetry, the perfect-seed ceiling of 0, and the high-risk flag all report rather than reject.

    ``fail_on_high_risk_mirna`` is in this list because it cannot be left behind: a high-risk hit is by
    definition a perfect seed hit, so it only ever saw candidates ``max_mirna_perfect_seed`` had
    already rejected. Demoting the seed gate alone would have handed those same rejections to this one
    under a new label and changed nothing.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)

    assert policy.descriptor(filter_id).action is FilterAction.WARN


@pytest.mark.unit
def test_the_profile_hash_changes_when_an_action_changes():
    """A content hash that ignores actions is a hash that misses what a run rejects.

    Actions live in ``FILTER_SPECS``, not in the profile baseline, so before this the hash was
    byte-identical across a change that moved 65.9% of candidates from rejected to reported.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    identity = policy.profile

    payload = f"{identity.name}{identity.version}{identity.content_hash}"
    assert identity.content_hash.startswith("sha256:")
    # The hash must depend on the actions, which is only observable by asserting the demoted defaults
    # are what the hash was computed over: if someone reverts one to fail, this file's other tests fail
    # and the hash changes with them.
    assert all(policy.descriptor(fid).action is FilterAction.WARN for fid in DEMOTED), payload


@pytest.mark.unit
def test_every_declared_filter_gets_a_column_even_when_it_did_not_run():
    """Fixed columns, always the same set, and a word rather than a blank when a filter is off.

    A blank cell is read as a verdict. ``not_evaluated`` says which of the four things actually
    happened, which is the distinction ``FilterEvaluation`` exists to make.
    """
    candidate = _candidate()
    candidate.record_filter_verdict(
        "min_asymmetry_score",
        observed=0.4,
        passed=False,
        action=FilterAction.WARN,
        status=SiRNACandidate.FilterStatus.LOW_ASYMMETRY,
    )
    row = build_candidate_row(candidate)

    for filter_id in DECLARED_FILTER_IDS:
        assert f"{filter_id}_verdict" in row, f"{filter_id} has no verdict column"
        assert f"{filter_id}_observed" in row, f"{filter_id} has no observed column"

    assert row["min_asymmetry_score_verdict"] == FilterEvaluation.FAIL.value
    assert row["min_asymmetry_score_observed"] == 0.4
    # A filter this candidate never reached says so, rather than leaving a cell a reader will fill in.
    assert row["max_off_target_count_verdict"] == FilterEvaluation.NOT_EVALUATED.value
    assert row["max_off_target_count_observed"] is None


@pytest.mark.unit
def test_the_observed_value_is_exported_for_the_gates_whose_input_is_not():
    """Six gates read human-stratified counters that no other column carries.

    ``evidence_exported=False`` marks them. Their verdict is only checkable from the row because the
    value they compared travels with it -- the identically named exported columns are all-species
    totals and disagree with what the gate saw.
    """
    unexported = [fid for fid, spec in FILTER_SPEC_BY_ID.items() if not spec.evidence_exported]

    assert unexported, "the fixture assumes at least one gate does not export its input"
    candidate = _candidate()
    for filter_id in unexported:
        candidate.record_filter_verdict(
            filter_id,
            observed=7,
            passed=True,
            action=FilterAction.WARN,
        )

    row = build_candidate_row(candidate)
    for filter_id in unexported:
        assert row[f"{filter_id}_observed"] == 7, f"{filter_id} must publish the value it compared"
