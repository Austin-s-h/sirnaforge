"""A screen that published nothing reaches every off-target gate, and none of them passes.

The residual half of #106. Its first half -- a *completed* screen that found no hits -- was fixed by
gating the no-hit branch on a measured zero (pinned in ``test_gate_verdict_completeness.py``). This
file covers the mirror case: the basic sequence-only fallback, ``nextflow_unavailable``,
``nextflow_failed`` and an output directory that never appeared all return without reaching
``_integrate_offtarget_results``, so no candidate reached ``_gate_offtarget_counts`` at all and every
channel-reading gate exported ``not_evaluated`` -- the same cell a run with no threshold configured
writes. A screen that never happened and a screen that came back clean must not export the same nine
cells.

Separate from ``test_gate_verdict_completeness.py`` (which owns #105's residual) so the two residuals
fail independently; the ``NO_SCREEN_PATHS`` driver is copied rather than imported for the same reason.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

import asyncio
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.models.policy import FilterAction, FilterEvaluation, RunMode, RunStatus, ScreeningChannel
from sirnaforge.models.sirna import (
    DesignResult,
    OffTargetFilterCriteria,
    SelectionState,
    SiRNACandidate,
    build_candidate_row,
)
from sirnaforge.reporting.payload import _evaluate as evaluate_descriptor
from sirnaforge.reporting.payload import observed_column
from sirnaforge.workflow import POST_SCREEN_FILTER_CHANNELS, SiRNAWorkflow, WorkflowConfig

GUIDE = "ACGTACGTACGTACGTACGTA"
NO_HITS: dict[str, object] = {"status": "completed", "results": {}}

#: ``sirnaforge.reporting.payload``'s own verdict code for "the run could not decide this gate".
REPORT_UNKNOWN = 2

#: The four ways out of screening on which nothing at all is published, with the exit each reports.
#: Named rather than stubbed at ``step5`` level, because the defect is in *which* paths reach the
#: gates: a route that skipped them would keep passing against a stub of itself.
NO_SCREEN_PATHS: dict[str, dict[str, str]] = {
    "nextflow_failed": {"status": "skipped", "reason": "nextflow_failed"},
    "nextflow_unavailable": {"status": "skipped", "reason": "nextflow_unavailable"},
    "basic_fallback": {"status": "partial", "method": "basic"},
    "missing_output": {"status": "partial", "method": "embedded_nextflow"},
}


def _candidate(candidate_id: str = "probe") -> SiRNACandidate:
    """A minimal passing candidate carrying enough component scores to be re-scored after screening."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE.translate(str.maketrans("ATCG", "TAGC"))[::-1],
        length=len(GUIDE),
        gc_content=47.6,
        asymmetry_score=0.5,
        design_score=50.0,
        component_scores={"target_accessibility": 0.5, "asymmetry": 0.5, "gc_content": 0.5},
    )


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    stated: dict[str, object] | None = None,
    run_mode: RunMode | None = None,
    screen_species: list[str] | None = None,
) -> SiRNAWorkflow:
    """A workflow whose policy is resolved through the one resolver, as every real entry point does."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        run_mode=run_mode,
        stated=stated,
        query_species="human",
        screen_species=screen_species or ["human"],
    )
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        resolved_policy=policy,
        screen_species=screen_species or ["human"],
        query_species="human",
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    return workflow


def _design_result(workflow: SiRNAWorkflow, candidates: list[SiRNACandidate]) -> DesignResult:
    return DesignResult(
        input_file="<test>",
        parameters=workflow.config.design_params,
        candidates=list(candidates),
        top_candidates=list(candidates),
        total_sequences=1,
        total_candidates=len(candidates),
        filtered_candidates=len(candidates),
        processing_time=0.0,
    )


def _run_step5_without_a_screen(
    workflow: SiRNAWorkflow, candidates: list[SiRNACandidate], path: str
) -> tuple[DesignResult, dict[str, object]]:
    """Drive ``step5_offtarget_analysis`` to one of ``NO_SCREEN_PATHS`` and hand back both answers.

    Reference resolution and repeat detection are stubbed on every path: both need a real cDNA
    reference, and neither is what is under test. Everything downstream of them is the product code.
    """
    design_result = _design_result(workflow, candidates)

    async def _references(*_args: object, **_kwargs: object) -> bool:
        return path != "basic_fallback"

    workflow._resolve_screening_references = _references  # type: ignore[method-assign]
    workflow._run_repeat_detection = lambda *_args, **_kwargs: {"status": "skipped"}  # type: ignore[method-assign]

    if path == "basic_fallback":
        # No index and no active species: the sequence-only fallback, which aligns nothing at all.
        workflow._resolve_active_screen_species = lambda *_args, **_kwargs: []  # type: ignore[method-assign]
    elif path == "nextflow_unavailable":
        workflow._validate_nextflow_environment = lambda _runner: False  # type: ignore[method-assign]
    elif path == "nextflow_failed":

        async def _explode(*_args: object, **_kwargs: object) -> dict[str, object]:
            raise RuntimeError("nextflow submit failed")

        workflow._run_nextflow_offtarget_analysis = _explode  # type: ignore[method-assign]
    elif path == "missing_output":

        async def _no_output(*_args: object, **_kwargs: object) -> dict[str, object]:
            # The runner reported success and wrote nothing: parsing returns "missing", so
            # _integrate_offtarget_results takes its no-completed-data exit.
            return await workflow._process_nextflow_results(
                candidates, workflow.config.output_dir / "absent", {"status": "completed"}
            )

        workflow._run_nextflow_offtarget_analysis = _no_output  # type: ignore[method-assign]
    else:  # pragma: no cover - a typo in the parametrisation, not a product path
        raise AssertionError(f"unknown path {path}")

    outcome = asyncio.run(workflow.step5_offtarget_analysis(design_result))
    for key, value in NO_SCREEN_PATHS[path].items():
        assert outcome[key] == value, f"{path} no longer exits the way this test drives it"
    return design_result, outcome


def _in_force_gates(workflow: SiRNAWorkflow) -> tuple[str, ...]:
    """The post-screen gates this run actually applies: a threshold, and an action that is not off.

    Derived from the resolved policy rather than listed, because which gates are in force is a
    registry decision (``max_mirna_1mm_seed`` and ``max_transcriptome_seed_perfect`` ship off in
    0.7.1) and a hard-coded list would assert against the registry instead of against the fix.
    """
    policy = workflow.config.resolved_policy
    in_force = []
    for filter_id in POST_SCREEN_FILTER_CHANNELS:
        descriptor = policy.descriptor(filter_id)
        if descriptor.action is not FilterAction.OFF and descriptor.threshold is not None:
            in_force.append(filter_id)
    return tuple(in_force)


@pytest.mark.unit
@pytest.mark.parametrize("path", sorted(NO_SCREEN_PATHS))
def test_every_in_force_gate_is_unknown_when_the_screen_published_nothing(tmp_path: Path, path: str) -> None:
    """#106 residual: a gate the run could not apply reports unknown, not ``not_evaluated``.

    None of these four exits reaches ``_integrate_offtarget_results``, so before the fix no candidate
    reached ``_gate_offtarget_counts`` and the row carried ``not_evaluated`` for every off-target gate
    -- indistinguishable from a run that configured no ceiling, and read by #103's report as a gate
    nobody asked for rather than as one the run failed to answer.

    The observed value stays empty on purpose: writing the unmeasured zero there would let the report
    re-derive a confident PASS the run never made. The exported cell is asserted alongside the
    in-memory verdict, because the cell is what a consumer actually reads.
    """
    workflow = _workflow(tmp_path, f"no_screen_{path}")
    assert workflow.config.resolved_policy.run_mode is RunMode.QUALIFIED, "the default mode claims completeness"
    candidate = _candidate()

    _run_step5_without_a_screen(workflow, [candidate], path)

    gates = _in_force_gates(workflow)
    assert len(gates) >= 5, f"the default policy must apply more than a token gate here: {gates}"
    row = build_candidate_row(candidate)
    populated = {key for key, value in row.items() if value is not None}
    for filter_id in gates:
        # The exported cell first: an absent verdict exports ``not_evaluated``, which is exactly the
        # reading this fix removes, so failing on the cell names the defect rather than a KeyError.
        assert row[f"{filter_id}_verdict"] == FilterEvaluation.UNKNOWN.value, filter_id
        assert candidate.filter_verdicts[filter_id] == FilterEvaluation.UNKNOWN.value, filter_id
        assert candidate.filter_observed[filter_id] is None, filter_id
        descriptor = workflow.config.resolved_policy.descriptor(filter_id)
        _, verdict_code, _ = evaluate_descriptor(descriptor, pd.Series(row), observed_column(descriptor, populated))
        assert verdict_code == REPORT_UNKNOWN, filter_id

    # An undecidable gate withholds; it never rejects. The two are different answers with different
    # fixes, and only the second belongs on ``passes_filters``.
    assert candidate.passes_filters is True
    summary = workflow._selection_summary
    assert summary["filter_excluded"] == 0
    assert summary["evidence_excluded"] == 1
    reasons = summary["evidence_shortfall_reasons"]
    # Disqualifying reasons are scoped to gates all of whose channels the run *requires* (#101), and
    # 0.7.1 calls every miRNA pair exploratory -- so the transcriptome gates are named and the miRNA
    # ones deliberately are not, even though both are UNKNOWN on the row above.
    for filter_id in gates:
        required_channels = POST_SCREEN_FILTER_CHANNELS[filter_id] == frozenset({ScreeningChannel.TRANSCRIPTOME})
        assert (reasons.get(f"unknown:{filter_id}") == 1) is required_channels, filter_id
    assert reasons["no_evidence:transcriptome:human"] == 1


@pytest.mark.unit
@pytest.mark.parametrize("path", sorted(NO_SCREEN_PATHS))
def test_a_gate_that_is_off_still_makes_no_claim_when_the_screen_published_nothing(tmp_path: Path, path: str) -> None:
    """The guard on the fix above: the gates are completed, not switched on.

    ``max_mirna_1mm_seed`` is read by no 0.7.1 gate and ``max_transcriptome_seed_perfect`` /
    ``max_total_offtarget_hits`` ship without a threshold, so all three make no claim at all. An
    ``UNKNOWN`` there would invent a shortfall out of a setting the user never made -- and in the
    default QUALIFIED mode that shortfall is what empties a shortlist.
    """
    workflow = _workflow(tmp_path, f"off_gates_{path}")
    candidate = _candidate()

    _run_step5_without_a_screen(workflow, [candidate], path)

    silent = tuple(filter_id for filter_id in POST_SCREEN_FILTER_CHANNELS if filter_id not in _in_force_gates(workflow))
    assert "max_mirna_1mm_seed" in silent, "the registry still ships this gate off; this test reads it"
    for filter_id in silent:
        assert candidate.filter_verdicts.get(filter_id, FilterEvaluation.NOT_EVALUATED.value) in {
            FilterEvaluation.NOT_EVALUATED.value,
            None,
        }, filter_id
        assert f"unknown:{filter_id}" not in workflow._selection_summary["evidence_shortfall_reasons"], filter_id


@pytest.mark.unit
def test_a_gate_scoped_to_every_screened_species_cannot_pass_on_a_zero_nobody_measured(tmp_path: Path) -> None:
    """``max_off_target_count`` declares no species scope, and that is what used to let it pass.

    A gate's evidence requirement is the channel x species pairs its own ``FilterScope`` names; an
    unrestricted scope means every species the run screened, which ``_integrate_offtarget_results``
    records and which a path that never integrated does not. With that set left empty the gate has no
    pair to be incomplete about, so a count of zero against a ceiling of 15 read as a decided PASS --
    on the default configuration, where this gate and ``max_transcriptome_seed_perfect`` are the two
    that declare no scope. The scope is therefore established from what was *requested* here too.
    """
    workflow = _workflow(
        tmp_path,
        "unscoped_gate",
        stated={"max_transcriptome_seed_perfect": 0},
        screen_species=["human", "mouse"],
    )
    assert workflow.config.resolved_policy.descriptor("max_off_target_count").scope.species == frozenset()
    candidate = _candidate()

    _run_step5_without_a_screen(workflow, [candidate], "nextflow_unavailable")

    assert candidate.filter_verdicts["max_off_target_count"] == FilterEvaluation.UNKNOWN.value
    # Opted in by the caller, so it is in force despite shipping off, and it declares no scope either.
    assert candidate.filter_verdicts["max_transcriptome_seed_perfect"] == FilterEvaluation.UNKNOWN.value
    assert workflow._screened_species_scope == frozenset({"human", "mouse"})


@pytest.mark.unit
def test_a_completed_clean_screen_still_records_a_measured_pass(tmp_path: Path) -> None:
    """The other side of the same distinction, so the fix cannot be made by blanket-unknowning.

    A screen that ran and found nothing is a *measurement* of zero and must keep its decided PASS with
    the observed value intact. This is the assertion that makes the test above mean something: if
    every path reported ``UNKNOWN`` the nine cells would agree again, just at the other extreme.
    """
    workflow = _workflow(tmp_path, "clean_screen")
    candidate = _candidate()

    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=["human"], mirna_screened=True)

    for filter_id in _in_force_gates(workflow):
        assert candidate.filter_verdicts[filter_id] == FilterEvaluation.PASS.value, filter_id
        assert candidate.filter_observed[filter_id] == 0, filter_id


@pytest.mark.unit
def test_a_verdict_decided_before_the_failure_is_not_overwritten(tmp_path: Path) -> None:
    """The call that held the counts owns the answer; the fallback only fills an empty cell.

    ``record_filter_verdict`` overwrites, and ``nextflow_failed`` is reachable from an exception raised
    *after* integration already gated some candidates -- a shortfall warning while persisting hit
    tables, say. Re-gating them would replace a real measurement with ``UNKNOWN`` and turn a completed
    screen into an unanswered one on the way out.
    """
    workflow = _workflow(tmp_path, "decided_survives")
    measured, untouched = _candidate("measured"), _candidate("untouched")

    async def _integrate_then_explode(*_args: object, **_kwargs: object) -> dict[str, object]:
        workflow._integrate_offtarget_results(
            [measured], dict(NO_HITS), OffTargetFilterCriteria(), screened_species=["human"], mirna_screened=True
        )
        raise RuntimeError("aggregation failed after the alignment was integrated")

    async def _references(*_args: object, **_kwargs: object) -> bool:
        return True

    workflow._resolve_screening_references = _references  # type: ignore[method-assign]
    workflow._run_repeat_detection = lambda *_args, **_kwargs: {"status": "skipped"}  # type: ignore[method-assign]
    workflow._run_nextflow_offtarget_analysis = _integrate_then_explode  # type: ignore[method-assign]

    outcome = asyncio.run(workflow.step5_offtarget_analysis(_design_result(workflow, [measured, untouched])))
    assert outcome["reason"] == "nextflow_failed"

    assert measured.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.PASS.value
    assert measured.filter_observed["max_transcriptome_hits_0mm"] == 0
    assert untouched.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.UNKNOWN.value


@pytest.mark.unit
def test_a_screen_the_user_switched_off_makes_no_claim(tmp_path: Path) -> None:
    """``--skip-off-targets`` is a channel that was never requested, and ``not_evaluated`` says so.

    The deliberate non-call. ``NOT_REQUESTED`` and ``FAILED`` are different plan statuses for a reason:
    a user who asked for no screen is not owed an undecidable verdict on every gate that would have
    read one, and reporting ``UNKNOWN`` there would make the exports of a design-only run and of a
    broken screen identical -- the exact collapse this slice exists to undo, in the other direction.
    """
    workflow = _workflow(tmp_path, "user_disabled", stated={"check_off_targets": False})
    assert workflow.config.design_params.check_off_targets is False
    candidate = _candidate()

    outcome = asyncio.run(workflow.step5_offtarget_analysis(_design_result(workflow, [candidate])))

    # ``run_status`` is W4's addition beside these two, and it agrees with the non-call below: a
    # channel nobody requested is NOT_REQUESTED, not an incomplete screen (#100).
    assert outcome == {
        "status": "skipped",
        "reason": "user_disabled",
        "run_status": RunStatus.NOT_REQUESTED.value,
    }
    for filter_id in POST_SCREEN_FILTER_CHANNELS:
        assert filter_id not in candidate.filter_verdicts, filter_id


@pytest.mark.unit
@pytest.mark.parametrize("path", sorted(NO_SCREEN_PATHS))
def test_an_exploratory_run_keeps_its_shortlist_and_labels_it_provisional(tmp_path: Path, path: str) -> None:
    """Nine new ``UNKNOWN`` verdicts must not empty an exploratory shortlist, only label it.

    ``EXPLORATORY`` requires nothing, so an undecided gate reorders and labels rather than withholds.
    Without this the fix would silently convert every unscreened exploratory run into an empty
    deliverable, which is the failure mode #100's own selection rules were written to avoid.
    """
    workflow = _workflow(tmp_path, f"exploratory_{path}", run_mode=RunMode.EXPLORATORY)
    candidate = _candidate()

    design_result, _ = _run_step5_without_a_screen(workflow, [candidate], path)

    assert candidate.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.UNKNOWN.value
    assert [c.id for c in design_result.top_candidates] == ["probe"]
    assert candidate.selection_state == SelectionState.PROVISIONAL.value
    assert workflow._selection_summary["evidence_excluded"] == 0


@pytest.mark.unit
@pytest.mark.parametrize("path", sorted(NO_SCREEN_PATHS))
def test_the_no_screen_exits_keep_their_status_and_reason(tmp_path: Path, path: str) -> None:
    """Run status is W4's, not this slice's: gating these paths must not restate their exit.

    ``_run_status_taxonomy`` is where ``{status, reason}`` becomes an enum and an exit code. Until
    then these four strings are what the CLI and ``workflow_summary.json`` read, and a gate verdict is
    not a reason to change them -- so they are pinned here rather than left to drift.
    """
    workflow = _workflow(tmp_path, f"exit_shape_{path}")

    _, outcome = _run_step5_without_a_screen(workflow, [_candidate()], path)

    for key, value in NO_SCREEN_PATHS[path].items():
        assert outcome[key] == value
