"""A design failure that drops a transcript is written down, so the target set cannot read as covered.

#100's design-stage half. ``step3_design_sirnas``' parallel branch catches a failed transcript batch
and carries on -- the right call for a 40-transcript gene, because one bad sequence must not lose the
other 39 -- but the loss was recorded nowhere. The batch's transcripts vanished from the candidate
pool while ``DesignResult.total_sequences`` kept counting them, so the two numbers diverged silently
and no artifact said which ids were lost or why. A ``qualified`` run could then publish a shortlist
over an incomplete target set.

The non-parallel branch has the inverse failure mode and needs no record: it has no ``try``/``except``
at all, so a design failure propagates and no result is published to be wrong. That asymmetry is
pinned here too, because it is the reason only one of the two paths is instrumented.

All fixtures are synthetic; the designer is stubbed on every test, so nothing here needs ViennaRNA,
a network or Nextflow. TP53 is the documented public example gene.
"""

import asyncio
from pathlib import Path

import pytest
import typer
from rich.progress import Progress

from sirnaforge.cli import _fail_if_nothing_could_qualify
from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.core.selection import missing_required_evidence
from sirnaforge.data.base import DatabaseType, TranscriptInfo
from sirnaforge.models.policy import ExitCode, RunMode
from sirnaforge.models.sirna import DesignResult, SelectionState, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

pytestmark = pytest.mark.unit

GUIDE = "ACGTACGTACGTACGTACGTA"

#: Two transcripts in one batch, so a batch-level failure has more than one id to name.
TRANSCRIPT_IDS = ("ENST00000000001", "ENST00000000002")


def _transcript(transcript_id: str) -> TranscriptInfo:
    return TranscriptInfo(
        transcript_id=transcript_id,
        gene_id="ENSG00000000001",
        gene_name="TP53",
        sequence="ACGT" * 40,
        database=DatabaseType.ENSEMBL,
    )


def _candidate(candidate_id: str, transcript_id: str = TRANSCRIPT_IDS[0]) -> SiRNACandidate:
    """A minimal passing candidate; the design-stage gates are not what these tests exercise."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id=transcript_id,
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE.translate(str.maketrans("ATCG", "TAGC"))[::-1],
        length=len(GUIDE),
        gc_content=47.6,
        asymmetry_score=0.5,
        design_score=50.0,
        component_scores={"target_accessibility": 0.5, "asymmetry": 0.5, "gc_content": 0.5},
    )


def _design_result(workflow: SiRNAWorkflow, candidates: list[SiRNACandidate]) -> DesignResult:
    return DesignResult(
        input_file="<test>",
        parameters=workflow.config.design_params,
        candidates=list(candidates),
        top_candidates=list(candidates),
        total_sequences=len(TRANSCRIPT_IDS),
        total_candidates=len(candidates),
        filtered_candidates=len(candidates),
        processing_time=0.0,
    )


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    run_mode: RunMode | None = None,
    input_fasta: Path | None = None,
) -> SiRNAWorkflow:
    """A workflow whose policy is resolved through the one resolver, as every real entry point does."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        run_mode=run_mode,
        query_species="human",
        screen_species=["human"],
    )
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        resolved_policy=policy,
        screen_species=["human"],
        query_species="human",
        input_fasta=input_fasta,
    )
    return SiRNAWorkflow(config)


def _run_design(workflow: SiRNAWorkflow, transcripts: list[TranscriptInfo]) -> DesignResult:
    """Drive the parallel design branch to completion and hand back its combined result."""
    with Progress() as progress:
        return asyncio.run(workflow.step3_design_sirnas(transcripts, progress))


def _designer_returning(workflow: SiRNAWorkflow, *, failing: set[str]) -> None:
    """Stub the designer so named transcripts raise and the rest each yield one candidate."""

    def _design_from_sequence(sequence: str, transcript_id: str, **_kwargs: object) -> DesignResult:
        if transcript_id in failing:
            raise RuntimeError(f"ViennaRNA folded nothing for {transcript_id}")
        return _design_result(workflow, [_candidate(f"{transcript_id}_1", transcript_id)])

    assert workflow.sirnaforgeer is not None
    workflow.sirnaforgeer.design_from_sequence = _design_from_sequence  # type: ignore[method-assign]


class TestABatchFailureIsRecorded:
    """The parallel branch's ``except``: what the run lost, and where it says so."""

    def test_a_failed_batch_names_every_transcript_it_dropped(self, tmp_path: Path) -> None:
        """A batch whose future raised contributes no candidate, so every id in it is unrepresented."""
        workflow = _workflow(tmp_path, "failed-batch")

        def _explode(batch: list[TranscriptInfo]) -> tuple[list[DesignResult], dict[str, set[str]], dict[str, str]]:
            raise RuntimeError("worker crashed")

        workflow._process_transcript_batch = _explode  # type: ignore[method-assign]
        result = _run_design(workflow, [_transcript(tid) for tid in TRANSCRIPT_IDS])

        assert set(workflow._design_input_shortfalls) == set(TRANSCRIPT_IDS)
        for reason in workflow._design_input_shortfalls.values():
            # The exception's own text, not a generic "design failed": the fix for a crashed worker is
            # not the fix for an unfoldable sequence.
            assert "worker crashed" in reason
        assert result.candidates == []

    def test_the_design_summary_publishes_the_dropped_ids(self, tmp_path: Path) -> None:
        """``design_summary`` is where the divergence it explains is visible.

        ``input_sequences`` still counts the dropped transcripts -- it is the number handed to design,
        and narrowing it would erase the discrepancy instead of naming it.
        """
        workflow = _workflow(tmp_path, "summary")

        def _explode(batch: list[TranscriptInfo]) -> tuple[list[DesignResult], dict[str, set[str]], dict[str, str]]:
            raise RuntimeError("worker crashed")

        workflow._process_transcript_batch = _explode  # type: ignore[method-assign]
        result = _run_design(workflow, [_transcript(tid) for tid in TRANSCRIPT_IDS])

        summary = workflow._summarize_design_results(result)
        assert set(summary["design_input_shortfalls"]) == set(TRANSCRIPT_IDS)
        assert summary["design_input_dropped_count"] == len(TRANSCRIPT_IDS)
        assert summary["input_sequences"] == len(TRANSCRIPT_IDS)

    def test_a_transcript_that_failed_inside_a_surviving_batch_is_recorded(self, tmp_path: Path) -> None:
        """The same loss by the other route: one transcript raises, its batch mates still design."""
        workflow = _workflow(tmp_path, "one-transcript")
        _designer_returning(workflow, failing={TRANSCRIPT_IDS[1]})

        result = _run_design(workflow, [_transcript(tid) for tid in TRANSCRIPT_IDS])

        assert set(workflow._design_input_shortfalls) == {TRANSCRIPT_IDS[1]}
        assert "folded nothing" in workflow._design_input_shortfalls[TRANSCRIPT_IDS[1]]
        # The surviving transcript's candidate is still published: dropping the batch's whole output
        # over one bad sequence is the failure mode this branch exists to avoid.
        assert any(candidate.transcript_id == TRANSCRIPT_IDS[0] for candidate in result.candidates)

    def test_the_worker_returns_its_shortfalls_rather_than_writing_them(self, tmp_path: Path) -> None:
        """``_process_transcript_batch`` runs in a thread pool, so only its caller mutates run state."""
        workflow = _workflow(tmp_path, "worker-returns")
        _designer_returning(workflow, failing={TRANSCRIPT_IDS[1]})

        _results, _mapping, shortfalls = workflow._process_transcript_batch(
            [_transcript(tid) for tid in TRANSCRIPT_IDS]
        )

        assert set(shortfalls) == {TRANSCRIPT_IDS[1]}
        assert workflow._design_input_shortfalls == {}

    def test_a_reason_is_never_blank(self, tmp_path: Path) -> None:
        """``str(exc)`` is empty for a bare ``RuntimeError()``; a shortfall with no reason is useless."""
        workflow = _workflow(tmp_path, "blank-reason")

        def _explode(batch: list[TranscriptInfo]) -> tuple[list[DesignResult], dict[str, set[str]], dict[str, str]]:
            raise RuntimeError

        workflow._process_transcript_batch = _explode  # type: ignore[method-assign]
        _run_design(workflow, [_transcript(TRANSCRIPT_IDS[0])])

        assert workflow._design_input_shortfalls[TRANSCRIPT_IDS[0]] == "RuntimeError"

    def test_a_clean_design_records_nothing(self, tmp_path: Path) -> None:
        """The record is written only by a failure, so a whole design still claims a whole target set."""
        workflow = _workflow(tmp_path, "clean")
        _designer_returning(workflow, failing=set())

        result = _run_design(workflow, [_transcript(tid) for tid in TRANSCRIPT_IDS])

        assert workflow._design_input_shortfalls == {}
        summary = workflow._summarize_design_results(result)
        assert summary["design_input_shortfalls"] == {}
        assert summary["design_input_dropped_count"] == 0

    def test_the_non_parallel_branch_hard_crashes_instead(self, tmp_path: Path) -> None:
        """The inverse failure mode, and why only the parallel branch carries the record.

        The single-call branch has no ``except``: the failure propagates, the run stops and nothing is
        published, so there is no summary that could claim an uncovered target set.
        """
        workflow = _workflow(tmp_path, "single-call", input_fasta=tmp_path / "guides.fasta")

        def _explode(_path: str) -> DesignResult:
            raise RuntimeError("designer died")

        assert workflow.sirnaforgeer is not None
        workflow.sirnaforgeer.design_from_file = _explode  # type: ignore[method-assign]

        with pytest.raises(RuntimeError, match="designer died"):
            _run_design(workflow, [_transcript(TRANSCRIPT_IDS[0])])
        assert workflow._design_input_shortfalls == {}


class TestSelectionCannotCallTheTargetSetComplete:
    """The record reaches the eligibility decision through ``_selection_inputs``."""

    def test_the_selection_inputs_carry_the_shortfall(self, tmp_path: Path) -> None:
        """One builder feeds every entry point, so the record travels with the run's policy."""
        workflow = _workflow(tmp_path, "inputs", run_mode=RunMode.QUALIFIED)
        workflow._design_input_shortfalls = {TRANSCRIPT_IDS[0]: "worker crashed"}

        inputs = workflow._selection_inputs()

        assert inputs.design_input_shortfalls == {TRANSCRIPT_IDS[0]: "worker crashed"}
        assert f"design_input:{TRANSCRIPT_IDS[0]}" in missing_required_evidence(inputs)

    def test_the_inputs_are_a_snapshot(self, tmp_path: Path) -> None:
        """A later write must not retroactively change a decision already made."""
        workflow = _workflow(tmp_path, "snapshot", run_mode=RunMode.QUALIFIED)
        inputs = workflow._selection_inputs()

        workflow._design_input_shortfalls[TRANSCRIPT_IDS[0]] = "worker crashed"

        assert inputs.design_input_shortfalls == {}

    def test_a_qualified_run_reports_the_dropped_target(self, tmp_path: Path) -> None:
        """``required_evidence_missing`` names it, and ``design_input_excluded`` counts it."""
        workflow = _workflow(tmp_path, "qualified", run_mode=RunMode.QUALIFIED)
        workflow._design_input_shortfalls = {TRANSCRIPT_IDS[1]: "worker crashed"}
        workflow._completed_evidence_pairs = frozenset({("transcriptome", "human")})
        candidate = _candidate("survivor")
        candidate.off_target_screened = True

        workflow._apply_selection([candidate])

        summary = workflow._selection_summary
        assert f"design_input:{TRANSCRIPT_IDS[1]}" in summary["required_evidence_missing"]
        assert summary["design_input_excluded"] == 1

    def test_the_guides_that_were_designed_are_still_eligible(self, tmp_path: Path) -> None:
        """A batch failure denies a claim about coverage, not the evidence behind a designed guide.

        The surviving guides were screened exactly as well as before, so withholding them would punish
        them for a transcript they never came from -- and would empty the shortlist of a run whose
        only fault is a target it must not claim to have covered.
        """
        workflow = _workflow(tmp_path, "still-eligible", run_mode=RunMode.QUALIFIED)
        workflow._design_input_shortfalls = {TRANSCRIPT_IDS[1]: "worker crashed"}
        workflow._completed_evidence_pairs = frozenset({("transcriptome", "human")})
        candidate = _candidate("survivor")
        candidate.off_target_screened = True

        shortlist = workflow._apply_selection([candidate])

        assert candidate.selection_state == SelectionState.ELIGIBLE.value
        assert [c.id for c in shortlist] == ["survivor"]
        assert workflow._selection_summary["evidence_excluded"] == 0

    def test_a_lost_target_alone_reaches_the_incomplete_evidence_exit(self, tmp_path: Path) -> None:
        """Exit 2 through the check W4 already wrote, with no CLI change of its own.

        Every candidate here was excluded by a gate before the evidence test, so the run reports no
        per-candidate reason at all -- the dropped transcript is the only named cause, and without it
        an empty shortlist over an uncovered target set would exit 0 as a legitimate result.
        """
        workflow = _workflow(tmp_path, "exit-code", run_mode=RunMode.QUALIFIED)
        workflow._design_input_shortfalls = {TRANSCRIPT_IDS[1]: "worker crashed"}
        workflow._completed_evidence_pairs = frozenset({("transcriptome", "human")})
        candidate = _candidate("gate-failer")
        candidate.off_target_screened = True
        candidate.passes_filters = False

        workflow._apply_selection([candidate])
        results = {"selection_summary": workflow._selection_summary}
        assert results["selection_summary"]["evidence_shortfall_reasons"] == {}

        with pytest.raises(typer.Exit) as raised:
            _fail_if_nothing_could_qualify(results, json_summary=True)
        assert int(raised.value.exit_code) == int(ExitCode.INCOMPLETE_EVIDENCE)

    def test_a_run_with_no_design_stage_reports_no_shortfall(self, tmp_path: Path) -> None:
        """``run_offtarget_only_workflow`` designs nothing, so it can lose no design input."""
        workflow = _workflow(tmp_path, "no-design", run_mode=RunMode.QUALIFIED)
        workflow._completed_evidence_pairs = frozenset({("transcriptome", "human")})
        candidate = _candidate("pre-designed")
        candidate.off_target_screened = True

        workflow._apply_selection([candidate])

        assert workflow._selection_summary["required_evidence_missing"] == ()
        assert workflow._selection_summary["design_input_excluded"] == 0
