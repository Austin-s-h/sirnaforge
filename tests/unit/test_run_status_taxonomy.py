"""One run-status vocabulary, and exit codes a caller can depend on (#100, W4).

Screening published seven ``{status, reason}`` string pairs invented at each of its seven exits and
backed by no enum: ``skipped`` meant both "the user switched the channel off" and "the pipeline
aborted", ``partial`` meant both "one species produced no evidence" and "no reference alignment ran
at all", and a step whose aggregate self-reported ``completed`` said so even when the reconciled
evidence knew a *required* unit had published nothing. The CLI collapsed all of it into one generic
"⚠️  Partial" cell, and collapsed every failure into one bare ``Exit(1)``.

Two rules this file exists to hold:

1. ``run_status`` is published BESIDE ``status``/``reason``, never instead of them. Every consumer
   comparing ``status == "completed"`` -- the table, #103's report, result directories already on
   disk -- must keep working, so each exit's own keys are asserted unchanged next to the new one.
2. The exit code is decided from *selection*, never from the screen's ``run_status``. An incomplete
   screen whose surviving candidates still qualified produced a usable deliverable and exits 0;
   deriving the code from the screen would fail every run on a machine with no Nextflow, and a
   complete run that legitimately qualified nobody has always exited 0.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any

import pytest
import typer

from sirnaforge.cli import (
    _fail_if_nothing_could_qualify,
    _fail_if_nothing_was_eligible,
    _offtarget_status_cell,
    app,
)
from sirnaforge.models.policy import ExitCode, RunMode, RunStatus
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig, screening_run_status

GUIDE = "ACGTACGTACGTACGTACGTA"

#: Offline orthologue evidence, so nothing here touches Ensembl Compara.
_ORTHOLOG_MAPPING_FIXTURE = Path(__file__).resolve().parent / "data" / "ortholog_mapping_synthetic.json"

#: Every exit out of screening, as ``{status, reason, method}`` and the run status it means. The
#: pairs are the exits' own literals: this table is the defect, restated as a mapping, and each row
#: is driven through the product code below rather than only through the classifier.
EXIT_TAXONOMY: dict[str, tuple[dict[str, str], RunStatus]] = {
    "no_candidates": ({"status": "skipped", "reason": "no_candidates"}, RunStatus.NO_ELIGIBLE),
    "user_disabled": ({"status": "skipped", "reason": "user_disabled"}, RunStatus.NOT_REQUESTED),
    "nextflow_failed": ({"status": "skipped", "reason": "nextflow_failed"}, RunStatus.EXECUTION_ERROR),
    "nextflow_unavailable": ({"status": "skipped", "reason": "nextflow_unavailable"}, RunStatus.EXECUTION_ERROR),
    "basic_fallback": ({"status": "partial", "method": "basic"}, RunStatus.INCOMPLETE),
    # ``_parse_nextflow_results``' own word for an output directory that never appeared. The route
    # that hits it publishes ``partial`` and carries ``missing_output`` as the reason, so both
    # spellings have to land on the same answer.
    "missing_results": ({"status": "missing", "method": "nextflow"}, RunStatus.EXECUTION_ERROR),
    "missing_output": ({"status": "partial", "reason": "missing_output"}, RunStatus.EXECUTION_ERROR),
    "completed": ({"status": "completed", "method": "embedded_nextflow"}, RunStatus.COMPLETED),
}

#: The four exits ``step5_offtarget_analysis`` can be driven to without a screen, with the keys each
#: publishes. Copied from ``test_zero_hit_screen_gating.py`` rather than imported, per the contract's
#: separate-files rule: the two slices' drivers must be able to fail independently.
NO_SCREEN_PATHS: dict[str, dict[str, str]] = {
    "nextflow_failed": {"status": "skipped", "reason": "nextflow_failed"},
    "nextflow_unavailable": {"status": "skipped", "reason": "nextflow_unavailable"},
    "basic_fallback": {"status": "partial", "method": "basic"},
    "missing_output": {"status": "partial", "method": "embedded_nextflow"},
}


def _candidate(candidate_id: str = "cand_1") -> SiRNACandidate:
    """A minimal candidate carrying one guide."""
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


def _reference_fasta(tmp_path: Path) -> Path:
    """A cDNA file standing in for a configured screening reference; nothing here resolves it."""
    path = tmp_path / "reference_cdna.fa"
    if not path.exists():
        path.write_text(">ENST00000000001 cdna gene:ENSG00000000001 gene_symbol:TP53\n" + GUIDE * 3 + "\n")
    return path


def _guides(tmp_path: Path) -> Path:
    """The staged guide FASTA every digest in a run is keyed on."""
    path = tmp_path / "input_candidates.fasta"
    path.write_text(f">cand_1\n{GUIDE}\n")
    return path


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    species: list[str] | None = None,
    check_off_targets: bool = True,
) -> SiRNAWorkflow:
    """A TP53 workflow with a configured reference, so its screen has planned units to report on.

    A run with no reference configured requests no transcriptome screen at all, and a channel nobody
    asked for has no plan entry by design -- so there would be nothing for a run status to disagree
    with.
    """
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        screen_species=species or ["human"],
        query_species="human",
        design_params=DesignParameters(check_off_targets=check_off_targets),
        transcriptome_fasta=str(_reference_fasta(tmp_path)),
        ortholog_mapping_file=_ORTHOLOG_MAPPING_FIXTURE,
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


def _run_step5_without_a_screen(workflow: SiRNAWorkflow, candidates: list[SiRNACandidate], path: str) -> dict[str, Any]:
    """Drive ``step5_offtarget_analysis`` to one of :data:`NO_SCREEN_PATHS` and return its summary.

    Reference resolution and repeat detection are stubbed on every path: both need a real cDNA
    reference, and neither is what is under test. Everything downstream of them is product code.
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
            # The runner reported success and wrote nothing: parsing returns "missing", so the run
            # never had an output directory to read.
            return await workflow._process_nextflow_results(
                candidates, workflow.config.output_dir / "absent", {"status": "completed"}
            )

        workflow._run_nextflow_offtarget_analysis = _no_output  # type: ignore[method-assign]
    else:  # pragma: no cover - a typo in the parametrisation, not a product path
        raise AssertionError(f"unknown path {path}")

    return asyncio.run(workflow.step5_offtarget_analysis(design_result))


def _legacy_results(workflow: SiRNAWorkflow, *, species_screened: list[str]) -> Path:
    """A results directory in the pre-#100 shape: an aggregate summary, and no evidence anywhere.

    The summary names no missing species, which is exactly why the run's own word stays ``completed``:
    an aggregate cannot report the species it never heard about.
    """
    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    (aggregated / "combined_summary.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "species_screened": species_screened,
                "species_analyzed": species_screened,
                "hits_per_species": dict.fromkeys(species_screened, 0),
                "total_candidates": 1,
            }
        )
    )
    return results_dir


def _selection(**overrides: Any) -> dict[str, Any]:
    """A qualified selection summary that could qualify nobody for want of evidence."""
    summary: dict[str, Any] = {
        "run_mode": RunMode.QUALIFIED.value,
        "eligible_candidates": 0,
        "evidence_excluded": 3,
        "filter_excluded": 0,
        "repeat_excluded": 0,
        "unscored_excluded": 0,
        "evidence_shortfall_reasons": {"no_evidence:transcriptome:human": 3},
        "required_evidence_missing": ["transcriptome:human"],
    }
    summary.update(overrides)
    return summary


# ---------------------------------------------------------------------------
# 1. The vocabulary
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_vocabulary_and_the_exit_codes_are_the_declared_wire_values() -> None:
    """Both enums are a published contract: a renamed member is a breaking change for a caller."""
    assert [status.value for status in RunStatus] == [
        "completed",
        "incomplete",
        "not_requested",
        "no_eligible_candidates",
        "execution_error",
    ]
    assert (ExitCode.SUCCESS, ExitCode.EXECUTION_ERROR) == (0, 1)
    assert (ExitCode.INCOMPLETE_EVIDENCE, ExitCode.NO_ELIGIBLE_CANDIDATES) == (2, 3)


# ---------------------------------------------------------------------------
# 2. One classifier for all seven exits
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("name", sorted(EXIT_TAXONOMY))
def test_every_exit_string_pair_has_one_run_status(name: str) -> None:
    """The seven pairs were uncoordinated strings; this is the single definition of what they mean.

    ``skipped`` is the pair that carried two opposite meanings, and the reason is what separates
    them: a channel the user switched off was never requested, an engine that could not be executed
    is an error the caller must act on.
    """
    summary, expected = EXIT_TAXONOMY[name]
    assert screening_run_status(status=summary["status"], reason=summary.get("reason", "")) is expected, (
        f"{name} must read as {expected.value}"
    )


@pytest.mark.unit
def test_a_completed_step_reads_incomplete_when_a_required_unit_published_nothing() -> None:
    """``status`` speaks for the aggregate's self-report; the run status speaks for the evidence.

    An aggregate that never heard of a species cannot report it missing, so ``completed`` is exactly
    the word a run whose aggregation failed publishes. Only the reconciled plan knows better.
    """
    assert screening_run_status(status="completed", required_evidence_missing=True) is RunStatus.INCOMPLETE
    assert screening_run_status(status="completed", required_evidence_missing=False) is RunStatus.COMPLETED
    # A reason that already names the outcome is not overridden by an evidence shortfall: a channel
    # nobody requested has no required unit to miss, and an execution error is the more useful word.
    assert (
        screening_run_status(status="skipped", reason="user_disabled", required_evidence_missing=True)
        is RunStatus.NOT_REQUESTED
    )


@pytest.mark.unit
def test_an_unrecognised_exit_claims_no_completeness() -> None:
    """A new exit must not arrive silently reading ``completed``: the default is the honest one."""
    assert screening_run_status(status="", reason="") is RunStatus.INCOMPLETE
    assert screening_run_status(status="skipped", reason="something_new") is RunStatus.INCOMPLETE


# ---------------------------------------------------------------------------
# 3. The real exits publish it, beside the keys they already published
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("path", sorted(NO_SCREEN_PATHS))
def test_each_no_screen_exit_publishes_its_run_status_beside_its_own_keys(tmp_path: Path, path: str) -> None:
    """The four exits that publish nothing now say which of them happened, in one word.

    ``nextflow_failed``/``nextflow_unavailable``/``missing_output`` are execution errors -- nothing
    ran, or ran and staged nothing -- while the sequence-only fallback did produce a (thin) answer
    and is merely incomplete. Before this slice all four read as ``skipped`` or ``partial`` and the
    CLI printed one "⚠️  Partial" cell for every one of them.
    """
    workflow = _workflow(tmp_path, f"exit_{path}")

    outcome = _run_step5_without_a_screen(workflow, [_candidate()], path)

    # Beside, never instead: every key these exits already published is untouched.
    for key, value in NO_SCREEN_PATHS[path].items():
        assert outcome[key] == value, "an existing consumer of status/reason/method must keep working"
    expected = RunStatus.INCOMPLETE if path == "basic_fallback" else RunStatus.EXECUTION_ERROR
    assert outcome["run_status"] == expected.value


@pytest.mark.unit
def test_a_disabled_channel_is_not_requested_rather_than_incomplete(tmp_path: Path) -> None:
    """``--skip-off-targets`` asked for no screen, so there is nothing for it to be missing.

    The distinction matters most for a design-only run: reporting it as incomplete would make the
    one mode that claims no screening evidence look like a broken screen.
    """
    workflow = _workflow(tmp_path, "disabled", check_off_targets=False)

    outcome = asyncio.run(workflow.step5_offtarget_analysis(_design_result(workflow, [_candidate()])))

    assert outcome["status"] == "skipped"
    assert outcome["reason"] == "user_disabled"
    assert outcome["run_status"] == RunStatus.NOT_REQUESTED.value


@pytest.mark.unit
def test_nothing_to_screen_is_reported_as_no_eligible_candidates(tmp_path: Path) -> None:
    """Design produced nothing, so screening has no work: a result, not a failure of the screen."""
    workflow = _workflow(tmp_path, "no_candidates")

    outcome = asyncio.run(workflow.step5_offtarget_analysis(_design_result(workflow, [])))

    assert outcome["status"] == "skipped"
    assert outcome["reason"] == "no_candidates"
    assert outcome["run_status"] == RunStatus.NO_ELIGIBLE.value


# ---------------------------------------------------------------------------
# 4. The evidence is what a "completed" screen is measured against
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_completed_screen_missing_its_required_species_publishes_incomplete(tmp_path: Path) -> None:
    """The run's own word and its evidence disagreed, and nothing published the disagreement (W1).

    The aggregate here says it screened rat and names no missing species -- so ``status`` stays
    ``completed`` and the alignment-evidence heuristic is satisfied -- while the *required* human
    transcriptome unit reconciles as FAILED. That is the shape of a run whose required index build
    crashed: every published field said "complete", and the only record of the gap was one log line.
    """
    workflow = _workflow(tmp_path, "required_missing", species=["human", "rat"])
    workflow._record_screening_plan(_guides(tmp_path), {})
    results_dir = _legacy_results(workflow, species_screened=["rat"])

    outcome = asyncio.run(workflow._process_nextflow_results([_candidate()], results_dir, {"status": "completed"}))

    assert workflow.config.resolved_policy.evidence_requirements.required_pairs == frozenset(
        {("transcriptome", "human")}
    )
    assert workflow._missing_required_evidence_units() == (("transcriptome", "human"),)
    # Beside, never instead: W1's pinned word is still what a legacy consumer reads.
    assert outcome["status"] == "completed"
    assert outcome["run_status"] == RunStatus.INCOMPLETE.value


@pytest.mark.unit
def test_a_screen_that_covered_every_required_unit_is_completed(tmp_path: Path) -> None:
    """The negative control: an exploratory gap does not make a run incomplete.

    Only ``EvidenceRequirements`` decides requiredness, and in 0.7.1 every miRNA pair and every
    non-query species is exploratory -- so a run that screened its required species reports
    ``completed`` even though the miRNA channel published nothing at all.
    """
    workflow = _workflow(tmp_path, "required_present", species=["human"])
    workflow._record_screening_plan(_guides(tmp_path), {})
    results_dir = _legacy_results(workflow, species_screened=["human"])

    outcome = asyncio.run(workflow._process_nextflow_results([_candidate()], results_dir, {"status": "completed"}))

    assert workflow._missing_required_evidence_units() == ()
    assert outcome["status"] == "completed"
    assert outcome["run_status"] == RunStatus.COMPLETED.value


@pytest.mark.unit
def test_a_run_that_reconciled_no_plan_claims_no_shortfall(tmp_path: Path) -> None:
    """A direct call into result processing recorded no plan, so it cannot name a missing unit.

    Absence of evidence is not evidence of a shortfall here: without a plan there is no record that
    anything was expected, and inventing one would report every unit tests drive directly as failed.
    """
    workflow = _workflow(tmp_path, "no_plan", species=["human"])

    assert workflow._screening_evidence is None
    assert workflow._missing_required_evidence_units() == ()
    assert workflow._publish_run_status({"status": "completed"})["run_status"] == RunStatus.COMPLETED.value, (
        "with nothing reconciled, the aggregate's own word stands"
    )


# ---------------------------------------------------------------------------
# 5. Exit codes
# ---------------------------------------------------------------------------


def _exit_code(results: dict[str, Any], *, fail_on_no_eligible: bool = False) -> int:
    """The code the ``workflow`` command would exit with, by running its two checks in its order.

    Both calls, in sequence, because the order is part of the taxonomy: a shortlist that is empty for
    want of evidence must report the shortfall (2), not merely its emptiness (3).
    """
    try:
        _fail_if_nothing_could_qualify(results, json_summary=True)
        _fail_if_nothing_was_eligible(results, enabled=fail_on_no_eligible)
    except typer.Exit as exc:
        return int(exc.exit_code)
    return int(ExitCode.SUCCESS)


@pytest.mark.unit
def test_incomplete_required_evidence_exits_two() -> None:
    """The one incomplete-evidence code, separated from the execution error it used to share.

    Both bare ``Exit(1)`` sites meant "something went wrong" and a caller could not tell a crashed
    run from one whose screen simply never produced the evidence it required.
    """
    assert _exit_code({"selection_summary": _selection()}) == int(ExitCode.INCOMPLETE_EVIDENCE)
    # And it still wins when the shortlist is also empty, which it always is here: the shortfall is
    # the more specific and more actionable of the two answers.
    assert _exit_code({"selection_summary": _selection()}, fail_on_no_eligible=True) == int(
        ExitCode.INCOMPLETE_EVIDENCE
    )


@pytest.mark.unit
def test_a_complete_run_that_qualified_nobody_still_exits_zero() -> None:
    """Every over-tight threshold would read as a tool failure otherwise.

    Nothing is missing here: the candidates were considered and rejected on their merits, so neither
    check fires. Turning this into a non-zero exit by default would break any CI already invoking
    the command, which is why 3 exists but is opt-in.
    """
    complete = _selection(
        evidence_excluded=0, filter_excluded=4, evidence_shortfall_reasons={}, required_evidence_missing=[]
    )

    assert _exit_code({"selection_summary": complete}) == int(ExitCode.SUCCESS)
    assert _exit_code({"selection_summary": complete}, fail_on_no_eligible=True) == int(ExitCode.NO_ELIGIBLE_CANDIDATES)


@pytest.mark.unit
@pytest.mark.parametrize(
    "selection",
    [
        # Something qualified: the shortlist is not empty.
        {"run_mode": RunMode.QUALIFIED.value, "eligible_candidates": 2, "filter_excluded": 4},
        # Nothing reached selection at all: no candidate was judged, so none was found ineligible.
        {"run_mode": RunMode.QUALIFIED.value, "eligible_candidates": 0},
        # ZFN publishes no selection summary.
        {},
    ],
)
def test_the_opt_in_flag_still_needs_an_empty_shortlist_to_fire(selection: dict[str, Any]) -> None:
    """Opt-in does not mean unconditional: the flag reports a verdict, never the absence of one."""
    assert _exit_code({"selection_summary": selection}, fail_on_no_eligible=True) == int(ExitCode.SUCCESS)


@pytest.mark.unit
def test_the_opt_in_flag_is_reachable_from_the_command_line() -> None:
    """A reserved exit code nobody can ask for is not an exit code: pins the flag's name and default.

    Read off the parsed command rather than out of ``--help``, whose rich rendering wraps long option
    names and would make this pass or fail on terminal width.
    """
    workflow_command = typer.main.get_command(app).commands["workflow"]  # type: ignore[attr-defined]
    flag = next(param for param in workflow_command.params if param.name == "fail_on_no_eligible")

    assert flag.opts == ["--fail-on-no-eligible"], "no paired --no- form: it must be opt-in only"
    assert flag.default is False, "an empty shortlist stays a zero-exit result unless asked otherwise"


@pytest.mark.unit
def test_a_failed_screen_does_not_by_itself_change_the_exit_code() -> None:
    """The exit code is decided from selection, never from the screen's own run status.

    A machine with no Nextflow reports ``execution_error`` for the screen on every run. If that
    decided the exit code, every such run would fail even though the design shortlist was produced
    and its candidates qualified on the evidence the run does have.
    """
    results = {
        "offtarget_summary": {
            "status": "skipped",
            "reason": "nextflow_unavailable",
            "run_status": RunStatus.EXECUTION_ERROR.value,
        },
        "selection_summary": {
            "run_mode": RunMode.EXPLORATORY.value,
            "eligible_candidates": 5,
            "filter_excluded": 1,
        },
    }

    assert _exit_code(results, fail_on_no_eligible=True) == int(ExitCode.SUCCESS)


# ---------------------------------------------------------------------------
# 6. What the user reads
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_summary_table_tells_the_five_outcomes_apart() -> None:
    """One generic "Partial" for an aborted pipeline and for a run that lacked one species is a lie.

    Both spellings the two tables used are pinned: the cell is the only place most users ever see
    this taxonomy.
    """
    cells = {status: _offtarget_status_cell({"status": "partial", "run_status": status.value}) for status in RunStatus}

    assert len(set(cells.values())) == len(RunStatus), f"two outcomes read alike: {cells}"
    assert cells[RunStatus.COMPLETED] == "✅ Complete"
    assert "Incomplete" in cells[RunStatus.INCOMPLETE]
    assert "Not requested" in cells[RunStatus.NOT_REQUESTED]
    assert "Did not run" in cells[RunStatus.EXECUTION_ERROR]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("summary", "expected"),
    [
        ({"status": "completed"}, "✅ Complete"),
        ({"status": "partial"}, "⚠️  Partial"),
        ({"status": "completed", "run_status": "a_word_from_the_future"}, "✅ Complete"),
    ],
)
def test_a_summary_with_no_usable_run_status_reads_exactly_as_it_did_before(
    summary: dict[str, Any], expected: str
) -> None:
    """A result dictionary from an older version still renders: the new key is additive on both sides."""
    assert _offtarget_status_cell(summary) == expected
