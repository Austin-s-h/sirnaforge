"""Regression tests for reference resolution on input-FASTA and skip-off-target runs.

Review finding 15: ``allow_transcriptome_with_input_fasta`` was flipped to True (commit
4b86fdf) and ``run_sirna_workflow`` hardcoded ``design_only=False``, so two runs that
should touch no remote reference resolved all four DEFAULT_TRANSCRIPTOME_SOURCES and
downloaded/indexed multi-gigabyte cDNA files:

* ``run_sirna_workflow(input_fasta=...)`` / ``sirnaforge workflow --input-fasta ...``
* ``check_off_targets=False`` / ``--skip-off-targets``

Every test here asserts on *resolution*, so none of them needs the network.
"""

import asyncio
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.config import ReferenceSelection
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig, run_sirna_workflow

GUIDE = "ATGCGATGCGATGCGATGCGC"


def _write_fasta(path: Path) -> Path:
    """Write a two-transcript FASTA usable as workflow input."""
    path.write_text(">trans1\n" + "ATG" + "A" * 300 + "TAA\n>trans2\n" + "ATG" + "C" * 300 + "TGA\n")
    return path


async def _capture_selection(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """Patch run_complete_workflow so the workflow config can be inspected without running it."""
    captured: dict[str, Any] = {}

    async def _fake_run(self: SiRNAWorkflow) -> dict[str, Any]:
        captured["selection"] = self.config.transcriptome_selection
        return {}

    monkeypatch.setattr(SiRNAWorkflow, "run_complete_workflow", _fake_run)
    return captured


@pytest.mark.asyncio
@pytest.mark.unit
async def test_library_input_fasta_run_resolves_no_remote_transcriptome(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A caller who supplies their own FASTA must not trigger a reference download."""
    captured = await _capture_selection(monkeypatch)

    await run_sirna_workflow(
        gene_query="TOY",
        output_dir=str(tmp_path / "out"),
        input_fasta=str(_write_fasta(tmp_path / "my.fa")),
    )

    selection: ReferenceSelection = captured["selection"]
    assert selection.enabled is False, "input FASTA must default to design-only, not four cDNA downloads"
    assert selection.choices == ()
    assert "design-only" in (selection.disabled_reason or "")


@pytest.mark.asyncio
@pytest.mark.unit
async def test_library_input_fasta_can_opt_in_to_defaults(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The opt-in is still available for callers who do want default screening."""
    captured = await _capture_selection(monkeypatch)

    await run_sirna_workflow(
        gene_query="TOY",
        output_dir=str(tmp_path / "out"),
        input_fasta=str(_write_fasta(tmp_path / "my.fa")),
        allow_transcriptome_with_input_fasta=True,
        default_transcriptome_sources=("ensembl_human_cdna",),
    )

    selection: ReferenceSelection = captured["selection"]
    assert selection.enabled is True
    assert [choice.value for choice in selection.choices] == ["ensembl_human_cdna"]


@pytest.mark.asyncio
@pytest.mark.unit
async def test_library_input_fasta_honours_explicit_transcriptome(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An explicit transcriptome argument still wins over the design-only default."""
    captured = await _capture_selection(monkeypatch)

    await run_sirna_workflow(
        gene_query="TOY",
        output_dir=str(tmp_path / "out"),
        input_fasta=str(_write_fasta(tmp_path / "my.fa")),
        transcriptome_fasta="ensembl_human_cdna",
    )

    selection: ReferenceSelection = captured["selection"]
    assert selection.enabled is True
    assert [choice.value for choice in selection.choices] == ["ensembl_human_cdna"]


@pytest.mark.asyncio
@pytest.mark.unit
async def test_library_skip_off_targets_resolves_no_remote_transcriptome(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """check_off_targets=False must not resolve references it will never screen against."""
    captured = await _capture_selection(monkeypatch)

    await run_sirna_workflow(
        gene_query="TOY",
        output_dir=str(tmp_path / "out"),
        check_off_targets=False,
    )

    selection: ReferenceSelection = captured["selection"]
    assert selection.enabled is False, "check_off_targets=False must not download/index references"
    assert "design-only" in (selection.disabled_reason or "")


def _cli_resolved_selection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra_args: str) -> ReferenceSelection:
    """Run `sirnaforge workflow ... --input-fasta` and return the selection the CLI handed down."""
    captured: dict[str, Any] = {}

    async def _fake_workflow(**kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return {"transcript_summary": {}, "design_summary": {}}

    monkeypatch.setattr("sirnaforge.cli.run_sirna_workflow", _fake_workflow)

    result = CliRunner().invoke(
        app,
        [
            "workflow",
            "TOY",
            "--input-fasta",
            str(_write_fasta(tmp_path / "toy.fa")),
            "--output-dir",
            str(tmp_path / "out"),
            "--species",
            "human",
            *extra_args,
        ],
    )

    assert result.exit_code == 0, result.output
    selection: ReferenceSelection = captured["transcriptome_selection"]
    return selection


@pytest.mark.unit
def test_cli_input_fasta_workflow_resolves_no_remote_transcriptome(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The exact container scenario: `sirnaforge workflow TOY --input-fasta ... --species human`."""
    selection = _cli_resolved_selection(tmp_path, monkeypatch)

    assert selection.enabled is False, "the toy container workflow must not download four cDNA references"
    assert "design-only" in (selection.disabled_reason or "")


@pytest.mark.unit
def test_cli_transcriptome_fasta_still_enables_screening(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The documented escape hatch keeps working: name a reference and screening is back on."""
    selection = _cli_resolved_selection(tmp_path, monkeypatch, "--transcriptome-fasta", "ensembl_human_cdna")

    assert selection.enabled is True
    assert [choice.value for choice in selection.choices] == ["ensembl_human_cdna"]


def _candidate(candidate_id: str) -> SiRNACandidate:
    """Minimal candidate, enough for step 5 to accept it as screenable."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE.translate(str.maketrans("ATCG", "TAGC")),
        gc_content=57.1,
        length=len(GUIDE),
        asymmetry_score=0.7,
        composite_score=50.0,
        component_scores={"empirical": 1.0},
    )


@pytest.mark.asyncio
@pytest.mark.unit
async def test_skip_off_targets_materializes_no_reference_and_runs_no_repeat_scan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Even with an enabled selection, check_off_targets=False must do no reference work.

    Pre-fix, step5 called _configure_transcriptome_inputs and _run_repeat_detection *above*
    the check_off_targets guard, so a run that asked to skip off-target analysis still
    downloaded, indexed and scanned before printing "skipped by user request".
    """
    params = DesignParameters(check_off_targets=False, apply_modifications=False)
    config = WorkflowConfig(
        output_dir=tmp_path / "skip_out",
        gene_query="TOY",
        design_params=params,
        # An explicitly enabled selection: the guard, not the resolver, must stop the work here.
        transcriptome_fasta="ensembl_human_cdna",
    )
    workflow = SiRNAWorkflow(config)
    assert workflow.config.transcriptome_selection.enabled is True

    calls: list[str] = []

    async def _no_materialize(self: SiRNAWorkflow, choice: object) -> None:
        calls.append("materialize")

    def _no_repeat_scan(self: SiRNAWorkflow, candidates: object) -> dict[str, Any]:
        calls.append("repeat_scan")
        return {"status": "completed"}

    monkeypatch.setattr(SiRNAWorkflow, "_materialize_transcriptome_reference", _no_materialize)
    monkeypatch.setattr(SiRNAWorkflow, "_run_repeat_detection", _no_repeat_scan)

    candidate = _candidate("cand_skip")
    design_result = DesignResult(
        input_file="<test>",
        parameters=params,
        candidates=[candidate],
        top_candidates=[candidate],
        total_sequences=1,
        total_candidates=1,
        filtered_candidates=1,
        processing_time=0.1,
    )

    result = await workflow.step5_offtarget_analysis(design_result)

    assert result == {"status": "skipped", "reason": "user_disabled"}
    assert calls == [], f"skip request must precede all reference work, but ran: {calls}"
    assert workflow._repeat_summary["reason"] == "user_disabled"
    assert design_result.top_candidates == [candidate], "ranking must still be rebuilt on the skip path"


@pytest.mark.unit
def test_input_fasta_step5_materializes_no_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Step 5 of the container scenario must fetch nothing, end to end.

    The policy decision happens at config time but the download happens in step 5, so this
    drives step 5 with the selection the real CLI hands down for the toy invocation rather
    than a hand-built one -- pre-fix that selection named four Ensembl cDNA sources and step 5
    materialized every one of them.
    """
    config = WorkflowConfig(
        output_dir=tmp_path / "toy_out",
        gene_query="TOY",
        input_fasta=_write_fasta(tmp_path / "toy.fa"),
        design_params=DesignParameters(apply_modifications=False),  # check_off_targets stays True
        transcriptome_selection=_cli_resolved_selection(tmp_path, monkeypatch),
    )
    workflow = SiRNAWorkflow(config)

    calls: list[str] = []

    async def _no_materialize(self: SiRNAWorkflow, choice: object) -> None:
        calls.append("materialize")

    monkeypatch.setattr(SiRNAWorkflow, "_materialize_transcriptome_reference", _no_materialize)

    candidate = _candidate("cand_toy")
    design_result = DesignResult(
        input_file="<test>",
        parameters=config.design_params,
        candidates=[candidate],
        top_candidates=[candidate],
        total_sequences=1,
        total_candidates=1,
        filtered_candidates=1,
        processing_time=0.1,
    )

    # Sync test on purpose: _cli_resolved_selection uses CliRunner, and the CLI's own
    # asyncio.run() cannot be nested inside a pytest-asyncio event loop.
    asyncio.run(workflow.step5_offtarget_analysis(design_result))

    assert calls == [], f"an input-FASTA run must fetch no transcriptome reference, but ran: {calls}"
    assert workflow._repeat_summary["reason"] == "reference_unavailable"
