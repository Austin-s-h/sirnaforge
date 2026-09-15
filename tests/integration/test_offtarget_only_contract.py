"""``sirnaforge offtarget`` publishes the same selection and the same tables as a design run (#100).

W3 of the WORKFLOW-LOCK chain. ``run_offtarget_only_workflow`` bypasses step5 and step6, and so used
to bypass everything those two steps added over #100:

- no gate reading a screening channel was reached on a screen that published nothing, so a screen that
  never happened exported the same cells as one that came back clean;
- selection never ran, so every guide the command screened kept ``SelectionState.NOT_SELECTED`` -- the
  model default, meaning "selection never looked at this candidate";
- no candidate table was written at all, so ``sirnaforge report <output_dir>`` -- which this command
  itself prints as the next step -- raised ``ReportInputError("no candidates_all.csv ...")``.

What the entry point must NOT gain is a failure. A pre-designed guide has no design score and no
transcript context, so it can hold a qualified *screen* result with no potency composite at all; that
is score-unavailable, and a run that qualified nobody is a result rather than an error exit.

The screen itself is stubbed at ``SiRNAWorkflow.run_nextflow_offtarget_analysis``: these tests are
about what the entry point does with a screen's outcome, and no Nextflow, Docker or network is
available here. Guides are synthetic.
"""

import asyncio
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.models.policy import FilterEvaluation, RunMode
from sirnaforge.models.sirna import SelectionState, SiRNACandidate
from sirnaforge.reporting import build_payload
from sirnaforge.workflow import _ZERO_OFFTARGET_COUNTS, SiRNAWorkflow, run_offtarget_only_workflow

#: Also in the dev tier, not only the release tier ``integration`` auto-assigns: the whole file runs in
#: under a second with no network, Nextflow or Docker, and it guards an entry point three later slices
#: build on.
pytestmark = pytest.mark.dev

GUIDES = {
    "HD-001": "ACGTACGTACGTACGTACGTA",
    "HD-002": "GGCATTACGGCATTACGGCAT",
    "HD-003": "TTGCAAGGCTTACGGATCCAA",
}


def _input_fasta(tmp_path: Path) -> Path:
    """A pre-designed guide FASTA, the only input this entry point takes."""
    fasta = tmp_path / "predesigned.fa"
    fasta.write_text("".join(f">{name} pre-designed guide\n{seq}\n" for name, seq in GUIDES.items()))
    return fasta


async def _screen_published_nothing(
    self: SiRNAWorkflow, candidates: list[SiRNACandidate], input_fasta: Path
) -> dict[str, Any]:
    """A screen whose aggregate claims completion while no unit published evidence.

    The honest worst case, and the one the run must not paper over: ``status: completed`` from the
    runner only means Nextflow exited 0.
    """
    return {"status": "completed", "method": "nextflow", "results": {}, "aggregated": {}}


async def _screen_completed_clean(
    self: SiRNAWorkflow, candidates: list[SiRNACandidate], input_fasta: Path
) -> dict[str, Any]:
    """A completed screen that found nothing, gating every candidate on a measured zero.

    Mirrors what ``_integrate_offtarget_results`` does on the no-hit branch: the evidence the policy
    requires is complete, and every off-target gate is decided against ``_ZERO_OFFTARGET_COUNTS``.
    """
    self._active_screen_species = ["human"]
    self._screened_species_scope = self._requested_species_scope()
    self._completed_evidence_pairs = self.config.resolved_policy.evidence_requirements.required_pairs
    stats: dict[str, Any] = {}
    for candidate in candidates:
        candidate.off_target_screened = True
        self._gate_offtarget_counts(
            candidate,
            counts=_ZERO_OFFTARGET_COUNTS,
            filter_criteria=self.config.design_params.offtarget_filters,
            complete_pairs=self._completed_evidence_pairs,
            stats=stats,
        )
    return {"status": "completed", "method": "nextflow", "results": {}, "aggregated": {}}


def _run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, screen: Any, **kwargs: Any) -> dict[str, Any]:
    """Run the off-target-only entry point with the screen stubbed out."""
    monkeypatch.setattr(SiRNAWorkflow, "run_nextflow_offtarget_analysis", screen)
    return asyncio.run(
        run_offtarget_only_workflow(
            input_candidates_fasta=str(_input_fasta(tmp_path)),
            output_dir=str(tmp_path / "out"),
            screen_species=["human"],
            query_species="human",
            **kwargs,
        )
    )


def _tables(tmp_path: Path) -> Path:
    return tmp_path / "out" / "sirnaforge"


def _rows(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


# ---------------------------------------------------------------------------------------------
# The candidate tables, and the report that reads them
# ---------------------------------------------------------------------------------------------


@pytest.mark.integration
def test_the_offtarget_only_run_writes_the_candidate_table_its_own_advice_needs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The command tells the user to run ``sirnaforge report <output_dir>``; that needs a table.

    Before W3 this path wrote ``input_candidates.fasta`` and the off-target results and nothing else,
    so the printed follow-up instruction failed on the run it was printed for.
    """
    results = _run(tmp_path, monkeypatch, _screen_completed_clean)
    base = _tables(tmp_path)

    assert (base / "candidates_all.csv").exists()
    assert (base / "candidates_qualified.csv").exists()
    assert (base / "candidates_provisional.csv").exists()
    assert (base / "manifest.json").exists()
    assert len(_rows(base / "candidates_all.csv")) == len(GUIDES)
    # The result dictionary names what it wrote, so a caller need not guess the layout.
    assert results["written_files"]["candidates_all_csv"] == str(base / "candidates_all.csv")


@pytest.mark.integration
def test_the_report_builds_from_an_offtarget_only_directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``build_payload`` must accept this directory shape, and must not invent target context.

    A pre-designed guide carries the placeholder transcript id and position 1, which is not an
    enumeration on a real transcript: the report reports no isoform table and no design map for it
    rather than drawing one guide per fabricated isoform.
    """
    _run(tmp_path, monkeypatch, _screen_completed_clean)

    payload = build_payload(tmp_path / "out")

    assert len(payload.guides) == len(GUIDES)
    assert payload.run["transcripts"] == []
    assert all(guide.isoforms == [] for guide in payload.guides)
    assert any("no transcript context" in caveat for caveat in payload.caveats)
    # The policy came from this run's own manifest, not from library defaults.
    assert payload.provenance["policy_source"] == "the run's own manifest"
    assert payload.provenance["run_mode"] == RunMode.QUALIFIED.value


@pytest.mark.integration
def test_a_qualified_screen_result_needs_no_potency_composite(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The contract's central case: qualified on screening evidence, score-unavailable on potency.

    The guides were designed elsewhere, so this tool has no design score for them and no transcript to
    measure accessibility against. Filling either in from the guide sequence would be a fabricated
    potency claim -- so both stay null and the candidate is still eligible.
    """
    results = _run(tmp_path, monkeypatch, _screen_completed_clean)
    qualified = _rows(_tables(tmp_path) / "candidates_qualified.csv")

    assert len(qualified) == len(GUIDES)
    assert qualified["composite_score"].isna().all()
    assert qualified["design_score"].isna().all()
    assert qualified["target_accessibility_p"].isna().all()
    assert not qualified["scored_after_screening"].any()
    assert results["selection_summary"]["scored_after_screening"] == 0
    assert results["selection_summary"]["eligible_candidates"] == len(GUIDES)
    # The order list exists, and its headers carry no score they do not have.
    fasta = (_tables(tmp_path) / "candidates_qualified.fasta").read_text()
    assert fasta.count(">") == len(GUIDES)
    assert "score=" not in fasta


# ---------------------------------------------------------------------------------------------
# Selection reaches this entry point at all
# ---------------------------------------------------------------------------------------------


@pytest.mark.integration
def test_no_candidate_is_left_at_the_model_default(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``not_selected`` means selection never looked at the candidate, and it must not survive a run."""
    _run(tmp_path, monkeypatch, _screen_completed_clean)
    states = set(_rows(_tables(tmp_path) / "candidates_all.csv")["selection_state"])

    assert SelectionState.NOT_SELECTED.value not in states
    assert states == {SelectionState.ELIGIBLE.value}


@pytest.mark.integration
def test_the_selection_summary_reaches_the_returned_results(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``selection_summary`` is additive beside the seven keys this entry point always returned."""
    results = _run(tmp_path, monkeypatch, _screen_completed_clean)

    assert set(results) == {
        "workflow_type",
        "input_candidates",
        "candidate_count",
        "output_dir",
        "processing_time",
        "offtarget_summary",
        "reference_summary",
        "selection_summary",
        "written_files",
    }
    assert results["selection_summary"]["run_mode"] == RunMode.QUALIFIED.value
    assert results["candidate_count"] == len(GUIDES)


# ---------------------------------------------------------------------------------------------
# A screen that published nothing
# ---------------------------------------------------------------------------------------------


@pytest.mark.integration
def test_a_screen_that_published_nothing_qualifies_nobody_and_says_why(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A completed-looking screen with no evidence withholds every guide, under a named reason.

    ``candidates_qualified.csv`` is written with a header and no rows: "nothing qualified" is an
    answer, and an absent file would read as "this version does not report it".
    """
    results = _run(tmp_path, monkeypatch, _screen_published_nothing)
    base = _tables(tmp_path)
    everything = _rows(base / "candidates_all.csv")

    assert list(_rows(base / "candidates_qualified.csv")["id"]) == []
    assert not (base / "candidates_qualified.fasta").exists()
    assert set(everything["selection_state"]) == {SelectionState.WITHHELD.value}
    assert results["selection_summary"]["evidence_excluded"] == len(GUIDES)
    reasons = results["selection_summary"]["evidence_shortfall_reasons"]
    assert "no_evidence:transcriptome:human" in reasons


@pytest.mark.integration
def test_a_screen_that_published_nothing_leaves_its_gates_undecided(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#106's rule on this entry point: an unrun gate reports UNKNOWN, never a fabricated pass.

    Without the gate call this path reached no gate at all and exported ``not_evaluated`` -- the same
    cell a run with no threshold configured writes.
    """
    _run(tmp_path, monkeypatch, _screen_published_nothing)
    everything = _rows(_tables(tmp_path) / "candidates_all.csv")

    assert set(everything["max_off_target_count_verdict"]) == {FilterEvaluation.UNKNOWN.value}
    assert set(everything["max_transcriptome_hits_0mm_verdict"]) == {FilterEvaluation.UNKNOWN.value}


@pytest.mark.integration
def test_an_exploratory_run_keeps_its_guides_and_labels_them_provisional(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exploratory claims less, so it retains the guides -- named, not silently reported as qualified."""
    results = _run(tmp_path, monkeypatch, _screen_published_nothing, run_mode=RunMode.EXPLORATORY.value)
    base = _tables(tmp_path)

    assert len(_rows(base / "candidates_provisional.csv")) == len(GUIDES)
    assert list(_rows(base / "candidates_qualified.csv")["id"]) == []
    assert results["selection_summary"]["provisional_candidates"] == len(GUIDES)
    assert results["selection_summary"]["run_mode"] == RunMode.EXPLORATORY.value


# ---------------------------------------------------------------------------------------------
# The manifest, and the exit code
# ---------------------------------------------------------------------------------------------


@pytest.mark.integration
def test_the_manifest_records_the_selection_exports_and_no_orf_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """This command validates no ORF, so the manifest claims none rather than reporting one missing."""
    _run(tmp_path, monkeypatch, _screen_completed_clean)
    manifest = json.loads((_tables(tmp_path) / "manifest.json").read_text())

    files = manifest["files"]
    assert files["candidates_all_csv"]["rows"] == len(GUIDES)
    assert files["candidates_qualified_csv"]["rows"] == len(GUIDES)
    assert "orf_validation_report" not in files
    # The gates the run applied, so a report built later reads this run's thresholds.
    assert manifest["run_policy"]["run_mode"] == RunMode.QUALIFIED.value


@pytest.mark.integration
def test_the_command_exits_zero_when_it_qualified_nobody(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A qualified run that could qualify no guide is a result, and must not become an error exit.

    The design command deliberately exits non-zero on this (``_fail_if_nothing_could_qualify``); this
    command does not, and gaining a selection summary must not change that. Deciding what exit code an
    incomplete run deserves is W4's, and this test is what would notice it happening here by accident.

    The real entry point runs first, and the CLI is then driven over its actual result dictionary --
    the repo's existing pattern for CliRunner tests, which stub the workflow function rather than
    running a workflow inside the runner's captured streams.
    """
    results = _run(tmp_path, monkeypatch, _screen_published_nothing)
    assert results["selection_summary"]["eligible_candidates"] == 0
    assert results["selection_summary"]["evidence_shortfall_reasons"]

    async def _canned(**kwargs: Any) -> dict[str, Any]:
        return results

    monkeypatch.setattr("sirnaforge.cli.run_offtarget_only_workflow", _canned)
    invoked = CliRunner().invoke(
        app,
        [
            "offtarget",
            "-i",
            str(_input_fasta(tmp_path)),
            "-o",
            str(tmp_path / "cli_out"),
            "--species",
            "human",
        ],
    )

    assert invoked.exit_code == 0, invoked.output
    assert (_tables(tmp_path) / "candidates_all.csv").exists()
