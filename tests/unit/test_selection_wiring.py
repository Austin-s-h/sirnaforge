"""Selection is one decision, and the run publishes it under its own names (#100).

W2 routes every gate through :mod:`sirnaforge.core.filtering` and the eligibility decision through
:mod:`sirnaforge.core.selection`, then names the result in files instead of leaving it in a column a
reader has to know about. Three defects that survived until it did:

- ``candidates_pass.csv`` was rebuilt by re-filtering on ``passes_filters == "PASS"``, so a qualified
  run that could qualify nobody still shipped every withheld candidate as a "pass" list, and there was
  no file anywhere saying which candidates the run could actually stand behind.
- an exploratory run labelled an evidence-incomplete candidate ``eligible`` -- the same word a fully
  evidenced one gets -- so a provisional result was indistinguishable from a qualified one.
- ``_score_candidate_post_screen``'s failure path set ``scored_after_screening = False`` and left every
  post-screen field from an earlier successful call, so a row could publish "not scored after
  screening" beside a full post-screen composite.

The behaviour-preserving half of the extraction is guarded by the tests that were already green:
tests/unit/test_gate_verdict_completeness.py, test_filter_verdicts.py,
test_post_screen_scoring_integrity.py and test_issue80_regressions.py. What is asserted here is only
what W2 changes. All fixtures are synthetic; TP53 is the documented public example gene.
"""

import asyncio
import json
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.core.hit_classification import HitClassCounts
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.policy import FilterEvaluation, RunMode, ScreeningChannel
from sirnaforge.models.sirna import (
    DesignParameters,
    DesignResult,
    SelectionState,
    SiRNACandidate,
)
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

GUIDE = "ACGTACGTACGTACGTACGTA"


def _candidate(candidate_id: str, *, design_score: float = 50.0) -> SiRNACandidate:
    """A passing candidate carrying enough component scores to be re-scored after screening.

    The thermodynamic fields are populated because ``SiRNACandidateSchema`` types them float64, and an
    all-null column arrives as ``object`` -- which fails validation before step6 writes any CSV.
    """
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE.translate(str.maketrans("ATCG", "TAGC"))[::-1],
        length=len(GUIDE),
        gc_content=47.6,
        asymmetry_score=0.5,
        mfe=-4.2,
        structure="." * len(GUIDE),
        duplex_stability=-39.0,
        paired_fraction=0.1,
        design_score=design_score,
        component_scores={
            "target_accessibility": 0.5,
            "asymmetry": 0.5,
            "gc_content": 0.5,
            "duplex_stability_score": 0.6,
            "melting_temp_c": 70.0,
            "dg_5p": -8.0,
            "dg_3p": -9.0,
            "delta_dg_end": 1.0,
        },
    )


def _workflow(tmp_path: Path, name: str, *, run_mode: RunMode | None = None) -> SiRNAWorkflow:
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


def _select_and_report(workflow: SiRNAWorkflow, candidates: list[SiRNACandidate]) -> Path:
    """Run the two steps a real run runs in this order, and return the output directory."""
    result = _design_result(workflow, candidates)
    workflow._apply_post_screen_ranking(result)
    asyncio.run(workflow.step6_generate_reports(result))
    return workflow.config.output_dir / "sirnaforge"


def _rows(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


# ---------------------------------------------------------------------------------------------
# The qualified export: named, and empty when nothing qualified
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_a_qualified_run_with_no_evidence_publishes_an_empty_qualified_export(tmp_path: Path) -> None:
    """#100: a run that could qualify nobody must publish that answer, not a full "pass" list.

    Every candidate passes its gates, so ``candidates_pass.csv`` keeps all three rows -- that file is
    the compat deliverable and narrowing it once deleted the whole output of an unscreened run. The
    answer to "which of these can I order?" therefore has to live somewhere else, and before W2 it
    lived nowhere: no ``candidates_qualified.csv`` was written at all.
    """
    workflow = _workflow(tmp_path, "qualified_nothing", run_mode=RunMode.QUALIFIED)
    assert ("transcriptome", "human") in workflow.config.resolved_policy.evidence_requirements.required_pairs
    candidates = [_candidate(f"cand_{i}") for i in range(3)]
    for candidate in candidates:
        candidate.off_target_screened = False  # nothing was screened, so nothing is evidenced

    base = _select_and_report(workflow, candidates)

    assert len(_rows(base / "candidates_all.csv")) == 3
    assert len(_rows(base / "candidates_pass.csv")) == 3, "the legacy pass predicate must not move"
    qualified = base / "candidates_qualified.csv"
    assert qualified.exists(), "a qualified run must publish its qualified set even when it is empty"
    assert len(_rows(qualified)) == 0
    # No order list either: an empty FASTA is not one, and a stale file would be worse.
    assert not (base / "candidates_qualified.fasta").exists()
    assert set(_rows(base / "candidates_all.csv")["selection_state"]) == {SelectionState.WITHHELD.value}


@pytest.mark.unit
def test_the_qualified_export_holds_the_eligible_rows_in_selection_order(tmp_path: Path) -> None:
    """The qualified file is a slice of the resolved selection, not a re-filter of something else.

    Order matters: it is the ranking selection produced, which is the order ``candidates_all.csv`` is
    written in. A file rebuilt from an independent predicate could disagree with ``top_candidates``,
    which is how a withheld candidate came to lead a "pass" list.
    """
    workflow = _workflow(tmp_path, "qualified_order", run_mode=RunMode.QUALIFIED)
    strong = _candidate("strong", design_score=90.0)
    weak = _candidate("weak", design_score=10.0)
    rejected = _candidate("rejected", design_score=99.0)
    rejected.passes_filters = SiRNACandidate.FilterStatus.GC_OUT_OF_RANGE
    for candidate in (strong, weak, rejected):
        candidate.off_target_screened = True

    base = _select_and_report(workflow, [weak, rejected, strong])

    qualified = _rows(base / "candidates_qualified.csv")
    assert list(qualified["id"]) == ["strong", "weak"]
    assert set(qualified["selection_state"]) == {SelectionState.ELIGIBLE.value}
    # The gate-rejected row is published, and says why, in candidates_all.csv only.
    all_rows = _rows(base / "candidates_all.csv")
    assert len(all_rows) == 3
    assert (base / "candidates_qualified.fasta").read_text().count(">") == 2


# ---------------------------------------------------------------------------------------------
# The provisional export: an exploratory result that cannot be mistaken for a qualified one
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_an_exploratory_shortfall_is_provisional_and_gets_its_own_export(tmp_path: Path) -> None:
    """#100: exploratory retains incomplete evidence, and must say that it is what it retained.

    Before W2 these rows read ``eligible`` -- byte-identical to a fully evidenced candidate from a
    qualified run -- so the scope of the claim survived nowhere on the row or in the run's outputs.
    They still rank and still enter ``top_candidates``: the state is a label, not an exclusion.
    """
    workflow = _workflow(tmp_path, "exploratory_provisional", run_mode=RunMode.EXPLORATORY)
    requirements = workflow.config.resolved_policy.evidence_requirements
    assert requirements.required_pairs == frozenset(), "nothing is required in exploratory mode"
    assert any(entry.channel is ScreeningChannel.TRANSCRIPTOME for entry in requirements.channel_requirements)

    candidates = [_candidate(f"cand_{i}") for i in range(2)]
    for candidate in candidates:
        candidate.off_target_screened = False

    result = _design_result(workflow, candidates)
    workflow._apply_post_screen_ranking(result)
    assert [c.id for c in result.top_candidates] == [c.id for c in candidates], "provisional still ranks"
    assert workflow._selection_summary["provisional_candidates"] == 2
    assert workflow._selection_summary["eligible_candidates"] == 0

    asyncio.run(workflow.step6_generate_reports(result))
    base = workflow.config.output_dir / "sirnaforge"
    all_rows = _rows(base / "candidates_all.csv")
    assert set(all_rows["selection_state"]) == {SelectionState.PROVISIONAL.value}
    provisional = base / "candidates_provisional.csv"
    assert provisional.exists()
    assert len(_rows(provisional)) == 2
    # Nothing is qualified, so the qualified export is empty rather than absent or full.
    assert len(_rows(base / "candidates_qualified.csv")) == 0
    # The order list names the scope on every header, so a sequence cannot leave the tool unlabelled.
    assert (base / "candidates_pass.fasta").read_text().count(f"selection={SelectionState.PROVISIONAL.value}") == 2


@pytest.mark.unit
def test_a_fully_evidenced_exploratory_candidate_is_eligible_not_provisional(tmp_path: Path) -> None:
    """The split has to cut: with the evidence present, exploratory returns the same answer as before.

    Otherwise ``provisional`` would just be exploratory mode's new name for everything, and the
    distinction it exists to draw would carry no information.
    """
    workflow = _workflow(tmp_path, "exploratory_complete", run_mode=RunMode.EXPLORATORY)
    candidate = _candidate("evidenced")
    workflow._integrate_offtarget_results(
        [candidate], {"status": "completed", "results": {}}, screened_species=["human"], mirna_screened=True
    )
    assert candidate.off_target_screened is True

    base = _select_and_report(workflow, [candidate])

    assert _rows(base / "candidates_all.csv")["selection_state"].tolist() == [SelectionState.ELIGIBLE.value]
    assert len(_rows(base / "candidates_qualified.csv")) == 1
    assert len(_rows(base / "candidates_provisional.csv")) == 0


# ---------------------------------------------------------------------------------------------
# A failed re-score cannot leave the score it replaced
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_a_failed_rescore_clears_the_post_screen_score_it_had_written(tmp_path: Path) -> None:
    """#100: ``scored_after_screening=False`` beside a post-screen composite is a contradiction.

    Selection demotes an unscored candidate because its number is not comparable with its neighbours';
    ``ranking_score`` nonetheless ranks on ``composite_score`` whenever one is present. So a retry that
    failed after an earlier success used to leave exactly the number the run had just disowned -- and
    every per-term contribution behind it -- on the row and in the ranking.
    """
    workflow = _workflow(tmp_path, "rescore", run_mode=RunMode.QUALIFIED)
    candidate = _candidate("rescored")
    counts = HitClassCounts(off_target=1)

    assert workflow._score_candidate_post_screen(candidate, counts, frozenset()) is True
    assert candidate.composite_score is not None
    assert candidate.score_off_target is not None

    # An out-of-range sub-score is the one thing compute_composite refuses outright.
    candidate.component_scores["asymmetry"] = 1.5
    assert workflow._score_candidate_post_screen(candidate, counts, frozenset()) is False

    assert candidate.scored_after_screening is False
    assert candidate.composite_score is None, "a disowned score must not survive on the row"
    assert candidate.score_off_target is None
    assert candidate.weight_vector == ""


@pytest.mark.unit
def test_a_first_scoring_failure_keeps_the_score_the_candidate_arrived_with(tmp_path: Path) -> None:
    """The clearing is scoped to a retry, because a design-time score is not this method's to delete.

    A candidate whose very first post-screen attempt fails keeps the design-stage number it arrived
    with: that is what the mixed-scale ranking rule reads, and #80 F10 pins that it is reported.
    """
    workflow = _workflow(tmp_path, "first_failure", run_mode=RunMode.QUALIFIED)
    candidate = _candidate("never_scored")
    candidate.composite_score = 99.0  # a design-stage composite, as MiRNADesigner writes
    candidate.component_scores["asymmetry"] = 1.5

    assert workflow._score_candidate_post_screen(candidate, HitClassCounts(), frozenset()) is False

    assert candidate.scored_after_screening is False
    assert candidate.composite_score == 99.0


# ---------------------------------------------------------------------------------------------
# The extraction itself: one evaluator, one selection decision
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_every_off_target_gate_is_evaluated_by_the_shared_evaluator(tmp_path: Path) -> None:
    """The gates read their rules from one place, and every one of them still records an outcome.

    A regression here means either the spec assembly dropped a gate (its column would go missing) or
    the shared evaluator stopped agreeing with what ``_check_offtarget_filters`` used to do inline.
    ``max_transcriptome_hits_0mm`` has a default ceiling of 1, so 5 perfect hits is a decided FAIL,
    while ``max_transcriptome_seed_perfect`` ships with no threshold and must stay ``not_evaluated``.
    """
    workflow = _workflow(tmp_path, "gates", run_mode=RunMode.QUALIFIED)
    candidate = _candidate("gated")
    criteria = workflow.config.design_params.offtarget_filters

    should_fail, status = workflow._check_offtarget_filters(
        5, 0, 0, 0, 0, 0, 5, 5, criteria, candidate, complete_pairs=frozenset({("transcriptome", "human")})
    )

    assert should_fail is True
    assert status is SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH
    assert candidate.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.FAIL.value
    assert candidate.filter_observed["max_transcriptome_hits_0mm"] == 5
    assert candidate.filter_verdicts["max_transcriptome_hits_1mm"] == FilterEvaluation.PASS.value
    assert candidate.filter_verdicts["max_transcriptome_seed_perfect"] == FilterEvaluation.NOT_EVALUATED.value
    # Every gate that read a channel whose evidence never completed is undecided, not a pass at zero.
    undecided = workflow._check_offtarget_filters(
        0, 0, 0, 0, 0, 0, 0, 0, criteria, _candidate("unscreened"), complete_pairs=frozenset()
    )
    assert undecided == (False, None)


@pytest.mark.unit
# 0.65 is the declared default floor, so these are above it, exactly on it, and just below it.
@pytest.mark.parametrize(("asymmetry", "expected"), [(0.9, True), (0.65, True), (0.6499, False)])
def test_the_asymmetry_gate_still_agrees_with_the_thermodynamic_helper(asymmetry: float, expected: bool) -> None:
    """The design gate compares through the shared evaluator now, so pin it against the old authority.

    ``_apply_score_filters`` used to call ``ThermodynamicCalculator.meets_asymmetry_threshold``
    directly. It now declares ``FilterComparator.GE`` instead, which is the same comparison -- and this
    is what would catch the two drifting apart, boundary included.
    """
    designer = SiRNADesigner(DesignParameters())
    threshold = designer.parameters.filters.min_asymmetry_score
    candidate = _candidate("asymmetry_probe")

    designer._apply_score_filters(candidate, asymmetry, empirical_score=0.5)

    verdict = candidate.filter_verdicts["min_asymmetry_score"]
    assert (verdict == FilterEvaluation.PASS.value) is expected
    assert ThermodynamicCalculator.meets_asymmetry_threshold(asymmetry, threshold) is expected


@pytest.mark.unit
def test_selection_and_the_exports_agree_on_one_answer(tmp_path: Path) -> None:
    """``top_candidates``, ``selection_state`` and the qualified export are one decision, not three.

    They were derived separately: ``top_candidates`` from the ranking loop, the state from a set of
    ``id()`` values, and the pass exports from a dataframe predicate. Anything that re-derives one of
    them independently can disagree with the other two, and a reader has no way to tell which is right.
    """
    workflow = _workflow(tmp_path, "agreement", run_mode=RunMode.QUALIFIED)
    evidenced = _candidate("evidenced")
    unevidenced = _candidate("unevidenced")
    evidenced.off_target_screened = True
    unevidenced.off_target_screened = False

    result = _design_result(workflow, [evidenced, unevidenced])
    workflow._apply_post_screen_ranking(result)
    asyncio.run(workflow.step6_generate_reports(result))
    base = workflow.config.output_dir / "sirnaforge"

    assert [c.id for c in result.top_candidates] == ["evidenced"]
    assert workflow._eligible_candidate_ids == frozenset({"evidenced"})
    assert list(_rows(base / "candidates_qualified.csv")["id"]) == ["evidenced"]
    all_rows = _rows(base / "candidates_all.csv")
    states = dict(zip(all_rows["id"], all_rows["selection_state"], strict=True))
    assert states == {"evidenced": SelectionState.ELIGIBLE.value, "unevidenced": SelectionState.WITHHELD.value}


@pytest.mark.unit
def test_the_selection_exports_are_recorded_in_the_fair_manifest(tmp_path: Path) -> None:
    """A published file absent from the manifest is a file no downstream consumer can find.

    The row counts are in the manifest too, which is what makes "0 qualified" a recorded outcome of the
    run rather than something a reader has to open the CSV to discover.
    """
    workflow = _workflow(tmp_path, "manifest", run_mode=RunMode.QUALIFIED)
    candidate = _candidate("evidenced")
    candidate.off_target_screened = True

    base = _select_and_report(workflow, [candidate])
    manifest = json.loads((base / "manifest.json").read_text())

    files = manifest["files"]
    assert files["candidates_qualified_csv"]["rows"] == 1
    assert files["candidates_provisional_csv"]["rows"] == 0
    assert files["candidates_qualified_fasta"]["sequences"] == 1
