"""Every gate records its own verdict, and a zero only means clean when it was measured.

Four defects found reviewing PR #104 and reproduced against its tip, plus the #100 eligibility rule
they share a root with. What links them is a single question: when a run publishes a verdict, can a
consumer tell whether the gate looked, and at what?

- #105 ``min_isoform_coverage`` was the one gate that assigned its rejection label by hand instead of
  going through the verdict recorder, so a run resolving it to ``warn`` still threw the candidate out
  and recorded nothing about having done so.
- #106 the clean-screen path skipped the off-target gates entirely, so a completed screen that found
  nothing exported every one of them as ``not_evaluated`` with a blank observed value -- which
  #103's report reads as *unknown* rather than as the measured pass it is.
- #108 the normal aggregate ingest never recorded its own path, so ``detail_files.mirna`` was empty
  beside non-zero miRNA counters, and the summary points at those files instead of carrying them.
- #100 selection excluded unscored candidates only on a *mixed* batch, so a wholly unscreened batch
  went out whole -- in a run mode whose entire claim is completed evidence.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

import asyncio
import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.models.policy import (
    FilterAction,
    FilterEvaluation,
    FilterStage,
    RunMode,
    ScreeningChannel,
)
from sirnaforge.models.sirna import (
    DesignParameters,
    DesignResult,
    FilterCriteria,
    OffTargetFilterCriteria,
    SiRNACandidate,
    build_candidate_row,
)
from sirnaforge.reporting.payload import REASON_OK, observed_column
from sirnaforge.reporting.payload import _evaluate as evaluate_descriptor
from sirnaforge.workflow import POST_SCREEN_FILTER_CHANNELS, SiRNAWorkflow, WorkflowConfig

GUIDE = "ACGTACGTACGTACGTACGTA"
NO_HITS: dict[str, object] = {"status": "completed", "results": {}}

#: Verdict codes as :mod:`sirnaforge.reporting.payload` emits them, so a test can say which one it
#: means without restating the mapping.
REPORT_PASS = 0
REPORT_UNKNOWN = 2

MIRNA_GATES = tuple(
    filter_id
    for filter_id, channels in POST_SCREEN_FILTER_CHANNELS.items()
    if channels == frozenset({ScreeningChannel.MIRNA_SEED})
)
TRANSCRIPTOME_GATES = tuple(
    filter_id
    for filter_id, channels in POST_SCREEN_FILTER_CHANNELS.items()
    if channels == frozenset({ScreeningChannel.TRANSCRIPTOME})
)


def _candidate(candidate_id: str = "probe", *, isoform_coverage: float | None = None) -> SiRNACandidate:
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
        isoform_coverage=isoform_coverage,
        design_score=50.0,
        component_scores={"target_accessibility": 0.5, "asymmetry": 0.5, "gc_content": 0.5},
    )


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    stated: dict[str, object] | None = None,
    filter_actions: dict[str, str] | None = None,
    run_mode: RunMode | None = None,
    entry_point: EntryPoint = EntryPoint.SCREENING_WORKFLOW,
    screen_species: list[str] | None = None,
) -> SiRNAWorkflow:
    """A workflow whose policy is resolved through the one resolver, as every real entry point does."""
    policy = resolve_run_policy(
        entry_point=entry_point,
        run_mode=run_mode,
        stated=stated,
        filter_actions=filter_actions,
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


# ---------------------------------------------------------------------------------------------
# #105: the isoform coverage gate obeys its resolved action
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    ("action", "expect_rejected"),
    [("warn", False), ("fail", True)],
)
def test_isoform_coverage_obeys_its_resolved_action(tmp_path: Path, action: str, expect_rejected: bool) -> None:
    """#105: a coverage floor resolved to ``warn`` records the failure without rejecting.

    The pre-fix gate assigned ``passes_filters = LOW_ISOFORM_COVERAGE`` directly, so the action was
    unreachable: warn and fail behaved identically and ``filter_verdicts`` stayed empty either way.
    Both halves are asserted, because recording the verdict without honouring the action, or
    honouring the action without recording the verdict, each looks like a fix on its own.
    """
    workflow = _workflow(
        tmp_path,
        f"iso_{action}",
        stated={"min_isoform_coverage": 0.8},
        filter_actions={"min_isoform_coverage": action},
    )
    assert workflow.config.resolved_policy.descriptor("min_isoform_coverage").action.value == action

    candidate = _candidate(isoform_coverage=0.2)
    rejected = workflow._apply_isoform_coverage_gate(candidate)

    assert rejected is expect_rejected
    # The verdict is the same either way -- the observation does not depend on what is done with it.
    assert candidate.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.FAIL.value
    assert candidate.filter_observed["min_isoform_coverage"] == pytest.approx(0.2)
    if expect_rejected:
        assert candidate.passes_filters == SiRNACandidate.FilterStatus.LOW_ISOFORM_COVERAGE
    else:
        assert candidate.passes_filters is True, "warn records the finding; it does not reject"


@pytest.mark.unit
def test_isoform_coverage_records_a_pass_and_an_unknown(tmp_path: Path) -> None:
    """#105: the two non-rejecting outcomes are recorded rather than inferred from silence.

    Coverage at or above the floor is a PASS a report can show. Coverage that could not be computed
    is UNKNOWN with no observed value: an annotation gap is not evidence of poor coverage, and it is
    also not evidence of good coverage, which is what a missing verdict would have implied.
    """
    workflow = _workflow(tmp_path, "iso_pass", stated={"min_isoform_coverage": 0.5})

    passing = _candidate("passing", isoform_coverage=0.9)
    assert workflow._apply_isoform_coverage_gate(passing) is False
    assert passing.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.PASS.value
    assert passing.filter_observed["min_isoform_coverage"] == pytest.approx(0.9)
    assert passing.passes_filters is True

    unknown = _candidate("unknown", isoform_coverage=None)
    assert workflow._apply_isoform_coverage_gate(unknown) is False
    assert unknown.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.UNKNOWN.value
    assert unknown.filter_observed["min_isoform_coverage"] is None
    assert unknown.passes_filters is True, "a missing annotation must not be read as a low coverage"


@pytest.mark.unit
def test_isoform_coverage_off_is_not_evaluated_but_still_reports_its_value(tmp_path: Path) -> None:
    """#105: with no floor configured the gate writes a word, not an empty cell.

    The exported verdict column is fixed per filter, so a blank there reads as a verdict. The
    measured coverage is still published -- the gate did not run, but the number exists.
    """
    workflow = _workflow(tmp_path, "iso_off")
    assert workflow.config.design_params.filters.min_isoform_coverage is None

    candidate = _candidate(isoform_coverage=0.1)
    assert workflow._apply_isoform_coverage_gate(candidate) is False
    assert candidate.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.NOT_EVALUATED.value
    assert candidate.filter_observed["min_isoform_coverage"] == pytest.approx(0.1)
    assert candidate.passes_filters is True


# ---------------------------------------------------------------------------------------------
# #106: a completed zero-hit screen is a measurement
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_a_completed_zero_hit_screen_records_measured_passes(tmp_path: Path) -> None:
    """#106: a clean screen must record PASS with an observed zero, not ``not_evaluated``.

    Pre-fix the no-hit branch scored the candidate and continued, reaching no off-target gate at all.
    The row then exported ``max_transcriptome_hits_0mm_verdict=not_evaluated`` with a blank observed
    value, and #103's report -- which re-derives each verdict from the observed column -- read a
    measured clean screen as *unknown*. The report evaluator is asserted here alongside the row, so
    the two cannot drift apart again: agreeing with the pipeline is the point of the column.
    """
    workflow = _workflow(tmp_path, "zero_hit")
    candidate = _candidate()

    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=["human"])

    assert candidate.off_target_screened is True
    assert candidate.scored_after_screening is True

    row = build_candidate_row(candidate)
    populated = {key for key, value in row.items() if value is not None}
    for filter_id in TRANSCRIPTOME_GATES + MIRNA_GATES:
        descriptor = workflow.config.resolved_policy.descriptor(filter_id)
        if descriptor.action is FilterAction.OFF or descriptor.threshold is None:
            continue
        assert candidate.filter_verdicts[filter_id] == FilterEvaluation.PASS.value, filter_id
        assert candidate.filter_observed[filter_id] == 0, filter_id
        assert row[f"{filter_id}_verdict"] == FilterEvaluation.PASS.value, filter_id
        value, verdict_code, reason = evaluate_descriptor(
            descriptor, pd.Series(row), observed_column(descriptor, populated)
        )
        assert (value, verdict_code, reason) == (0, REPORT_PASS, REASON_OK), filter_id


@pytest.mark.unit
def test_a_missing_mirna_channel_is_unknown_while_the_alignment_still_passes(tmp_path: Path) -> None:
    """#106/#100: an absent channel yields UNKNOWN; a completed one is unaffected by it.

    The two channels are independent evidence, so the miRNA scan not running cannot turn the
    alignment's measured zero into an unknown, and it must not let the miRNA gates pass on a count
    of zero nobody looked for. The observed value stays empty deliberately: writing the lower bound
    there would make the report re-derive a confident PASS and contradict the run.
    """
    workflow = _workflow(tmp_path, "no_mirna")
    candidate = _candidate()

    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=["human"], mirna_screened=False)

    for filter_id in MIRNA_GATES:
        descriptor = workflow.config.resolved_policy.descriptor(filter_id)
        if descriptor.action is FilterAction.OFF or descriptor.threshold is None:
            continue
        assert candidate.filter_verdicts[filter_id] == FilterEvaluation.UNKNOWN.value, filter_id
        assert candidate.filter_observed[filter_id] is None, filter_id
        row = build_candidate_row(candidate)
        populated = {key for key, value in row.items() if value is not None}
        _, verdict_code, _ = evaluate_descriptor(descriptor, pd.Series(row), observed_column(descriptor, populated))
        assert verdict_code == REPORT_UNKNOWN, filter_id

    assert candidate.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.PASS.value
    # Two channels, so this one is decided by neither alone.
    assert candidate.filter_verdicts["max_total_offtarget_hits"] in {
        FilterEvaluation.UNKNOWN.value,
        FilterEvaluation.NOT_EVALUATED.value,
    }


@pytest.mark.unit
def test_an_incomplete_channel_still_fails_a_count_over_its_ceiling(tmp_path: Path) -> None:
    """#100: a lower bound above the threshold is proof enough, so a known failure still fails.

    Only the *pass* direction needs a completed channel. Hits that were found are real evidence
    whether or not the search finished, so unknown must not become a get-out for a candidate whose
    measured count already breaks the ceiling.
    """
    workflow = _workflow(tmp_path, "lower_bound")
    candidate = _candidate()

    should_fail, status = workflow._check_offtarget_filters(
        5,  # transcriptome_0mm, over its ceiling of 1
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        OffTargetFilterCriteria(),
        candidate,
        complete_channels=frozenset(),  # nothing completed
    )

    assert should_fail is True
    assert status == SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH
    assert candidate.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.FAIL.value
    assert candidate.filter_observed["max_transcriptome_hits_0mm"] == 5
    # The gates that measured zero on the same incomplete evidence cannot claim a pass.
    assert candidate.filter_verdicts["max_transcriptome_hits_1mm"] == FilterEvaluation.UNKNOWN.value


@pytest.mark.unit
@pytest.mark.parametrize(
    ("summary", "expected"),
    [
        (None, False),
        ({}, False),
        ({"total_candidates": 0, "analysis_files_processed": 1}, False),
        ({"total_candidates": 1, "analysis_files_processed": 1}, True),
        ({"analysis_files_processed": 2}, True),
    ],
)
def test_mirna_channel_completion_prefers_parsed_files_over_globbed_ones(
    summary: dict[str, int] | None, expected: bool
) -> None:
    """#106: ``total_candidates`` is the evidence, because it counts files that actually parsed.

    ``analysis_files_processed`` counts what the aggregator globbed, including a 0-byte file it then
    failed to read -- the same weak signal that once let a 0-byte ``*_analysis.tsv`` publish a
    completed transcriptome screen. The third case is exactly that shape and must read as incomplete.
    """
    assert SiRNAWorkflow._mirna_channel_completed(summary) is expected


# ---------------------------------------------------------------------------------------------
# #108: the miRNA detail pointer names the file that was actually read
# ---------------------------------------------------------------------------------------------

_MIRNA_COLUMNS = ["qname", "mirna_id", "database", "species", "seed_mismatches", "offtarget_score"]


def _mirna_rows(qname: str = "probe") -> list[dict[str, str]]:
    return [
        {
            "qname": qname,
            "mirna_id": "hsa-miR-0000",
            "database": "mirgenedb",
            "species": "human",
            "seed_mismatches": "0",
            "offtarget_score": "1.0",
        }
    ]


def _write_mirna_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_MIRNA_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.unit
def test_the_mirna_aggregate_is_published_as_a_detail_pointer(tmp_path: Path) -> None:
    """#108: the normal successful path records the aggregate it ingested.

    Pre-fix only the fallback appended a path, so a run that read ``combined_mirna_hits.tsv``
    published ``detail_files.mirna: []`` beside non-zero miRNA counters. Since the summary points at
    these files rather than carrying their rows, that empty list was the only route a consumer had
    to the hits it was being told about.
    """
    workflow = _workflow(tmp_path, "mirna_pointer")
    results_dir = tmp_path / "mirna_pointer" / "off_target" / "results"
    aggregate = results_dir / "aggregated" / "combined_mirna_hits.tsv"
    _write_mirna_tsv(aggregate, _mirna_rows())

    parsed = asyncio.run(workflow._parse_nextflow_results(results_dir))

    assert parsed["mirna_hit_files"] == [aggregate]
    detail = SiRNAWorkflow._offtarget_detail_files(parsed, results_dir)
    assert detail["mirna"] == [str(aggregate)]


@pytest.mark.unit
def test_the_mirna_json_fallback_is_published_too(tmp_path: Path) -> None:
    """#108: the JSON aggregate is a route to the same rows, so it is named on the same footing."""
    workflow = _workflow(tmp_path, "mirna_json")
    results_dir = tmp_path / "mirna_json" / "off_target" / "results"
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True)
    json_path = aggregated / "combined_mirna_hits.json"
    json_path.write_text(json.dumps(_mirna_rows()))

    parsed = asyncio.run(workflow._parse_nextflow_results(results_dir))

    assert parsed["mirna_hit_files"] == [json_path]


@pytest.mark.unit
def test_a_per_file_fallback_pointer_is_recorded_exactly_once(tmp_path: Path) -> None:
    """#108: the fallback still publishes, and no path is named twice.

    The pointer list is a set of locations. Naming one file twice would read as two independent
    sources of evidence for the same hits, which is how the miRNA counters came to be exactly 2x.
    """
    workflow = _workflow(tmp_path, "mirna_fallback")
    results_dir = tmp_path / "mirna_fallback" / "off_target" / "results"
    per_file = results_dir / "mirna" / "mirna_analysis.tsv"
    _write_mirna_tsv(per_file, _mirna_rows())
    # An empty aggregate is what sends the parser to the fallback in the first place.
    (results_dir / "aggregated").mkdir(parents=True, exist_ok=True)
    (results_dir / "aggregated" / "combined_mirna_hits.tsv").write_text("")

    parsed = asyncio.run(workflow._parse_nextflow_results(results_dir))

    assert parsed["mirna_hit_files"] == [per_file]
    assert len(parsed["mirna_hit_files"]) == len(set(parsed["mirna_hit_files"]))


# ---------------------------------------------------------------------------------------------
# #100: run-mode eligibility, not batch mixedness
# ---------------------------------------------------------------------------------------------


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


@pytest.mark.unit
def test_a_wholly_unscreened_batch_cannot_qualify(tmp_path: Path) -> None:
    """#100: the reproduced defect -- an entirely unscreened passing batch entered top_candidates.

    The exclusion used to require a *mixed* batch, and a wholly unscored list was called internally
    consistent and shipped whole. That is precisely the shape a failed or absent screen produces, so
    the one case the rule had to catch was the one it waved through. The shortfall reason is asserted
    too: an empty shortlist that does not say why is indistinguishable from a run that found nothing.
    """
    workflow = _workflow(tmp_path, "qualified", run_mode=RunMode.QUALIFIED)
    candidate = _candidate()

    # screened_species=[] is a run with no alignment evidence for the query species at all.
    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=[])
    assert candidate.off_target_screened is False

    result = _design_result(workflow, [candidate])
    workflow._apply_post_screen_ranking(result)

    assert result.top_candidates == []
    summary = workflow._selection_summary
    assert summary["run_mode"] == RunMode.QUALIFIED.value
    assert summary["evidence_excluded"] == 1
    reasons = summary["evidence_shortfall_reasons"]
    assert reasons["no_evidence:transcriptome:human"] == 1
    # The same absent channel also leaves its gates undecided, and those are reported by name too --
    # but only gates reading the required channel, never one reading an exploratory one.
    assert set(reasons) - {"no_evidence:transcriptome:human"} <= {
        f"unknown:{filter_id}" for filter_id in TRANSCRIPTOME_GATES + ("max_total_offtarget_hits",)
    }


@pytest.mark.unit
@pytest.mark.parametrize("run_mode", [RunMode.EXPLORATORY, RunMode.DESIGN_ONLY])
def test_only_qualified_mode_withholds_the_shortlist(tmp_path: Path, run_mode: RunMode) -> None:
    """#100: exploratory retains incomplete evidence, and a design-only run keeps its design answer.

    The rule is about what the run *claimed*, so the two modes that claim less are unaffected. Were
    this not asserted, tightening qualified mode could silently delete the shortlist from
    ``--skip-off-targets``, which never promised screening evidence in the first place.
    """
    workflow = _workflow(tmp_path, f"mode_{run_mode.value}", run_mode=run_mode)
    candidate = _candidate()
    candidate.off_target_screened = False

    result = _design_result(workflow, [candidate])
    workflow._apply_post_screen_ranking(result)

    assert result.top_candidates == [candidate]
    assert workflow._selection_summary["evidence_excluded"] == 0


@pytest.mark.unit
def test_a_qualified_candidate_with_complete_evidence_still_qualifies(tmp_path: Path) -> None:
    """#100: the rule must not empty the shortlist on the runs it is meant to certify.

    A screen that completed and found nothing is the common case, and it has to come out the other
    side eligible -- otherwise the fix for the unscreened batch would have replaced one wrong answer
    with another.
    """
    workflow = _workflow(tmp_path, "qualified_clean", run_mode=RunMode.QUALIFIED)
    candidate = _candidate()

    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=["human"])
    result = _design_result(workflow, [candidate])
    workflow._apply_post_screen_ranking(result)

    assert result.top_candidates == [candidate]
    assert workflow._selection_summary["evidence_excluded"] == 0
    assert workflow._selection_summary["eligible_candidates"] == 1


@pytest.mark.unit
def test_an_unknown_required_coverage_blocks_qualification(tmp_path: Path) -> None:
    """#100: unknown required coverage blocks qualification even with no annotation to blame.

    ``min_isoform_coverage`` reads annotation rather than a screening channel, so an incomplete
    channel never excuses it: if the gate is in force and its input is missing, the candidate's
    coverage is unknown and a qualified run cannot certify it. The same candidate qualifies once the
    coverage is known and above the floor, which is what shows the block is the unknown and not the
    gate merely being configured.
    """
    workflow = _workflow(
        tmp_path,
        "unknown_coverage",
        stated={"min_isoform_coverage": 0.5},
        run_mode=RunMode.QUALIFIED,
    )
    assert workflow.config.resolved_policy.evidence_requirements.unknown_evidence_action.value == "fail"

    unknown = _candidate("unknown_cov", isoform_coverage=None)
    known = _candidate("known_cov", isoform_coverage=0.9)
    for candidate in (unknown, known):
        candidate.off_target_screened = True
        workflow._apply_isoform_coverage_gate(candidate)

    result = _design_result(workflow, [unknown, known])
    workflow._apply_post_screen_ranking(result)

    assert unknown not in result.top_candidates
    assert known in result.top_candidates
    assert workflow._selection_summary["evidence_shortfall_reasons"] == {"unknown:min_isoform_coverage": 1}


@pytest.mark.unit
def test_an_exploratory_channels_unknown_does_not_disqualify(tmp_path: Path) -> None:
    """#101/#100: dropping an exploratory channel cannot change a required gate's effect.

    Every miRNA pair is exploratory in 0.7.1, so a run whose miRNA aggregate is simply absent must
    not disqualify every candidate it screened. Without this scope the mirror defect appears: an
    exploratory channel silently deciding a qualified shortlist.
    """
    workflow = _workflow(tmp_path, "exploratory_unknown", run_mode=RunMode.QUALIFIED)
    candidate = _candidate()

    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=["human"], mirna_screened=False)
    assert candidate.filter_verdicts["max_mirna_perfect_seed"] == FilterEvaluation.UNKNOWN.value

    result = _design_result(workflow, [candidate])
    workflow._apply_post_screen_ranking(result)

    assert result.top_candidates == [candidate]
    assert workflow._selection_summary["evidence_excluded"] == 0


@pytest.mark.unit
def test_a_failed_screen_reselects_on_the_way_out(tmp_path: Path) -> None:
    """#100: final selection runs on every step5 exit, including the exception path.

    Without the re-selection the shortlist stays as the pre-screen call left it, so a qualified run
    whose screen threw published a shortlist built from design evidence alone -- the exact claim the
    mode exists to refuse. Asserted through step5 rather than the private method, because the defect
    was in which paths called it.
    """
    workflow = _workflow(tmp_path, "failed_screen", run_mode=RunMode.QUALIFIED)
    candidate = _candidate()
    candidate.passes_filters = True
    result = _design_result(workflow, [candidate])

    async def _explode(*_args: object, **_kwargs: object) -> dict[str, object]:
        raise RuntimeError("nextflow unavailable")

    workflow._run_repeat_detection = lambda *_args, **_kwargs: {"status": "skipped"}  # type: ignore[method-assign]
    workflow._run_nextflow_offtarget_analysis = _explode  # type: ignore[method-assign]

    outcome = asyncio.run(workflow.step5_offtarget_analysis(result))

    assert outcome["status"] == "skipped"
    assert outcome["reason"] == "nextflow_failed"
    assert result.top_candidates == [], "a screen that never produced evidence cannot leave a qualified shortlist"
    assert workflow._selection_summary["evidence_excluded"] == 1


@pytest.mark.unit
def test_a_refusal_before_screening_still_reselects_and_still_refuses(tmp_path: Path) -> None:
    """#100: an exception that propagates out of step5 also gets a final selection.

    ``_prepare_offtarget_input`` refuses duplicate candidate ids outright, because fanning one
    candidate's results onto another's sequence is worse than stopping. That refusal must stay a
    refusal -- not be caught and downgraded to "skipped" -- while the shortlist is still brought into
    line with the evidence on the way out.
    """
    workflow = _workflow(tmp_path, "duplicate_ids", run_mode=RunMode.QUALIFIED)
    first, second = _candidate("clash"), _candidate("clash")
    result = _design_result(workflow, [first, second])

    with pytest.raises(ValueError, match="Duplicate candidate id"):
        asyncio.run(workflow.step5_offtarget_analysis(result))

    assert result.top_candidates == []
    assert workflow._selection_summary["evidence_excluded"] == 2


@pytest.mark.unit
def test_every_declared_gate_names_its_channel_or_none(tmp_path: Path) -> None:
    """One channel map, and eligibility depends on it, so an unmapped post-screen gate is a defect.

    A gate added to ``_check_offtarget_filters`` without an entry here would raise; a *declared*
    post-screen gate missing from it would silently be treated as reading no channel, and so as
    always in scope for the unknown-blocks-qualification rule. Pinned against the resolver's own
    registry rather than a hand-written list.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human")
    post_screen = {
        resolved.filter_id for resolved in policy.filters if resolved.descriptor.stage is FilterStage.POST_SCREEN
    }
    # isoform coverage is post-screen but reads annotation, not a screening channel: it is
    # deliberately absent, and its absence is what puts it always in scope.
    assert post_screen - set(POST_SCREEN_FILTER_CHANNELS) == {"min_isoform_coverage"}
    assert set(POST_SCREEN_FILTER_CHANNELS) <= post_screen


@pytest.mark.unit
def test_the_run_summary_publishes_why_the_shortlist_is_empty(tmp_path: Path) -> None:
    """#100: success, no-eligible-candidates and incomplete must not read alike.

    All three leave ``top_candidates`` empty, so the distinguishing information has to be published
    rather than inferred. These counters are what a caller reads to tell a clean run that found
    nothing good from a run that never produced the evidence to judge.
    """
    workflow = _workflow(tmp_path, "summary", run_mode=RunMode.QUALIFIED)
    rejected = _candidate("rejected")
    rejected.passes_filters = SiRNACandidate.FilterStatus.GC_OUT_OF_RANGE
    rejected.off_target_screened = True

    result = _design_result(workflow, [rejected])
    workflow._apply_post_screen_ranking(result)

    assert result.top_candidates == []
    assert workflow._selection_summary["filter_excluded"] == 1
    assert workflow._selection_summary["evidence_excluded"] == 0
    assert workflow._selection_summary["evidence_shortfall_reasons"] == {}


@pytest.mark.unit
def test_isoform_coverage_gate_action_reaches_the_integration_path(tmp_path: Path) -> None:
    """#105 end to end: the action survives the real integration call, not just the private gate.

    The gate is reached from ``_score_and_gate``, so a fix applied only where the unit test calls it
    would leave the product path unchanged. ``failed_isoform_coverage`` counts rejections, so a warn
    action must leave it at zero while still recording the verdict.
    """
    workflow = _workflow(
        tmp_path,
        "iso_integration",
        stated={"min_isoform_coverage": 0.9},
        filter_actions={"min_isoform_coverage": "warn"},
    )
    candidate = _candidate()
    # One of three protein-coding transcripts carries this guide: coverage 1/3, below the floor.
    workflow._protein_coding_transcript_ids = {"ENST1", "ENST2", "ENST3"}
    workflow._protein_coding_transcript_count = 3
    workflow._guide_to_transcripts = {GUIDE: {"ENST1"}}

    _, stats = workflow._integrate_offtarget_results(
        [candidate], dict(NO_HITS), OffTargetFilterCriteria(), screened_species=["human"]
    )

    assert candidate.isoform_coverage == pytest.approx(1 / 3)
    assert candidate.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.FAIL.value
    assert candidate.passes_filters is True, "warn does not reject"
    assert stats["failed_isoform_coverage"] == 0, "the counter counts rejections, and warn is not one"


@pytest.mark.unit
def test_a_fail_action_isoform_gate_is_still_counted(tmp_path: Path) -> None:
    """#105: the default action is unchanged, so the run still counts and labels the rejection."""
    workflow = _workflow(tmp_path, "iso_integration_fail", stated={"min_isoform_coverage": 0.9})
    candidate = _candidate()
    workflow._protein_coding_transcript_ids = {"ENST1", "ENST2", "ENST3"}
    workflow._protein_coding_transcript_count = 3
    workflow._guide_to_transcripts = {GUIDE: {"ENST1"}}

    _, stats = workflow._integrate_offtarget_results(
        [candidate], dict(NO_HITS), OffTargetFilterCriteria(), screened_species=["human"]
    )

    assert candidate.passes_filters == SiRNACandidate.FilterStatus.LOW_ISOFORM_COVERAGE
    assert stats["failed_isoform_coverage"] == 1


@pytest.mark.unit
def test_design_parameters_still_own_the_coverage_floor(tmp_path: Path) -> None:
    """The resolver, not the gate, owns the number: a stated floor must reach FilterCriteria.

    Pinned because the gate reads ``design_params.filters.min_isoform_coverage`` while the action
    comes from the descriptor. If those two ever resolved from different places, a run could gate on
    one threshold and report another.
    """
    workflow = _workflow(tmp_path, "floor_owner", stated={"min_isoform_coverage": 0.42})
    assert workflow.config.design_params.filters.min_isoform_coverage == pytest.approx(0.42)
    assert workflow.config.resolved_policy.descriptor("min_isoform_coverage").threshold == pytest.approx(0.42)
    assert FilterCriteria().min_isoform_coverage is None, "still off by default"


@pytest.mark.unit
def test_default_design_parameters_are_unchanged_by_these_fixes() -> None:
    """None of this changes what a default run gates on. A behaviour fix that moves a default is two changes."""
    defaults = DesignParameters()
    assert defaults.filters.min_isoform_coverage is None
    assert defaults.offtarget_filters.max_transcriptome_hits_0mm == 1
    assert defaults.offtarget_filters.max_mirna_perfect_seed == 0


@pytest.mark.unit
def test_a_withheld_candidate_is_labelled_not_deleted(tmp_path: Path) -> None:
    """#100: the qualified shortlist narrows; the published rows record why, and are never dropped.

    An earlier pass at this filtered ``candidates_pass.csv``/``.fasta`` down to the eligible set, on the
    reasoning that they must not disagree with ``top_candidates``. Measured against the previous tip
    that deleted a whole deliverable: ``--input-fasta`` with no ``--transcriptome-fasta`` resolves no
    reference, screens nothing, and 152 guides became 0 with the FASTA removed outright. #100 asks for
    the opposite -- "legacy pass outputs record selection semantics" -- so the state is a column and
    only ``top_candidates`` is narrowed.
    """
    workflow = _workflow(tmp_path, "withheld_label", run_mode=RunMode.QUALIFIED)
    eligible, withheld = _candidate("eligible"), _candidate("withheld")
    eligible.off_target_screened = True
    eligible.scored_after_screening = True
    withheld.off_target_screened = False
    withheld.scored_after_screening = True

    result = _design_result(workflow, [eligible, withheld])
    workflow._apply_post_screen_ranking(result)

    assert [c.id for c in result.top_candidates] == ["eligible"]
    assert eligible.selection_state == "eligible"
    assert withheld.selection_state == "withheld_incomplete_evidence"
    # Both still pass their gates, so both still reach the published pass rows.
    assert withheld.passes_filters is True
    assert build_candidate_row(withheld)["selection_state"] == "withheld_incomplete_evidence"


@pytest.mark.unit
def test_selection_state_is_not_selected_before_any_selection_runs() -> None:
    """A row from a run that never selected must not claim to have been judged.

    ``sirnaforge design`` writes candidate rows through the same builder and runs no selection, so the
    default has to be a third state rather than either verdict.
    """
    assert build_candidate_row(_candidate())["selection_state"] == "not_selected"


@pytest.mark.unit
def test_unresolved_orthology_reports_no_conservation_rather_than_zero(tmp_path: Path) -> None:
    """#101: an unresolvable orthologue lookup publishes null conservation, never a uniform 0.0.

    ``conservation_score`` is (non-query species with an ortholog hit) / (non-query species screened).
    A species whose Compara lookup never completed cannot contribute to that numerator, so the ratio
    reported "conserved in 0 of N" for a question nobody managed to ask -- a 0.0 that reads as a
    measurement, and strictly worse than the null it replaced. The denominator is deliberately NOT
    shrunk to the resolvable species: that was tried and it let a degraded run outscore the complete
    run it degraded from.

    The same run with orthology resolved still reports the measured 0.0, which is what shows the null
    is about the lookup and not about the absence of hits.
    """
    workflow = _workflow(tmp_path, "conservation", screen_species=["human", "mouse"])
    workflow._active_screen_species = ["human", "mouse"]

    unresolved = _candidate("unresolved")
    workflow._integrate_offtarget_results(
        [unresolved],
        dict(NO_HITS),
        screened_species=["human", "mouse"],
        unresolved_orthology_species=frozenset({"mouse"}),
    )
    assert unresolved.scored_after_screening is True, "conservation is reported, not scored, so scoring still succeeds"
    assert unresolved.conservation_score is None

    resolved = _candidate("resolved")
    workflow._integrate_offtarget_results([resolved], dict(NO_HITS), screened_species=["human", "mouse"])
    assert resolved.conservation_score == pytest.approx(0.0), "the lookup completed and found nothing: a real zero"


@pytest.mark.unit
def test_an_unresolved_species_outside_the_denominator_does_not_null_conservation(tmp_path: Path) -> None:
    """Only the species the ratio actually counts can make it unknown.

    The denominator is the non-query species handed to the aligner. An unresolved lookup for a species
    this run never screened says nothing about the ratio, and nulling on it would delete conservation
    from runs whose own evidence was complete.
    """
    workflow = _workflow(tmp_path, "conservation_scope", screen_species=["human", "mouse"])
    workflow._active_screen_species = ["human", "mouse"]

    candidate = _candidate()
    workflow._integrate_offtarget_results(
        [candidate],
        dict(NO_HITS),
        screened_species=["human", "mouse"],
        unresolved_orthology_species=frozenset({"rat"}),
    )
    assert candidate.conservation_score == pytest.approx(0.0)


@pytest.mark.unit
def test_a_two_channel_gate_does_not_disqualify_on_an_exploratory_absence(tmp_path: Path) -> None:
    """#101/#100: a gate reading one required and one exploratory channel cannot disqualify.

    ``max_total_offtarget_hits`` sums the transcriptome and miRNA channels, so an absent miRNA
    aggregate alone makes it ``UNKNOWN``. Scoping the unknown-blocks rule by *intersection* with the
    required channels put that gate in scope anyway and emptied the shortlist of a complete, clean
    transcriptome screen -- the exact mirror of the defect the scoping exists to prevent. The test is a
    subset: every channel a gate reads must be required before its unknown can disqualify.
    """
    workflow = _workflow(tmp_path, "two_channel", run_mode=RunMode.QUALIFIED)
    criteria = workflow.config.design_params.offtarget_filters.model_copy(update={"max_total_offtarget_hits": 200})
    candidate = _candidate()

    workflow._integrate_offtarget_results(
        [candidate], dict(NO_HITS), criteria, screened_species=["human"], mirna_screened=False
    )
    assert candidate.filter_verdicts["max_total_offtarget_hits"] == FilterEvaluation.UNKNOWN.value
    assert candidate.filter_verdicts["max_transcriptome_hits_0mm"] == FilterEvaluation.PASS.value

    result = _design_result(workflow, [candidate])
    workflow._apply_post_screen_ranking(result)

    assert result.top_candidates == [candidate]
    assert workflow._selection_summary["evidence_shortfall_reasons"] == {}


@pytest.mark.unit
def test_exploratory_sorts_incomplete_evidence_below_complete(tmp_path: Path) -> None:
    """#100: a provisional result must not outrank a qualified one in the same list.

    Exploratory mode retains incomplete evidence, so ordering is the only thing left that can keep it
    from leading the shortlist. It cannot be expressed with the *required* shortfall, because nothing
    is required in exploratory mode and that set is always empty -- which is why the ordering question
    and the disqualifying question are asked separately. Here the incomplete candidate has the higher
    score, so a sort that ignored evidence would put it first.
    """
    workflow = _workflow(
        tmp_path, "exploratory_order", run_mode=RunMode.EXPLORATORY, stated={"min_isoform_coverage": 0.5}
    )
    unknown_coverage = _candidate("high_score_unknown_coverage")
    known_coverage = _candidate("low_score_known_coverage")

    workflow._integrate_offtarget_results(
        [unknown_coverage, known_coverage], dict(NO_HITS), screened_species=["human"], mirna_screened=True
    )
    unknown_coverage.isoform_coverage, known_coverage.isoform_coverage = None, 0.9
    unknown_coverage.design_score, known_coverage.design_score = 90.0, 10.0
    for candidate in (unknown_coverage, known_coverage):
        candidate.composite_score = None
        candidate.filter_verdicts.pop("min_isoform_coverage", None)
        workflow._apply_isoform_coverage_gate(candidate)

    result = _design_result(workflow, [unknown_coverage, known_coverage])
    workflow._apply_post_screen_ranking(result)

    assert [c.id for c in result.top_candidates] == ["low_score_known_coverage", "high_score_unknown_coverage"]
    assert unknown_coverage.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.UNKNOWN.value


@pytest.mark.unit
def test_the_report_cannot_re_derive_a_pass_from_a_verdict_the_run_called_unknown(tmp_path: Path) -> None:
    """#103/#106: the run's UNKNOWN is authoritative; re-thresholding cannot conjure the measurement.

    The gate records ``UNKNOWN`` with an *empty* observed value so nothing re-derives a pass from it.
    But ``observed_column`` falls back to the descriptor's own column when the observed one is empty for
    every row, and three of those fallbacks are exported and default to 0 -- so the report called an
    unscreened guide's ``max_off_target_count`` a confident pass at 0, from a run that said it did not
    know. Reading the recorded verdict first is what closes that.
    """
    workflow = _workflow(tmp_path, "report_unknown", run_mode=RunMode.QUALIFIED)
    candidate = _candidate()
    workflow._integrate_offtarget_results([candidate], dict(NO_HITS), screened_species=[], mirna_screened=False)
    assert candidate.filter_verdicts["max_off_target_count"] == FilterEvaluation.UNKNOWN.value

    row = build_candidate_row(candidate)
    # off_target_count is exported and defaults to 0, so it is the fallback column that used to lie.
    assert row["off_target_count"] == 0
    descriptor = workflow.config.resolved_policy.descriptor("max_off_target_count")
    populated = {key for key, value in row.items() if value is not None}
    value, verdict_code, _ = evaluate_descriptor(descriptor, pd.Series(row), observed_column(descriptor, populated))

    assert verdict_code == REPORT_UNKNOWN
    assert value is None
