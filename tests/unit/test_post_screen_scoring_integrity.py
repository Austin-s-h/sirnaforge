"""Regression tests for post-screen scoring integrity (issue #80 findings 8-10).

The reorder made screening run before final scoring, so ``_integrate_offtarget_results`` now owns
the composite score every candidate is ranked by. Ways that ownership shipped wrong rankings
silently:

8. miRNA mode: the biogenesis bonuses were folded into composite_score only, so rebuilding the
   composite after screening dropped them and ``--design-mode mirna`` stopped affecting the ranking.
   The bonus normaliser also has to reach rows that never went through MiRNADesigner (dirty
   controls), or those rows land on a different scale in the same CSV.
9. A run whose query-species alignment never happened scored every candidate as perfectly specific.
   The evidence for "it happened" has to be POSITIVE: the aggregate's self-reported missing_species
   is empty precisely when aggregation itself failed.
10. Conservation was scored against the CLI species list rather than the set handed to the aligner,
    so an extra species could make the numerator exceed the denominator and abandon scoring
    altogether. Subtracting the *failed* species from that denominator then let a partial screen
    outscore the complete screen it degraded from.

Each half of a fix is pinned by a test that changes answer when only that half is reverted; the
two conservation halves mask each other under weaker assertions.
"""

import asyncio
import json
from pathlib import Path

import pytest

from sirnaforge.core.design import (
    MIRNA_BONUS_KEY,
    MIRNA_BONUS_MAX_KEY,
    MiRNADesigner,
    mirna_max_biogenesis_bonus,
)
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    OffTargetFilterCriteria,
    SiRNACandidate,
)
from sirnaforge.utils.control_candidates import DIRTY_CONTROL_LABEL, DIRTY_CONTROL_SUFFIX
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

# Design-time composites put HIGH_BASE_GUIDE above BONUS_GUIDE, and the miRNA biogenesis bonuses
# (A/U at guide position 1, pos1 pairing, 3' supplementary) reverse that. The pair therefore
# distinguishes a miRNA-aware ranking from a plain siRNA one.
BONUS_GUIDE = "TACCTCCGGGTAGCGTTGGCA"
HIGH_BASE_GUIDE = "CCCCCGATTACCGGCCGCCTT"

NO_HITS = {"status": "completed", "results": {}}


def _candidate(candidate_id: str, guide: str) -> SiRNACandidate:
    """Build an unscored candidate; the designers fill in every score field."""
    passenger = guide.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    gc = (guide.count("G") + guide.count("C")) / len(guide) * 100
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=guide,
        passenger_sequence=passenger,
        gc_content=gc,
        length=len(guide),
        asymmetry_score=0.7,
        composite_score=0.0,
        mfe=-4.2,
        duplex_stability=-39.0,
        structure="." * len(guide),
    )


def _workflow(
    tmp_path: Path,
    out_name: str,
    design_params: DesignParameters | None = None,
    genome_species: list[str] | None = None,
    nextflow_config: dict[str, str] | None = None,
) -> SiRNAWorkflow:
    """Build a workflow instance around a query gene, without running it."""
    config = WorkflowConfig(
        output_dir=tmp_path / out_name,
        gene_query="TP53",
        genome_species=genome_species,
        nextflow_config=nextflow_config,
        design_params=design_params or DesignParameters(),
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    return workflow


def _build_species_indices(workflow: SiRNAWorkflow, tmp_path: Path) -> None:
    """Give the classifier a TP53 ortholog in human, mouse and rat so hits classify as ORTHOLOG."""
    for species, transcript, symbol in (
        ("human", "ENST00000000001", "TP53"),
        ("mouse", "ENSMUST00000000002", "Tp53"),
        ("rat", "ENSRNOT00000000003", "Tp53"),
    ):
        reference = tmp_path / f"{species}_cdna.fasta"
        reference.write_text(
            f">{transcript} gene:ENSG_{species} gene_symbol:{symbol} transcript_biotype:protein_coding\nACGT\n"
        )
        workflow._transcript_index.build(species, reference)


@pytest.mark.unit
def test_mirna_bonus_survives_off_target_screening(tmp_path: Path) -> None:
    """F8: post-screen scoring must not discard the miRNA-specific score components.

    MiRNADesigner writes the ago-start / pos1 / supplementary bonuses into composite_score only,
    never into component_scores' term set. Pre-fix, _score_candidate_post_screen rebuilt the
    composite from that term set and overwrote composite_score, so every screened miRNA run was
    ranked by the plain siRNA composite and --design-mode mirna had no effect on the output.
    """
    mirna_params = DesignParameters(design_mode=DesignMode.MIRNA)
    designer = MiRNADesigner(mirna_params)

    def _scored_pair() -> tuple[SiRNACandidate, SiRNACandidate]:
        bonus, high_base = _candidate("cand_bonus", BONUS_GUIDE), _candidate("cand_high_base", HIGH_BASE_GUIDE)
        designer._score_candidates([bonus, high_base])
        return bonus, high_base

    bonus, high_base = _scored_pair()
    assert bonus.component_scores[MIRNA_BONUS_KEY] > high_base.component_scores[MIRNA_BONUS_KEY]
    assert bonus.composite_score > high_base.composite_score, "miRNA design ranks the bonus guide first"

    mirna_workflow = _workflow(tmp_path, "mirna_out", mirna_params)
    mirna_workflow._integrate_offtarget_results([bonus, high_base], NO_HITS, OffTargetFilterCriteria())

    # Control: the same candidates, screened identically, but in plain siRNA mode.
    sirna_bonus, sirna_high_base = _scored_pair()
    sirna_workflow = _workflow(tmp_path, "sirna_out", DesignParameters())
    sirna_workflow._integrate_offtarget_results([sirna_bonus, sirna_high_base], NO_HITS, OffTargetFilterCriteria())

    assert bonus.scored_after_screening is True
    assert high_base.scored_after_screening is True
    assert sirna_high_base.composite_score > sirna_bonus.composite_score, (
        "without the miRNA bonuses the higher base score wins; otherwise this pair proves nothing"
    )
    assert bonus.composite_score > high_base.composite_score, (
        "miRNA mode must keep ranking the bonus guide first after screening"
    )


@pytest.mark.unit
def test_partial_run_missing_query_species_does_not_score_perfect_specificity(tmp_path: Path) -> None:
    """F9: an alignment that never ran must not be scored as zero off-targets.

    _process_nextflow_results already detects that a species produced no alignment files and
    downgrades the run to "partial". Pre-fix it still integrated the (empty) results, so every
    candidate took the no-hits branch, received off_target_sub_score(0) = 1.0 -- perfect
    specificity -- and was marked off_target_screened=True.
    """
    workflow = _workflow(tmp_path, "partial_out", genome_species=["human", "mouse"])
    candidate = _candidate("cand_partial", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    design_time_score = candidate.composite_score

    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated_dir = results_dir / "aggregated"
    aggregated_dir.mkdir(parents=True, exist_ok=True)
    # The human BWA-MEM2 index OOMed: no alignment files, so no hits for any candidate.
    (aggregated_dir / "combined_summary.json").write_text(
        json.dumps({"total_results": 0, "human_hits": 0, "other_species_hits": 0, "missing_species": ["human"]})
    )

    outcome = asyncio.run(workflow._process_nextflow_results([candidate], results_dir, {}))

    assert outcome["status"] == "partial"
    assert candidate.off_target_screened is False, "no alignment ran, so the counts are unknown"
    assert candidate.scored_after_screening is False
    assert candidate.composite_score == design_time_score, "the design-time score must be kept, not overwritten"
    assert candidate.score_off_target is None, "no off-target term may be claimed"
    assert outcome["filtering_stats"]["candidates_not_scored_after_screening"] == 1


@pytest.mark.unit
def test_candidate_never_submitted_is_not_scored_as_clean(tmp_path: Path) -> None:
    """F9: absence from the results is only 'clean' for a candidate that reached the aligner.

    Pre-fix, off_target_screened was set True before any check that the candidate had been
    submitted, and the no-hits branch then awarded it perfect specificity -- the same output a
    genuinely clean candidate gets.
    """
    workflow = _workflow(tmp_path, "unsubmitted_out")
    submitted = _candidate("cand_submitted", BONUS_GUIDE)
    unsubmitted = _candidate("cand_unsubmitted", HIGH_BASE_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([submitted, unsubmitted])
    design_time_score = unsubmitted.composite_score

    # Only one of the two ever reaches the aligner.
    asyncio.run(workflow._prepare_offtarget_input([submitted]))

    _, stats = workflow._integrate_offtarget_results([submitted, unsubmitted], NO_HITS, OffTargetFilterCriteria())

    assert submitted.off_target_screened is True
    assert submitted.scored_after_screening is True
    assert unsubmitted.off_target_screened is False, "never aligned, so the counts are unknown"
    assert unsubmitted.scored_after_screening is False
    assert unsubmitted.composite_score == design_time_score
    assert unsubmitted.score_off_target is None
    assert stats["candidates_not_scored_after_screening"] == 1


@pytest.mark.unit
def test_conservation_denominator_counts_species_screened_only_via_indices(tmp_path: Path) -> None:
    """F10 (denominator half): the denominator is the species handed to the aligner, not the CLI list.

    Species supplied only through genome_fastas/genome_indices are appended to the ACTIVE list
    handed to Nextflow but never to mirna_genome_species. The ortholog hit here is in mouse ONLY,
    so the two candidate denominators give different answers and the test can tell them apart:
    the screened set {mouse, rat} scores 1/2, while the CLI list {mouse} scores a perfect 1/1.
    Pinning the *value* rather than merely "scoring did not abort" is deliberate -- the numerator
    intersection added alongside this fix silently repairs the ratio when only the denominator is
    wrong, so a weaker assertion passes with the denominator reverted.
    """
    workflow = _workflow(
        tmp_path,
        "conservation_out",
        genome_species=["human", "mouse"],
        nextflow_config={"genome_fastas": "human:human.fa,mouse:mouse.fa,rat:rat.fa"},
    )
    candidate = _candidate("cand_conservation", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    design_time_score = candidate.composite_score
    _build_species_indices(workflow, tmp_path)

    assert workflow._resolve_active_genome_species({}) == ["human", "mouse", "rat"], (
        "rat is screened even though it was never named on --genome-species"
    )

    hits = [{"rname": "ENSMUST00000000002", "species": "mouse", "nm": 0, "seed_mismatches": 0}]
    _, stats = workflow._integrate_offtarget_results(
        [candidate],
        {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}},
        OffTargetFilterCriteria(),
    )

    assert candidate.ortholog_hits == 1
    assert candidate.scored_after_screening is True, "scoring must not be abandoned"
    assert candidate.conservation_score == pytest.approx(0.5), "one ortholog out of the two screened non-query species"
    assert candidate.composite_score != design_time_score
    assert stats["candidates_not_scored_after_screening"] == 0


@pytest.mark.unit
def test_ortholog_hits_outside_the_requested_species_do_not_break_conservation(tmp_path: Path) -> None:
    """F10 (numerator half): the numerator is intersected with the denominator.

    A hit table can name a species this run never asked the aligner for (a stale cached result, a
    hand-edited TSV). Counting it would push the conservation numerator above its denominator,
    conservation_sub_score raises, and the blanket except abandons scoring for the whole candidate
    -- which then keeps its optimistic design-time score while its neighbours carry post-screen
    ones. Reverting only the intersection leaves 3 ortholog species over a denominator of 2 here.
    """
    workflow = _workflow(
        tmp_path,
        "unrequested_out",
        genome_species=["human", "mouse"],
        nextflow_config={"genome_fastas": "human:human.fa,mouse:mouse.fa,rat:rat.fa"},
    )
    candidate = _candidate("cand_unrequested", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    _build_species_indices(workflow, tmp_path)
    # dog was never handed to the aligner, but a dog hit is in the table anyway.
    dog_reference = tmp_path / "dog_cdna.fasta"
    dog_reference.write_text(
        ">ENSCAFT00000000004 gene:ENSG_dog gene_symbol:Tp53 transcript_biotype:protein_coding\nACGT\n"
    )
    workflow._transcript_index.build("dog", dog_reference)

    assert workflow._resolve_active_genome_species({}) == ["human", "mouse", "rat"]

    hits = [
        {"rname": "ENSMUST00000000002", "species": "mouse", "nm": 0, "seed_mismatches": 0},
        {"rname": "ENSRNOT00000000003", "species": "rat", "nm": 0, "seed_mismatches": 0},
        {"rname": "ENSCAFT00000000004", "species": "dog", "nm": 0, "seed_mismatches": 0},
    ]
    _, stats = workflow._integrate_offtarget_results(
        [candidate],
        {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}},
        OffTargetFilterCriteria(),
    )

    assert candidate.ortholog_hits == 3, "the dog hit is still reported as an ortholog hit"
    assert candidate.scored_after_screening is True, "an unrequested species must not abort scoring"
    assert candidate.conservation_score == pytest.approx(1.0), "mouse and rat, the two species actually requested"
    assert stats["candidates_not_scored_after_screening"] == 0


@pytest.mark.unit
def test_partial_screen_does_not_outscore_the_equivalent_complete_screen(tmp_path: Path) -> None:
    """A degraded run must never rank above the complete run it degraded from.

    Dropping an unscreened species from the conservation denominator looks conservative and is the
    opposite: when the missing species is the only non-query one the denominator empties, the term
    goes inactive, and compute_composite renormalises its 0.10 weight onto the surviving terms. The
    same candidate then scored 85.2 on the broken screen against 76.7 on the good one. A species
    whose alignment produced nothing stays in the denominator instead: it cannot supply an ortholog
    hit, so conservation becomes a lower bound rather than a bonus.
    """

    def _run(out_name: str, screened: list[str]) -> SiRNACandidate:
        workflow = _workflow(tmp_path, out_name, genome_species=["human", "mouse"])
        candidate = _candidate(f"cand_{out_name}", BONUS_GUIDE)
        MiRNADesigner(DesignParameters())._score_candidates([candidate])
        workflow._integrate_offtarget_results(
            [candidate], NO_HITS, OffTargetFilterCriteria(), screened_species=screened
        )
        return candidate

    complete = _run("complete_screen", ["human", "mouse"])
    partial = _run("partial_screen", ["human"])  # the mouse alignment produced no files

    assert complete.scored_after_screening is True
    assert partial.scored_after_screening is True, "the query species aligned, so scoring still runs"
    assert partial.composite_score <= complete.composite_score, (
        f"a partial screen ({partial.composite_score:.2f}) must never score above the equivalent "
        f"complete screen ({complete.composite_score:.2f})"
    )
    assert complete.conservation_score == pytest.approx(0.0)
    assert partial.conservation_score == pytest.approx(0.0), "mouse stays in the denominator, unhit"


@pytest.mark.unit
def test_transcriptome_index_species_stay_in_the_conservation_denominator(tmp_path: Path) -> None:
    """main.nf screens --transcriptome-indices species, so conservation must count them.

    main.nf mixes transcriptome_indices into the same ch_genomes alignment channel as
    genome_indices/genome_fastas, and the subworkflow derives the aggregation species list from
    that channel. _resolve_active_genome_species scanned only the two genome_* keys, so a
    transcriptome-only species was dropped from the ACTIVE list -- and with mouse dropped the
    denominator here empties, conservation goes inactive, and its weight is redistributed to the
    remaining terms even though mouse was aligned and its ortholog hit is right there in the table.
    """
    workflow = _workflow(tmp_path, "tx_indices_out", genome_species=["human", "mouse"])
    candidate = _candidate("cand_tx_indices", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    _build_species_indices(workflow, tmp_path)

    params = {"genome_indices": "human:human.idx", "transcriptome_indices": "mouse:mouse.idx"}
    assert workflow._resolve_active_genome_species(params) == ["human", "mouse"]

    hits = [{"rname": "ENSMUST00000000002", "species": "mouse", "nm": 0, "seed_mismatches": 0}]
    workflow._integrate_offtarget_results(
        [candidate],
        {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}},
        OffTargetFilterCriteria(),
    )

    assert candidate.conservation_score == pytest.approx(1.0), "mouse was screened, so it is in the denominator"
    assert candidate.score_conservation is not None, "the conservation term must be active"


@pytest.mark.unit
def test_unscored_candidates_do_not_outrank_post_screen_scores(tmp_path: Path) -> None:
    """F10 (ranking half): a candidate that failed to score must not compete on its old score.

    A design-time composite lacks the off_target/isoform_coverage/conservation terms and is
    systematically the more optimistic number, so sorting it against post-screen neighbours put
    exactly the candidates whose evidence is missing at the top of candidates_pass.csv.
    """
    params = DesignParameters(apply_modifications=False)
    workflow = _workflow(tmp_path, "ranking_out", params)

    scored = _candidate("cand_scored", BONUS_GUIDE)
    unscored = _candidate("cand_unscored", HIGH_BASE_GUIDE)
    MiRNADesigner(params)._score_candidates([scored, unscored])
    for candidate in (scored, unscored):
        candidate.passes_filters = True
    # An out-of-range sub-score is the one thing compute_composite refuses outright, so this
    # candidate reaches ranking with a design-time score and scored_after_screening False.
    unscored.component_scores["asymmetry"] = 1.5
    unscored.composite_score = 99.0

    _, stats = workflow._integrate_offtarget_results([scored, unscored], NO_HITS, OffTargetFilterCriteria())

    assert scored.scored_after_screening is True
    assert unscored.scored_after_screening is False

    design_result = DesignResult(
        input_file="<test>",
        parameters=params,
        candidates=[scored, unscored],
        top_candidates=[scored, unscored],
        total_sequences=2,
        total_candidates=2,
        filtered_candidates=2,
        processing_time=0.1,
    )
    workflow._apply_post_screen_ranking(design_result)

    assert unscored.composite_score == 99.0, "the design-time score is kept for reporting"
    assert design_result.candidates[0] is scored, "post-screen scores rank above design-time ones"
    assert unscored not in design_result.top_candidates
    assert scored in design_result.top_candidates
    assert stats["candidates_not_scored_after_screening"] == 1, "a degraded run must be countable"


@pytest.mark.unit
def test_missing_aggregate_is_not_treated_as_a_clean_screen(tmp_path: Path) -> None:
    """Aggregation failing is the case where the aggregate's own missing_species list lies.

    _parse_nextflow_results reports "completed" whenever the output directory merely exists, and
    _load_offtarget_aggregates returns {} when combined_summary.json is absent -- so keying the
    guard on the aggregate's self-reported missing_species saw an empty list exactly when nothing
    had been aggregated at all, and every candidate was marked off_target_screened=True with
    off_target_sub_score(0) = 1.0. Positive evidence (a published alignment file per species) is
    the only reading that survives its own absence.
    """
    workflow = _workflow(tmp_path, "no_aggregate_out", genome_species=["human", "mouse"])
    candidate = _candidate("cand_no_aggregate", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    design_time_score = candidate.composite_score

    # The directory exists (Nextflow made it) but nothing was ever aggregated into it.
    results_dir = workflow.config.output_dir / "off_target" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    outcome = asyncio.run(workflow._process_nextflow_results([candidate], results_dir, {}))

    assert candidate.score_off_target is None, "no off-target term may be claimed with nothing aligned"
    assert candidate.off_target_screened is False
    assert candidate.scored_after_screening is False
    assert candidate.composite_score == design_time_score
    assert outcome["status"] == "partial", "no alignment evidence is not a completed run"
    assert outcome["filtering_stats"]["candidates_not_scored_after_screening"] == 1
    assert any("no aggregated summary" in warning for warning in outcome["warnings"])


@pytest.mark.unit
def test_mirna_only_mode_does_not_claim_transcriptome_specificity(tmp_path: Path) -> None:
    """miRNA-only mode aligns against no transcriptome, so it cannot earn the off-target term.

    sirna_offtarget_analysis.nf derives the aggregation species list from ch_genome_indices with
    `.ifEmpty { '' }`, so a miRNA-only run publishes a summary naming no species and reporting no
    missing ones. That shape used to read as a complete, perfectly specific screen.
    """
    workflow = _workflow(tmp_path, "mirna_only_out", genome_species=["human", "mouse"])
    candidate = _candidate("cand_mirna_only", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    design_time_score = candidate.composite_score

    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated_dir = results_dir / "aggregated"
    aggregated_dir.mkdir(parents=True, exist_ok=True)
    (aggregated_dir / "combined_summary.json").write_text(
        json.dumps({"total_results": 0, "species_analyzed": [], "missing_species": [], "species_file_counts": {}})
    )

    outcome = asyncio.run(workflow._process_nextflow_results([candidate], results_dir, {}))

    assert candidate.score_off_target is None, "a miRNA-only run has no transcriptome evidence"
    assert candidate.off_target_screened is False
    assert candidate.scored_after_screening is False
    assert candidate.composite_score == design_time_score
    assert outcome["status"] == "partial"


@pytest.mark.unit
def test_partial_screen_reports_its_hits_as_a_lower_bound(tmp_path: Path) -> None:
    """A row with off_target_screened=False may still carry hits: they are a floor, not a total.

    Contract test for the semantics documented in models/sirna.py and docs/models_and_scoring.md.
    Hits found in the species that did align are real evidence -- they are reported and still fed
    to the filters -- so the flag means "this screen was incomplete", not "these counts are zero".
    Zeroing them would produce the worse contradiction of a filter failure with no hits behind it.
    """
    workflow = _workflow(tmp_path, "lower_bound_out", genome_species=["human", "mouse"])
    candidate = _candidate("cand_lower_bound", BONUS_GUIDE)
    MiRNADesigner(DesignParameters())._score_candidates([candidate])
    design_time_score = candidate.composite_score
    _build_species_indices(workflow, tmp_path)

    # Mouse aligned and found an ortholog; the human alignment never ran.
    hits = [{"rname": "ENSMUST00000000002", "species": "mouse", "nm": 0, "seed_mismatches": 0}]
    workflow._integrate_offtarget_results(
        [candidate],
        {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}},
        OffTargetFilterCriteria(),
        screened_species=["mouse"],
    )

    assert candidate.off_target_screened is False, "the query species never aligned"
    assert candidate.ortholog_hits == 1, "evidence that was found is still reported"
    assert candidate.ortholog_species == "mouse"
    assert candidate.scored_after_screening is False
    assert candidate.composite_score == design_time_score


@pytest.mark.unit
def test_dirty_controls_share_the_mirna_score_scale(tmp_path: Path) -> None:
    """Every row of a miRNA run must be divided by the same maximum-bonus normaliser.

    Dirty controls are deep copies of *rejected* candidates, which never reach
    MiRNADesigner._score_candidates and so carry no recorded bonus maximum. Defaulting that
    maximum to 0.0 left those rows undivided while every real candidate was divided by 1.25,
    putting the controls ~25% high in the same candidates_all.csv.
    """
    params = DesignParameters(design_mode=DesignMode.MIRNA)
    workflow = _workflow(tmp_path, "controls_out", params)

    real = _candidate("cand_real", BONUS_GUIDE)
    MiRNADesigner(params)._score_candidates([real])
    assert real.component_scores[MIRNA_BONUS_KEY] > 0.0, "this guide must earn a bonus or the test proves nothing"

    control = real.model_copy(deep=True)
    control.id = f"{real.id}{DIRTY_CONTROL_SUFFIX}_1"
    control.quality_issues = [DIRTY_CONTROL_LABEL]
    for key in (MIRNA_BONUS_KEY, MIRNA_BONUS_MAX_KEY):
        control.component_scores.pop(key)

    workflow._integrate_offtarget_results([real, control], NO_HITS, OffTargetFilterCriteria())

    contributions = sum(
        value
        for value in (
            control.score_asymmetry,
            control.score_gc_content,
            control.score_accessibility,
            control.score_empirical,
            control.score_off_target,
            control.score_isoform_coverage,
            control.score_conservation,
        )
        if value is not None
    )
    assert control.composite_score == pytest.approx(contributions / (1.0 + mirna_max_biogenesis_bonus())), (
        "a control with no recorded bonus is still divided by the mode's maximum bonus"
    )
    assert control.composite_score < real.composite_score, (
        "a control cloned from a rejected candidate must not outscore the real candidate it copied"
    )
