"""Regression tests for post-screen scoring integrity (issue #80 findings 8-10).

The reorder made screening run before final scoring, so ``_integrate_offtarget_results`` now owns
the composite score every candidate is ranked by. Three ways that ownership shipped wrong rankings
silently:

8. miRNA mode: the biogenesis bonuses were folded into composite_score only, so rebuilding the
   composite after screening dropped them and ``--design-mode mirna`` stopped affecting the ranking.
9. A run whose query-species alignment never happened scored every candidate as perfectly specific.
10. Conservation was scored against the CLI species list rather than the screened one, so an extra
    species could make the numerator exceed the denominator and abandon scoring altogether.
"""

import asyncio
import json
from pathlib import Path

import pytest

from sirnaforge.core.design import MIRNA_BONUS_KEY, MiRNADesigner
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    OffTargetFilterCriteria,
    SiRNACandidate,
)
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
def test_conservation_denominator_covers_species_screened_only_via_indices(tmp_path: Path) -> None:
    """F10: the conservation denominator must be the screened species, not the CLI species list.

    Species supplied only through genome_fastas/genome_indices are appended to the ACTIVE list
    handed to Nextflow but never to mirna_genome_species. Pre-fix, an ortholog hit in such a
    species made the conservation numerator exceed its denominator, conservation_sub_score raised,
    and the blanket except abandoned scoring for the whole candidate -- which then kept its
    optimistic design-time score while its neighbours carried post-screen ones.
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

    assert workflow._resolve_active_genome_species({}) == ["human", "mouse", "rat"], (
        "rat is screened even though it was never named on --genome-species"
    )

    hits = [
        {"rname": "ENSMUST00000000002", "species": "mouse", "nm": 0, "seed_mismatches": 0},
        {"rname": "ENSRNOT00000000003", "species": "rat", "nm": 0, "seed_mismatches": 0},
    ]
    _, stats = workflow._integrate_offtarget_results(
        [candidate],
        {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}},
        OffTargetFilterCriteria(),
    )

    assert candidate.ortholog_hits == 2
    assert candidate.scored_after_screening is True, "scoring must not be abandoned"
    assert candidate.conservation_score == pytest.approx(1.0), "orthologs in both screened non-query species"
    assert candidate.composite_score != design_time_score
    assert stats["candidates_not_scored_after_screening"] == 0


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
