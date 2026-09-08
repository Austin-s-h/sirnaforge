"""Regression tests for issue #80 adversarial-review findings F1-F6.

Each test targets one finding from the round-2 review of feat/issue-80-hit-classification
(see .dev/ISSUE_80_CONTRACT.md) and is written to fail against the pre-fix behaviour, not
merely to exercise the fixed code path.
"""

import asyncio
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    OffTargetFilterCriteria,
    SiRNACandidate,
)
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig


def _full_candidate(
    candidate_id: str,
    guide_sequence: str = "ATGCGATGCGATGCGATGCGC",
    transcript_id: str = "ENST00000000001",
) -> SiRNACandidate:
    """Build a SiRNACandidate with every issue #80 output field populated to a real value."""
    passenger = guide_sequence.translate(str.maketrans("ATCG", "TAGC"))
    gc = (guide_sequence.count("G") + guide_sequence.count("C")) / len(guide_sequence) * 100
    return SiRNACandidate(
        id=candidate_id,
        transcript_id=transcript_id,
        position=1,
        guide_sequence=guide_sequence,
        passenger_sequence=passenger,
        gc_content=gc,
        length=len(guide_sequence),
        asymmetry_score=0.7,
        composite_score=82.0,
        mfe=-4.2,
        duplex_stability=-39.0,
        structure="." * len(guide_sequence),
        component_scores={
            "duplex_stability_score": 0.6,
            "dg_5p": -8.0,
            "dg_3p": -9.0,
            "delta_dg_end": 1.0,
            "melting_temp_c": 70.0,
        },
        on_target_hits=2,
        ortholog_hits=1,
        repeat_hits=1,
        ortholog_species="mouse",
        repeat_flagged=False,
        repeat_transcript_fraction=0.0005,
        isoform_coverage=0.66,
        conservation_score=0.5,
        score_asymmetry=8.4,
        score_gc_content=5.7,
        score_accessibility=9.1,
        score_empirical=12.0,
        score_off_target=20.0,
        score_isoform_coverage=9.9,
        score_conservation=5.0,
        scored_after_screening=True,
        weight_set_version="2.0.0",
    )


def _minimal_workflow(tmp_path: Path, out_name: str, design_params: DesignParameters | None = None) -> SiRNAWorkflow:
    """Build a workflow instance around a query gene, without running it."""
    config = WorkflowConfig(
        output_dir=tmp_path / out_name,
        gene_query="TP53",
        design_params=design_params or DesignParameters(),
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    return workflow


@pytest.mark.unit
def test_workflow_csv_emits_every_issue80_column_and_matches_save_csv(tmp_path: Path) -> None:
    """F2: the workflow's CSV writer must not drift from DesignResult.save_csv's columns.

    Pre-fix, step6_generate_reports built its own hand-rolled row dict that omitted
    on_target_hits, ortholog_hits, repeat_hits, ortholog_species, repeat_flagged,
    repeat_transcript_fraction, isoform_coverage, conservation_score, every score_* term,
    scored_after_screening and weight_set_version.
    """
    candidates = [
        _full_candidate("cand_a"),
        _full_candidate("cand_b", guide_sequence="GGCATTACGCATTACGCATTA", transcript_id="ENST00000000002"),
    ]
    params = DesignParameters(design_mode=DesignMode.MIRNA, apply_modifications=False)
    design_result = DesignResult(
        input_file="<test>",
        parameters=params,
        candidates=candidates,
        top_candidates=candidates,
        total_sequences=2,
        total_candidates=2,
        filtered_candidates=2,
        processing_time=0.1,
    )

    # Path 1: DesignResult.save_csv -- the `sirnaforge design` path.
    save_csv_path = tmp_path / "design_candidates.csv"
    design_result.save_csv(str(save_csv_path))
    save_csv_columns = set(pd.read_csv(save_csv_path, nrows=0).columns)

    # Path 2: SiRNAWorkflow.step6_generate_reports -- the `sirnaforge workflow` path.
    config = WorkflowConfig(output_dir=tmp_path / "wf_out", gene_query="TP53", design_params=params)
    workflow = SiRNAWorkflow(config)
    asyncio.run(workflow.step6_generate_reports(design_result))
    workflow_csv_path = config.output_dir / "sirnaforge" / "candidates_all.csv"
    assert workflow_csv_path.exists(), "step6_generate_reports must write candidates_all.csv"
    workflow_columns = set(pd.read_csv(workflow_csv_path, nrows=0).columns)

    new_issue80_columns = {
        "on_target_hits",
        "ortholog_hits",
        "repeat_hits",
        "ortholog_species",
        "repeat_flagged",
        "repeat_transcript_fraction",
        "isoform_coverage",
        "conservation_score",
        "score_asymmetry",
        "score_gc_content",
        "score_accessibility",
        "score_empirical",
        "score_off_target",
        "score_isoform_coverage",
        "score_conservation",
        "scored_after_screening",
        "weight_set_version",
    }
    missing = new_issue80_columns - workflow_columns
    assert not missing, f"workflow CSV is missing issue #80 columns: {sorted(missing)}"

    # The two writers must never again drift on which columns they emit.
    assert workflow_columns == save_csv_columns


@pytest.mark.unit
def test_enumerate_candidates_ids_do_not_collide_on_shared_prefix() -> None:
    """F1(i): truncating sanitized transcript ids at 24 chars used to collide.

    TRINITY_DN1234567_c0_g1_i1 and ..._i2 share their first 24 sanitized characters.
    Pre-fix, bare ``[:24]`` truncation made every enumerated candidate id identical
    between the two transcripts at the same sliding-window position.
    """
    designer = SiRNADesigner(DesignParameters())
    sequence = "ACGT" * 15  # 60 nt: several 21-nt sliding windows

    candidates_1, rejected_1 = designer._enumerate_candidates(sequence, "TRINITY_DN1234567_c0_g1_i1")
    candidates_2, rejected_2 = designer._enumerate_candidates(sequence, "TRINITY_DN1234567_c0_g1_i2")

    ids_1 = [c.id for c in candidates_1 + rejected_1]
    ids_2 = [c.id for c in candidates_2 + rejected_2]
    all_ids = ids_1 + ids_2

    assert set(ids_1).isdisjoint(ids_2), "the two transcripts' candidate ids must never collide"
    assert len(all_ids) == len(set(all_ids)), "every enumerated candidate id must be globally distinct"


@pytest.mark.unit
def test_prepare_offtarget_input_rejects_duplicate_candidate_ids(tmp_path: Path) -> None:
    """F1(ii): a duplicate candidate id must raise, not silently clobber the fan-out map.

    Pre-fix there was no uniqueness check: two candidates sharing an id but carrying
    different sequences would each become their own "representative", and the second
    would silently overwrite the first's entry in representative_to_candidates.
    """
    workflow = _minimal_workflow(tmp_path, "dup_out")

    first = _full_candidate("dup_id", guide_sequence="ATGCGATGCGATGCGATGCGC")
    second = _full_candidate("dup_id", guide_sequence="GGCATTACGCATTACGCATTA")

    with pytest.raises(ValueError, match="dup_id"):
        asyncio.run(workflow._prepare_offtarget_input([first, second]))


@pytest.mark.unit
def test_nm3_off_target_hit_keeps_transcriptome_hits_total_in_agreement(tmp_path: Path) -> None:
    """F3: hits beyond nm=2 must still count toward off_target_count and transcriptome_hits_total.

    Pre-fix, transcriptome_hits_total was ``sum(transcriptome_totals.values())``, and
    transcriptome_totals only has buckets for nm in {0, 1, 2}, so an nm=3 hit inflated
    off_target_count without ever reaching the total -- the two fields silently disagreed.
    """
    workflow = _minimal_workflow(tmp_path, "nm3_out")
    candidate = _full_candidate("cand_nm3")
    hits = [
        {"rname": "ENST00000000098", "species": "human", "nm": 0, "seed_mismatches": 0},
        {"rname": "ENST00000000099", "species": "human", "nm": 3, "seed_mismatches": 3},
    ]
    offtarget_data = {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}}

    workflow._integrate_offtarget_results([candidate], offtarget_data, OffTargetFilterCriteria())

    assert candidate.off_target_count == 2, "both hits are unrelated human transcripts: genuine off-targets"
    assert candidate.transcriptome_hits_total == candidate.off_target_count
    bucketed = candidate.transcriptome_hits_0mm + candidate.transcriptome_hits_1mm + candidate.transcriptome_hits_2mm
    assert bucketed <= candidate.transcriptome_hits_total, "0/1/2mm buckets are a subset of the total, not an addend"
    assert candidate.transcriptome_hits_0mm == 1


@pytest.mark.unit
def test_repeat_flagged_candidate_excluded_when_check_off_targets_disabled(tmp_path: Path) -> None:
    """F4/F5: the repeat exclusion must hold even when the user disabled screening entirely.

    Pre-fix, _apply_post_screen_ranking only ran after a successful Nextflow screening call,
    so a repeat-flagged candidate kept its design-time rank -- and could still land in
    top_candidates -- on every path that never reaches that point.
    """
    params = DesignParameters(check_off_targets=False, apply_modifications=False)
    workflow = _minimal_workflow(tmp_path, "disabled_out", params)

    repeat_candidate = _full_candidate("cand_repeat", guide_sequence="ATGCGATGCGATGCGATGCGC")
    repeat_candidate.composite_score = 99.0
    repeat_candidate.repeat_flagged = True
    clean_candidate = _full_candidate("cand_clean", guide_sequence="GGCATTACGCATTACGCATTA")
    clean_candidate.composite_score = 50.0

    design_result = DesignResult(
        input_file="<test>",
        parameters=params,
        candidates=[repeat_candidate, clean_candidate],
        top_candidates=[repeat_candidate, clean_candidate],
        total_sequences=2,
        total_candidates=2,
        filtered_candidates=2,
        processing_time=0.1,
    )

    result = asyncio.run(workflow.step5_offtarget_analysis(design_result))

    assert result == {"status": "skipped", "reason": "user_disabled"}
    assert repeat_candidate not in design_result.top_candidates
    assert repeat_candidate in design_result.candidates, "excluded from ranking, not dropped from reporting"

    asyncio.run(workflow.step6_generate_reports(design_result))
    all_csv = pd.read_csv(workflow.config.output_dir / "sirnaforge" / "candidates_all.csv")
    assert "cand_repeat" in set(all_csv["id"]), "repeat-flagged candidates must still be reported"


@pytest.mark.unit
def test_repeat_flagged_candidate_excluded_when_nextflow_screening_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """F4/F5: the repeat exclusion must also hold when Nextflow screening itself fails.

    Pre-fix, the only ranking call lived after a successful Nextflow run, so a failed run
    left the repeat-flagged candidate's design-time (pre-exclusion) rank untouched.
    """
    params = DesignParameters(apply_modifications=False)  # check_off_targets defaults True
    workflow = _minimal_workflow(tmp_path, "failed_out", params)

    async def _boom(self: SiRNAWorkflow, *args: object, **kwargs: object) -> dict[str, object]:
        raise RuntimeError("nextflow unavailable in this test")

    monkeypatch.setattr(SiRNAWorkflow, "_run_nextflow_offtarget_analysis", _boom)

    repeat_candidate = _full_candidate("cand_repeat2", guide_sequence="ATGCGATGCGATGCGATGCGC")
    repeat_candidate.composite_score = 99.0
    repeat_candidate.repeat_flagged = True
    clean_candidate = _full_candidate("cand_clean2", guide_sequence="GGCATTACGCATTACGCATTA")
    clean_candidate.composite_score = 50.0

    design_result = DesignResult(
        input_file="<test>",
        parameters=params,
        candidates=[repeat_candidate, clean_candidate],
        top_candidates=[repeat_candidate, clean_candidate],
        total_sequences=2,
        total_candidates=2,
        filtered_candidates=2,
        processing_time=0.1,
    )

    result = asyncio.run(workflow.step5_offtarget_analysis(design_result))

    assert result["status"] == "skipped"
    assert result["reason"] == "nextflow_failed"
    assert repeat_candidate not in design_result.top_candidates
    assert repeat_candidate in design_result.candidates


@pytest.mark.unit
def test_per_species_breakdown_is_populated_for_every_screened_species(tmp_path: Path) -> None:
    """F6: per_species must report real per-species class counts, not a permanently empty dict.

    Pre-fix, stats['per_species'] was hard-coded to {} and never written to, so any
    downstream reader (e.g. offtarget_summary) always saw an empty breakdown regardless
    of how many species were actually screened.
    """
    config = WorkflowConfig(
        output_dir=tmp_path / "species_out",
        gene_query="TP53",
        genome_species=["human", "mouse"],
        design_params=DesignParameters(),
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}

    human_reference = tmp_path / "human_cdna.fasta"
    human_reference.write_text(
        ">ENST00000000001 gene:ENSG00000000001 gene_symbol:TP53 transcript_biotype:protein_coding\n"
        "ACGT\n"
        ">ENST00000000009 gene:ENSG00000000009 gene_symbol:OTHER transcript_biotype:protein_coding\n"
        "ACGT\n"
    )
    # Mouse spells the ortholog in Title Case; the index uppercases symbols on ingest.
    mouse_reference = tmp_path / "mouse_cdna.fasta"
    mouse_reference.write_text(
        ">ENSMUST00000000002 gene:ENSMUSG00000000002 gene_symbol:Tp53 transcript_biotype:protein_coding\nACGT\n"
    )
    workflow._transcript_index.build("human", human_reference)
    workflow._transcript_index.build("mouse", mouse_reference)

    candidate = _full_candidate("cand_species")
    hits = [
        {"rname": "ENST00000000001", "species": "human", "nm": 0, "seed_mismatches": 0},
        {"rname": "ENST00000000009", "species": "human", "nm": 0, "seed_mismatches": 0},
        {"rname": "ENSMUST00000000002", "species": "mouse", "nm": 0, "seed_mismatches": 0},
    ]
    offtarget_data = {"status": "completed", "results": {candidate.id: {"hits": hits, "off_target_score": 0.0}}}

    _, stats = workflow._integrate_offtarget_results([candidate], offtarget_data, OffTargetFilterCriteria())

    empty_bucket = {
        "on_target": 0,
        "ortholog": 0,
        "repeat": 0,
        "off_target": 0,
        "symbol_lookup_missing": 0,
        "species_index_missing": 0,
    }
    assert stats["per_species"]["human"] == {**empty_bucket, "on_target": 1, "off_target": 1}
    assert stats["per_species"]["mouse"] == {**empty_bucket, "ortholog": 1}
