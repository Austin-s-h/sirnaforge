"""Tests for the off-target screening flag.

Hit counts default to 0, so without an explicit flag a candidate that was never
submitted to off-target analysis is indistinguishable from one that was screened
and came back clean (issue #78 §3).
"""

import pytest

from sirnaforge.models.sirna import DesignParameters, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

GUIDE = "ATGCGATGCGATGCGATGCGC"
PASSENGER = "GCGCATCGCATCGCATCGCAT"


def _candidate(candidate_id: str) -> SiRNACandidate:
    """Build a candidate with default (unscreened) off-target fields."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=PASSENGER,
        gc_content=57.1,
        length=21,
        asymmetry_score=0.7,
        composite_score=50.0,
        # Populated so save_csv sees float columns rather than all-null ones.
        mfe=-4.2,
        duplex_stability=-39.0,
        structure="." * 21,
        component_scores={
            "duplex_stability_score": 0.6,
            "dg_5p": -8.0,
            "dg_3p": -9.0,
            "delta_dg_end": 1.0,
            "melting_temp_c": 70.0,
        },
    )


def _workflow(tmp_path) -> SiRNAWorkflow:
    """Build a workflow instance without running it."""
    config = WorkflowConfig(output_dir=tmp_path / "out", gene_query="tp53", design_params=DesignParameters())
    return SiRNAWorkflow(config)


@pytest.mark.unit
def test_candidates_start_unscreened():
    """A freshly designed candidate has not been screened."""
    assert _candidate("c1").off_target_screened is False


@pytest.mark.unit
def test_integration_marks_screened_candidates_including_clean_ones(tmp_path):
    """A submitted candidate is screened even when the pipeline reports no hits."""
    workflow = _workflow(tmp_path)
    with_hits, clean = _candidate("c_hits"), _candidate("c_clean")
    offtarget_data = {
        "status": "completed",
        "results": {"c_hits": {"hits": [{"nm": 1, "seed_mismatches": 2, "offtarget_score": 20.0}]}},
    }

    updated, _stats = workflow._integrate_offtarget_results([with_hits, clean], offtarget_data)

    assert all(c.off_target_screened for c in updated)
    assert clean.transcriptome_hits_total == 0, "clean candidate should keep zero hits"
    assert with_hits.transcriptome_hits_1mm == 1


@pytest.mark.unit
def test_failed_run_leaves_candidates_unscreened(tmp_path):
    """A skipped or failed off-target run must not claim the candidates were screened."""
    workflow = _workflow(tmp_path)
    candidate = _candidate("c1")

    workflow._integrate_offtarget_results([candidate], {"status": "skipped", "reason": "nextflow_failed"})

    assert candidate.off_target_screened is False


@pytest.mark.unit
def test_screening_flag_reaches_the_candidate_csv(tmp_path):
    """The flag is exported so downstream readers can tell 0 from unknown."""
    from sirnaforge.models.sirna import DesignResult  # noqa: PLC0415 - keeps the module import list minimal

    screened, unscreened = _candidate("c_screened"), _candidate("c_unscreened")
    screened.off_target_screened = True
    result = DesignResult(
        input_file="<test>",
        parameters=DesignParameters(),
        candidates=[screened, unscreened],
        top_candidates=[screened],
        total_sequences=1,
        total_candidates=2,
        filtered_candidates=2,
        processing_time=0.0,
    )

    frame = result.save_csv(str(tmp_path / "candidates.csv"))

    assert list(frame["off_target_screened"]) == [True, False]
