"""Tests for injecting dirty control candidates."""

from __future__ import annotations

from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.utils.control_candidates import DIRTY_CONTROL_LABEL, inject_dirty_controls


def _candidate(idx: int, score: float, sequence: str | None = None) -> SiRNACandidate:
    base_seq = sequence or ("A" * 21)
    return SiRNACandidate(
        id=f"cand{idx}",
        transcript_id="tx1",
        position=idx,
        guide_sequence=base_seq,
        passenger_sequence=base_seq,
        gc_content=40.0,
        length=21,
        asymmetry_score=0.5,
        composite_score=score,
    )


def _rejected(idx: int, sequence: str, reason: SiRNACandidate.FilterStatus) -> SiRNACandidate:
    cand = _candidate(idx, 0.0, sequence)
    cand.passes_filters = reason
    cand.quality_issues.append(reason.value)
    return cand


def _design_result() -> DesignResult:
    candidates = [_candidate(1, 95.0), _candidate(2, 12.5), _candidate(3, 5.0)]
    rejected = [
        _rejected(101, "GGGGGGGGGGGGGGGGGGGGG", SiRNACandidate.FilterStatus.GC_OUT_OF_RANGE),
        _rejected(102, "CCCCCCCCCCCCCCCCCCCCC", SiRNACandidate.FilterStatus.POLY_RUNS),
    ]
    return DesignResult(
        input_file="test.fasta",
        parameters=DesignParameters(top_n=1),
        candidates=candidates,
        top_candidates=[candidates[0]],
        total_sequences=1,
        total_candidates=len(candidates),
        filtered_candidates=1,
        processing_time=0.1,
        tool_versions={},
        rejected_candidates=rejected,
    )


def test_inject_dirty_controls_appends_controls():
    """Dirty controls come directly from rejected candidates and remain labelled."""
    result = _design_result()
    controls = inject_dirty_controls(result)

    assert len(controls) == 2
    assert result.total_candidates == 5
    assert len(result.top_candidates) == 3  # one top + two controls

    control_ids = [c.id for c in controls]
    assert all("__DIRTY_CONTROL_" in cid for cid in control_ids)
    assert all(c.passes_filters.value == DIRTY_CONTROL_LABEL for c in controls)
    rejected_sequences = {rc.guide_sequence for rc in result.rejected_candidates}
    assert {c.guide_sequence for c in controls}.issubset(rejected_sequences)


def test_inject_dirty_controls_is_idempotent():
    """Running the injector twice should not create duplicate control entries."""
    result = _design_result()
    first = inject_dirty_controls(result)
    second = inject_dirty_controls(result)

    assert len(first) == 2
    assert second == []
    # Ensure no duplicate controls were added
    assert sum(1 for c in result.candidates if "__DIRTY_CONTROL_" in c.id) == 2
