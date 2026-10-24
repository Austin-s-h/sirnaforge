"""Strict genomic-to-transcript mapping and allele-aware candidate pool behavior."""

from types import SimpleNamespace

import pytest

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.data.base import DatabaseType, TranscriptInfo
from sirnaforge.data.variant_helpers import apply_variant_workflow, map_variant_to_transcript
from sirnaforge.models.sirna import DesignParameters, DesignResult, FilterCriteria, SiRNACandidate
from sirnaforge.models.transcript_annotation import Interval, TranscriptAnnotation
from sirnaforge.models.variant import VariantMode, VariantRecord


def _annotation(strand: int = 1) -> TranscriptAnnotation:
    exons = [
        Interval(seq_region_name="1", start=100, end=104, strand=strand),
        Interval(seq_region_name="1", start=200, end=204, strand=strand),
    ]
    return TranscriptAnnotation(
        transcript_id="ENST1.2",
        gene_id="ENSG1",
        seq_region_name="1",
        start=100,
        end=204,
        strand=strand,
        exons=exons,
        provider="fixture",
        reference_choice="GRCh38",
    )


def _candidate(position: int, sequence: str, suffix: str = "") -> SiRNACandidate:
    return SiRNACandidate(
        id=f"guide_{position}{suffix}",
        transcript_id="ENST1.2",
        position=position,
        guide_sequence="A" * 21,
        passenger_sequence=sequence[position - 1 : position + 20],
        gc_content=40.0,
        length=21,
        asymmetry_score=0.5,
        design_score=50.0,
    )


def _result(candidates: list[SiRNACandidate]) -> DesignResult:
    return DesignResult(
        input_file="fixture",
        parameters=DesignParameters(),
        candidates=candidates,
        top_candidates=list(candidates),
        total_sequences=1,
        total_candidates=len(candidates),
        filtered_candidates=len(candidates),
        processing_time=0.0,
    )


@pytest.mark.parametrize(("strand", "position", "expected"), [(1, 201, 6), (-1, 103, 6)])
def test_maps_spliced_variant_by_exon_order_and_strand(strand: int, position: int, expected: int) -> None:
    """Map genomic variants through spliced exons on forward and reverse strands."""
    # The variant is the second cDNA base of exon two after the five-base first exon.
    variant = VariantRecord(chr="chr1", pos=position, ref="C", alt="T")
    assert map_variant_to_transcript(variant, _annotation(strand)) == expected


def test_intronic_or_unmapped_variant_fails_explicitly() -> None:
    """Reject variants without an annotated exon mapping."""
    annotation = _annotation()
    variant = VariantRecord(chr="chr1", pos=150, ref="C", alt="T")
    assert map_variant_to_transcript(variant, annotation) is None

    transcript = SimpleNamespace(transcript_id="ENST1.2", sequence="A" * 30)
    with pytest.raises(ValueError, match="mapping is unknown"):
        apply_variant_workflow(_result([]), [variant], [transcript], {"ENST1.2": annotation}, "avoid", lambda *_: None)


def test_avoid_removes_overlapping_reference_candidate() -> None:
    """Remove guides whose transcript interval overlaps a mapped SNV."""
    variant = VariantRecord(chr="chr1", pos=201, ref="A", alt="T")
    sequence = "A" * 30
    transcript = SimpleNamespace(transcript_id="ENST1.2", sequence=sequence)
    overlapping = _candidate(6, sequence)
    clear = _candidate(10, sequence)

    result = apply_variant_workflow(
        _result([overlapping, clear]),
        [variant],
        [transcript],
        {"ENST1.2": _annotation()},
        VariantMode.AVOID.value,
        lambda *_: pytest.fail("avoid mode must not design alternate sequences"),
    )

    assert result.candidates == [clear]


def test_both_keeps_reference_and_scores_alternate_full_context() -> None:
    """Keep reference guides and pass a full edited transcript to the designer."""
    variant = VariantRecord(chr="chr1", pos=201, ref="A", alt="T")
    sequence = "A" * 30
    transcript = SimpleNamespace(transcript_id="ENST1.2", sequence=sequence)
    reference = _candidate(6, sequence)
    designed_sequences: list[str] = []

    def design_alt(edited_sequence: str, transcript_id: str) -> SimpleNamespace:
        designed_sequences.append(edited_sequence)
        return SimpleNamespace(candidates=[_candidate(6, edited_sequence, "_generated")])

    result = apply_variant_workflow(
        _result([reference]),
        [variant],
        [transcript],
        {"ENST1": _annotation()},
        VariantMode.BOTH.value,
        design_alt,
    )

    assert designed_sequences == ["A" * 6 + "T" + "A" * 23]
    assert len(result.candidates) == 2
    assert result.candidates[0].targeted_alleles == ["ref"]
    assert result.candidates[1].targeted_alleles == ["alt"]
    assert result.candidates[1].passenger_sequence[1] == "T"


def test_actual_designer_generates_target_allele_guides_and_avoid_excludes_site() -> None:
    """Run alternate full-context design, then check target and avoid pools."""
    sequence = "ACGT" * 50
    transcript = TranscriptInfo(
        transcript_id="ENST1.2",
        gene_id="ENSG1",
        gene_name="FIXTURE",
        sequence=sequence,
        database=DatabaseType.ENSEMBL,
    )
    annotation = TranscriptAnnotation(
        transcript_id="ENST1.2",
        gene_id="ENSG1",
        seq_region_name="1",
        start=100,
        end=299,
        strand=1,
        exons=[Interval(seq_region_name="1", start=100, end=299, strand=1)],
        provider="offline_fixture",
        reference_choice="GRCh38",
    )
    variant = VariantRecord(chr="chr1", pos=120, ref="A", alt="G")
    params = DesignParameters(
        filters=FilterCriteria(
            gc_min=0,
            gc_max=100,
            max_poly_runs=100,
            max_repeat_transcript_fraction=1,
            max_paired_fraction=1,
            min_asymmetry_score=0,
            min_empirical_score=0.4,
        )
    )
    designer = SiRNADesigner(params)
    reference_result = designer.design_from_sequence(sequence, transcript.transcript_id)
    reference_at_site = {
        candidate.position: candidate
        for candidate in reference_result.candidates
        if candidate.position - 1 <= 20 < candidate.position - 1 + candidate.length
    }

    target_result = apply_variant_workflow(
        designer.design_from_sequence(sequence, transcript.transcript_id),
        [variant],
        [transcript],
        {transcript.transcript_id: annotation},
        VariantMode.TARGET,
        designer.design_from_sequence,
    )
    assert target_result.candidates
    assert all(candidate.targeted_alleles == ["alt"] for candidate in target_result.candidates)
    assert all(candidate.variant_mode == "target" for candidate in target_result.candidates)
    assert all(
        candidate.position - 1 <= 20 < candidate.position - 1 + candidate.length
        for candidate in target_result.candidates
    )
    assert all(
        candidate.passenger_sequence[20 - (candidate.position - 1)] == "G" for candidate in target_result.candidates
    )
    for candidate in target_result.candidates:
        calculated_gc = 100 * sum(base in "GC" for base in candidate.passenger_sequence) / candidate.length
        assert candidate.gc_content == pytest.approx(calculated_gc)
    paired_positions = set(reference_at_site).intersection(c.position for c in target_result.candidates)
    assert any(
        target.gc_content != reference_at_site[target.position].gc_content
        and target.score_gc_content != reference_at_site[target.position].score_gc_content
        for target in target_result.candidates
        if target.position in paired_positions
    )

    both_result = apply_variant_workflow(
        designer.design_from_sequence(sequence, transcript.transcript_id),
        [variant],
        [transcript],
        {transcript.transcript_id: annotation},
        VariantMode.BOTH,
        designer.design_from_sequence,
    )
    both_by_position: dict[int, dict[str, SiRNACandidate]] = {}
    for candidate in both_result.candidates:
        if candidate.position - 1 <= 20 < candidate.position - 1 + candidate.length:
            allele = candidate.targeted_alleles[0]
            both_by_position.setdefault(candidate.position, {})[allele] = candidate
    assert both_by_position
    for alleles in both_by_position.values():
        assert set(alleles) == {"ref", "alt"}
        relative_position = 20 - (alleles["ref"].position - 1)
        assert alleles["ref"].passenger_sequence[relative_position] == "A"
        assert alleles["alt"].passenger_sequence[relative_position] == "G"
        assert alleles["ref"].guide_sequence != alleles["alt"].guide_sequence
        assert alleles["ref"].gc_content != alleles["alt"].gc_content
        assert alleles["ref"].score_gc_content != alleles["alt"].score_gc_content

    avoid_result = apply_variant_workflow(
        designer.design_from_sequence(sequence, transcript.transcript_id),
        [variant],
        [transcript],
        {transcript.transcript_id: annotation},
        VariantMode.AVOID,
        designer.design_from_sequence,
    )
    assert avoid_result.candidates
    assert all(
        not (candidate.position - 1 <= 20 < candidate.position - 1 + candidate.length)
        for candidate in avoid_result.candidates
    )
