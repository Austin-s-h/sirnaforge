"""Synthetic fixtures for frozen target intent and exact transcript coverage."""

import pytest
from pydantic import ValidationError

from sirnaforge.core.target_intent import (
    TargetIntentRequest,
    TargetIntentVerdict,
    evaluate_target_intent_guide,
    resolve_target_intent,
)
from sirnaforge.data.base import DatabaseType, TranscriptInfo
from sirnaforge.models.policy import TargetSelectivity

GUIDE = "AACCTGTAAG"
TARGET = GUIDE.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _transcript(transcript_id: str, sequence: str | None) -> TranscriptInfo:
    """Build one producer transcript record."""
    return TranscriptInfo(
        transcript_id=transcript_id,
        gene_id="ENSG1",
        sequence=sequence,
        database=DatabaseType.ENSEMBL,
    )


def _resolve(
    selectivity: TargetSelectivity,
    records: list[TranscriptInfo],
    *,
    annotation_ids: list[str],
    enumeration_ids: list[str],
    required: frozenset[str] = frozenset(),
    excluded: frozenset[str] = frozenset(),
):
    """Resolve a compact synthetic annotation fixture."""
    return resolve_target_intent(
        TargetIntentRequest(
            selectivity=selectivity,
            target_species=frozenset({"human", "mouse"}),
            offtarget_screen_species=frozenset({"rat"}),
            required_transcript_ids=required,
            excluded_transcript_ids=excluded,
        ),
        records,
        annotation_transcript_ids=annotation_ids,
        annotation_provenance="test-annotation/v1",
        enumeration_input_ids=enumeration_ids,
    )


def test_pan_isoform_freezes_full_annotation_denominator_before_enumeration_subset() -> None:
    """Pan-isoform requires every annotated transcript, including unenumerated isoforms."""
    resolved = _resolve(
        TargetSelectivity.PAN_ISOFORM,
        [_transcript("ENST1.2", f"AAA{TARGET}AAA"), _transcript("ENST2.1", "CCCGGGCCC")],
        annotation_ids=["ENST1.2", "ENST2.1"],
        enumeration_ids=["ENST1.2"],
    )
    assert resolved.intent.annotation_universe == frozenset({"ENST1", "ENST2"})
    assert resolved.intent.coverage_denominator == frozenset({"ENST1", "ENST2"})
    assert resolved.intent.enumeration_inputs == frozenset({"ENST1"})
    result = evaluate_target_intent_guide(GUIDE, resolved)
    assert result.verdict is TargetIntentVerdict.FAIL
    assert result.uncovered_required_transcript_ids == frozenset({"ENST2"})
    assert resolved.intent.target_species == frozenset({"human", "mouse"})
    assert resolved.intent.offtarget_screen_species == frozenset({"rat"})


def test_isoform_selective_can_pass_when_pan_isoform_fails_for_other_isoforms() -> None:
    """A selective guide may cover its required isoform without covering other isoforms."""
    resolved = _resolve(
        TargetSelectivity.ISOFORM_SELECTIVE,
        [_transcript("ENST1.4", f"AAA{TARGET}AAA"), _transcript("ENST2.1", "CCCGGGCCC")],
        annotation_ids=["ENST1.4", "ENST2.1"],
        enumeration_ids=["ENST1.4"],
        required=frozenset({"ENST1"}),
    )
    result = evaluate_target_intent_guide(GUIDE, resolved)
    assert result.verdict is TargetIntentVerdict.PASS
    assert resolved.intent.coverage_denominator == frozenset({"ENST1"})


def test_missing_required_sequence_yields_unknown_not_pass() -> None:
    """A denominator member without a sequence leaves pan coverage unknown."""
    resolved = _resolve(
        TargetSelectivity.PAN_ISOFORM,
        [_transcript("ENST1", f"{TARGET}AAA")],
        annotation_ids=["ENST1", "ENST2"],
        enumeration_ids=["ENST1"],
    )
    result = evaluate_target_intent_guide(GUIDE, resolved)
    assert result.verdict is TargetIntentVerdict.UNKNOWN
    assert result.missing_sequence_transcript_ids == frozenset({"ENST2"})


def test_excluded_hit_is_detected_even_when_transcript_was_not_enumerated() -> None:
    """Eligibility searches the frozen universe, not only enumeration inputs."""
    resolved = _resolve(
        TargetSelectivity.ISOFORM_SELECTIVE,
        [_transcript("ENST1", f"{TARGET}AAA"), _transcript("ENST2", f"CCC{TARGET}CCC")],
        annotation_ids=["ENST1", "ENST2"],
        enumeration_ids=["ENST1"],
        required=frozenset({"ENST1"}),
        excluded=frozenset({"ENST2"}),
    )
    result = evaluate_target_intent_guide(GUIDE, resolved)
    assert result.verdict is TargetIntentVerdict.FAIL
    assert result.excluded_hit_transcript_ids == frozenset({"ENST2"})


def test_version_conflicts_and_caller_fabricated_denominator_are_rejected() -> None:
    """Conflicting transcript identities fail, and request config cannot set resolved fields."""
    with pytest.raises(ValueError, match="duplicate/version conflict"):
        _resolve(
            TargetSelectivity.PAN_ISOFORM,
            [_transcript("ENST1.1", f"{TARGET}"), _transcript("ENST1.2", f"{TARGET}")],
            annotation_ids=["ENST1.1", "ENST1.2"],
            enumeration_ids=["ENST1.1"],
        )
    with pytest.raises(ValidationError):
        TargetIntentRequest(
            selectivity=TargetSelectivity.PAN_ISOFORM,
            annotation_universe=frozenset({"caller-fabricated"}),
        )
    with pytest.raises(ValueError, match="duplicate/version conflict"):
        _resolve(
            TargetSelectivity.ISOFORM_SELECTIVE,
            [_transcript("ENST1.1", f"{TARGET}")],
            annotation_ids=["ENST1.1"],
            enumeration_ids=["ENST1.1"],
            required=frozenset({"ENST1.1", "ENST1.2"}),
        )


def test_rna_and_dna_sequences_normalize_to_same_exact_coverage() -> None:
    """Guide and cDNA T/U spellings do not change full-length complement detection."""
    resolved = _resolve(
        TargetSelectivity.PAN_ISOFORM,
        [_transcript("ENST1", f"{TARGET.replace('T', 'U')}AAA")],
        annotation_ids=["ENST1"],
        enumeration_ids=["ENST1"],
    )
    assert evaluate_target_intent_guide(GUIDE.replace("T", "U"), resolved).verdict is TargetIntentVerdict.PASS
