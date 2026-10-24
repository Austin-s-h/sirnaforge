"""Resolve target intent against annotation and evaluate exact transcript coverage."""

import re
from collections.abc import Collection, Sequence
from enum import Enum

from pydantic import BaseModel, ConfigDict, model_validator

from sirnaforge.data.base import TranscriptInfo
from sirnaforge.models.policy import FilterEvaluation, TargetIntent, TargetSelectivity


class TargetIntentRequest(BaseModel):
    """Caller-declared intent; annotation universe and denominator are resolver-owned."""

    selectivity: TargetSelectivity
    target_species: frozenset[str] = frozenset({"human"})
    offtarget_screen_species: frozenset[str] = frozenset()
    required_transcript_ids: frozenset[str] = frozenset()
    excluded_transcript_ids: frozenset[str] = frozenset()

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def declared_ids_do_not_overlap(self) -> "TargetIntentRequest":
        """Reject contradictory user declarations before annotation resolution."""
        if self.required_transcript_ids & self.excluded_transcript_ids:
            raise ValueError("transcript IDs cannot be both required and excluded")
        return self


class ResolvedTranscript(BaseModel):
    """One annotation-universe member and its optional full-length cDNA sequence."""

    transcript_id: str
    gene_id: str | None
    sequence: str | None

    model_config = ConfigDict(frozen=True, extra="forbid")


class ResolvedTargetIntent(BaseModel):
    """Frozen intent and complete transcript universe used for guide eligibility."""

    intent: TargetIntent
    transcripts: tuple[ResolvedTranscript, ...]

    model_config = ConfigDict(frozen=True, extra="forbid")


class TargetIntentVerdict(str, Enum):
    """Three-valued guide eligibility result."""

    PASS = "pass"
    FAIL = "fail"
    UNKNOWN = "unknown"


class TargetIntentEvaluation(BaseModel):
    """Coverage decision and per-transcript reasons for one guide sequence."""

    verdict: TargetIntentVerdict
    covered_transcript_ids: frozenset[str]
    uncovered_required_transcript_ids: frozenset[str]
    excluded_hit_transcript_ids: frozenset[str]
    missing_sequence_transcript_ids: frozenset[str]

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def filter_evaluation(self) -> FilterEvaluation:
        """Adapt to the shared filter verdict vocabulary."""
        return FilterEvaluation(self.verdict.value)


def _versionless(identifier: str) -> str:
    """Drop a terminal numeric transcript version while preserving stable-ID spelling."""
    return re.sub(r"\.\d+$", "", identifier.strip())


def _normalize_sequence(sequence: str, *, identifier: str) -> str:
    """Normalize RNA/DNA cDNA to DNA alphabet and reject ambiguity."""
    normalized = sequence.upper().replace("U", "T")
    invalid = sorted(set(normalized) - set("ACGT"))
    if not normalized or invalid:
        reason = f"unsupported nucleotide(s): {', '.join(invalid)}" if invalid else "empty sequence"
        raise ValueError(f"transcript {identifier} has {reason}")
    return normalized


def _index_annotation(
    annotation_transcript_ids: Collection[str],
    transcript_records: Sequence[TranscriptInfo],
) -> tuple[dict[str, str], dict[str, TranscriptInfo]]:
    """Index annotation members and producer records, rejecting collapsed identities."""
    universe: dict[str, str] = {}
    for raw_id in annotation_transcript_ids:
        identifier = _versionless(raw_id)
        if not identifier:
            raise ValueError("annotation transcript IDs must be non-empty")
        if identifier in universe:
            raise ValueError(f"duplicate/version conflict in annotation transcript IDs for {identifier}")
        universe[identifier] = raw_id
    if not universe:
        raise ValueError("annotation transcript universe is empty")

    records_by_id: dict[str, TranscriptInfo] = {}
    for record in transcript_records:
        identifier = _versionless(record.transcript_id)
        if identifier in records_by_id:
            raise ValueError(f"duplicate/version conflict in transcript records for {identifier}")
        if identifier not in universe:
            raise ValueError(f"transcript record {record.transcript_id} is outside the annotation universe")
        records_by_id[identifier] = record
    return universe, records_by_id


def resolve_target_intent(
    request: TargetIntentRequest,
    transcript_records: Sequence[TranscriptInfo],
    *,
    annotation_transcript_ids: Collection[str],
    annotation_provenance: str,
    enumeration_input_ids: Collection[str],
) -> ResolvedTargetIntent:
    """Build an intent using the complete annotation set before design narrows enumeration.

    All IDs are stored versionless. Every annotation ID remains in ``transcripts`` and in the
    frozen denominator even when it was omitted from ``enumeration_input_ids`` or has no sequence.
    """
    if not annotation_provenance.strip():
        raise ValueError("annotation_provenance must identify the annotation source")
    normalized_universe, records_by_id = _index_annotation(annotation_transcript_ids, transcript_records)

    def resolve_declared_ids(values: Collection[str], label: str) -> frozenset[str]:
        resolved = frozenset(_versionless(value) for value in values)
        if len(resolved) != len(values):
            raise ValueError(f"duplicate/version conflict in {label} transcript IDs")
        unknown = resolved - normalized_universe.keys()
        if unknown:
            raise ValueError(f"{label} transcript IDs are absent from annotation universe: {', '.join(sorted(unknown))}")
        return resolved

    required_declared = resolve_declared_ids(request.required_transcript_ids, "required")
    excluded = resolve_declared_ids(request.excluded_transcript_ids, "excluded")
    if required_declared & excluded:
        raise ValueError("transcript IDs cannot be both required and excluded after version normalization")
    enumeration_inputs = resolve_declared_ids(enumeration_input_ids, "enumeration")
    if not enumeration_inputs:
        raise ValueError("at least one enumeration input transcript is required")
    required = frozenset(normalized_universe) - excluded if request.selectivity is TargetSelectivity.PAN_ISOFORM else required_declared
    if not required:
        raise ValueError("target intent resolves to no required transcripts")
    if request.selectivity is TargetSelectivity.ISOFORM_SELECTIVE and not required:
        raise ValueError("isoform-selective intent requires at least one required transcript")

    transcripts: list[ResolvedTranscript] = []
    for identifier in sorted(normalized_universe):
        transcript_record = records_by_id.get(identifier)
        sequence = (
            None
            if transcript_record is None or transcript_record.sequence is None
            else _normalize_sequence(transcript_record.sequence, identifier=identifier)
        )
        transcripts.append(
            ResolvedTranscript(
                transcript_id=identifier,
                gene_id=None if transcript_record is None else transcript_record.gene_id,
                sequence=sequence,
            )
        )
    denominator = required
    intent = TargetIntent(
        selectivity=request.selectivity,
        target_species=request.target_species,
        offtarget_screen_species=request.offtarget_screen_species,
        required_transcript_ids=required,
        excluded_transcript_ids=excluded,
        annotation_universe=frozenset(normalized_universe),
        coverage_denominator=denominator,
        enumeration_inputs=enumeration_inputs,
        annotation_provenance=annotation_provenance.strip(),
    )
    return ResolvedTargetIntent(intent=intent, transcripts=tuple(transcripts))


def evaluate_target_intent_guide(guide_sequence: str, resolved: ResolvedTargetIntent) -> TargetIntentEvaluation:
    """Evaluate exact full-length guide reverse-complement coverage in the frozen universe."""
    guide = _normalize_sequence(guide_sequence, identifier="guide")
    reverse_complement = guide.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sequences = {transcript.transcript_id: transcript.sequence for transcript in resolved.transcripts}
    covered = frozenset(
        identifier
        for identifier, sequence in sequences.items()
        if sequence is not None and reverse_complement in sequence
    )
    required = resolved.intent.required_transcript_ids
    excluded = resolved.intent.excluded_transcript_ids
    selective_liabilities = (
        resolved.intent.annotation_universe - required
        if TargetSelectivity(resolved.intent.selectivity) is TargetSelectivity.ISOFORM_SELECTIVE
        else frozenset()
    )
    forbidden = excluded | selective_liabilities
    uncovered_required = required - covered
    excluded_hits = forbidden & covered
    missing = frozenset(
        identifier
        for identifier in required | forbidden
        if sequences.get(identifier) is None
    )
    if excluded_hits or uncovered_required - missing:
        verdict = TargetIntentVerdict.FAIL
    elif missing:
        verdict = TargetIntentVerdict.UNKNOWN
    else:
        verdict = TargetIntentVerdict.PASS
    return TargetIntentEvaluation(
        verdict=verdict,
        covered_transcript_ids=covered,
        uncovered_required_transcript_ids=uncovered_required,
        excluded_hit_transcript_ids=excluded_hits,
        missing_sequence_transcript_ids=missing,
    )
