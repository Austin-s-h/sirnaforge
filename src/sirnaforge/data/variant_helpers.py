"""Helpers for applying genomic variants to transcript-oriented siRNA design."""

from collections.abc import Callable, Mapping, Sequence

from Bio.Seq import Seq

from sirnaforge.data.base import TranscriptInfo
from sirnaforge.data.variant_resolution_helpers import normalize_chromosome
from sirnaforge.models.sirna import (
    DesignResult,
    SiRNACandidate,
    ranking_order_key,
    unique_guide_representatives,
)
from sirnaforge.models.transcript_annotation import TranscriptAnnotation
from sirnaforge.models.variant import VariantMode, VariantRecord
from sirnaforge.utils.ensembl_ids import strip_version
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


def map_variant_to_transcript(
    variant: VariantRecord,
    annotation: TranscriptAnnotation,
) -> int | None:
    """Map a genomic SNV to a zero-based cDNA index using transcript exons.

    Returns ``None`` when chromosome or exon membership does not match. A genomic
    coordinate is never used directly as an index into transcript sequence.
    """
    if len(variant.ref) != 1 or len(variant.alt) != 1:
        raise ValueError(f"Only SNVs are supported for allele-specific design: {variant.to_vcf_style()}")
    if variant.assembly != "GRCh38":
        raise ValueError(f"Unsupported variant assembly {variant.assembly!r}: {variant.to_vcf_style()}")
    if annotation.strand not in {-1, 1}:
        raise ValueError(f"Unknown transcript strand for {annotation.transcript_id}")
    reference = (annotation.reference_choice or "").upper()
    if not reference.startswith("GRCH38"):
        raise ValueError(
            f"Unknown or unsupported annotation assembly {annotation.reference_choice!r} for {annotation.transcript_id}"
        )

    variant_chromosome = normalize_chromosome(variant.chr)
    annotation_chromosome = normalize_chromosome(annotation.chr)
    if variant_chromosome is None or annotation_chromosome is None or variant_chromosome != annotation_chromosome:
        return None

    exons = [exon for exon in annotation.exons if exon.strand in {None, annotation.strand}]
    if not exons:
        raise ValueError(f"Transcript {annotation.transcript_id} has no usable exon annotation")
    exons.sort(key=lambda exon: exon.start, reverse=annotation.strand == -1)

    offset = 0
    for exon in exons:
        if exon.start <= variant.pos <= exon.end:
            within_exon = variant.pos - exon.start if annotation.strand == 1 else exon.end - variant.pos
            return offset + within_exon
        offset += exon.end - exon.start + 1
    return None


def _transcript_allele(allele: str, strand: int) -> str:
    """Convert a genomic allele to transcript orientation."""
    return allele.upper() if strand == 1 else str(Seq(allele.upper()).reverse_complement())


def apply_variant_to_transcript_sequence(
    sequence: str,
    variant: VariantRecord,
    annotation: TranscriptAnnotation,
) -> tuple[str, int | None]:
    """Apply a strictly mapped GRCh38 variant to cDNA, returning its cDNA position.

    A non-overlapping transcript is returned unchanged with ``None``. Mapped variants must match
    the transcript-oriented reference base; this function never treats a genomic coordinate as a
    sequence offset.
    """
    position = map_variant_to_transcript(variant, annotation)
    if position is None:
        return sequence.upper().replace("U", "T"), None
    normalized = sequence.upper().replace("U", "T")
    if position >= len(normalized):
        raise ValueError(
            f"Mapped variant {variant.to_vcf_style()} exceeds transcript sequence "
            f"{annotation.transcript_id} length ({len(normalized)})"
        )
    ref = _transcript_allele(variant.ref, annotation.strand)
    actual = normalized[position : position + 1]
    if actual != ref:
        raise ValueError(
            f"Reference mismatch for {variant.to_vcf_style()} on {annotation.transcript_id} at cDNA index "
            f"{position}: annotation maps {ref}, sequence contains {actual}"
        )
    edited = normalized[:position] + _transcript_allele(variant.alt, annotation.strand) + normalized[position + 1 :]
    return edited, position


def apply_variant_workflow(  # noqa: PLR0912
    design_result: DesignResult,
    variants: Sequence[VariantRecord],
    transcripts: Sequence[TranscriptInfo],
    annotations: Mapping[str, TranscriptAnnotation],
    mode: str | VariantMode,
    design_sequence: Callable[[str, str], DesignResult],
) -> DesignResult:
    """Filter reference guides and add allele-specific SNV guides to a design result.

    ``avoid`` removes candidates overlapping mapped SNVs. ``target`` replaces
    overlapping reference candidates with alternate-allele designs. ``both`` keeps
    overlapping reference designs and adds alternate-allele designs. A variant that
    cannot be mapped to any annotated exon fails with an explicit unknown-mapping
    error instead of silently interpreting its genomic position as cDNA.
    """
    mode = mode.value if isinstance(mode, VariantMode) else mode.lower()
    if mode not in {"avoid", "target", "both"}:
        raise ValueError(f"Invalid variant mode: {mode}")

    transcript_by_id = {str(t.transcript_id): t for t in transcripts}
    mapped: dict[str, list[tuple[VariantRecord, int]]] = {}
    for variant in variants:
        count = 0
        for transcript_id, transcript in transcript_by_id.items():
            annotation = annotations.get(transcript_id) or annotations.get(strip_version(transcript_id))
            if annotation is None:
                annotation = next(
                    (a for key, a in annotations.items() if strip_version(key) == strip_version(transcript_id)),
                    None,
                )
            if annotation is None:
                raise ValueError(
                    f"Variant mapping is unknown for {variant.to_vcf_style()}: "
                    f"transcript {transcript_id} has no genomic annotation"
                )
            position = map_variant_to_transcript(variant, annotation)
            if position is None:
                continue
            if not transcript.sequence:
                raise ValueError(f"Transcript sequence is unknown for variant mapping: {transcript_id}")
            sequence = transcript.sequence.upper().replace("U", "T")
            if position >= len(sequence):
                raise ValueError(
                    f"Mapped variant {variant.to_vcf_style()} exceeds transcript sequence "
                    f"{transcript_id} length ({len(sequence)})"
                )
            ref = _transcript_allele(variant.ref, annotation.strand)
            actual = sequence[position : position + 1]
            if actual != ref:
                raise ValueError(
                    f"Reference mismatch for {variant.to_vcf_style()} on {transcript_id} at cDNA index "
                    f"{position}: annotation maps {ref}, sequence contains {actual}"
                )
            mapped.setdefault(transcript_id, []).append((variant, position))
            count += 1
        if count == 0:
            raise ValueError(
                f"Variant mapping is unknown for {variant.to_vcf_style()}: no annotated transcript exon overlaps it"
            )

    def overlaps(candidate: SiRNACandidate, position: int) -> bool:
        start = candidate.position - 1
        return start <= position < start + candidate.length

    original = list(design_result.candidates)
    if mode == "avoid":
        retained = [
            candidate
            for candidate in original
            if not any(overlaps(candidate, position) for _, position in mapped.get(candidate.transcript_id, []))
        ]
        design_result.candidates = retained
    elif mode == "target":
        # Target mode is variant-specific: return only alternate-allele guides
        # spanning a mapped variant, never generic guides from the reference pool.
        design_result.candidates = []
    else:
        for candidate in original:
            for variant, position in mapped.get(candidate.transcript_id, []):
                if overlaps(candidate, position):
                    annotate_candidate_with_variant(candidate, variant, "ref", mode)
        design_result.candidates = original

    if mode in {"target", "both"}:
        for transcript_id, transcript in transcript_by_id.items():
            variant_positions = mapped.get(transcript_id, [])
            if not variant_positions:
                continue
            assert transcript.sequence is not None
            sequence = transcript.sequence.upper().replace("U", "T")
            for variant, position in variant_positions:
                annotation = annotations.get(transcript_id) or next(
                    a for key, a in annotations.items() if strip_version(key) == strip_version(transcript_id)
                )
                edited, applied_position = apply_variant_to_transcript_sequence(sequence, variant, annotation)
                if applied_position != position:
                    raise RuntimeError("variant position changed during strict transcript edit")
                # Score this allele on its own full transcript context. Combining nearby
                # variants here would invent a haplotype the input did not specify.
                alt_result = design_sequence(edited, transcript_id)
                for candidate in alt_result.candidates:
                    if not overlaps(candidate, position):
                        continue
                    candidate.id = f"{candidate.id}_ALT_{variant.chr}_{variant.pos}"
                    annotate_candidate_with_variant(candidate, variant, "alt", mode)
                    design_result.candidates.append(candidate)

    # Re-rank the actual output pool and keep its summary counters coherent.
    design_result.candidates.sort(key=ranking_order_key)
    passing = [
        candidate
        for candidate in design_result.candidates
        if candidate.passes_filters is True or candidate.passes_filters == SiRNACandidate.FilterStatus.PASS
    ]
    design_result.top_candidates = unique_guide_representatives(passing)[: design_result.parameters.top_n]
    design_result.total_candidates = len(design_result.candidates)
    design_result.filtered_candidates = len(passing)
    return design_result


def generate_contexts_for_variant(
    variant: VariantRecord,
    reference_sequence: str,
    transcript_start: int,
    transcript_id: str,
    flank_size: int = 50,
) -> dict[str, tuple[str, int, int]]:
    """Generate reference and alternate sequence contexts for a variant.

    Args:
        variant: Variant record with position and alleles
        reference_sequence: Full reference transcript sequence
        transcript_start: Genomic start position of the transcript (1-based)
        transcript_id: Transcript identifier
        flank_size: Number of nucleotides to include on each side of variant (default: 50)

    Returns:
        Dictionary with 'ref' and 'alt' keys, each containing:
        - sequence context (str)
        - relative position of variant in context (int, 0-based)
        - length of the allele (int)

    Raises:
        ValueError: If variant position is outside transcript boundaries
    """
    # Calculate variant position within transcript (0-based)
    # variant.pos is 1-based genomic position
    # transcript_start is 1-based genomic start of transcript
    variant_transcript_pos = variant.pos - transcript_start

    if variant_transcript_pos < 0 or variant_transcript_pos >= len(reference_sequence):
        raise ValueError(
            f"Variant position {variant.pos} is outside transcript {transcript_id} "
            f"boundaries (start: {transcript_start}, length: {len(reference_sequence)})"
        )

    # Extract reference context
    start = max(0, variant_transcript_pos - flank_size)
    end = min(len(reference_sequence), variant_transcript_pos + len(variant.ref) + flank_size)

    ref_context = reference_sequence[start:end]
    ref_relative_pos = variant_transcript_pos - start

    # Generate alternate context by replacing ref allele with alt allele
    alt_context = (
        reference_sequence[start:variant_transcript_pos]
        + variant.alt
        + reference_sequence[variant_transcript_pos + len(variant.ref) : end]
    )
    alt_relative_pos = variant_transcript_pos - start

    logger.debug(
        f"Generated contexts for variant {variant.to_vcf_style()}: "
        f"ref_context_len={len(ref_context)}, alt_context_len={len(alt_context)}, "
        f"ref_pos={ref_relative_pos}, alt_pos={alt_relative_pos}"
    )

    return {
        "ref": (ref_context, ref_relative_pos, len(variant.ref)),
        "alt": (alt_context, alt_relative_pos, len(variant.alt)),
    }


def check_candidate_overlaps_variant(
    candidate_pos: int,
    candidate_length: int,
    variant: VariantRecord,
    transcript_start: int,
) -> bool:
    """Check if a siRNA candidate overlaps with a variant position.

    Args:
        candidate_pos: 1-based start position of candidate in transcript
        candidate_length: Length of the siRNA candidate
        variant: Variant record to check
        transcript_start: 1-based genomic start position of transcript

    Returns:
        True if candidate overlaps the variant, False otherwise
    """
    # Convert candidate position to genomic coordinates
    candidate_genomic_start = transcript_start + candidate_pos - 1
    candidate_genomic_end = candidate_genomic_start + candidate_length - 1

    # Check if variant position falls within candidate region
    # variant.pos is 1-based genomic position
    variant_genomic_end = variant.pos + len(variant.ref) - 1

    overlaps = not (candidate_genomic_end < variant.pos or candidate_genomic_start > variant_genomic_end)

    if overlaps:
        logger.debug(
            f"Candidate at {candidate_genomic_start}-{candidate_genomic_end} overlaps "
            f"variant at {variant.pos}-{variant_genomic_end}"
        )

    return overlaps


def annotate_candidate_with_variant(
    candidate: SiRNACandidate,
    variant: VariantRecord,
    allele: str,
    variant_mode: str,
) -> None:
    """Annotate a siRNA candidate with variant information.

    Modifies the candidate in place to add variant-specific metadata.

    Args:
        candidate: SiRNACandidate to annotate
        variant: Variant record that overlaps the candidate
        allele: Which allele the candidate targets ('ref' or 'alt')
        variant_mode: Variant handling mode ('target', 'avoid', 'both')
    """
    # Serialize VariantRecord to dict for storage
    variant_dict = variant.model_dump()

    # Add variant to overlapped_variants list
    if variant_dict not in candidate.overlapped_variants:
        candidate.overlapped_variants.append(variant_dict)

    # Mark as allele-specific
    candidate.allele_specific = True

    # Track targeted alleles
    if allele not in candidate.targeted_alleles:
        candidate.targeted_alleles.append(allele)

    # Set variant mode
    candidate.variant_mode = variant_mode

    logger.debug(f"Annotated candidate {candidate.id} with variant {variant.to_vcf_style()} (allele: {allele})")


def apply_variant_to_sequence(
    sequence: str,
    variant: VariantRecord,
    transcript_start: int,
    allele: str = "alt",
) -> str:
    """Apply a variant to a reference sequence to generate an alternate sequence.

    Args:
        sequence: Reference sequence
        variant: Variant record with position and alleles
        transcript_start: 1-based genomic start position of transcript
        allele: Which allele to apply ('ref' returns unchanged sequence, 'alt' applies variant)

    Returns:
        Modified sequence with variant applied

    Raises:
        ValueError: If variant position is outside sequence boundaries or allele is invalid
    """
    if allele == "ref":
        return sequence

    if allele != "alt":
        raise ValueError(f"Invalid allele: {allele}. Must be 'ref' or 'alt'")

    # Calculate variant position within transcript (0-based)
    variant_transcript_pos = variant.pos - transcript_start

    if variant_transcript_pos < 0 or variant_transcript_pos >= len(sequence):
        raise ValueError(
            f"Variant position {variant.pos} is outside sequence boundaries "
            f"(transcript_start: {transcript_start}, sequence_length: {len(sequence)})"
        )

    # Verify reference allele matches sequence
    ref_in_seq = sequence[variant_transcript_pos : variant_transcript_pos + len(variant.ref)]
    if ref_in_seq.upper() != variant.ref.upper():
        logger.warning(
            f"Reference allele mismatch at position {variant.pos}: "
            f"expected {variant.ref}, found {ref_in_seq} in sequence"
        )

    # Apply the variant by replacing ref with alt
    modified_sequence = (
        sequence[:variant_transcript_pos] + variant.alt + sequence[variant_transcript_pos + len(variant.ref) :]
    )

    logger.debug(f"Applied variant {variant.to_vcf_style()} to sequence (allele: {allele})")

    return modified_sequence


def get_variant_position_in_transcript(
    variant: VariantRecord,
    transcript_start: int,
) -> int:
    """Calculate the 0-based position of a variant within a transcript.

    Args:
        variant: Variant record
        transcript_start: 1-based genomic start position of transcript

    Returns:
        0-based position of variant within transcript
    """
    return variant.pos - transcript_start
