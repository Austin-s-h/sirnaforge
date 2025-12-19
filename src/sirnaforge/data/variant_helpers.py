"""Helper functions for generating sequence contexts with variant alleles."""

from sirnaforge.models.variant import VariantRecord
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


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
    candidate: "SiRNACandidate",  # type: ignore # noqa: F821
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
