"""Unit tests for variant helper functions."""

import pytest

from sirnaforge.data.variant_helpers import (
    apply_variant_to_sequence,
    check_candidate_overlaps_variant,
    generate_contexts_for_variant,
    get_variant_position_in_transcript,
)
from sirnaforge.models.variant import VariantRecord


class TestGenerateContextsForVariant:
    """Tests for generate_contexts_for_variant."""

    def test_generate_contexts_simple_snp(self):
        """Test generating contexts for a simple SNP."""
        # Create a simple variant: G>A at position 50 in transcript starting at genomic position 1
        variant = VariantRecord(
            chr="chr1",
            pos=51,  # 1-based genomic position (position 50 in 0-based transcript)
            ref="G",
            alt="A",
        )

        # Reference sequence with G at position 50 (0-based)
        reference_sequence = "A" * 50 + "G" + "C" * 50  # 101 bp total

        contexts = generate_contexts_for_variant(
            variant=variant,
            reference_sequence=reference_sequence,
            transcript_start=1,  # Transcript starts at genomic position 1
            transcript_id="TEST_TRANSCRIPT",
            flank_size=20,
        )

        # Check ref context
        assert "ref" in contexts
        ref_seq, ref_pos, ref_len = contexts["ref"]
        assert "G" in ref_seq  # Contains the reference allele
        assert ref_len == 1
        assert ref_seq[ref_pos : ref_pos + ref_len] == "G"

        # Check alt context
        assert "alt" in contexts
        alt_seq, alt_pos, alt_len = contexts["alt"]
        assert "A" in alt_seq  # Contains the alternate allele
        assert alt_len == 1
        assert alt_seq[alt_pos : alt_pos + alt_len] == "A"

    def test_generate_contexts_deletion(self):
        """Test generating contexts for a deletion variant."""
        variant = VariantRecord(
            chr="chr1",
            pos=21,  # Delete 3 bases starting at position 21
            ref="ATG",
            alt="A",
        )

        reference_sequence = "A" * 20 + "ATG" + "C" * 20  # 43 bp total

        contexts = generate_contexts_for_variant(
            variant=variant,
            reference_sequence=reference_sequence,
            transcript_start=1,
            transcript_id="TEST",
            flank_size=10,
        )

        ref_seq, ref_pos, ref_len = contexts["ref"]
        assert ref_len == 3
        assert ref_seq[ref_pos : ref_pos + ref_len] == "ATG"

        alt_seq, alt_pos, alt_len = contexts["alt"]
        assert alt_len == 1
        assert alt_seq[alt_pos : alt_pos + alt_len] == "A"
        assert len(alt_seq) < len(ref_seq)  # Alt sequence is shorter due to deletion

    def test_generate_contexts_insertion(self):
        """Test generating contexts for an insertion variant."""
        variant = VariantRecord(
            chr="chr1",
            pos=21,
            ref="A",
            alt="ATGC",  # Insert TGC after A
        )

        reference_sequence = "A" * 20 + "A" + "C" * 20  # 41 bp total

        contexts = generate_contexts_for_variant(
            variant=variant,
            reference_sequence=reference_sequence,
            transcript_start=1,
            transcript_id="TEST",
            flank_size=10,
        )

        ref_seq, ref_pos, ref_len = contexts["ref"]
        assert ref_len == 1

        alt_seq, alt_pos, alt_len = contexts["alt"]
        assert alt_len == 4
        assert alt_seq[alt_pos : alt_pos + alt_len] == "ATGC"
        assert len(alt_seq) > len(ref_seq)  # Alt sequence is longer due to insertion

    def test_generate_contexts_out_of_bounds(self):
        """Test that variant outside transcript raises ValueError."""
        variant = VariantRecord(
            chr="chr1",
            pos=200,  # Beyond transcript length
            ref="G",
            alt="A",
        )

        reference_sequence = "A" * 100

        with pytest.raises(ValueError, match="outside transcript"):
            generate_contexts_for_variant(
                variant=variant,
                reference_sequence=reference_sequence,
                transcript_start=1,
                transcript_id="TEST",
            )

    def test_generate_contexts_at_boundaries(self):
        """Test generating contexts at transcript boundaries."""
        # Variant at the very start
        variant_start = VariantRecord(chr="chr1", pos=1, ref="A", alt="G")
        reference_sequence = "A" * 100

        contexts_start = generate_contexts_for_variant(
            variant=variant_start,
            reference_sequence=reference_sequence,
            transcript_start=1,
            transcript_id="TEST",
            flank_size=10,
        )

        assert "ref" in contexts_start
        assert "alt" in contexts_start

        # Variant at the end
        variant_end = VariantRecord(chr="chr1", pos=100, ref="A", alt="G")

        contexts_end = generate_contexts_for_variant(
            variant=variant_end,
            reference_sequence=reference_sequence,
            transcript_start=1,
            transcript_id="TEST",
            flank_size=10,
        )

        assert "ref" in contexts_end
        assert "alt" in contexts_end


class TestCheckCandidateOverlapsVariant:
    """Tests for check_candidate_overlaps_variant."""

    def test_candidate_overlaps_variant(self):
        """Test detecting overlap between candidate and variant."""
        variant = VariantRecord(chr="chr1", pos=50, ref="G", alt="A")

        # Candidate that overlaps: positions 40-60 in transcript (genomic 41-61)
        assert (
            check_candidate_overlaps_variant(
                candidate_pos=40,
                candidate_length=21,
                variant=variant,
                transcript_start=1,
            )
            is True
        )

    def test_candidate_does_not_overlap_variant(self):
        """Test detecting no overlap between candidate and variant."""
        variant = VariantRecord(chr="chr1", pos=100, ref="G", alt="A")

        # Candidate that doesn't overlap: positions 10-30 in transcript (genomic 11-31)
        assert (
            check_candidate_overlaps_variant(
                candidate_pos=10,
                candidate_length=21,
                variant=variant,
                transcript_start=1,
            )
            is False
        )

    def test_candidate_overlaps_deletion(self):
        """Test overlap detection with deletion variant."""
        variant = VariantRecord(chr="chr1", pos=50, ref="ATG", alt="A")

        # Candidate overlapping the deletion
        assert (
            check_candidate_overlaps_variant(
                candidate_pos=45,
                candidate_length=21,
                variant=variant,
                transcript_start=1,
            )
            is True
        )

    def test_candidate_adjacent_to_variant(self):
        """Test that adjacent candidates don't overlap."""
        variant = VariantRecord(chr="chr1", pos=50, ref="G", alt="A")

        # Candidate ending just before variant (positions 28-48 in transcript, genomic 29-49)
        assert (
            check_candidate_overlaps_variant(
                candidate_pos=28,
                candidate_length=21,
                variant=variant,
                transcript_start=1,
            )
            is False
        )

        # Candidate starting just after variant (positions 51-71 in transcript, genomic 52-72)
        assert (
            check_candidate_overlaps_variant(
                candidate_pos=51,
                candidate_length=21,
                variant=variant,
                transcript_start=1,
            )
            is False
        )


class TestApplyVariantToSequence:
    """Tests for apply_variant_to_sequence."""

    def test_apply_variant_ref_allele(self):
        """Test applying ref allele returns unchanged sequence."""
        variant = VariantRecord(chr="chr1", pos=51, ref="G", alt="A")
        sequence = "A" * 50 + "G" + "C" * 50

        result = apply_variant_to_sequence(sequence, variant, transcript_start=1, allele="ref")

        assert result == sequence

    def test_apply_variant_alt_allele_snp(self):
        """Test applying alt allele for SNP."""
        variant = VariantRecord(chr="chr1", pos=51, ref="G", alt="A")
        sequence = "A" * 50 + "G" + "C" * 50

        result = apply_variant_to_sequence(sequence, variant, transcript_start=1, allele="alt")

        expected = "A" * 50 + "A" + "C" * 50
        assert result == expected

    def test_apply_variant_deletion(self):
        """Test applying deletion variant."""
        variant = VariantRecord(chr="chr1", pos=21, ref="ATG", alt="A")
        sequence = "A" * 20 + "ATG" + "C" * 20

        result = apply_variant_to_sequence(sequence, variant, transcript_start=1, allele="alt")

        expected = "A" * 20 + "A" + "C" * 20
        assert result == expected
        assert len(result) < len(sequence)

    def test_apply_variant_insertion(self):
        """Test applying insertion variant."""
        variant = VariantRecord(chr="chr1", pos=21, ref="A", alt="ATGC")
        sequence = "A" * 20 + "A" + "C" * 20

        result = apply_variant_to_sequence(sequence, variant, transcript_start=1, allele="alt")

        expected = "A" * 20 + "ATGC" + "C" * 20
        assert result == expected
        assert len(result) > len(sequence)

    def test_apply_variant_invalid_allele(self):
        """Test that invalid allele raises ValueError."""
        variant = VariantRecord(chr="chr1", pos=51, ref="G", alt="A")
        sequence = "A" * 100

        with pytest.raises(ValueError, match="Invalid allele"):
            apply_variant_to_sequence(sequence, variant, transcript_start=1, allele="invalid")

    def test_apply_variant_out_of_bounds(self):
        """Test that variant outside sequence raises ValueError."""
        variant = VariantRecord(chr="chr1", pos=200, ref="G", alt="A")
        sequence = "A" * 100

        with pytest.raises(ValueError, match="outside sequence boundaries"):
            apply_variant_to_sequence(sequence, variant, transcript_start=1, allele="alt")


class TestGetVariantPositionInTranscript:
    """Tests for get_variant_position_in_transcript."""

    def test_get_variant_position_simple(self):
        """Test getting variant position in transcript."""
        variant = VariantRecord(chr="chr1", pos=51, ref="G", alt="A")

        # Transcript starts at genomic position 1
        pos = get_variant_position_in_transcript(variant, transcript_start=1)

        assert pos == 50  # 0-based position in transcript

    def test_get_variant_position_offset_transcript(self):
        """Test variant position with transcript not starting at 1."""
        variant = VariantRecord(chr="chr1", pos=1051, ref="G", alt="A")

        # Transcript starts at genomic position 1001
        pos = get_variant_position_in_transcript(variant, transcript_start=1001)

        assert pos == 50  # 0-based position in transcript
