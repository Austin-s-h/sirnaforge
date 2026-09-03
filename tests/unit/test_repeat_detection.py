"""Unit tests for reference-relative repeat detection."""

from __future__ import annotations

import logging
import time
from pathlib import Path

import pytest

from sirnaforge.core import repeat_detection
from sirnaforge.core.repeat_detection import (
    DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
    RepeatDetector,
    _compute_kmer_code,
    _normalize_sequence,
    _reverse_complement,
    _rolling_kmer_codes,
    _sequence_to_codes,
)


@pytest.mark.unit
def test_normalize_sequence() -> None:
    """Normalization converts to uppercase and U -> T."""
    assert _normalize_sequence("acguACGU") == "ACGTACGT"
    assert _normalize_sequence("uuuu") == "TTTT"
    assert _normalize_sequence("ACGT") == "ACGT"


@pytest.mark.unit
def test_reverse_complement() -> None:
    """Reverse complement is computed correctly."""
    assert _reverse_complement("ACGT") == "ACGT"
    assert _reverse_complement("AAAA") == "TTTT"
    assert _reverse_complement("GCTA") == "TAGC"
    assert _reverse_complement("ATCGATCG") == "CGATCGAT"


@pytest.mark.unit
def test_planted_repeat_is_flagged(tmp_path: Path) -> None:
    """A sequence planted in multiple transcripts above threshold is flagged."""
    # Create a reference with 100 transcripts
    # Plant the sequence in 20 of them (20% = well above 0.1% threshold)
    ref = tmp_path / "ref.fa"
    repeat_seq = "ACGTACGTACGTACGTACG"  # 19 nt
    unique_seq = "TTTTTTTTTTTTTTTTTTT"

    with ref.open("w") as f:
        for i in range(20):
            # These contain the repeat
            f.write(f">TX{i:03d}\n")
            f.write(f"AAAAAAAAAA{repeat_seq}CCCCCCCCCC\n")
        for i in range(20, 100):
            # These don't
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    detector = RepeatDetector()
    result = detector.scan([repeat_seq, unique_seq], ref)

    assert result.reference_transcript_count == 100
    assert result.threshold_fraction == DEFAULT_REPEAT_TRANSCRIPT_FRACTION

    repeat_obs = result.observations[repeat_seq]
    assert repeat_obs.transcript_count == 20
    assert repeat_obs.transcript_fraction == 0.20
    assert repeat_obs.is_repeat is True

    unique_obs = result.observations[unique_seq]
    assert unique_obs.transcript_count == 0
    assert unique_obs.transcript_fraction == 0.0
    assert unique_obs.is_repeat is False


@pytest.mark.unit
def test_fractional_threshold_scales_with_reference_size(tmp_path: Path) -> None:
    """The same planted repeat at the same RATE is caught in both small and large references.

    This is the key test proving that fractional semantics work correctly.
    """
    repeat_seq = "ACGTACGTACGTACGTACG"  # 19 nt

    # Small reference: 10 transcripts, 4 contain the repeat (40%)
    small_ref = tmp_path / "small.fa"
    with small_ref.open("w") as f:
        for i in range(4):
            f.write(f">SMALL_TX{i:03d}\n")
            f.write(f"AAAAAAAAAA{repeat_seq}CCCCCCCCCC\n")
        for i in range(4, 10):
            f.write(f">SMALL_TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    # Large reference: 400 transcripts, 160 contain the repeat (40%)
    large_ref = tmp_path / "large.fa"
    with large_ref.open("w") as f:
        for i in range(160):
            f.write(f">LARGE_TX{i:03d}\n")
            f.write(f"AAAAAAAAAA{repeat_seq}CCCCCCCCCC\n")
        for i in range(160, 400):
            f.write(f">LARGE_TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    # Use the same threshold for both
    detector = RepeatDetector(threshold_fraction=0.001)

    small_result = detector.scan([repeat_seq], small_ref)
    large_result = detector.scan([repeat_seq], large_ref)

    # Both should flag the repeat
    assert small_result.observations[repeat_seq].is_repeat is True
    assert large_result.observations[repeat_seq].is_repeat is True

    # Check the fractions are the same (40%)
    assert small_result.observations[repeat_seq].transcript_fraction == pytest.approx(0.4)
    assert large_result.observations[repeat_seq].transcript_fraction == pytest.approx(0.4)

    # Prove that an absolute-count rule calibrated to the large reference
    # would NOT have caught it in the small reference
    absolute_threshold = 50  # e.g., "> 50 transcripts"

    small_count = small_result.observations[repeat_seq].transcript_count
    large_count = large_result.observations[repeat_seq].transcript_count

    assert small_count == 4
    assert large_count == 160

    # The absolute rule would have missed the small reference
    assert small_count <= absolute_threshold, "Absolute rule should NOT flag small reference"
    assert large_count > absolute_threshold, "Absolute rule SHOULD flag large reference"


@pytest.mark.unit
def test_reverse_complement_orientation_detected(tmp_path: Path) -> None:
    """A sequence's reverse complement is also counted."""
    ref = tmp_path / "ref.fa"
    # Forward sequence
    fwd = "ACGTACGTACGTACGTACG"
    # Reverse complement
    rc = _reverse_complement(fwd)

    # Plant the RC in transcripts
    with ref.open("w") as f:
        for i in range(20):
            f.write(f">TX{i:03d}\n")
            f.write(f"AAAAAAAAAA{rc}CCCCCCCCCC\n")
        for i in range(20, 100):
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    detector = RepeatDetector()
    # Query with the forward sequence
    result = detector.scan([fwd], ref)

    # Should find it via RC matching
    obs = result.observations[fwd]
    assert obs.transcript_count == 20
    assert obs.transcript_fraction == 0.20
    assert obs.is_repeat is True


@pytest.mark.unit
def test_transcript_counted_once_despite_multiple_occurrences(tmp_path: Path) -> None:
    """A transcript containing the sequence three times counts once."""
    ref = tmp_path / "ref.fa"
    seq = "ACGTACGTACGTACGTACG"

    # One transcript with the sequence repeated 3 times
    with ref.open("w") as f:
        f.write(">TX001\n")
        f.write(f"{seq}AAAAA{seq}CCCCC{seq}\n")
        for i in range(1, 100):
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    detector = RepeatDetector()
    result = detector.scan([seq], ref)

    obs = result.observations[seq]
    # Should count the transcript once, not three times
    assert obs.transcript_count == 1
    assert obs.transcript_fraction == 0.01
    assert obs.is_repeat is True  # 0.01 > 0.001


@pytest.mark.unit
def test_n_containing_windows_do_not_match(tmp_path: Path) -> None:
    """Windows containing N or other non-ACGT bases do not match."""
    ref = tmp_path / "ref.fa"
    seq = "ACGTACGTACGTACGTACG"

    # Plant the sequence with an N in the middle in one transcript
    with ref.open("w") as f:
        f.write(">TX001\n")
        f.write(f"{seq[:10]}N{seq[10:]}\n")  # Break the sequence with N
        f.write(">TX002\n")
        f.write(f"{seq}AAAA\n")  # Clean copy
        for i in range(2, 100):
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    detector = RepeatDetector()
    result = detector.scan([seq], ref)

    obs = result.observations[seq]
    # Should only match TX002, not TX001 (which has the N-broken version)
    assert obs.transcript_count == 1
    assert obs.transcript_fraction == 0.01


@pytest.mark.unit
def test_mixed_needle_lengths_in_one_scan(tmp_path: Path) -> None:
    """Multiple distinct needle lengths (19, 21, 23) all resolve correctly."""
    ref = tmp_path / "ref.fa"

    # Use random-looking distinct sequences with no repeats or overlap
    seq_19 = "AGCTACGTAGCATCGATCA"  # 19 nt
    seq_21 = "TGACGATCGATCGACGATCGA"  # 21 nt
    seq_23 = "CAGTACGTACGTACGTACGTACG"  # 23 nt

    with ref.open("w") as f:
        # TX001 contains seq_19
        f.write(">TX001\n")
        f.write(f"AAATTTCCC{seq_19}GGGCCCAAA\n")
        # TX002 contains seq_21
        f.write(">TX002\n")
        f.write(f"TTTAAAGGG{seq_21}CCCGGGTTT\n")
        # TX003 contains seq_23
        f.write(">TX003\n")
        f.write(f"GGGAAACCC{seq_23}TTTGGGCCC\n")
        # TX004-100 contain nothing overlapping
        for i in range(3, 100):
            f.write(f">TX{i:03d}\n")
            f.write("ATATATATATATATATATATAT\n")

    detector = RepeatDetector()
    result = detector.scan([seq_19, seq_21, seq_23], ref)

    assert result.observations[seq_19].transcript_count == 1
    assert result.observations[seq_21].transcript_count == 1
    assert result.observations[seq_23].transcript_count == 1


@pytest.mark.unit
def test_brute_force_cross_check(tmp_path: Path) -> None:
    """Vectorised counts match a naive string-based implementation."""
    ref = tmp_path / "ref.fa"

    # Small deterministic synthetic reference
    transcripts = {
        "TX001": "ACGTACGTACGTACGTACGTACGTACGT",
        "TX002": "TTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        "TX003": "ACGTTTTTACGTACGTACGTGGGGGGGG",
        "TX004": "CCCCCCCCACGTACGTACGTACGTCCCC",
        "TX005": "GGGGGGGGGGGGGGGGGGGGGGGGGGGG",
    }

    with ref.open("w") as f:
        for tid, seq in transcripts.items():
            f.write(f">{tid}\n{seq}\n")

    needles = [
        "ACGTACGTACGTACG",  # 15 nt - appears in multiple
        "TTTTTTTTTTTTTT",  # 14 nt - appears in TX002
        "GGGGGGGGGGGGGG",  # 14 nt - appears in TX003 and TX005
    ]

    # Run the vectorized implementation
    detector = RepeatDetector()
    result = detector.scan(needles, ref)

    # Naive string-based counting for cross-check
    def naive_count(needle: str, transcripts: dict[str, str]) -> int:
        """Count transcripts containing needle or its RC."""
        norm_needle = _normalize_sequence(needle)
        rc = _reverse_complement(norm_needle)
        count = 0
        for seq in transcripts.values():
            norm_seq = _normalize_sequence(seq)
            if norm_needle in norm_seq or rc in norm_seq:
                count += 1
        return count

    for needle in needles:
        norm = _normalize_sequence(needle)
        expected_count = naive_count(needle, transcripts)
        actual_count = result.observations[norm].transcript_count
        assert actual_count == expected_count, f"Mismatch for {norm}: expected {expected_count}, got {actual_count}"


@pytest.mark.unit
def test_missing_reference_file_returns_empty_result(tmp_path: Path) -> None:
    """Missing reference file returns empty result with is_repeat=False, no raise."""
    missing = tmp_path / "does_not_exist.fa"
    detector = RepeatDetector()

    result = detector.scan(["ACGTACGTACGTACGTACG"], missing)

    assert result.reference_transcript_count == 0
    assert len(result.observations) == 1
    obs = result.observations["ACGTACGTACGTACGTACG"]
    assert obs.transcript_count == 0
    assert obs.transcript_fraction == 0.0
    assert obs.is_repeat is False


@pytest.mark.unit
def test_empty_needle_set_returns_empty_result(tmp_path: Path) -> None:
    """Empty needle set returns empty observations."""
    ref = tmp_path / "ref.fa"
    with ref.open("w") as f:
        f.write(">TX001\nACGTACGTACGTACGTACGT\n")

    detector = RepeatDetector()
    result = detector.scan([], ref)

    assert result.reference_transcript_count == 0
    assert len(result.observations) == 0


@pytest.mark.unit
def test_empty_reference_returns_empty_result(tmp_path: Path) -> None:
    """Empty reference (no sequences) returns reference_transcript_count=0."""
    ref = tmp_path / "empty.fa"
    ref.touch()  # Create empty file

    detector = RepeatDetector()
    result = detector.scan(["ACGTACGTACGTACGTACG"], ref)

    assert result.reference_transcript_count == 0
    obs = result.observations["ACGTACGTACGTACGTACG"]
    assert obs.is_repeat is False


@pytest.mark.unit
def test_sequence_longer_than_all_transcripts(tmp_path: Path) -> None:
    """A needle longer than every transcript in the reference is handled gracefully."""
    ref = tmp_path / "ref.fa"
    with ref.open("w") as f:
        f.write(">TX001\nACGT\n")  # 4 nt
        f.write(">TX002\nTTTT\n")  # 4 nt

    long_needle = "A" * 100  # 100 nt

    detector = RepeatDetector()
    result = detector.scan([long_needle], ref)

    assert result.reference_transcript_count == 2
    obs = result.observations[long_needle]
    assert obs.transcript_count == 0
    assert obs.is_repeat is False


@pytest.mark.unit
def test_repeat_sequences_property(tmp_path: Path) -> None:
    """repeat_sequences property returns only the flagged sequences."""
    ref = tmp_path / "ref.fa"
    repeat_seq = "AGCTACGTAGCATCGATCA"  # 19 nt, appears in many transcripts
    unique_seq = "TGACGATCGATCGACGACT"  # 19 nt, appears in only one

    # Use 2000 transcripts so 1/2000 = 0.0005 < 0.001 threshold
    with ref.open("w") as f:
        for i in range(100):
            f.write(f">TX{i:04d}\n")
            f.write(f"TTTCCCAAA{repeat_seq}GGGAAATTT\n")
        # One transcript with the unique sequence
        f.write(">TX0100\n")
        f.write(f"AAATTTGGG{unique_seq}CCCGGGTTT\n")
        for i in range(101, 2000):
            f.write(f">TX{i:04d}\n")
            f.write("ATATATATATATATATATATAT\n")

    detector = RepeatDetector()
    result = detector.scan([repeat_seq, unique_seq], ref)

    repeats = result.repeat_sequences
    assert repeat_seq in repeats  # 100/2000 = 0.05 > 0.001
    assert unique_seq not in repeats  # 1/2000 = 0.0005 < 0.001
    assert len(repeats) == 1


@pytest.mark.unit
def test_rna_input_sequences_normalized(tmp_path: Path) -> None:
    """RNA input sequences (with U) are normalized to DNA (T) before matching."""
    ref = tmp_path / "ref.fa"
    dna_seq = "ACGTACGTACGTACGTACG"  # DNA

    with ref.open("w") as f:
        for i in range(20):
            f.write(f">TX{i:03d}\n")
            f.write(f"AAA{dna_seq}CCC\n")
        for i in range(20, 100):
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    # Query with RNA (U instead of T)
    rna_seq = "ACGUACGUACGUACGUACG"

    detector = RepeatDetector()
    result = detector.scan([rna_seq], ref)

    # Should find matches because U -> T normalization
    norm = _normalize_sequence(rna_seq)
    obs = result.observations[norm]
    assert obs.transcript_count == 20
    assert obs.is_repeat is True


@pytest.mark.unit
def test_threshold_at_boundary(tmp_path: Path) -> None:
    """is_repeat verdict is correct at the threshold boundary."""
    ref = tmp_path / "ref.fa"
    seq = "ACGTACGTACGTACGTACG"  # 19 nt

    # Create reference with exactly 1000 transcripts
    # Plant the sequence in exactly 1 transcript (0.001 = 0.1%)
    with ref.open("w") as f:
        f.write(">TX000\n")
        f.write(f"CCC{seq}GGG\n")
        for i in range(1, 1000):
            f.write(f">TX{i:03d}\n")
            f.write("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")

    detector = RepeatDetector(threshold_fraction=0.001)
    result = detector.scan([seq], ref)

    obs = result.observations[seq]
    assert obs.transcript_count == 1
    assert obs.transcript_fraction == 0.001

    # At exactly the threshold, should NOT be flagged (> not >=)
    assert obs.is_repeat is False

    # Just above threshold should be flagged
    detector2 = RepeatDetector(threshold_fraction=0.0009)
    result2 = detector2.scan([seq], ref)
    obs2 = result2.observations[seq]
    assert obs2.is_repeat is True


@pytest.mark.unit
def test_performance_guard_2mb_reference(tmp_path: Path) -> None:
    """Performance guard: 2MB synthetic reference scans in under 2.0 seconds.

    Pre-rewrite performance was 6.57s for this size, which extrapolates to ~55 minutes
    for a 1GB reference. This test discriminates between the old and new implementations.
    """
    ref = tmp_path / "ref.fa"

    # Generate a deterministic 2MB reference: 400 transcripts x 5000 nt
    # Use a repeating base pattern to avoid random module
    base = "ACGTTGCAAGGCTTAACCGGATCGATCGTTACG" * 200  # 6800 nt pattern

    with ref.open("w") as f:
        for i in range(400):
            f.write(f">T{i}\n")
            # Offset by i % 100 for variation
            start = i % 100
            f.write(base[start : start + 5000] + "\n")

    # Plant 50 needles of length 21
    needles = [base[i * 7 : i * 7 + 21] for i in range(50)]

    detector = RepeatDetector()
    start = time.perf_counter()
    detector.scan(needles, ref)
    elapsed = time.perf_counter() - start

    assert elapsed < 2.0, f"Scan took {elapsed:.3f}s, should be under 2.0s"


@pytest.mark.unit
def test_vectorized_codes_match_scalar() -> None:
    """Vectorized rolling codes equal the scalar implementation for every window."""
    # Hand-written test sequence with a mix of bases and an N
    seq = "ACGTACGTNACGTACGTACGT"
    seq_codes = _sequence_to_codes(seq)

    k = 5
    codes, valid = _rolling_kmer_codes(seq_codes, k)

    # Check every window against the scalar implementation
    for i in range(len(seq) - k + 1):
        window = seq[i : i + k]
        expected_code = _compute_kmer_code(window, normalize=False)

        if expected_code is None:
            # Invalid window should be marked as invalid
            assert not valid[i], f"Window {i} ('{window}') should be invalid"
        else:
            # Valid window should match
            assert valid[i], f"Window {i} ('{window}') should be valid"
            assert codes[i] == expected_code, f"Window {i} ('{window}') code mismatch"


@pytest.mark.unit
def test_needle_spanning_transcript_boundary_not_counted(tmp_path: Path) -> None:
    """A needle that would span a transcript boundary is not counted.

    This tests the specific bug the sentinel padding exists to prevent.
    """
    ref = tmp_path / "ref.fa"

    # Two transcripts that, if concatenated, would contain the needle
    # but neither one alone contains it
    needle = "ACGTACGTACGTACGTACG"  # 19 nt

    # tx1 ends with the first 12 bases of the needle
    tx1 = "TTTTTTTT" + needle[:12]  # "TTTTTTTTACGTACGTACGT"

    # tx2 starts with the last 7 bases of the needle
    tx2 = needle[12:] + "GGGGGGGG"  # "ACGTACGGGGGGGGG"

    # If concatenated without sentinel: tx1 + tx2 would contain the needle
    # But neither alone should match

    with ref.open("w") as f:
        f.write(">TX001\n")
        f.write(f"{tx1}\n")
        f.write(">TX002\n")
        f.write(f"{tx2}\n")

    detector = RepeatDetector()
    result = detector.scan([needle], ref)

    obs = result.observations[needle]
    # Neither transcript alone contains the full needle
    assert obs.transcript_count == 0


@pytest.mark.unit
def test_lowercase_reference_still_matches(tmp_path: Path) -> None:
    """Lowercase reference sequences still match uppercase needles."""
    ref = tmp_path / "ref.fa"
    needle = "ACGTACGTACGTACGTACG"  # Uppercase

    # Write reference in lowercase
    with ref.open("w") as f:
        for i in range(20):
            f.write(f">TX{i:03d}\n")
            f.write(f"aaa{needle.lower()}ccc\n")  # lowercase
        for i in range(20, 100):
            f.write(f">TX{i:03d}\n")
            f.write("gggggggggggggggggggggggggggggg\n")

    detector = RepeatDetector()
    result = detector.scan([needle], ref)

    obs = result.observations[needle]
    assert obs.transcript_count == 20
    assert obs.is_repeat is True


@pytest.mark.unit
def test_single_transcript_exceeding_chunk_size(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A single transcript longer than _CHUNK_SIZE is handled correctly."""
    # Temporarily shrink the chunk size for testing
    original_chunk_size = repeat_detection._CHUNK_SIZE
    monkeypatch.setattr(repeat_detection, "_CHUNK_SIZE", 500)

    ref = tmp_path / "ref.fa"
    needle = "ACGTACGTACGTACGTACG"  # 19 nt

    # Build a transcript longer than the patched chunk size (500 bytes)
    # Plant the needle in the middle
    long_seq = "A" * 300 + needle + "T" * 300  # ~620 nt

    with ref.open("w") as f:
        f.write(">TX001\n")
        f.write(f"{long_seq}\n")
        for i in range(1, 100):
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    detector = RepeatDetector()
    result = detector.scan([needle], ref)

    obs = result.observations[needle]
    assert obs.transcript_count == 1

    # Restore original chunk size
    monkeypatch.setattr(repeat_detection, "_CHUNK_SIZE", original_chunk_size)


@pytest.mark.unit
def test_needle_longer_than_31_bases_skipped_with_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Needles longer than 31 bases are skipped with a warning (code width guard)."""
    ref = tmp_path / "ref.fa"

    # Create a needle longer than 31 bases
    long_needle = "A" * 35

    with ref.open("w") as f:
        f.write(">TX001\n")
        f.write(f"{long_needle}CCCCCCCCC\n")
        for i in range(1, 100):
            f.write(f">TX{i:03d}\n")
            f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

    detector = RepeatDetector()
    with caplog.at_level(logging.WARNING):
        result = detector.scan([long_needle], ref)

    # Should have logged a warning
    assert any("code width limit" in record.message for record in caplog.records)

    # Needle should be in observations but with count 0 (no match attempted)
    obs = result.observations.get(long_needle)
    # Actually, it won't be in observations because it was skipped entirely
    # Let me check the implementation logic
    # Looking at the code, if a needle is >31, it's skipped in _prepare_needle_index
    # and won't be in length_data, so it won't be processed at all
    # But the scan method initializes transcript_hits for all needles in needle_set
    # So it should still be in the observations with count 0
    assert obs is not None
    assert obs.transcript_count == 0
    assert obs.is_repeat is False
