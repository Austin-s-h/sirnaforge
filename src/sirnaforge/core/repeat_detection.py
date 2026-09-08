"""Reference-relative repeat detection for siRNA candidate sequences.

This module identifies guide sequences that occur frequently across a reference transcriptome,
flagging them as potential repeat-element overlaps. The detection threshold is expressed as a
**fraction of reference transcripts** rather than an absolute count, ensuring the verdict behaves
consistently across gene families of different sizes (user story 7).

Threshold derivation
--------------------
A clean guide occurs only in its target gene's isoforms — typically 1-20 of ~190,000 human cDNA
entries (~0.01%). A guide overlapping a dispersed repeat element (e.g., Alu-like) occurs in
hundreds of transcripts (>450 observed, ~0.24%). The default threshold of 0.001 (0.1%) sits an
order of magnitude clear of both populations.

**This is a starting point, not a calibrated value.** Empirical tuning on real candidate sets
may refine the threshold.

Implementation notes
--------------------
- Normalisation: uppercase, U -> T (guide sequences may be RNA).
- Orientation: both the sequence and its reverse complement are counted. A guide is antisense to
  its target, so the sense k-mer sits in the cDNA. BWA reports either strand, and this makes the
  design-time flag comparable to the post-alignment cross-check.
- Transcript counting: a transcript containing the sequence once or multiple times counts once.
- Invalid windows: any k-mer overlapping a non-ACGT base matches nothing.
- Performance: one streaming pass over the reference using vectorised numpy operations. Peak memory
  is O(needles + max_transcript_length), independent of total reference size.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from sirnaforge.utils.fasta import _iter_fasta_records

logger = logging.getLogger(__name__)

DEFAULT_REPEAT_TRANSCRIPT_FRACTION = 0.001  # 0.1% of reference transcripts
_CHUNK_SIZE = 10_000_000  # 10 MB buffer; 21-mers => ~80 MB int64 code array
_INVALID = np.uint8(255)  # Marker for non-ACGT bases

# Membership prefilter: a resident bitmap on the code's low bits. One gather into 4 MB is far
# cheaper per window than searchsorted's ~10 random probes, and only needles plus ~0.05% false
# positives reach the exact check.
_PREFILTER_BITS = 22
_PREFILTER_SIZE = 1 << _PREFILTER_BITS
_PREFILTER_MASK = _PREFILTER_SIZE - 1

# Base encoding lookup table (256-entry)
_LUT = np.full(256, _INVALID, dtype=np.uint8)
_LUT[ord("A")] = 0
_LUT[ord("C")] = 1
_LUT[ord("G")] = 2
_LUT[ord("T")] = 3
_LUT[ord("a")] = 0
_LUT[ord("c")] = 1
_LUT[ord("g")] = 2
_LUT[ord("t")] = 3


def normalize_guide_sequence(seq: str) -> str:
    """Normalize a nucleotide sequence for matching: uppercase, U -> T.

    Public because the hit classifier must key repeat-flagged guides exactly the way this
    module counted them; two normalisers that drift would silently stop the REPEAT class firing.

    Args:
        seq: Nucleotide sequence, RNA or DNA, any case.

    Returns:
        Upper-case DNA form of the sequence.
    """
    return seq.upper().replace("U", "T")


# Internal alias retained so this module's existing call sites stay terse.
_normalize_sequence = normalize_guide_sequence


def _reverse_complement(seq: str) -> str:
    """Compute the reverse complement of a DNA sequence.

    Assumes input is already normalized (uppercase, T not U).
    """
    complement_table = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement_table)[::-1]


def _encode_base(base: str) -> int:
    """Encode a single nucleotide to 2-bit integer: A=0, C=1, G=2, T=3.

    Returns -1 for non-ACGT bases.
    """
    return {"A": 0, "C": 1, "G": 2, "T": 3}.get(base, -1)


def _sequence_to_codes(seq: str) -> np.ndarray:
    """Encode a normalized DNA sequence to a uint8 array of 2-bit codes.

    Uses a 256-entry lookup table for fast vectorized encoding.
    Returns an array where ACGT -> 0/1/2/3 and non-ACGT -> 255 (marker).
    """
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def _compute_kmer_code(seq: str, normalize: bool = True) -> int | None:
    """Encode a k-mer to a single integer using 2-bit encoding.

    Returns None if the sequence contains non-ACGT bases or is too long (>31 bases).
    """
    if normalize:
        seq = _normalize_sequence(seq)

    if len(seq) > 31:
        return None

    code = 0
    for base in seq:
        base_code = _encode_base(base)
        if base_code == -1:
            return None
        code = (code << 2) | base_code

    return code


def _rolling_kmer_codes(seq_codes: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute rolling k-mer codes from a sequence's 2-bit code array.

    Fully vectorized: k passes over the array, each a whole-array operation.

    Args:
        seq_codes: uint8 array where ACGT=0/1/2/3, invalid=255
        k: k-mer length

    Returns:
        codes: int64 array of k-mer codes (length n-k+1)
        valid: bool array marking valid k-mers (no invalid bases in window)
    """
    n = len(seq_codes)
    if n < k:
        return np.array([], dtype=np.int64), np.array([], dtype=bool)

    m = n - k + 1

    # Replace invalid markers with 0 for arithmetic (but track them for validity)
    safe = np.where(seq_codes == _INVALID, 0, seq_codes).astype(np.int64)

    # Rolling codes via k vectorized passes
    codes = np.zeros(m, dtype=np.int64)
    for j in range(k):
        codes = (codes << 2) | safe[j : j + m]

    return codes, _window_validity(seq_codes, k)


def _window_validity(seq_codes: np.ndarray, k: int) -> np.ndarray:
    """Mark length-k windows containing no invalid base, via a prefix sum.

    Args:
        seq_codes: uint8 array where ACGT=0/1/2/3, invalid=255
        k: window length

    Returns:
        Boolean array of length n-k+1.
    """
    inv = (seq_codes == _INVALID).astype(np.int32)
    csum = np.concatenate(([0], np.cumsum(inv)))
    valid: np.ndarray = (csum[k:] - csum[:-k]) == 0
    return valid


def _rolling_kmer_codes_multi(
    seq_codes: np.ndarray, lengths: Sequence[int]
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Compute rolling codes for several k at once, deriving short k from the longest.

    Codes are big-endian prefix-aligned, so the length-k code of a window is the top 2k bits
    of that window's length-kmax code. One accumulation pass therefore serves every length,
    instead of one pass per length.

    Args:
        seq_codes: uint8 array where ACGT=0/1/2/3, invalid=255
        lengths: distinct window lengths to compute

    Returns:
        Mapping from length to (codes, valid) as returned by :func:`_rolling_kmer_codes`.
    """
    usable = sorted({k for k in lengths if k <= len(seq_codes)})
    if not usable:
        return {}

    kmax = usable[-1]
    max_codes, _ = _rolling_kmer_codes(seq_codes, kmax)

    result: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for k in usable:
        # Windows of length k start at every position a kmax window does, plus the kmax-k
        # trailing starts that no kmax window covers -- those need their own accumulation.
        shifted = max_codes >> (2 * (kmax - k))
        if k == kmax:
            codes = shifted
        else:
            tail, _ = _rolling_kmer_codes(seq_codes[len(seq_codes) - (kmax - k) - k + 1 :], k)
            codes = np.concatenate((shifted, tail))
        result[k] = (codes, _window_validity(seq_codes, k))
    return result


@dataclass(frozen=True)
class RepeatObservation:
    """Observation of a sequence's occurrence rate in a reference transcriptome.

    Attributes:
        sequence: Normalized sequence (uppercase, T not U).
        transcript_count: Number of distinct reference transcripts containing the sequence.
        transcript_fraction: Fraction of reference transcripts containing the sequence.
        is_repeat: Verdict: True if transcript_fraction > threshold.
    """

    sequence: str
    transcript_count: int
    transcript_fraction: float
    is_repeat: bool


@dataclass(frozen=True)
class RepeatScanResult:
    """Result of scanning a set of sequences against a reference transcriptome.

    Attributes:
        reference_transcript_count: Total number of transcripts in the reference.
        threshold_fraction: Threshold used for the is_repeat verdict.
        observations: Mapping from normalized sequence to its observation.
    """

    reference_transcript_count: int
    threshold_fraction: float
    observations: dict[str, RepeatObservation]

    @property
    def repeat_sequences(self) -> frozenset[str]:
        """Set of normalized sequences flagged as repeats."""
        return frozenset(seq for seq, obs in self.observations.items() if obs.is_repeat)


class RepeatDetector:
    """Detector for repeat elements in siRNA candidate sequences.

    Identifies guide sequences that occur frequently across a reference transcriptome,
    using a fraction-of-transcripts threshold to ensure consistent behavior across
    references of different sizes.

    Attributes:
        threshold_fraction: Fraction of reference transcripts above which a sequence
            is flagged as a repeat.
    """

    def __init__(self, threshold_fraction: float = DEFAULT_REPEAT_TRANSCRIPT_FRACTION) -> None:
        """Initialize the repeat detector.

        Args:
            threshold_fraction: Threshold for repeat flagging. Default is 0.001 (0.1%).
        """
        self.threshold_fraction = threshold_fraction

    def scan(
        self,
        sequences: Iterable[str],
        reference_fasta: Path,
    ) -> RepeatScanResult:
        """Scan a set of sequences against a reference transcriptome.

        Args:
            sequences: Guide sequences to scan (may be RNA or DNA).
            reference_fasta: Path to reference transcriptome FASTA file.

        Returns:
            RepeatScanResult with observations for each sequence.
        """
        needle_set = {_normalize_sequence(seq) for seq in sequences}

        if not needle_set:
            return RepeatScanResult(
                reference_transcript_count=0,
                threshold_fraction=self.threshold_fraction,
                observations={},
            )

        if not reference_fasta.exists():
            logger.warning(f"Reference FASTA not found: {reference_fasta}")
            return self._empty_result(needle_set)

        length_data = self._prepare_needle_index(needle_set)
        transcript_hits: dict[str, set[str]] = {needle: set() for needle in needle_set}
        transcript_count = 0

        try:
            # Batch transcripts for vectorized processing
            max_needle_len = max((k for k in length_data), default=0)
            buffer = bytearray()
            transcript_map: list[tuple[str, int]] = []  # (transcript_id, start_offset)

            for transcript_id, sequence in _iter_fasta_records(reference_fasta):
                transcript_count += 1
                norm_seq = _normalize_sequence(sequence)

                # Add separator if buffer is not empty
                if buffer:
                    # Add (max_needle_len - 1) N's as sentinel to prevent cross-transcript matches
                    buffer.extend(b"N" * (max_needle_len - 1))

                transcript_map.append((transcript_id, len(buffer)))
                buffer.extend(norm_seq.encode("ascii"))

                # Flush buffer if it exceeds chunk size
                if len(buffer) >= _CHUNK_SIZE:
                    self._process_batch(buffer, transcript_map, length_data, transcript_hits)
                    buffer.clear()
                    transcript_map.clear()

            # Flush remaining buffer
            if buffer:
                self._process_batch(buffer, transcript_map, length_data, transcript_hits)

        except Exception as e:
            logger.warning(f"Error reading reference FASTA {reference_fasta}: {e}")
            return self._empty_result(needle_set)

        return self._build_result(needle_set, transcript_hits, transcript_count)

    def _prepare_needle_index(
        self, needle_set: set[str]
    ) -> dict[int, tuple[np.ndarray, dict[int, set[str]], np.ndarray]]:
        """Prepare sorted code arrays and mappings for each needle length.

        Needles longer than 31 bases are skipped with a warning (code width guard).
        """
        needles_by_length: dict[int, set[str]] = {}
        for needle in needle_set:
            if len(needle) > 31:
                logger.warning(f"Needle longer than 31 bases skipped (code width limit): {needle[:40]}...")
                continue
            needles_by_length.setdefault(len(needle), set()).add(needle)

        length_data: dict[int, tuple[np.ndarray, dict[int, set[str]], np.ndarray]] = {}
        for length, needles in needles_by_length.items():
            codes_to_needles: dict[int, set[str]] = {}
            all_codes = []

            for needle in needles:
                for seq in [needle, _reverse_complement(needle)]:
                    code = _compute_kmer_code(seq, normalize=False)
                    if code is not None:
                        codes_to_needles.setdefault(code, set()).add(needle)
                        all_codes.append(code)

            if all_codes:
                sorted_codes = np.array(sorted(set(all_codes)), dtype=np.int64)
                prefilter = np.zeros(_PREFILTER_SIZE, dtype=bool)
                prefilter[sorted_codes & _PREFILTER_MASK] = True
                length_data[length] = (sorted_codes, codes_to_needles, prefilter)

        return length_data

    def _process_batch(
        self,
        buffer: bytearray,
        transcript_map: list[tuple[str, int]],
        length_data: dict[int, tuple[np.ndarray, dict[int, set[str]], np.ndarray]],
        transcript_hits: dict[str, set[str]],
    ) -> None:
        """Process a batch of transcripts in one vectorized pass."""
        if not buffer or not transcript_map:
            return

        # Encode the entire buffer
        seq_codes = _sequence_to_codes(buffer.decode("ascii"))

        # Build transcript offset array for searchsorted
        transcript_starts = np.array([offset for _, offset in transcript_map], dtype=np.int64)

        # One accumulation pass serves every needle length (codes are prefix-aligned)
        rolled = _rolling_kmer_codes_multi(seq_codes, list(length_data))

        for length, (sorted_codes, codes_to_needles, prefilter) in length_data.items():
            rolled_entry = rolled.get(length)
            if rolled_entry is None:
                continue
            kmer_codes, valid_mask = rolled_entry
            if len(kmer_codes) == 0:
                continue

            positions, codes = self._find_match_positions(kmer_codes, valid_mask, sorted_codes, prefilter)
            if positions.size == 0:
                continue

            # Attribute every match to its transcript in one searchsorted, then collapse to
            # distinct (transcript, code) pairs -- a transcript counts once however many windows
            # matched, so the Python loop below is over pairs, not over windows.
            owners = np.searchsorted(transcript_starts, positions, side="right") - 1
            pairs = np.unique(np.stack((owners, codes), axis=1), axis=0)
            for owner_idx, code in pairs.tolist():
                if 0 <= owner_idx < len(transcript_map):
                    transcript_id = transcript_map[owner_idx][0]
                    for needle in codes_to_needles[code]:
                        transcript_hits[needle].add(transcript_id)

    def _find_match_positions(
        self,
        kmer_codes: np.ndarray,
        valid_mask: np.ndarray,
        sorted_codes: np.ndarray,
        prefilter: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Find matching window positions and their codes using binary search.

        Args:
            kmer_codes: All k-mer codes (including invalid)
            valid_mask: Boolean mask of valid k-mers
            sorted_codes: Sorted array of needle codes
            prefilter: Optional low-bit membership bitmap for the same needle codes

        Returns:
            Parallel arrays of (window positions, matched codes). Stays in numpy so the caller
            can attribute matches to transcripts without a per-window Python call.
        """
        empty = np.empty(0, dtype=np.int64)
        valid_codes = kmer_codes[valid_mask]
        if valid_codes.size == 0:
            return empty, empty

        valid_positions = np.flatnonzero(valid_mask)

        # Prefilter through a small resident bitmap keyed on the code's low bits. searchsorted
        # over every window costs ~10 random probes each; the bitmap is one gather into a few MB
        # and lets through only the needles plus a fraction of a percent of false positives.
        if prefilter is not None:
            survivors = prefilter[valid_codes & _PREFILTER_MASK]
            valid_codes = valid_codes[survivors]
            valid_positions = valid_positions[survivors]
            if valid_codes.size == 0:
                return empty, empty

        indices = np.searchsorted(sorted_codes, valid_codes)
        # Clip rather than mask so the equality test below does the filtering in one pass.
        clipped = np.minimum(indices, len(sorted_codes) - 1)
        matches = sorted_codes[clipped] == valid_codes

        return valid_positions[matches], valid_codes[matches]

    def _build_result(
        self, needle_set: set[str], transcript_hits: dict[str, set[str]], transcript_count: int
    ) -> RepeatScanResult:
        """Build observations from hit counts."""
        observations = {}
        for needle in needle_set:
            count = len(transcript_hits[needle])
            fraction = count / transcript_count if transcript_count > 0 else 0.0
            observations[needle] = RepeatObservation(
                sequence=needle,
                transcript_count=count,
                transcript_fraction=fraction,
                is_repeat=fraction > self.threshold_fraction,
            )

        return RepeatScanResult(
            reference_transcript_count=transcript_count,
            threshold_fraction=self.threshold_fraction,
            observations=observations,
        )

    def _empty_result(self, needles: set[str]) -> RepeatScanResult:
        """Build an empty result with all is_repeat=False."""
        observations = {
            needle: RepeatObservation(
                sequence=needle,
                transcript_count=0,
                transcript_fraction=0.0,
                is_repeat=False,
            )
            for needle in needles
        }
        return RepeatScanResult(
            reference_transcript_count=0,
            threshold_fraction=self.threshold_fraction,
            observations=observations,
        )
