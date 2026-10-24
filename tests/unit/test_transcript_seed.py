"""Oracle fixtures for canonical transcript-seed site discovery."""

import pytest

from sirnaforge.core.transcript_seed import scan_transcript_seed_sites
from sirnaforge.models.transcript_seed import TranscriptSeedTranscript, TranscriptSequenceScope

GUIDE = "AACCTGTAAG"


@pytest.mark.parametrize(
    ("target", "expected_sites"),
    [
        ("ACAGGT", [("6mer", 1, 6)]),
        ("ACAGGTA", [("7mer-A1", 1, 7)]),
        ("TACAGGT", [("7mer-m8", 1, 7), ("6mer", 2, 7)]),
        ("TACAGGTA", [("8mer", 1, 8), ("7mer-A1", 2, 8)]),
    ],
)
def test_seed_class_and_antiparallel_guide_geometry(
    target: str,
    expected_sites: list[tuple[str, int, int]],
) -> None:
    """Report matches at each distinct transcript register."""
    result = scan_transcript_seed_sites(
        {"g": GUIDE},
        [TranscriptSeedTranscript(transcript_id="tx", gene_id="gene", cdna_sequence=target)],
    )
    assert [(site.site_class.value, site.start, site.end) for site in result.sites] == expected_sites


def test_rna_dna_normalization_and_reverse_complement_orientation() -> None:
    """Normalize U/T while preserving antiparallel target matching."""
    guide_rna = GUIDE.replace("T", "U")
    transcript_rna = "GGGUACAGGUAGGG"
    result = scan_transcript_seed_sites(
        {"g": guide_rna},
        [TranscriptSeedTranscript(transcript_id="tx", cdna_sequence=transcript_rna)],
    )
    assert [(site.site_class.value, site.start, site.end) for site in result.sites] == [
        ("8mer", 4, 11),
        ("7mer-A1", 5, 11),
    ]
    same_orientation = scan_transcript_seed_sites(
        {"g": GUIDE}, [TranscriptSeedTranscript(transcript_id="sense", cdna_sequence=GUIDE)]
    )
    assert same_orientation.site_count == 0


def test_deduplicates_same_site_class_but_preserves_distinct_positions_and_guides() -> None:
    """Keep each guide/site identity once while counting repeated positions."""
    result = scan_transcript_seed_sites(
        {"g1": GUIDE, "g2": GUIDE},
        [TranscriptSeedTranscript(transcript_id="tx", cdna_sequence="ACAGGTCCACAGGT")],
    )
    assert result.site_count == 4
    assert {(site.guide_id, site.start, site.site_class.value) for site in result.sites} == {
        (guide, start, "6mer") for guide in ("g1", "g2") for start in (1, 9)
    }
    assert (result.transcript_count, result.gene_count, result.unresolved_gene_count) == (1, 0, 1)


def test_three_prime_utr_scope_restores_transcript_coordinates_and_tracks_missing_annotation() -> None:
    """Scope matches to annotated UTR bases and report absent annotations as unavailable."""
    result = scan_transcript_seed_sites(
        {"g": GUIDE},
        [
            TranscriptSeedTranscript(
                transcript_id="annotated",
                gene_id="gene",
                cdna_sequence="GGGGTACAGGTAGGG",
                utr3_sequence="TACAGGTA",
                utr3_start=5,
            ),
            TranscriptSeedTranscript(transcript_id="missing", cdna_sequence="TACAGGTA"),
            TranscriptSeedTranscript(transcript_id="known_empty", cdna_sequence="GGGG", utr3_sequence="", utr3_start=5),
        ],
        scope=TranscriptSequenceScope.THREE_PRIME_UTR,
    )
    assert [(site.transcript_id, site.start, site.end, site.site_class.value) for site in result.sites] == [
        ("annotated", 5, 12, "8mer"),
        ("annotated", 6, 12, "7mer-A1"),
    ]
    assert result.transcripts_scanned == 2
    assert result.transcripts_unavailable == ("missing",)
    assert (result.site_count, result.transcript_count, result.gene_count) == (2, 1, 1)


def test_malformed_input_and_mismatched_utr_interval_fail_explicitly() -> None:
    """Reject ambiguous bases and UTR sequences that disagree with the cDNA interval."""
    with pytest.raises(ValueError, match="unsupported nucleotide"):
        scan_transcript_seed_sites(
            {"g": "ACGTACGN"}, [TranscriptSeedTranscript(transcript_id="tx", cdna_sequence="ACGT")]
        )
    with pytest.raises(ValueError, match="does not match its cDNA interval"):
        scan_transcript_seed_sites(
            {"g": GUIDE},
            [TranscriptSeedTranscript(transcript_id="tx", cdna_sequence="AAAA", utr3_sequence="TTTT", utr3_start=1)],
            scope=TranscriptSequenceScope.THREE_PRIME_UTR,
        )
