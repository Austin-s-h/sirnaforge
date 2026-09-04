"""Tests for SAM record parsing in :class:`BwaAnalyzer`.

``design.py`` builds every guide as ``reverse_complement(target_seq)``, so BWA reports
the majority of genuine transcriptome hits on the minus strand. For those records SAM
SEQ/MD are written in the reverse-complemented read frame, while the seed window
(``seed_start``..``seed_end``) is defined in guide coordinates -- if the two frames are
not reconciled the seed-mismatch logic is inverted (issue #80 follow-up, finding 3).
"""

import pytest

from sirnaforge.core.off_target import BwaAnalyzer, _parse_cigar_query_layout

GUIDE = "ATGCGATGCGATGCGATGCGC"
GUIDE_RC = "GCGCATCGCATCGCATCGCAT"


def _sam_line(
    *,
    flag: int,
    cigar: str = "21M",
    md: str | None = None,
    nm: int = 0,
    as_score: int | None = 42,
    qname: str = "cand_1",
    rname: str = "ENST00000000001",
    pos: int = 500,
    seq: str | None = None,
) -> str:
    """Build a single-record SAM line for the 21nt GUIDE query."""
    if seq is None:
        seq = GUIDE_RC if flag & 16 else GUIDE
    fields = [qname, str(flag), rname, str(pos), "60", cigar, "*", "0", "0", seq, "*"]
    fields.append(f"NM:i:{nm}")
    if as_score is not None:
        fields.append(f"AS:i:{as_score}")
    if md is not None:
        fields.append(f"MD:Z:{md}")
    return "\t".join(fields)


def _parse(line: str) -> dict:
    """Parse one SAM line through a transcriptome-mode analyzer."""
    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    hits = analyzer._parse_sam_output(line, {"cand_1": GUIDE})
    assert len(hits) == 1
    return hits[0]


@pytest.mark.unit
def test_plus_strand_md_position_is_guide_position():
    """On the plus strand the MD frame already is the guide frame."""
    # MD "2A18": one mismatch at read position 3, inside the 2-8 seed.
    hit = _parse(_sam_line(flag=0, md="2A18", nm=1))

    assert hit["strand"] == "+"
    assert hit["mismatch_positions"] == [3]
    assert hit["seed_mismatches"] == 1
    assert hit["offtarget_score"] == 15.0


@pytest.mark.unit
def test_minus_strand_md_position_is_mirrored_onto_guide():
    """A minus-strand mismatch at read position 19 is guide position 3 (in the seed)."""
    # 21nt guide: read position p maps to guide position 22 - p, so MD "18A2" -> guide 3.
    hit = _parse(_sam_line(flag=16, md="18A2", nm=1))

    assert hit["strand"] == "-"
    assert hit["mismatch_positions"] == [3]
    assert hit["seed_mismatches"] == 1
    assert hit["offtarget_score"] == 15.0


@pytest.mark.unit
def test_minus_strand_mismatch_outside_seed_stays_outside():
    """The mirror is symmetric: a read-frame seed mismatch is a guide 3' mismatch."""
    # MD "2A18" on the minus strand is read position 3 -> guide position 19.
    hit = _parse(_sam_line(flag=16, md="2A18", nm=1))

    assert hit["mismatch_positions"] == [19]
    assert hit["seed_mismatches"] == 0
    assert hit["offtarget_score"] == 11.0


@pytest.mark.unit
def test_soft_clipped_bases_count_as_mismatches():
    """A 15M6S partial hit must not look like a full-length perfect match."""
    hit = _parse(_sam_line(flag=0, cigar="15M6S", md="15", nm=0))

    # The 6 clipped guide bases never paired with the target.
    assert hit["nm"] == 6
    assert hit["mismatch_positions"] == [16, 17, 18, 19, 20, 21]
    assert hit["seed_mismatches"] == 0
    assert hit["offtarget_score"] > 0.0


@pytest.mark.unit
def test_leading_soft_clip_shifts_md_positions():
    """MD positions are aligned-block relative, so a leading clip offsets them."""
    # 3S18M with MD "2A15": read position 3 within the aligned block is guide position 6.
    hit = _parse(_sam_line(flag=0, cigar="3S18M", md="2A15", nm=1))

    assert hit["nm"] == 4  # 1 substitution + 3 clipped guide bases
    assert hit["mismatch_positions"] == [1, 2, 3, 6]
    assert hit["seed_mismatches"] == 3  # guide positions 2, 3 and 6 fall in the 2-8 seed


@pytest.mark.unit
def test_minus_strand_soft_clip_maps_to_guide_three_prime_end():
    """A leading clip on a minus-strand read is a clip at the guide 3' end."""
    hit = _parse(_sam_line(flag=16, cigar="6S15M", md="15", nm=0))

    assert hit["nm"] == 6
    assert hit["mismatch_positions"] == [16, 17, 18, 19, 20, 21]


@pytest.mark.unit
def test_filter_and_rank_tolerates_missing_as_tag():
    """Records without an AS tag carry as_score=None and must still sort."""
    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    sam = "\n".join(
        [
            _sam_line(flag=0, md="21", nm=0, as_score=None),
            _sam_line(flag=16, md="18A2", nm=1, as_score=40),
        ]
    )

    ranked = analyzer._filter_and_rank(analyzer._parse_sam_output(sam, {"cand_1": GUIDE}))

    assert [hit["offtarget_score"] for hit in ranked] == [0.0, 15.0]
    assert ranked[0]["as_score"] is None


# ---------------------------------------------------------------------------
# Missing CIGAR: the mirror must still fire. Round 1 guarded it on the
# CIGAR-derived query length, which a `*` CIGAR left at 0.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_minus_strand_mirror_survives_an_absent_cigar():
    """A `*` CIGAR must not silently reinstate the read-frame seed positions."""
    hit = _parse(_sam_line(flag=16, cigar="*", md="18A2", nm=1))

    # Without the fallback query length this returned position 19 / seed_mismatches 0 / 11.0,
    # i.e. byte-for-byte the numbers the original finding reported.
    assert hit["mismatch_positions"] == [3]
    assert hit["seed_mismatches"] == 1
    assert hit["offtarget_score"] == 15.0


@pytest.mark.unit
def test_minus_strand_mirror_degrades_loudly_when_the_query_length_is_unknown(caplog):
    """With no CIGAR, no SEQ and an unrecognised qname there is no frame to map into."""
    line = _sam_line(flag=16, cigar="*", md="18A2", nm=1, qname="unknown", seq="*")

    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    with caplog.at_level("WARNING"):
        hit = analyzer._parse_sam_output(line, {})[0]

    assert "mismatch positions stay in the read frame" in caplog.text
    assert hit["mismatch_positions"] == []
    # The aligner's own edit distance is still honoured, and the score floor keeps the hit from
    # sorting ahead of genuine perfect matches despite having no placed positions.
    assert hit["nm"] == 1
    assert hit["offtarget_score"] == 10.0


# ---------------------------------------------------------------------------
# Gapped records. `nm == len(mismatch_positions)` does NOT hold for deletions,
# and a 0.0 score sorts FIRST because _filter_and_rank sorts ascending.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_deletion_is_not_scored_as_zero_risk():
    """A 2bp deletion must not look like a perfect match."""
    hit = _parse(_sam_line(flag=0, cigar="10M2D11M", md="10^AC11", nm=2))

    assert hit["nm"] == 2
    # Deleted *reference* bases have no guide position of their own, so this is the one case
    # where len(mismatch_positions) < nm; the penalty is anchored on the flanking guide base.
    assert hit["mismatch_positions"] == []
    assert hit["offtarget_score"] == 26.0


@pytest.mark.unit
def test_insertion_counts_the_unpaired_guide_base():
    """An inserted guide base has no target base opposite it, exactly like a clip."""
    hit = _parse(_sam_line(flag=0, cigar="10M1I10M", md="20", nm=1))

    assert hit["nm"] == 1
    assert hit["mismatch_positions"] == [11]
    assert hit["offtarget_score"] == 11.0


@pytest.mark.unit
def test_gapped_hit_sorts_behind_a_genuine_perfect_match():
    """max_hits truncation keeps the head of the ascending sort, so order is load-bearing."""
    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    sam = "\n".join(
        [
            _sam_line(flag=0, cigar="10M2D11M", md="10^AC11", nm=2, rname="ENST_GAPPED"),
            _sam_line(flag=0, cigar="21M", md="21", nm=0, rname="ENST_PERFECT"),
        ]
    )

    ranked = analyzer._filter_and_rank(analyzer._parse_sam_output(sam, {"cand_1": GUIDE}))

    assert [hit["rname"] for hit in ranked] == ["ENST_PERFECT", "ENST_GAPPED"]


# ---------------------------------------------------------------------------
# Minus strand + TRAILING clip: the clip lands in the guide 5' seed, which is
# the highest-impact combination and the one round 1 left unpinned.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_minus_strand_trailing_clip_lands_in_the_guide_seed():
    """A trailing clip on a minus-strand read strips the guide 5' end, seed included."""
    hit = _parse(_sam_line(flag=16, cigar="15M6S", md="15", nm=0))

    # Read positions 16-21 mirror onto guide positions 6..1, so 5 of the 7 seed positions (2-8)
    # never paired: this is emphatically not a seed match, hence the large penalty.
    assert hit["mismatch_positions"] == [1, 2, 3, 4, 5, 6]
    assert hit["seed_mismatches"] == 5
    assert hit["nm"] == 6
    assert hit["offtarget_score"] == 98.0


# ---------------------------------------------------------------------------
# mirna_seed mode hands BWA only the extracted seed window, so read positions
# are window positions and need the same seed_start-1 shift the in-process
# backends apply.
# ---------------------------------------------------------------------------


def _parse_mirna_seed(line: str) -> dict:
    """Parse one SAM line through a mirna_seed-mode analyzer (7nt window query)."""
    analyzer = BwaAnalyzer("/nonexistent/index", mode="mirna_seed")
    hits = analyzer._parse_sam_output(line, {"cand_1": GUIDE})
    assert len(hits) == 1
    return hits[0]


@pytest.mark.unit
def test_mirna_seed_mode_shifts_window_positions_onto_guide_coordinates():
    """Window position 2 of the 7nt seed window is guide position 3, not 2."""
    hit = _parse_mirna_seed(_sam_line(flag=0, cigar="7M", md="1A5", nm=1, seq=GUIDE[1:8]))

    assert hit["mismatch_positions"] == [3]
    assert hit["seed_mismatches"] == 1
    assert hit["offtarget_score"] == 15.0


@pytest.mark.unit
def test_mirna_seed_mode_mirrors_within_the_window_not_the_guide():
    """The minus-strand mirror uses the 7nt window length, so it stays inside the seed."""
    hit = _parse_mirna_seed(_sam_line(flag=16, cigar="7M", md="1A5", nm=1, seq=GUIDE[1:8]))

    # Window position 2 of 7 mirrors to window position 6 -> guide position 7. Every position of
    # a seed-window query is inside the seed on either strand; mirroring with the 21nt *guide*
    # length put this at 15 -- outside the seed -- inventing a plus/minus asymmetry.
    assert hit["mismatch_positions"] == [7]
    assert hit["seed_mismatches"] == 1
    assert hit["offtarget_score"] == 15.0


# ---------------------------------------------------------------------------
# _parse_cigar_query_layout branches.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_layout_counts_hard_clips_as_unpaired_read_bases():
    """Hard clips keep read offsets comparable to the original query, and never pair."""
    layout = _parse_cigar_query_layout("3H18M")

    assert layout.read_length == 21
    assert layout.clipped == (1, 2, 3)
    assert layout.md_to_read[1] == 4  # the first MD offset is read position 4
    assert layout.inserted == ()
    assert layout.deletions == ()


@pytest.mark.unit
def test_layout_handles_clipping_at_both_ends():
    """Both-ends clipping: the aligned block sits in the middle of the read."""
    layout = _parse_cigar_query_layout("2S17M2S")

    assert layout.read_length == 21
    assert layout.clipped == (1, 2, 20, 21)
    assert layout.md_to_read[1] == 3


@pytest.mark.unit
@pytest.mark.parametrize("cigar", ["", "*"])
def test_layout_falls_back_to_the_query_length_without_a_cigar(cigar):
    """An absent CIGAR means "assume the whole read aligned ungapped"."""
    layout = _parse_cigar_query_layout(cigar, read_length=21)

    assert layout.read_length == 21
    assert layout.clipped == ()
    assert layout.deletions == ()
    # MD offsets are already read positions in that case.
    assert layout.md_to_read[1] == 1
    assert layout.md_to_read[21] == 21


@pytest.mark.unit
def test_layout_records_deletions_with_a_flanking_anchor():
    """D/N ops consume no read bases, so the penalty anchors on the base 5' of the gap."""
    layout = _parse_cigar_query_layout("10M2D11M")

    assert layout.read_length == 21
    assert layout.deletions == ((10, 2),)
    assert layout.clipped == ()
