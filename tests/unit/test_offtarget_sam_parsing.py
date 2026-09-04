"""Tests for SAM record parsing in :class:`BwaAnalyzer`.

``design.py`` builds every guide as ``reverse_complement(target_seq)``, so BWA reports
the majority of genuine transcriptome hits on the minus strand. For those records SAM
SEQ/MD are written in the reverse-complemented read frame, while the seed window
(``seed_start``..``seed_end``) is defined in guide coordinates -- if the two frames are
not reconciled the seed-mismatch logic is inverted (issue #80 follow-up, finding 3).
"""

import pytest

from sirnaforge.core.off_target import BwaAnalyzer

GUIDE = "ATGCGATGCGATGCGATGCGC"
GUIDE_RC = "GCGCATCGCATCGCATCGCAT"


def _sam_line(
    *,
    flag: int,
    cigar: str = "21M",
    md: str | None = None,
    nm: int = 0,
    as_score: int | None = 42,
) -> str:
    """Build a single-record SAM line for the 21nt GUIDE query."""
    seq = GUIDE_RC if flag & 16 else GUIDE
    fields = ["cand_1", str(flag), "ENST00000000001", "500", "60", cigar, "*", "0", "0", seq, "*"]
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
