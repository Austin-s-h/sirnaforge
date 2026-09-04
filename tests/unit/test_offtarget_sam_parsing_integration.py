"""End-to-end test of the SAM frame fix: BWA record -> mismatch strata -> filter verdict.

The unit tests in ``test_offtarget_sam_parsing.py`` pin ``_parse_sam_output`` in isolation, but
the defect that mattered lives at the module boundary: ``nm`` is what
``SiRNAWorkflow._integrate_offtarget_results`` stratifies into
``transcriptome_hits_{0,1,2}mm``, and those counters are what
``OffTargetFilterCriteria.max_transcriptome_hits_*`` gates on. Redefining ``nm`` therefore flips
pass/fail verdicts, and nothing crossed that boundary in a test before.
"""

import statistics

import pytest

from sirnaforge.core.off_target import BwaAnalyzer
from sirnaforge.models.off_target import AnalysisSummary
from sirnaforge.models.sirna import DesignParameters, OffTargetFilterCriteria, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

GUIDE = "ATGCGATGCGATGCGATGCGC"
GUIDE_RC = "GCGCATCGCATCGCATCGCAT"


def _candidate() -> SiRNACandidate:
    """A minimal, unscreened, non-repeat-flagged candidate for GUIDE."""
    return SiRNACandidate(
        id="cand_1",
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE_RC,
        gc_content=57.1,
        length=21,
        asymmetry_score=0.7,
        composite_score=50.0,
    )


def _sam(*, cigar: str, md: str, nm: int, rnames: tuple[str, ...]) -> str:
    """Build one minus-strand SAM record per reference name.

    Minus strand because ``core/design.py`` builds every guide as
    ``reverse_complement(target_seq)``, which makes ``flag & 16`` the common case.
    """
    return "\n".join(
        "\t".join(
            [
                "cand_1",
                "16",
                rname,
                "500",
                "60",
                cigar,
                "*",
                "0",
                "0",
                GUIDE_RC,
                "*",
                f"NM:i:{nm}",
                "AS:i:30",
                f"MD:Z:{md}",
            ]
        )
        for rname in rnames
    )


def _strata_and_verdict(
    tmp_path, sam: str, criteria: OffTargetFilterCriteria | None = None
) -> tuple[SiRNACandidate, dict]:
    """Push a SAM blob through BWA parsing and off-target integration onto one candidate.

    ``off_target_score`` is set the way ``_parse_offtarget_results`` sets it -- the MAX
    ``offtarget_score`` over the candidate's hits -- so ``off_target_penalty`` carries a realistic
    value rather than the 0.0 default.
    """
    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    hits = analyzer._filter_and_rank(analyzer._parse_sam_output(sam, {"cand_1": GUIDE}))

    workflow = SiRNAWorkflow(
        WorkflowConfig(output_dir=tmp_path / "out", gene_query="tp53", design_params=DesignParameters())
    )
    candidate = _candidate()
    _, stats = workflow._integrate_offtarget_results(
        [candidate],
        {
            "status": "completed",
            "results": {
                "cand_1": {
                    "hits": hits,
                    "off_target_score": max((hit["offtarget_score"] for hit in hits), default=0.0),
                }
            },
        },
        criteria,
    )
    return candidate, stats


@pytest.mark.unit
def test_full_length_perfect_hits_still_fail_as_perfect_matches(tmp_path):
    """The control: two genuine 21/21 off-target hits must still trip the 0mm gate."""
    candidate, _ = _strata_and_verdict(tmp_path, _sam(cigar="21M", md="21", nm=0, rnames=("ENST_A", "ENST_B")))

    assert candidate.transcriptome_hits_0mm == 2
    assert candidate.passes_filters == SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH


@pytest.mark.unit
def test_clipped_minus_strand_partial_hits_leave_the_perfect_match_stratum(tmp_path):
    """A 15M6S/NM:i:0 hit is a 15/21 partial match, not a perfect match.

    ``bwa mem -T 15`` emits these routinely. Before the frame fix they arrived with ``nm=0`` and
    landed in ``transcriptome_hits_0mm``, failing the candidate as TRANSCRIPTOME_PERFECT_MATCH on
    the strength of a partial hit whose entire seed never paired. They now carry the guide-level
    distance (6), so they fall outside the nm<=2 strata -- but they are still counted in full by
    ``transcriptome_hits_total`` / ``off_target_count``, which ``max_off_target_count`` gates.
    """
    candidate, _ = _strata_and_verdict(tmp_path, _sam(cigar="15M6S", md="15", nm=0, rnames=("ENST_A", "ENST_B")))

    assert candidate.transcriptome_hits_0mm == 0
    assert candidate.transcriptome_hits_1mm == 0
    assert candidate.transcriptome_hits_2mm == 0
    # Not "not an off-target": both hits are still counted, and the seed-perfect counter agrees
    # that a hit whose guide 5' end never paired is not a seed match.
    assert candidate.transcriptome_hits_total == 2
    assert candidate.off_target_count == 2
    assert candidate.transcriptome_hits_seed_0mm == 0
    assert candidate.passes_filters is True


@pytest.mark.unit
def test_single_base_clip_lands_in_the_one_mismatch_stratum(tmp_path):
    """One unpaired guide base is one mismatch-equivalent, so 20M1S is a 1mm hit."""
    candidate, _ = _strata_and_verdict(tmp_path, _sam(cigar="20M1S", md="20", nm=0, rnames=("ENST_A",)))

    assert candidate.transcriptome_hits_0mm == 0
    assert candidate.transcriptome_hits_1mm == 1
    assert candidate.passes_filters is True


@pytest.mark.unit
def test_deleted_bases_are_stratified_and_never_score_as_zero_risk(tmp_path):
    """A 2bp deletion is a 2mm-equivalent hit, not a perfect match at the head of the ranking."""
    candidate, _ = _strata_and_verdict(tmp_path, _sam(cigar="10M2D11M", md="10^AC11", nm=2, rnames=("ENST_A",)))

    assert candidate.transcriptome_hits_0mm == 0
    assert candidate.transcriptome_hits_2mm == 1
    assert candidate.passes_filters is True


# A leading clip on a minus-strand record lands on the guide's 3' end (read positions 1-6 mirror
# to guide positions 16-21), so guide positions 2-8 -- the seed -- pair perfectly. The guide-level
# nm is 6, which is outside every nm<=2 stratum: this is the shape that no
# max_transcriptome_hits_{0,1,2}mm threshold can see.
_SEED_PERFECT_PARTIAL = {"cigar": "6S15M", "md": "15", "nm": 0}


@pytest.mark.unit
def test_seed_perfect_partial_hit_is_counted_but_lands_in_no_mismatch_stratum(tmp_path):
    """Pin the documented gap: seed-perfect partials are invisible to the mismatch gates.

    This is the contract ``docs/models_and_scoring.md`` §1.4 and the CHANGELOG now describe. It is
    not a bug to fix here -- picking a default ceiling is a product decision -- but it must not
    change silently, because on stock defaults ``max_off_target_count`` (15) is the only thing
    standing between a run of these hits and a PASS.
    """
    candidate, _ = _strata_and_verdict(tmp_path, _sam(rnames=("ENST_A", "ENST_B"), **_SEED_PERFECT_PARTIAL))

    # The seed paired perfectly, and both hits are counted in the unstratified totals...
    assert candidate.transcriptome_hits_seed_0mm == 2
    assert candidate.transcriptome_hits_total == 2
    assert candidate.off_target_count == 2
    # ...yet no mismatch stratum contains them, so no max_transcriptome_hits_* gate applies.
    assert (candidate.transcriptome_hits_0mm, candidate.transcriptome_hits_1mm, candidate.transcriptome_hits_2mm) == (
        0,
        0,
        0,
    )
    assert candidate.passes_filters is True

    # Tightening the perfect-match gate all the way to zero still cannot see them -- proof that
    # max_transcriptome_hits_0mm is the wrong knob for this shape ...
    tightened, _ = _strata_and_verdict(
        tmp_path,
        _sam(rnames=("ENST_A", "ENST_B"), **_SEED_PERFECT_PARTIAL),
        OffTargetFilterCriteria(max_transcriptome_hits_0mm=0),
    )
    assert tightened.passes_filters is True

    # ... while the seed gate is the knob that does.
    gated, _ = _strata_and_verdict(
        tmp_path,
        _sam(rnames=("ENST_A", "ENST_B"), **_SEED_PERFECT_PARTIAL),
        OffTargetFilterCriteria(max_transcriptome_seed_perfect=0),
    )
    assert gated.passes_filters == SiRNACandidate.FilterStatus.TRANSCRIPTOME_SEED_PERFECT


@pytest.mark.unit
def test_max_transcriptome_seed_perfect_gates_the_seed_perfect_partial_hits(tmp_path):
    """Setting the opt-in threshold must actually fail the candidate.

    ``max_transcriptome_seed_perfect`` existed on ``OffTargetFilterCriteria`` but was never read by
    ``_check_offtarget_filters``, so a user who set it got silence -- against the two hits above
    even ``max_transcriptome_seed_perfect=0`` passed. Both boundary cases are asserted here so the
    gate cannot silently become inert again, nor turn into a ``>=`` comparison.
    """
    over, stats = _strata_and_verdict(
        tmp_path,
        _sam(rnames=("ENST_A", "ENST_B"), **_SEED_PERFECT_PARTIAL),
        OffTargetFilterCriteria(max_transcriptome_seed_perfect=1),
    )
    assert over.transcriptome_hits_seed_0mm == 2
    assert over.passes_filters == SiRNACandidate.FilterStatus.TRANSCRIPTOME_SEED_PERFECT
    assert stats["failed_transcriptome_seed_perfect"] == 1

    # It is an upper bound, so exactly at the threshold the candidate still passes ...
    at_threshold, _ = _strata_and_verdict(
        tmp_path,
        _sam(rnames=("ENST_A", "ENST_B"), **_SEED_PERFECT_PARTIAL),
        OffTargetFilterCriteria(max_transcriptome_seed_perfect=2),
    )
    assert at_threshold.passes_filters is True
    # ... and it ships off, so no existing run changes verdict.
    assert OffTargetFilterCriteria().max_transcriptome_seed_perfect is None


@pytest.mark.unit
def test_off_target_penalty_is_max_offtarget_score_so_higher_means_safer(tmp_path):
    """``off_target_penalty`` is not a penalty post-screen; pin the direction it really has.

    It stores the MAX ``offtarget_score``, where ``0.0`` is reserved for a full-length exact match
    (the highest-risk hit there is). Widening ``nm`` widened these values, so the field now reports
    a large number for the *safer* candidate -- the opposite of "lower is better".
    """
    perfect, _ = _strata_and_verdict(tmp_path, _sam(cigar="21M", md="21", nm=0, rnames=("ENST_A",)))
    partial, _ = _strata_and_verdict(tmp_path, _sam(rnames=("ENST_A",), **_SEED_PERFECT_PARTIAL))

    assert perfect.off_target_penalty == 0.0
    assert partial.off_target_penalty > 0.0
    # The riskiest hit reports the lowest "penalty".
    assert perfect.off_target_penalty < partial.off_target_penalty
    description = SiRNACandidate.model_fields["off_target_penalty"].description or ""
    assert "lower is better" not in description
    assert "higher = safer" in description


@pytest.mark.unit
def test_summary_mean_mismatches_is_the_mean_of_the_redefined_nm():
    """``*_summary.json``'s ``mean_mismatches`` changed meaning with ``nm``; say so in the model.

    ``run_bwa_alignment_analysis`` derives it as ``fmean(hit.nm for hit in all_hits)``, so once
    ``nm`` became a guide-level distance this stopped being the aligner's mean edit distance. A
    consumer diffing summaries across the 0.6.0 boundary sees the number move with no schema
    change, which is exactly why the field has to carry the warning itself.
    """
    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    hits = analyzer._parse_sam_output(
        _sam(cigar="15M6S", md="15", nm=0, rnames=("ENST_A", "ENST_B")), {"cand_1": GUIDE}
    )

    # Both records carry NM:i:0, so the aligner's mean edit distance is 0.0 ...
    assert statistics.fmean(0 for _ in hits) == 0.0
    # ... while the value the summary actually reports is the guide-level mean.
    assert statistics.fmean(hit["nm"] for hit in hits) == 6.0

    description = AnalysisSummary.model_fields["mean_mismatches"].description or ""
    assert description != "Mean number of mismatches"
    assert "nm" in description and "NM tag" in description
