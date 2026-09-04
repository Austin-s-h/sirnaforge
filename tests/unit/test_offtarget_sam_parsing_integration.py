"""End-to-end test of the SAM frame fix: BWA record -> mismatch strata -> filter verdict.

The unit tests in ``test_offtarget_sam_parsing.py`` pin ``_parse_sam_output`` in isolation, but
the defect that mattered lives at the module boundary: ``nm`` is what
``SiRNAWorkflow._integrate_offtarget_results`` stratifies into
``transcriptome_hits_{0,1,2}mm``, and those counters are what
``OffTargetFilterCriteria.max_transcriptome_hits_*`` gates on. Redefining ``nm`` therefore flips
pass/fail verdicts, and nothing crossed that boundary in a test before.
"""

import pytest

from sirnaforge.core.off_target import BwaAnalyzer
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate
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


def _strata_and_verdict(tmp_path, sam: str) -> SiRNACandidate:
    """Push a SAM blob through BWA parsing and off-target integration onto one candidate."""
    analyzer = BwaAnalyzer("/nonexistent/index", mode="transcriptome")
    hits = analyzer._filter_and_rank(analyzer._parse_sam_output(sam, {"cand_1": GUIDE}))

    workflow = SiRNAWorkflow(
        WorkflowConfig(output_dir=tmp_path / "out", gene_query="tp53", design_params=DesignParameters())
    )
    candidate = _candidate()
    workflow._integrate_offtarget_results(
        [candidate],
        {"status": "completed", "results": {"cand_1": {"hits": hits}}},
    )
    return candidate


@pytest.mark.unit
def test_full_length_perfect_hits_still_fail_as_perfect_matches(tmp_path):
    """The control: two genuine 21/21 off-target hits must still trip the 0mm gate."""
    candidate = _strata_and_verdict(tmp_path, _sam(cigar="21M", md="21", nm=0, rnames=("ENST_A", "ENST_B")))

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
    candidate = _strata_and_verdict(tmp_path, _sam(cigar="15M6S", md="15", nm=0, rnames=("ENST_A", "ENST_B")))

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
    candidate = _strata_and_verdict(tmp_path, _sam(cigar="20M1S", md="20", nm=0, rnames=("ENST_A",)))

    assert candidate.transcriptome_hits_0mm == 0
    assert candidate.transcriptome_hits_1mm == 1
    assert candidate.passes_filters is True


@pytest.mark.unit
def test_deleted_bases_are_stratified_and_never_score_as_zero_risk(tmp_path):
    """A 2bp deletion is a 2mm-equivalent hit, not a perfect match at the head of the ranking."""
    candidate = _strata_and_verdict(tmp_path, _sam(cigar="10M2D11M", md="10^AC11", nm=2, rnames=("ENST_A",)))

    assert candidate.transcriptome_hits_0mm == 0
    assert candidate.transcriptome_hits_2mm == 1
    assert candidate.passes_filters is True
