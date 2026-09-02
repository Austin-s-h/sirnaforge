"""Tests for the design-stage scoring components and the filters they drive.

Covers the gaps called out in issue #78: the duplex-stability normalisation, the
empirical (Reynolds) rule and its own threshold, the asymmetry gate, and the
miRNA-mode composite score ceiling.
"""

import importlib.util

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from sirnaforge.core.design import DUPLEX_DG_PER_NT_STRONG, DUPLEX_DG_PER_NT_WEAK, SiRNADesigner
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate

VIENNA_AVAILABLE = importlib.util.find_spec("RNA") is not None

pytestmark = pytest.mark.skipif(not VIENNA_AVAILABLE, reason="ViennaRNA not available")

GUIDE_21 = "CAAATTTCCTTCCACTCGGAT"


def _candidate(guide: str, candidate_id: str = "test_1_21") -> SiRNACandidate:
    """Build a minimal candidate whose passenger is the guide's reverse complement."""
    passenger = str(Seq(guide).reverse_complement())
    gc = (guide.count("G") + guide.count("C")) / len(guide) * 100
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="test",
        position=1,
        guide_sequence=guide,
        passenger_sequence=passenger,
        gc_content=gc,
        length=len(guide),
        asymmetry_score=0.0,
        composite_score=0.0,
    )


def _design_gaphd(fasta_path, **param_overrides):
    """Run the designer over the packaged GAPDH transcript."""
    record = next(SeqIO.parse(fasta_path, "fasta"))
    designer = SiRNADesigner(DesignParameters(**param_overrides))
    return designer.design_from_sequence(str(record.seq).upper(), record.id)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("dg_per_nt", "expected"),
    [
        (DUPLEX_DG_PER_NT_STRONG, 1.0),
        (DUPLEX_DG_PER_NT_WEAK, 0.0),
        ((DUPLEX_DG_PER_NT_STRONG + DUPLEX_DG_PER_NT_WEAK) / 2, 0.5),
        (-3.0, 1.0),
        (-0.5, 0.0),
    ],
)
def test_duplex_score_maps_per_nucleotide_dg(monkeypatch, dg_per_nt, expected):
    """The duplex score is a clamped linear map of ΔG per nucleotide."""
    monkeypatch.setattr(
        ThermodynamicCalculator,
        "calculate_duplex_stability",
        lambda _self, guide, _passenger: dg_per_nt * len(guide),
    )
    designer = SiRNADesigner(DesignParameters())

    score, dg = designer._calculate_duplex_score(_candidate(GUIDE_21))

    assert score == pytest.approx(expected)
    assert dg == pytest.approx(dg_per_nt * len(GUIDE_21))


@pytest.mark.unit
def test_duplex_score_is_length_normalised(monkeypatch):
    """Two designs with the same ΔG per nucleotide score the same at any length."""
    monkeypatch.setattr(
        ThermodynamicCalculator,
        "calculate_duplex_stability",
        lambda _self, guide, _passenger: -1.9 * len(guide),
    )
    designer = SiRNADesigner(DesignParameters())

    short_score, _ = designer._calculate_duplex_score(_candidate(GUIDE_21[:19]))
    long_score, _ = designer._calculate_duplex_score(_candidate(GUIDE_21 + "GC"))

    assert short_score == pytest.approx(long_score)


@pytest.mark.unit
def test_duplex_score_does_not_saturate(realistic_transcripts_fasta):
    """The duplex score must discriminate, not pin candidates at the window edge.

    The previous [-40, -5] kcal/mol window was calibrated against the self-fold
    ΔG; against the corrected duplex ΔG it put 209/505 candidates at exactly 1.0.
    """
    result = _design_gaphd(realistic_transcripts_fasta)
    scores = [c.component_scores["duplex_stability_score"] for c in result.candidates]

    assert len(scores) > 100, "expected a substantial candidate pool"
    assert not any(score == 1.0 for score in scores), "duplex score is clamped at the strong edge"
    assert max(scores) - min(scores) > 0.3, "duplex score does not discriminate"


@pytest.mark.unit
def test_delta_dg_end_sign_convention(realistic_transcripts_fasta):
    """delta_dg_end is dg_5p - dg_3p, so positive means a destabilised 5' end."""
    result = _design_gaphd(realistic_transcripts_fasta)

    for candidate in result.candidates[:25]:
        scores = candidate.component_scores
        assert scores["delta_dg_end"] == pytest.approx(scores["dg_5p"] - scores["dg_3p"])
