"""Tests for the design-stage scoring components and the filters they drive.

Covers the gaps called out in issue #78: the duplex-stability normalisation, the
empirical (Reynolds) rule and its own threshold, the asymmetry gate, and the
miRNA-mode composite score ceiling.
"""

import importlib.util
import inspect

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from pydantic import ValidationError

import sirnaforge.cli as cli_module
from sirnaforge.core.design import DUPLEX_DG_PER_NT_STRONG, DUPLEX_DG_PER_NT_WEAK, MiRNADesigner, SiRNADesigner
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.sirna import (
    DEFAULT_MIN_ASYMMETRY_SCORE,
    EMPIRICAL_SCORE_MAX,
    EMPIRICAL_SCORE_MIN,
    DesignMode,
    DesignParameters,
    FilterCriteria,
    SiRNACandidate,
)
from sirnaforge.workflow import run_sirna_workflow

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
def test_empirical_score_range_matches_the_declared_bounds():
    """The declared EMPIRICAL_SCORE_MIN/MAX must match what the rule can produce.

    `min_empirical_score` is bounded by those constants, so drift between them and
    the rule would reintroduce an unsatisfiable threshold.
    """
    designer = SiRNADesigner(DesignParameters())
    middle = "CGTACGTACGTACGTAC"  # 17 nt between position 1 and position 19
    scores = {
        designer._calculate_empirical_score(_candidate(f"{first}{middle}{last}")) for first in "ACGT" for last in "ACGT"
    }

    assert min(scores) == pytest.approx(EMPIRICAL_SCORE_MIN)
    assert max(scores) == pytest.approx(EMPIRICAL_SCORE_MAX)


@pytest.mark.unit
def test_empirical_score_credits_t_as_the_dna_spelling_of_u():
    """Guides are DNA, so T at position 19 must earn the A/U bonus."""
    designer = SiRNADesigner(DesignParameters())
    middle = "CGTACGTACGTACGTAC"

    with_t = designer._calculate_empirical_score(_candidate(f"G{middle}T"))
    with_a = designer._calculate_empirical_score(_candidate(f"G{middle}A"))
    with_c = designer._calculate_empirical_score(_candidate(f"G{middle}C"))

    assert with_t == pytest.approx(with_a)
    assert with_t > with_c


@pytest.mark.unit
def test_unsatisfiable_empirical_threshold_is_rejected():
    """A threshold above the attainable maximum must fail loudly, not reject everything."""
    with pytest.raises(ValidationError):
        FilterCriteria(min_empirical_score=EMPIRICAL_SCORE_MAX + 0.01)


@pytest.mark.unit
def test_asymmetry_threshold_default_is_shared():
    """FilterCriteria and the thermodynamics gate must not drift apart."""
    gate_default = inspect.signature(ThermodynamicCalculator.is_thermodynamically_favorable).parameters["threshold"]

    assert FilterCriteria().min_asymmetry_score == DEFAULT_MIN_ASYMMETRY_SCORE
    assert gate_default.default == DEFAULT_MIN_ASYMMETRY_SCORE


@pytest.mark.unit
def test_low_asymmetry_label_tracks_the_asymmetry_score(realistic_transcripts_fasta):
    """The LOW_ASYMMETRY label must be consistent with the column it is named after."""
    result = _design_gaphd(realistic_transcripts_fasta)
    threshold = result.parameters.filters.min_asymmetry_score
    labelled = [c for c in result.candidates if c.passes_filters == SiRNACandidate.FilterStatus.LOW_ASYMMETRY]
    passing = [c for c in result.candidates if c.passes_filters is True]

    assert labelled, "expected some candidates to fail the asymmetry gate"
    assert passing, "expected some candidates to pass"
    assert all(c.asymmetry_score < threshold for c in labelled)
    assert all(c.asymmetry_score >= threshold for c in passing)


@pytest.mark.unit
def test_low_empirical_label_tracks_the_empirical_score(realistic_transcripts_fasta):
    """The empirical rule gets its own status, gated on its own threshold.

    The threshold is set explicitly because the shipped default is EMPIRICAL_SCORE_MIN, i.e. the
    gate is inert: the rubric's positional rules sit at the guide 3' end, where measured knockdown
    shows no signal. This exercises the gate mechanism, not the default.
    """
    result = _design_gaphd(realistic_transcripts_fasta, filters=FilterCriteria(min_empirical_score=0.5))
    threshold = result.parameters.filters.min_empirical_score
    labelled = [c for c in result.candidates if c.passes_filters == SiRNACandidate.FilterStatus.LOW_EMPIRICAL_SCORE]

    assert labelled, "expected some candidates to fail the empirical gate"
    assert all(c.component_scores["empirical"] < threshold for c in labelled)
    assert all(c.component_scores["empirical"] >= threshold for c in result.candidates if c.passes_filters is True)


@pytest.mark.unit
def test_pass_is_no_longer_decided_by_two_nucleotides(realistic_transcripts_fasta):
    """Eligibility must not reduce to `guide[0] in {G,C} and guide[18] == 'A'`.

    That two-nucleotide rule reproduced the emitted labels exactly while
    min_asymmetry_score was applied to the empirical score.
    """
    result = _design_gaphd(realistic_transcripts_fasta)
    passing = [c for c in result.candidates if c.passes_filters is True]

    assert passing
    assert any(c.guide_sequence[18] != "A" for c in passing), "pass still requires A at guide position 19"
    assert any(c.guide_sequence[0] not in ("G", "C") for c in passing), "pass still requires G/C at position 1"


@pytest.mark.unit
def test_new_filter_label_survives_csv_schema_validation(realistic_transcripts_fasta, tmp_path):
    """LOW_EMPIRICAL_SCORE must be registered with the candidate CSV schema.

    Threshold set explicitly: the shipped default leaves this gate inert.
    """
    result = _design_gaphd(realistic_transcripts_fasta, filters=FilterCriteria(min_empirical_score=0.5))
    output = tmp_path / "candidates.csv"

    validated = result.save_csv(str(output))

    assert "LOW_EMPIRICAL_SCORE" in set(validated["passes_filters"])


@pytest.mark.unit
def test_design_scores_do_not_pile_up_at_the_ceiling(realistic_transcripts_fasta):
    """The design score must discriminate at the top, not park candidates at 100.

    An earlier miRNA implementation clamped the bonus-inflated score, parking the best candidates
    at exactly 100.0 -- nine of them on this transcript -- so the top of the ranking carried no
    ordering information. design_v4 sums to 1.0 and takes features in [0, 1], so 100 is reachable
    only by a candidate perfect on all three terms; nothing is clamped into it.
    """
    record = next(SeqIO.parse(realistic_transcripts_fasta, "fasta"))
    designer = MiRNADesigner(DesignParameters(design_mode=DesignMode.MIRNA))

    result = designer.design_from_sequence(str(record.seq).upper(), record.id)
    scores = [c.design_score for c in result.candidates if c.design_score is not None]

    assert scores, "no candidate received a design score"
    assert all(0.0 <= score <= 100.0 for score in scores)
    assert sum(1 for score in scores if score == 100.0) == 0, "candidates are still clamped at the ceiling"


@pytest.mark.unit
def test_mirna_mode_does_not_alter_the_design_score(realistic_transcripts_fasta):
    """Issue #96: the design stage scores one vector, so miRNA mode and siRNA mode agree exactly.

    Before the fix miRNA mode folded the biogenesis bonuses into the design-stage composite and
    divided the result by 1.25, so every declared weight was silently scaled by 0.80 and a
    candidate earning no bonus kept only 80% of its score. The biogenesis terms are now declared
    members of postscreen_mirna_v4 and enter only once off_target exists.
    """
    record = next(SeqIO.parse(realistic_transcripts_fasta, "fasta"))
    sequence = str(record.seq).upper()

    base = SiRNADesigner(DesignParameters()).design_from_sequence(sequence, record.id)
    mirna = MiRNADesigner(DesignParameters(design_mode=DesignMode.MIRNA)).design_from_sequence(sequence, record.id)
    base_scores = {c.id: c.design_score for c in base.candidates}

    assert len(mirna.candidates) == len(base.candidates)
    for candidate in mirna.candidates:
        assert candidate.design_score == pytest.approx(base_scores[candidate.id])
        assert candidate.weight_vector == "design_v4"
        # The biogenesis evidence is still recorded -- it is just not scored yet.
        assert candidate.supp_13_16_score is not None
        assert candidate.component_scores["ago_start"] in (0.0, 1.0)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("guide_base", "passenger_base", "expected"),
    [
        ("A", "T", "perfect"),  # DNA spelling of A:U
        ("T", "A", "perfect"),
        ("A", "U", "perfect"),
        ("G", "C", "perfect"),
        ("G", "T", "wobble"),  # DNA spelling of the G:U wobble
        ("G", "U", "wobble"),
        ("A", "G", "mismatch"),
        ("C", "T", "mismatch"),
    ],
)
def test_pos1_pairing_reads_dna_as_rna(guide_base, passenger_base, expected):
    """A:T is a Watson-Crick pair, not a mismatch."""
    designer = MiRNADesigner(DesignParameters(design_mode=DesignMode.MIRNA))

    assert designer._classify_pos1_pairing(guide_base, passenger_base) == expected


@pytest.mark.unit
def test_supplementary_score_counts_t_as_u():
    """Positions 13-16 are AU-rich whether spelled with T or U."""
    designer = MiRNADesigner(DesignParameters(design_mode=DesignMode.MIRNA))

    assert designer._calculate_supplementary_score("GCGCGCGCGCGCTTTTGCGCG") == pytest.approx(1.0)
    assert designer._calculate_supplementary_score("GCGCGCGCGCGCGCGCGCGCG") == pytest.approx(0.0)


@pytest.mark.unit
def test_designed_duplexes_are_perfectly_paired_at_position_1(realistic_transcripts_fasta):
    """The passenger is the exact reverse complement, so position 1 always pairs.

    Reading DNA as DNA classified all 257 A:T pairs on this transcript as
    mismatches, handing them an undeserved pos1_mismatch_bonus.
    """
    record = next(SeqIO.parse(realistic_transcripts_fasta, "fasta"))
    designer = MiRNADesigner(DesignParameters(design_mode=DesignMode.MIRNA))

    result = designer.design_from_sequence(str(record.seq).upper(), record.id)

    assert result.candidates
    assert {c.pos1_pairing_state for c in result.candidates} == {"perfect"}


@pytest.mark.unit
def test_delta_dg_end_sign_convention(realistic_transcripts_fasta):
    """delta_dg_end is dg_5p - dg_3p, so positive means a destabilised 5' end."""
    result = _design_gaphd(realistic_transcripts_fasta)

    for candidate in result.candidates[:25]:
        scores = candidate.component_scores
        assert scores["delta_dg_end"] == pytest.approx(scores["dg_5p"] - scores["dg_3p"])


@pytest.mark.unit
def test_threshold_filters_are_reachable_from_the_entry_points():
    """min_asymmetry_score, max_paired_fraction and min_empirical_score must plumb through.

    All three previously took their model defaults unconditionally: both entry points built
    FilterCriteria from gc_min/gc_max (plus max_poly_runs in the CLI) and never passed the
    thermodynamic or empirical thresholds, so setting them anywhere was a silent no-op -- the
    same defect class as max_off_target_count before it was exposed.
    """
    workflow_params = inspect.signature(run_sirna_workflow).parameters
    for name in ("min_asymmetry_score", "max_paired_fraction", "min_empirical_score"):
        assert name in workflow_params, f"run_sirna_workflow cannot receive {name}"
        assert workflow_params[name].default is None, f"{name} must default to None, not a pinned value"

    cli_source = inspect.getsource(cli_module)
    for flag in ("--min-asymmetry", "--max-paired-fraction", "--min-empirical"):
        assert flag in cli_source, f"{flag} is not exposed on the CLI"


@pytest.mark.unit
def test_threshold_overrides_take_effect_and_stay_validated():
    """An override must reach the model, and an out-of-range value must raise, not apply silently."""
    defaults = FilterCriteria(gc_min=30.0, gc_max=52.0)
    assert defaults.min_asymmetry_score == DEFAULT_MIN_ASYMMETRY_SCORE

    overridden = FilterCriteria(
        gc_min=30.0,
        gc_max=52.0,
        min_asymmetry_score=0.50,
        max_paired_fraction=0.75,
    )
    assert overridden.min_asymmetry_score == pytest.approx(0.50)
    assert overridden.max_paired_fraction == pytest.approx(0.75)

    # Constructed, not model_copy'd, so the declared bounds still hold.
    with pytest.raises(ValidationError):
        FilterCriteria(gc_min=30.0, gc_max=52.0, min_asymmetry_score=0.0)
