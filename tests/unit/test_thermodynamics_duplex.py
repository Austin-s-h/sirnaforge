"""Regression tests for guide:passenger duplex thermodynamics.

These tests pin the strand orientation used by :class:`ThermodynamicCalculator`.
Before the fix for issue #78 both ViennaRNA call sites reverse-complemented the
passenger strand -- which is *already* the guide's complement -- so every duplex
collapsed to a self-fold. At the default 21 nt design length the two end folds
became literally identical, which pinned ``asymmetry_score`` to exactly 0.5 for
every candidate.
"""

import importlib.util

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from sirnaforge.core.thermodynamics import END_WINDOW_NT, ThermodynamicCalculator
from sirnaforge.models.sirna import SiRNACandidate

VIENNA_AVAILABLE = importlib.util.find_spec("RNA") is not None

pytestmark = pytest.mark.skipif(not VIENNA_AVAILABLE, reason="ViennaRNA not available")

# Real human target windows (TP53 and GAPDH); passenger is the exact reverse complement.
GUIDE_21 = "CAAATTTCCTTCCACTCGGAT"
PASSENGER_21 = "ATCCGAGTGGAAGGAAATTTG"


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


def _cofold_mfe(strand_a: str, strand_b: str) -> float:
    """Independently fold two strands so tests do not reuse the implementation."""
    import RNA  # noqa: PLC0415 - optional dependency, imported after the skip guard

    model_details = RNA.md()
    model_details.temperature = 37.0
    _, mfe = RNA.fold_compound(f"{strand_a}&{strand_b}", model_details).mfe()
    return float(mfe)


def _guides_from(fasta_path, length: int = 21, limit: int = 40) -> list[str]:
    """Return distinct guide sequences tiled across the first record of a FASTA."""
    sequence = str(next(SeqIO.parse(fasta_path, "fasta")).seq).upper()
    guides = []
    for start in range(0, len(sequence) - length + 1, 11):
        guides.append(str(Seq(sequence[start : start + length]).reverse_complement()))
        if len(guides) == limit:
            break
    return guides


@pytest.mark.unit
def test_complementary_duplex_is_stable():
    """A fully complementary 21mer duplex must have a strongly negative ΔG."""
    calc = ThermodynamicCalculator(temperature=37.0)

    dg = calc.calculate_duplex_stability(GUIDE_21, PASSENGER_21)

    assert dg < 0.0, f"complementary duplex reported as non-binding: {dg}"
    # 21 canonical base pairs; the self-fold bug reported roughly -4 to -11 kcal/mol.
    assert dg < -20.0, f"ΔG={dg} is too weak for 21 paired bases"


@pytest.mark.unit
def test_duplex_stability_pairs_the_two_strands_not_the_guide_with_itself():
    """The full duplex must fold guide against passenger, not against a copy of itself."""
    calc = ThermodynamicCalculator(temperature=37.0)

    dg = calc.calculate_duplex_stability(GUIDE_21, PASSENGER_21)

    assert dg == pytest.approx(_cofold_mfe(GUIDE_21, PASSENGER_21))
    assert dg != pytest.approx(_cofold_mfe(GUIDE_21, GUIDE_21))


@pytest.mark.unit
def test_end_stabilities_are_negative_for_a_well_formed_duplex():
    """Both duplex ends pair, so both end ΔG values must be negative."""
    calc = ThermodynamicCalculator(temperature=37.0)

    dg_5p, dg_3p, _ = calc.calculate_asymmetry_score(_candidate(GUIDE_21))

    assert dg_5p < 0.0, f"5' end ΔG must be negative, got {dg_5p}"
    assert dg_3p < 0.0, f"3' end ΔG must be negative, got {dg_3p}"


@pytest.mark.unit
def test_end_windows_pair_opposite_ends_of_the_duplex():
    """The guide 5' end pairs with the passenger 3' end, and vice versa."""
    calc = ThermodynamicCalculator(temperature=37.0)
    candidate = _candidate(GUIDE_21)
    guide, passenger = candidate.guide_sequence, candidate.passenger_sequence

    dg_5p, dg_3p, _ = calc.calculate_asymmetry_score(candidate)

    assert dg_5p == pytest.approx(_cofold_mfe(guide[:END_WINDOW_NT], passenger[-END_WINDOW_NT:]))
    assert dg_3p == pytest.approx(_cofold_mfe(guide[-END_WINDOW_NT:], passenger[:END_WINDOW_NT]))


@pytest.mark.unit
@pytest.mark.parametrize("length", [19, 20, 21, 22, 23])
def test_end_windows_are_symmetric_for_every_design_length(length):
    """Both ends use the same window width, so their ΔG values stay comparable."""
    calc = ThermodynamicCalculator(temperature=37.0)
    candidate = _candidate(GUIDE_21[:length] if length <= 21 else GUIDE_21 + "GC"[: length - 21])
    guide, passenger = candidate.guide_sequence, candidate.passenger_sequence

    dg_5p, dg_3p, _ = calc.calculate_asymmetry_score(candidate)

    window = min(END_WINDOW_NT, length)
    assert dg_5p == pytest.approx(_cofold_mfe(guide[:window], passenger[-window:]))
    assert dg_3p == pytest.approx(_cofold_mfe(guide[-window:], passenger[:window]))


@pytest.mark.unit
def test_asymmetry_score_is_not_a_constant(realistic_transcripts_fasta):
    """asymmetry_score must discriminate between candidates.

    The self-fold bug made the 5' and 3' folds identical at the default 21 nt
    length, so every candidate scored exactly 0.5.
    """
    calc = ThermodynamicCalculator(temperature=37.0)
    guides = _guides_from(realistic_transcripts_fasta)
    assert len(guides) >= 20, "need a spread of real candidates to test discrimination"

    scores = [calc.calculate_asymmetry_score(_candidate(g, f"c_{i}"))[2] for i, g in enumerate(guides)]

    assert len(set(scores)) > 1, "asymmetry_score is constant across candidates"
    assert not all(score == 0.5 for score in scores), "asymmetry_score collapsed to the 0.5 midpoint"


@pytest.mark.unit
def test_melting_temperature_is_physically_plausible():
    """A 21mer siRNA duplex melts around 60-80 °C, not above 100 °C.

    Scaling the corrected duplex ΔG (about -39 kcal/mol) by the old empirical
    factor returned 101-124 °C.
    """
    calc = ThermodynamicCalculator(temperature=37.0)

    tm = calc.calculate_melting_temperature(GUIDE_21, PASSENGER_21)

    assert 50.0 < tm < 90.0, f"implausible melting temperature: {tm} °C"


@pytest.mark.unit
def test_melting_temperature_tracks_gc_content():
    """Tm must rise with GC content."""
    calc = ThermodynamicCalculator(temperature=37.0)
    gc_rich = _candidate("GCGCGCGCGGCCGGCCGCGCG")
    au_rich = _candidate("ATATATATAATTAATTATATA")

    tm_gc = calc.calculate_melting_temperature(gc_rich.guide_sequence, gc_rich.passenger_sequence)
    tm_au = calc.calculate_melting_temperature(au_rich.guide_sequence, au_rich.passenger_sequence)

    assert tm_gc > tm_au


@pytest.mark.unit
def test_melting_temperature_is_alphabet_agnostic():
    """The same duplex spelled with T or U must give the same Tm."""
    calc = ThermodynamicCalculator(temperature=37.0)

    tm_dna = calc.calculate_melting_temperature(GUIDE_21, PASSENGER_21)
    tm_rna = calc.calculate_melting_temperature(GUIDE_21.replace("T", "U"), PASSENGER_21.replace("T", "U"))

    assert tm_dna == pytest.approx(tm_rna)


@pytest.mark.unit
def test_melting_temperature_rejects_length_mismatch():
    """A passenger of a different length is a caller error, not a silent result."""
    calc = ThermodynamicCalculator(temperature=37.0)

    with pytest.raises(ValueError, match="same length"):
        calc.calculate_melting_temperature(GUIDE_21, PASSENGER_21[:-1])


@pytest.mark.unit
def test_asymmetry_score_follows_end_stability_difference():
    """An AU-rich guide 5' end (less stable) must score above the 0.5 midpoint."""
    calc = ThermodynamicCalculator(temperature=37.0)

    # 5' end AU-rich / 3' end GC-rich -> destabilised 5' end -> favoured for RISC loading.
    au_five_prime = _candidate("AAUAAUAAUCAUCGGCCGGCC".replace("U", "T"))
    dg_5p, dg_3p, favourable = calc.calculate_asymmetry_score(au_five_prime)
    assert dg_5p > dg_3p
    assert favourable > 0.5

    # Reversed composition -> stable 5' end -> penalised.
    gc_five_prime = _candidate("GGCCGGCCGACTAATAATAAT")
    dg_5p, dg_3p, unfavourable = calc.calculate_asymmetry_score(gc_five_prime)
    assert dg_5p < dg_3p
    assert unfavourable < 0.5
