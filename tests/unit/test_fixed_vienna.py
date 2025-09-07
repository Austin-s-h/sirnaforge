"""
Test the fixed ViennaRNA integration.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.sirna import SiRNACandidate

VIENNA_AVAILABLE = importlib.util.find_spec("RNA") is not None


def test_fixed_duplex_stability():
    """Test that duplex stability calculation works without segfault."""
    if not VIENNA_AVAILABLE:
        pytest.skip("ViennaRNA not available")

    calc = ThermodynamicCalculator(temperature=37.0)

    guide = "CAAATTTCCTTCCACTCGGAT"  # Real TP53 sequence
    passenger = "ATCCGAGTGGAAGGAAATTTG"  # Reverse complement

    print(f"Testing guide: {guide}")
    print(f"Testing passenger: {passenger}")

    # This should work now without segfault
    dg = calc.calculate_duplex_stability(guide, passenger)
    assert isinstance(dg, float)
    print(f"✓ Duplex stability: {dg}")


def test_fixed_asymmetry_calculation():
    """Test that asymmetry score calculation works without segfault."""
    if not VIENNA_AVAILABLE:
        pytest.skip("ViennaRNA not available")

    calc = ThermodynamicCalculator(temperature=37.0)

    candidate = SiRNACandidate(
        id="test_1_21",
        transcript_id="test",
        position=1,
        guide_sequence="CAAATTTCCTTCCACTCGGAT",
        passenger_sequence="ATCCGAGTGGAAGGAAATTTG",
        gc_content=52.4,
        length=21,
        asymmetry_score=0.0,
        composite_score=0.0
    )

    # This should work now without segfault
    dg_5p, dg_3p, asymmetry = calc.calculate_asymmetry_score(candidate)
    assert isinstance(dg_5p, float)
    assert isinstance(dg_3p, float)
    assert isinstance(asymmetry, float)
    assert 0.0 <= asymmetry <= 1.0
    print(f"✓ Asymmetry calculation: 5'={dg_5p}, 3'={dg_3p}, asym={asymmetry}")


def test_multiple_sequences():
    """Test multiple sequences to ensure stability."""
    if not VIENNA_AVAILABLE:
        pytest.skip("ViennaRNA not available")

    calc = ThermodynamicCalculator(temperature=37.0)

    # Multiple real sequences from TP53
    sequences = [
        ("CAAATTTCCTTCCACTCGGAT", "ATCCGAGTGGAAGGAAATTTG"),
        ("ATCAAATCATCCATTGCTTGG", "CCAAGCAATGGATGATTTGAT"),
        ("AAATCATCCATTGCTTGGGAC", "GTCCCAAGCAATGGATGATTT"),
        ("GAACCGGAGGGAGCCGCAGTC", "GACTGCGGCTCCCTCCGGTTC"),
        ("ACCGGAGGGAGCCGCAGTCAG", "CTGACTGCGGCTCCCTCCGGT"),
    ]

    for i, (guide, passenger) in enumerate(sequences):
        print(f"Testing sequence {i+1}: {guide}")

        # Test duplex stability
        dg = calc.calculate_duplex_stability(guide, passenger)
        assert isinstance(dg, float)
        print(f"  ✓ Duplex stability: {dg}")

        # Test end stability
        dg_5p = calc._calculate_end_stability(guide[:7], passenger[:7])
        dg_3p = calc._calculate_end_stability(guide[14:], passenger[14:])
        assert isinstance(dg_5p, float)
        assert isinstance(dg_3p, float)
        print(f"  ✓ End stabilities: 5'={dg_5p}, 3'={dg_3p}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
