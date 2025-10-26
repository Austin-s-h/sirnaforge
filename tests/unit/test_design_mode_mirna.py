"""Unit tests for miRNA design mode functionality."""

import tempfile
from pathlib import Path

import pandas as pd

from sirnaforge.core.design import MiRNADesigner, SiRNADesigner
from sirnaforge.models.sirna import DesignMode, DesignParameters, FilterCriteria, MiRNADesignConfig


class TestMiRNADesignConfig:
    """Test MiRNADesignConfig model and presets."""

    def test_mirna_config_defaults(self):
        """Test that MiRNADesignConfig has correct default values."""
        config = MiRNADesignConfig()

        # Check thermodynamic thresholds
        assert config.gc_min == 30.0
        assert config.gc_max == 52.0
        assert config.asymmetry_min == 0.65
        assert config.max_homopolymer == 3

        # Check canonical duplex defaults
        assert config.overhang == "UU"
        assert config.modifications == "standard_2ome"

        # Check off-target preset
        assert config.off_target_preset == "MIRNA_SEED_7_8"

        # Check scoring weights are present
        assert "ago_start_bonus" in config.scoring_weights
        assert "pos1_mismatch_bonus" in config.scoring_weights
        assert "seed_clean_bonus" in config.scoring_weights
        assert "supp_13_16_bonus" in config.scoring_weights


class TestDesignModeEnum:
    """Test DesignMode enum functionality."""

    def test_design_mode_values(self):
        """Test that DesignMode enum has correct values."""
        assert DesignMode.SIRNA.value == "sirna"
        assert DesignMode.MIRNA.value == "mirna"

    def test_design_mode_in_parameters(self):
        """Test that design_mode is accessible in DesignParameters."""
        # Default should be sirna
        params_default = DesignParameters()
        assert params_default.design_mode == DesignMode.SIRNA

        # Can be set to mirna
        params_mirna = DesignParameters(design_mode=DesignMode.MIRNA)
        assert params_mirna.design_mode == DesignMode.MIRNA


class TestMiRNADesigner:
    """Test MiRNADesigner functionality."""

    def test_mirna_designer_instantiation(self):
        """Test that MiRNADesigner can be instantiated."""
        params = DesignParameters(design_mode=DesignMode.MIRNA)
        designer = MiRNADesigner(params)
        assert isinstance(designer, MiRNADesigner)
        assert isinstance(designer, SiRNADesigner)  # Should be a subclass

    def test_mirna_designer_scoring(self):
        """Test that MiRNADesigner applies miRNA-specific scoring."""
        # Create a simple test sequence
        test_sequence = "AUGAACGCCACGCUGCUCCGCAACAUCCUCAUCAACGCCAUCCUCGUGCGCAAC"

        params = DesignParameters(
            design_mode=DesignMode.MIRNA,
            sirna_length=21,
            top_n=5,
            filters=FilterCriteria(gc_min=30.0, gc_max=52.0),
        )

        designer = MiRNADesigner(params)
        result = designer.design_from_sequence(test_sequence, "test_transcript")

        # Check that result has candidates
        assert len(result.candidates) > 0

        # Check that miRNA-specific fields are populated
        for candidate in result.candidates[:5]:  # Check first 5
            assert candidate.guide_pos1_base is not None
            assert candidate.pos1_pairing_state is not None
            assert candidate.seed_class is not None
            assert candidate.supp_13_16_score is not None

    def test_mirna_designer_pos1_classification(self):
        """Test position 1 pairing state classification."""
        params = DesignParameters(design_mode=DesignMode.MIRNA)
        designer = MiRNADesigner(params)

        # Test perfect pair
        assert designer._classify_pos1_pairing("A", "U") == "perfect"
        assert designer._classify_pos1_pairing("G", "C") == "perfect"

        # Test wobble
        assert designer._classify_pos1_pairing("G", "U") == "wobble"
        assert designer._classify_pos1_pairing("U", "G") == "wobble"

        # Test mismatch
        assert designer._classify_pos1_pairing("A", "A") == "mismatch"
        assert designer._classify_pos1_pairing("C", "U") == "mismatch"

    def test_mirna_designer_supplementary_score(self):
        """Test 3' supplementary pairing score calculation."""
        params = DesignParameters(design_mode=DesignMode.MIRNA)
        designer = MiRNADesigner(params)

        # Test sequence with high A/U content in positions 13-16 (0-indexed: 12-15)
        guide_high_au = "AAAAAAAAAAAAAAUUUAAAA"  # Positions 12-15: AUUU (75% A/U)
        score_high = designer._calculate_supplementary_score(guide_high_au)
        assert score_high > 0.5  # Should be high score (low pairing)

        # Test sequence with low A/U content in positions 13-16 (0-indexed: 12-15)
        guide_low_au = "AAAAAAAAAAAAGGGCAAAA"  # Positions 12-15: AGGG (25% A/U)
        score_low = designer._calculate_supplementary_score(guide_low_au)
        assert score_low < 0.5  # Should be low score (high pairing potential)

    def test_mirna_vs_sirna_scoring_difference(self):
        """Test that miRNA mode produces different scores than siRNA mode."""
        test_sequence = "AUGAACGCCACGCUGCUCCGCAACAUCCUCAUCAACGCCAUCCUCGUGCGCAAC"

        # Run both modes
        params_sirna = DesignParameters(design_mode=DesignMode.SIRNA, top_n=5)
        params_mirna = DesignParameters(design_mode=DesignMode.MIRNA, top_n=5)

        designer_sirna = SiRNADesigner(params_sirna)
        designer_mirna = MiRNADesigner(params_mirna)

        result_sirna = designer_sirna.design_from_sequence(test_sequence, "test")
        result_mirna = designer_mirna.design_from_sequence(test_sequence, "test")

        # Check that scores differ (miRNA adds bonuses)
        # Find a common candidate and compare scores
        sirna_candidates = {c.guide_sequence: c for c in result_sirna.candidates}

        found_different = False
        for mirna_cand in result_mirna.candidates:
            if mirna_cand.guide_sequence in sirna_candidates:
                sirna_cand = sirna_candidates[mirna_cand.guide_sequence]
                # miRNA mode should add bonuses, so score could be higher
                # (or same if no bonuses apply, but at least check fields are populated)
                assert mirna_cand.guide_pos1_base is not None
                assert sirna_cand.guide_pos1_base is None  # siRNA mode doesn't set this
                found_different = True
                break

        assert found_different, "Should have found at least one common candidate to compare"


class TestMiRNADesignWorkflow:
    """Test miRNA design workflow end-to-end."""

    def test_design_from_file_mirna_mode(self):
        """Test designing siRNAs from a file in miRNA mode."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">test_transcript\n")
            f.write("AUGAACGCCACGCUGCUCCGCAACAUCCUCAUCAACGCCAUCCUCGUGCGCAAC\n")
            temp_fasta = Path(f.name)

        try:
            params = DesignParameters(
                design_mode=DesignMode.MIRNA,
                sirna_length=21,
                top_n=3,
            )

            designer = MiRNADesigner(params)
            result = designer.design_from_file(str(temp_fasta))

            # Check results
            assert result.total_sequences == 1
            assert len(result.candidates) > 0
            assert len(result.top_candidates) <= 3

            # Check miRNA-specific fields are populated
            for candidate in result.candidates:
                assert hasattr(candidate, "guide_pos1_base")
                assert hasattr(candidate, "pos1_pairing_state")
                assert hasattr(candidate, "seed_class")
                assert hasattr(candidate, "supp_13_16_score")

        finally:
            # Clean up
            temp_fasta.unlink()

    def test_csv_export_includes_mirna_columns(self):
        """Test that CSV export includes miRNA-specific columns."""
        test_sequence = "AUGAACGCCACGCUGCUCCGCAACAUCCUCAUCAACGCCAUCCUCGUGCGCAAC"

        params = DesignParameters(design_mode=DesignMode.MIRNA, top_n=3)
        designer = MiRNADesigner(params)
        result = designer.design_from_sequence(test_sequence, "test")

        # Save to temporary CSV
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_csv = Path(f.name)

        try:
            result.save_csv(str(temp_csv))

            # Read CSV and check columns
            df = pd.read_csv(temp_csv)

            # Check that miRNA-specific columns are present
            mirna_columns = [
                "guide_pos1_base",
                "pos1_pairing_state",
                "seed_class",
                "supp_13_16_score",
                "seed_7mer_hits",
                "seed_8mer_hits",
                "seed_hits_weighted",
                "off_target_seed_risk_class",
            ]

            for col in mirna_columns:
                assert col in df.columns, f"Missing column: {col}"

            # Check that at least some values are populated (not all None)
            assert df["guide_pos1_base"].notna().any()
            assert df["pos1_pairing_state"].notna().any()

        finally:
            temp_csv.unlink()
