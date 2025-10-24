"""Test example modification pattern files."""

import json
from pathlib import Path

import pytest

from sirnaforge.models.modifications import ChemicalModification, StrandMetadata
from sirnaforge.modifications import load_metadata


@pytest.mark.unit
@pytest.mark.local_python
class TestExamplePatternFiles:
    """Validate example modification pattern files."""

    def test_example_files_exist(self):
        """Verify all expected example files exist."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
        assert examples_dir.exists(), f"Examples directory not found: {examples_dir}"

        expected_files = [
            "standard_2ome.json",
            "minimal_terminal.json",
            "maximal_stability.json",
            "fda_approved_onpattro.json",
            "README.md",
        ]

        for filename in expected_files:
            filepath = examples_dir / filename
            assert filepath.exists(), f"Expected file not found: {filename}"

    def test_pattern_files_valid_json(self):
        """Verify pattern files are valid JSON."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"

        json_files = [
            "standard_2ome.json",
            "minimal_terminal.json",
            "maximal_stability.json",
            "fda_approved_onpattro.json",
        ]

        for filename in json_files:
            filepath = examples_dir / filename
            with filepath.open() as f:
                data = json.load(f)
                assert isinstance(data, dict), f"{filename} should contain a JSON object"

    def test_fda_approved_onpattro_structure(self):
        """Test FDA-approved Onpattro example has correct structure."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
        filepath = examples_dir / "fda_approved_onpattro.json"

        with filepath.open() as f:
            data = json.load(f)

        # Check top-level fields
        assert "drug_name" in data
        assert data["drug_name"] == "Patisiran (Onpattro)"
        assert "target_gene" in data
        assert data["target_gene"] == "TTR"
        assert "sequences" in data

        # Check guide sequence
        sequences = data["sequences"]
        assert "patisiran_ttr_guide" in sequences
        guide = sequences["patisiran_ttr_guide"]

        # Validate guide can be loaded as StrandMetadata
        guide_metadata = StrandMetadata.model_validate(guide)
        assert guide_metadata.id == "patisiran_ttr_guide"
        assert guide_metadata.sequence == "AUGGAAUACUCUUGGUUAC"
        assert guide_metadata.overhang == "dTdT"
        assert len(guide_metadata.chem_mods) == 1
        assert guide_metadata.chem_mods[0].type == "2OMe"
        assert guide_metadata.chem_mods[0].positions == [1, 4, 6, 11, 13, 16, 19]

    def test_onpattro_sequences_loadable(self):
        """Test Onpattro sequences can be loaded with load_metadata."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
        filepath = examples_dir / "fda_approved_onpattro.json"

        with filepath.open() as f:
            data = json.load(f)

        # Extract just the sequences section
        sequences = data.get("sequences", {})
        assert len(sequences) >= 1

        # Create temporary file with sequences in expected format
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
            json.dump(sequences, tmp)
            tmp_path = Path(tmp.name)

        try:
            # Load using our metadata loader
            metadata = load_metadata(tmp_path)
            assert "patisiran_ttr_guide" in metadata

            guide = metadata["patisiran_ttr_guide"]
            assert isinstance(guide, StrandMetadata)
            assert guide.sequence == "AUGGAAUACUCUUGGUUAC"
            assert len(guide.chem_mods) == 1
        finally:
            tmp_path.unlink()

    def test_pattern_modification_positions_valid(self):
        """Test that all modification positions are valid for their sequences."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
        onpattro_file = examples_dir / "fda_approved_onpattro.json"

        with onpattro_file.open() as f:
            data = json.load(f)

        for seq_id, seq_data in data.get("sequences", {}).items():
            metadata = StrandMetadata.model_validate(seq_data)

            # Verify all modification positions are within sequence length
            for mod in metadata.chem_mods:
                if mod.positions:
                    max_pos = max(mod.positions)
                    seq_len = len(metadata.sequence)
                    assert max_pos <= seq_len, (
                        f"{seq_id}: Modification {mod.type} has position {max_pos} "
                        f"but sequence length is only {seq_len}"
                    )

    def test_standard_pattern_structure(self):
        """Test standard pattern file has expected structure."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
        pattern_files = ["standard_2ome.json", "minimal_terminal.json", "maximal_stability.json"]

        for filename in pattern_files:
            filepath = examples_dir / filename
            with filepath.open() as f:
                data = json.load(f)

            # Check for expected top-level fields
            assert "pattern_name" in data, f"{filename} should have pattern_name"
            assert "description" in data, f"{filename} should have description"
            assert "guide_modifications" in data, f"{filename} should have guide_modifications"
            assert "passenger_modifications" in data, f"{filename} should have passenger_modifications"

            # Check guide modifications structure
            guide_mods = data["guide_modifications"]
            assert isinstance(guide_mods, dict), f"{filename} guide_modifications should be dict"

            # Should have at least one modification type
            assert len(guide_mods) > 0, f"{filename} should have at least one guide modification type"

    def test_example_readme_exists(self):
        """Test that README.md exists and is not empty."""
        examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
        readme = examples_dir / "README.md"

        assert readme.exists(), "README.md should exist"
        content = readme.read_text()
        assert len(content) > 100, "README.md should have substantial content"
        assert "modification" in content.lower(), "README should mention modifications"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
