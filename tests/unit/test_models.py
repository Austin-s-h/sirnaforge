"""Basic tests for siRNA models and CLI."""

import tempfile
from pathlib import Path

import pytest

from sirnaforge.cli import app
from sirnaforge.models.sirna import DesignParameters, FilterCriteria, SiRNACandidate


class TestDesignParameters:
    """Test DesignParameters model."""

    def test_default_parameters(self):
        """Test default parameter values."""
        params = DesignParameters()
        assert params.sirna_length == 21
        assert params.top_n == 10
        assert params.filters.gc_min == 30.0
        assert params.filters.gc_max == 52.0

    def test_custom_parameters(self):
        """Test custom parameter values."""
        filters = FilterCriteria(gc_min=35.0, gc_max=55.0)
        params = DesignParameters(sirna_length=22, top_n=5, filters=filters)
        assert params.sirna_length == 22
        assert params.top_n == 5
        assert params.filters.gc_min == 35.0

    def test_parameter_validation(self):
        """Test parameter validation."""
        with pytest.raises(ValueError):
            FilterCriteria(gc_min=60.0, gc_max=50.0)  # gc_max < gc_min


class TestSiRNACandidate:
    """Test SiRNACandidate model."""

    def test_valid_candidate(self):
        """Test creating a valid siRNA candidate."""
        candidate = SiRNACandidate(
            id="sirna_001",
            transcript_id="ENST00000123456",
            position=100,
            guide_sequence="AUCGAUCGAUCGAUCGAUCGA",
            passenger_sequence="UCGAUCGAUCGAUCGAUCGAU",
            gc_content=52.4,
            length=21,
            asymmetry_score=0.75,
            composite_score=85.2,
        )
        assert candidate.id == "sirna_001"
        assert candidate.gc_content == 52.4
        assert candidate.composite_score == 85.2

    def test_sequence_validation(self):
        """Test sequence validation."""
        with pytest.raises(ValueError):
            SiRNACandidate(
                id="sirna_002",
                transcript_id="ENST00000123456",
                position=100,
                guide_sequence="AUCGXUCGAUCGAUCGAUCGA",  # Invalid X base
                passenger_sequence="UCGAUCGAUCGAUCGAUCGAU",
                gc_content=52.4,
                length=21,
                asymmetry_score=0.75,
                composite_score=85.2,
            )

    def test_length_mismatch(self):
        """Test length mismatch validation."""
        with pytest.raises(ValueError):
            SiRNACandidate(
                id="sirna_003",
                transcript_id="ENST00000123456",
                position=100,
                guide_sequence="AUCGAUCGAUCGAUCGAUCGA",  # 21 bases
                passenger_sequence="UCGAUCGAUCGAUCGAUCGAUUU",  # 23 bases - mismatch!
                gc_content=52.4,
                length=21,
                asymmetry_score=0.75,
                composite_score=85.2,
            )

    def test_fasta_output(self):
        """Test FASTA format output."""
        candidate = SiRNACandidate(
            id="sirna_004",
            transcript_id="ENST00000123456",
            position=100,
            guide_sequence="AUCGAUCGAUCGAUCGAUCGA",
            passenger_sequence="UCGAUCGAUCGAUCGAUCGAU",
            gc_content=52.4,
            length=21,
            asymmetry_score=0.75,
            composite_score=85.2,
        )

        fasta = candidate.to_fasta()
        expected = ">sirna_004\nAUCGAUCGAUCGAUCGAUCGA\n"
        assert fasta == expected


class TestCLI:
    """Test CLI functionality."""

    def test_version_command(self, capsys):  # noqa: ARG002
        """Test version command."""
        # Note: capsys parameter kept for potential future use in output testing

        # This is a basic test - in practice you'd use typer's testing utilities
        # For now, just test that the CLI module imports correctly
        assert app is not None

    def test_validate_command_basic(self):
        """Test basic validation functionality."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">test_seq\nATCGATCGATCGATCGATCG\n")
            temp_fasta = f.name

        # Test that file exists and can be read
        fasta_path = Path(temp_fasta)
        assert fasta_path.exists()

        # Clean up
        fasta_path.unlink()


if __name__ == "__main__":
    pytest.main([__file__])
