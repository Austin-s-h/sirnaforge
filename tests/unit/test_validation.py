"""Tests for the enhanced validation system."""

import json
import tempfile
from pathlib import Path

import pandas as pd
import pydantic_core
import pytest

from sirnaforge.models.sirna import DesignParameters, DesignResult, FilterCriteria, SiRNACandidate
from sirnaforge.validation import ValidationConfig, ValidationMiddleware, ValidationUtils
from sirnaforge.validation.config import ValidationLevel, ValidationStage


class TestValidationUtils:
    """Test validation utility functions."""

    def test_validate_nucleotide_sequence_valid(self):
        """Test nucleotide sequence validation with valid sequences."""
        result = ValidationUtils.validate_nucleotide_sequence("ATCGATCGATCG")
        assert result.is_valid
        assert len(result.errors) == 0
        assert result.metadata["length"] == 12
        assert "gc_content" in result.metadata

    def test_validate_nucleotide_sequence_invalid_chars(self):
        """Test nucleotide sequence validation with invalid characters."""
        result = ValidationUtils.validate_nucleotide_sequence("ATCGXYZ")
        assert not result.is_valid
        assert "Invalid nucleotides found" in result.errors[0]

    def test_validate_nucleotide_sequence_poly_runs(self):
        """Test nucleotide sequence validation detects poly-runs."""
        result = ValidationUtils.validate_nucleotide_sequence("ATCGAAAAGCTA")
        assert result.is_valid  # Should be valid but with warnings
        assert len(result.warnings) > 0
        assert "Poly-A run detected" in result.warnings[0]

    def test_validate_sirna_length_valid(self):
        """Test siRNA length validation with valid length."""
        result = ValidationUtils.validate_sirna_length("ATCGATCGATCGATCGATCGA")  # 21 nt
        assert result.is_valid
        assert result.metadata["length"] == 21

    def test_validate_sirna_length_invalid(self):
        """Test siRNA length validation with invalid length."""
        result = ValidationUtils.validate_sirna_length("ATCG")  # Too short
        assert not result.is_valid
        assert "outside valid range" in result.errors[0]

    def test_validate_parameter_consistency_valid(self):
        """Test parameter consistency validation with valid parameters."""
        params = DesignParameters()
        result = ValidationUtils.validate_parameter_consistency(params)
        assert result.is_valid

    def test_validate_parameter_consistency_invalid_gc_range(self):
        """Test parameter consistency validation with invalid GC range."""
        filters = FilterCriteria(gc_min=60.0, gc_max=40.0)  # Invalid: min > max
        params = DesignParameters(filters=filters)
        result = ValidationUtils.validate_parameter_consistency(params)
        assert not result.is_valid
        assert "gc_min cannot be greater than gc_max" in result.errors[0]

    def test_validate_candidate_consistency_valid(self):
        """Test candidate consistency validation with valid candidate."""
        candidate = SiRNACandidate(
            id="test_1",
            transcript_id="ENST123",
            position=100,
            guide_sequence="ATCGATCGATCGATCGATCGA",
            passenger_sequence="TCGATCGATCGATCGATCGAT",
            gc_content=50.0,
            length=21,
            asymmetry_score=0.5,
            composite_score=75.0,
        )
        result = ValidationUtils.validate_candidate_consistency(candidate)
        assert result.is_valid

    def test_validate_candidate_consistency_mismatched_lengths(self):
        """Test candidate consistency validation with mismatched sequence lengths."""
        with pytest.raises(pydantic_core._pydantic_core.ValidationError, match="at most 23 characters"):
            _ = SiRNACandidate(
                id="test_1",
                transcript_id="ENST123",
                position=100,
                guide_sequence="ATCGATCGATCGATCGATCGA",  # 21 nt
                passenger_sequence="TCGATCGATCGATCGATCGAT123",  # 24 nt
                gc_content=50.0,
                length=21,
                asymmetry_score=0.5,
                composite_score=75.0,
            )


class TestValidationConfig:
    """Test validation configuration."""

    def test_default_config(self):
        """Test default validation configuration."""
        config = ValidationConfig()
        assert config.default_level == ValidationLevel.STRICT
        assert config.validate_sequences is True
        assert config.validate_biology is True

    def test_stage_level_override(self):
        """Test stage-specific validation level overrides."""
        config = ValidationConfig(
            default_level=ValidationLevel.WARNING, stage_levels={ValidationStage.INPUT: ValidationLevel.STRICT}
        )
        assert config.get_level_for_stage(ValidationStage.INPUT) == ValidationLevel.STRICT
        assert config.get_level_for_stage(ValidationStage.DESIGN) == ValidationLevel.WARNING

    def test_validation_enabled_check(self):
        """Test validation enabled check for stages."""
        config = ValidationConfig(stage_levels={ValidationStage.OUTPUT: ValidationLevel.DISABLED})
        assert config.is_enabled_for_stage(ValidationStage.INPUT) is True
        assert config.is_enabled_for_stage(ValidationStage.OUTPUT) is False


class TestValidationMiddleware:
    """Test validation middleware integration."""

    def test_validate_input_parameters(self):
        """Test input parameter validation."""
        config = ValidationConfig()
        middleware = ValidationMiddleware(config)

        params = DesignParameters()
        report = middleware.validate_input_parameters(params)

        assert report.stage == ValidationStage.INPUT
        assert report.overall_result.is_valid

    def test_validate_design_results_with_valid_data(self):
        """Test design results validation with valid data."""
        config = ValidationConfig()
        middleware = ValidationMiddleware(config)

        candidate = SiRNACandidate(
            id="test_1",
            transcript_id="ENST123",
            position=100,
            guide_sequence="ATCGATCGATCGATCGATCGA",
            passenger_sequence="TCGATCGATCGATCGATCGAT",
            gc_content=50.0,
            length=21,
            asymmetry_score=0.5,
            composite_score=75.0,
        )

        design_result = DesignResult(
            input_file="test.fasta",
            parameters=DesignParameters(),
            candidates=[candidate],
            top_candidates=[candidate],
            total_sequences=1,
            total_candidates=1,
            filtered_candidates=1,
            processing_time=1.0,
        )

        report = middleware.validate_design_results(design_result)
        assert report.stage == ValidationStage.DESIGN
        assert report.overall_result.is_valid

    def test_validate_dataframe_schema(self):
        """Test DataFrame schema validation."""
        config = ValidationConfig()
        middleware = ValidationMiddleware(config)

        # Create valid siRNA candidates DataFrame
        df = pd.DataFrame(
            {
                "id": ["test_1"],
                "transcript_id": ["ENST123"],
                "position": [100],
                "guide_sequence": ["ATCGATCGATCGATCGATCGA"],
                "passenger_sequence": ["TCGATCGATCGATCGATCGAT"],
                "gc_content": [50.0],
                "asymmetry_score": [0.5],
                "paired_fraction": [0.3],
                "off_target_count": [0],
                "transcript_hit_count": [1],
                "transcript_hit_fraction": [1.0],
                "composite_score": [75.0],
                "passes_filters": [True],
            }
        )

        report = middleware.validate_dataframe_output(df, "sirna_candidates")
        assert report.stage == ValidationStage.OUTPUT
        assert report.overall_result.is_valid

    def test_validation_report_generation(self):
        """Test validation report generation."""
        config = ValidationConfig()
        middleware = ValidationMiddleware(config)

        # Run some validations
        params = DesignParameters()
        middleware.validate_input_parameters(params)

        # Save report
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "validation_report.json"
            middleware.save_validation_report(report_path)

            assert report_path.exists()

            # Check report content
            with report_path.open() as f:
                report_data = json.load(f)

            assert "validation_config" in report_data
            assert "stage_reports" in report_data
            assert "summary" in report_data
