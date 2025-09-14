"""Validation configuration and settings."""

from enum import Enum

from pydantic import BaseModel, Field


class ValidationLevel(str, Enum):
    """Validation strictness levels."""

    STRICT = "strict"  # Fail on any validation error
    WARNING = "warning"  # Log warnings but continue
    DISABLED = "disabled"  # Skip validation


class ValidationStage(str, Enum):
    """Pipeline stages where validation can be applied."""

    INPUT = "input"  # Input data validation
    TRANSCRIPT_RETRIEVAL = "transcript_retrieval"  # After transcript fetching
    ORF_ANALYSIS = "orf_analysis"  # After ORF analysis
    DESIGN = "design"  # After siRNA design
    FILTERING = "filtering"  # After filtering
    SCORING = "scoring"  # After scoring
    OFF_TARGET = "off_target"  # After off-target analysis
    OUTPUT = "output"  # Final output validation


class ValidationConfig(BaseModel):
    """Configuration for validation system."""

    # Global validation settings
    default_level: ValidationLevel = Field(
        default=ValidationLevel.STRICT, description="Default validation level for all stages"
    )

    # Stage-specific validation levels
    stage_levels: dict[ValidationStage, ValidationLevel] = Field(
        default_factory=dict, description="Per-stage validation level overrides"
    )

    # Validation options
    validate_sequences: bool = Field(default=True, description="Enable nucleotide sequence validation")

    validate_ranges: bool = Field(default=True, description="Enable numerical range validation")

    validate_consistency: bool = Field(default=True, description="Enable cross-field consistency validation")

    validate_biology: bool = Field(default=True, description="Enable bioinformatics-specific validation")

    # Error handling
    max_validation_errors: int = Field(
        default=100, ge=1, description="Maximum number of validation errors to collect before stopping"
    )

    collect_all_errors: bool = Field(default=False, description="Collect all validation errors instead of failing fast")

    # Performance settings
    batch_size: int = Field(default=1000, ge=1, description="Batch size for large dataset validation")

    enable_caching: bool = Field(default=True, description="Enable validation result caching")

    def get_level_for_stage(self, stage: ValidationStage) -> ValidationLevel:
        """Get validation level for a specific stage."""
        return self.stage_levels.get(stage, self.default_level)

    def is_enabled_for_stage(self, stage: ValidationStage) -> bool:
        """Check if validation is enabled for a stage."""
        return self.get_level_for_stage(stage) != ValidationLevel.DISABLED

    def should_fail_on_error(self, stage: ValidationStage) -> bool:
        """Check if validation errors should cause failures."""
        return self.get_level_for_stage(stage) == ValidationLevel.STRICT


# Default configurations for different use cases
class ValidationPresets:
    """Predefined validation configurations."""

    @staticmethod
    def development() -> ValidationConfig:
        """Configuration for development environment."""
        return ValidationConfig(
            default_level=ValidationLevel.WARNING, collect_all_errors=True, max_validation_errors=50
        )

    @staticmethod
    def production() -> ValidationConfig:
        """Configuration for production environment."""
        return ValidationConfig(
            default_level=ValidationLevel.STRICT, collect_all_errors=False, max_validation_errors=10
        )

    @staticmethod
    def testing() -> ValidationConfig:
        """Configuration for testing environment."""
        return ValidationConfig(default_level=ValidationLevel.STRICT, collect_all_errors=True, enable_caching=False)

    @staticmethod
    def performance() -> ValidationConfig:
        """Configuration optimized for performance."""
        return ValidationConfig(
            default_level=ValidationLevel.WARNING,
            stage_levels={
                ValidationStage.INPUT: ValidationLevel.STRICT,
                ValidationStage.OUTPUT: ValidationLevel.STRICT,
            },
            validate_biology=False,
            batch_size=5000,
        )
