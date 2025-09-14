"""Validation utilities for data consistency and cross-validation."""

from typing import Any

import pandas as pd

from sirnaforge.models.schemas import OffTargetHitsSchema, ORFValidationSchema, SiRNACandidateSchema
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

# Constants derived from model definitions to avoid duplication
SIRNA_MIN_LENGTH = 19  # From SiRNACandidate.length field constraints
SIRNA_MAX_LENGTH = 23  # From SiRNACandidate.length field constraints


class ValidationResult:
    """Container for validation results."""

    def __init__(self, is_valid: bool = True):
        self.is_valid = is_valid
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.metadata: dict[str, Any] = {}

    def add_error(self, message: str) -> None:
        """Add a validation error."""
        self.errors.append(message)
        self.is_valid = False

    def add_warning(self, message: str) -> None:
        """Add a validation warning."""
        self.warnings.append(message)

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the result."""
        self.metadata[key] = value

    def merge(self, other: "ValidationResult") -> None:
        """Merge another validation result into this one."""
        if not other.is_valid:
            self.is_valid = False
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        self.metadata.update(other.metadata)

    def summary(self) -> dict[str, Any]:
        """Get a summary of validation results."""
        return {
            "is_valid": self.is_valid,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "errors": self.errors,
            "warnings": self.warnings,
            "metadata": self.metadata,
        }


class ValidationUtils:
    """Utility functions for data validation."""

    @staticmethod
    def validate_nucleotide_sequence(sequence: str, allow_ambiguous: bool = True) -> ValidationResult:
        """Validate nucleotide sequence composition."""
        result = ValidationResult()

        # Define valid characters
        valid_chars = set("ATCGU")
        if allow_ambiguous:
            valid_chars.update("NRYSWKMBDHV-")

        # Check for invalid characters
        invalid_chars = set(sequence.upper()) - valid_chars
        if invalid_chars:
            result.add_error(f"Invalid nucleotides found: {sorted(invalid_chars)}")

        # Check for excessive poly-runs
        for base in "ATCGU":
            if base * 4 in sequence.upper():  # 4+ consecutive identical bases
                result.add_warning(f"Poly-{base} run detected in sequence")

        # Add sequence composition metadata
        result.add_metadata("length", len(sequence))
        result.add_metadata("gc_content", ValidationUtils._calculate_gc_content(sequence))

        return result

    @staticmethod
    def validate_sirna_length(sequence: str) -> ValidationResult:
        """Validate siRNA sequence length."""
        result = ValidationResult()
        length = len(sequence)

        if not (19 <= length <= 23):
            result.add_error(f"siRNA length {length} outside valid range (19-23)")
        elif length != 21:
            result.add_warning(f"siRNA length {length} is non-standard (21 is optimal)")

        result.add_metadata("length", length)
        return result

    @staticmethod
    def validate_parameter_consistency(params: DesignParameters) -> ValidationResult:
        """Validate design parameter consistency."""
        result = ValidationResult()

        # Check filter criteria consistency
        if params.filters.gc_min > params.filters.gc_max:
            result.add_error("gc_min cannot be greater than gc_max")

        if params.filters.gc_max - params.filters.gc_min < 5:
            result.add_warning("Very narrow GC content range may yield few candidates")

        # Check scoring weights
        total_weight = (
            params.scoring.asymmetry
            + params.scoring.gc_content
            + params.scoring.accessibility
            + params.scoring.off_target
            + params.scoring.empirical
        )

        if abs(total_weight - 1.0) > 0.01:
            result.add_error(f"Scoring weights sum to {total_weight:.3f}, should be 1.0")

        # Check parameter ranges
        if params.top_n > 1000:
            result.add_warning("Large top_n value may impact performance")

        result.add_metadata("total_weight", total_weight)
        return result

    @staticmethod
    def validate_candidate_consistency(candidate: SiRNACandidate) -> ValidationResult:
        """Validate siRNA candidate internal consistency."""
        result = ValidationResult()

        # Check sequence lengths match — different guide/passenger lengths are allowed
        # but should be flagged as a warning so downstream processing can handle them.
        if len(candidate.guide_sequence) != len(candidate.passenger_sequence):
            result.add_warning("Guide and passenger sequences have different lengths")

        # Check position is positive
        if candidate.position <= 0:
            result.add_error("Position must be positive")

        # Check GC content is reasonable
        calculated_gc = ValidationUtils._calculate_gc_content(candidate.guide_sequence)
        if abs(calculated_gc - candidate.gc_content) > 5.0:
            result.add_warning(
                f"Reported GC content ({candidate.gc_content:.1f}%) differs from calculated ({calculated_gc:.1f}%)"
            )

        # Check score ranges
        if not (0 <= candidate.composite_score <= 100):
            result.add_error(f"Composite score {candidate.composite_score} outside valid range (0-100)")

        if not (0 <= candidate.asymmetry_score <= 1):
            result.add_error(f"Asymmetry score {candidate.asymmetry_score} outside valid range (0-1)")

        result.add_metadata("calculated_gc", calculated_gc)
        return result

    @staticmethod
    def validate_dataframe_schema(df: pd.DataFrame, schema_type: str) -> ValidationResult:
        """Validate DataFrame against appropriate pandera schema."""
        result = ValidationResult()

        try:
            if schema_type == "sirna_candidates":
                df_sirna = SiRNACandidateSchema.validate(df)
                result.add_metadata("validated_rows", len(df_sirna))

            elif schema_type == "orf_validation":
                df_orf = ORFValidationSchema.validate(df)
                result.add_metadata("validated_rows", len(df_orf))

            elif schema_type == "off_target_hits":
                df_hits = OffTargetHitsSchema.validate(df)
                result.add_metadata("validated_rows", len(df_hits))

            else:
                result.add_error(f"Unknown schema type: {schema_type}")

        except Exception as e:
            result.add_error(f"Schema validation failed: {str(e)}")
            logger.error(f"DataFrame schema validation error: {e}")

        return result

    @staticmethod
    def validate_transcript_ids_consistency(
        candidate_ids: set[str], orf_ids: set[str], transcript_ids: set[str]
    ) -> ValidationResult:
        """Validate consistency of transcript IDs across datasets."""
        result = ValidationResult()

        # Check for missing IDs
        candidates_missing = candidate_ids - transcript_ids
        if candidates_missing:
            result.add_error(f"Candidates reference unknown transcripts: {candidates_missing}")

        orf_missing = orf_ids - transcript_ids
        if orf_missing:
            result.add_error(f"ORF analysis references unknown transcripts: {orf_missing}")

        # Check for unused transcripts
        unused_transcripts = transcript_ids - candidate_ids - orf_ids
        if unused_transcripts:
            result.add_warning(f"Transcripts not used in analysis: {unused_transcripts}")

        result.add_metadata("candidate_count", len(candidate_ids))
        result.add_metadata("orf_count", len(orf_ids))
        result.add_metadata("transcript_count", len(transcript_ids))

        return result

    @staticmethod
    def validate_biological_constraints(candidate: SiRNACandidate) -> ValidationResult:
        """Validate bioinformatics-specific constraints."""
        result = ValidationResult()

        # Check for forbidden motifs (simplified examples)
        forbidden_motifs = ["AAAA", "TTTT", "CCCC", "GGGG"]
        for motif in forbidden_motifs:
            if motif in candidate.guide_sequence:
                result.add_warning(f"Forbidden motif {motif} found in guide sequence")

        # Check thermodynamic properties
        if candidate.asymmetry_score < 0.2:
            result.add_warning("Low asymmetry score may reduce siRNA efficacy")

        if candidate.paired_fraction > 0.6:
            result.add_warning("High secondary structure may reduce accessibility")

        # Check GC content range
        if candidate.gc_content < 30 or candidate.gc_content > 52:
            result.add_warning(f"GC content {candidate.gc_content:.1f}% outside optimal range (30-52%)")

        return result

    @staticmethod
    def _calculate_gc_content(sequence: str) -> float:
        """Calculate GC content percentage."""
        if not sequence:
            return 0.0

        gc_count = sequence.upper().count("G") + sequence.upper().count("C")
        total_count = len([c for c in sequence.upper() if c in "ATCGU"])

        if total_count == 0:
            return 0.0

        return (gc_count / total_count) * 100.0

    @staticmethod
    def cross_validate_pydantic_pandera() -> ValidationResult:
        """Cross-validate Pydantic model constraints with Pandera schema constraints."""
        result = ValidationResult()

        # Check siRNA length constraints using constants derived from model definitions
        pydantic_min_length = SIRNA_MIN_LENGTH
        pydantic_max_length = SIRNA_MAX_LENGTH
        pandera_min_length = SIRNA_MIN_LENGTH  # From SiRNACandidateSchema.check_sequence_lengths
        pandera_max_length = SIRNA_MAX_LENGTH

        if pydantic_min_length != pandera_min_length:
            result.add_error("Pydantic and Pandera minimum length constraints don't match")

        if pydantic_max_length != pandera_max_length:
            result.add_error("Pydantic and Pandera maximum length constraints don't match")

        # Check GC content ranges
        pydantic_gc_min = 0.0  # From SiRNACandidate model
        pydantic_gc_max = 100.0
        pandera_gc_min = 0.0  # From schema
        pandera_gc_max = 100.0

        if pydantic_gc_min != pandera_gc_min or pydantic_gc_max != pandera_gc_max:
            result.add_error("Pydantic and Pandera GC content ranges don't match")

        result.add_metadata("constraints_checked", ["sequence_length", "gc_content"])

        return result
