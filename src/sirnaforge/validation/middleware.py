"""Validation middleware for integrating validation into the siRNA design workflow."""

import json
import time
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from sirnaforge.data.base import TranscriptInfo
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.validation.config import ValidationConfig, ValidationLevel, ValidationStage
from sirnaforge.validation.utils import ValidationResult, ValidationUtils

logger = get_logger(__name__)


class ValidationReport:
    """Comprehensive validation report for a workflow stage."""

    def __init__(self, stage: ValidationStage):
        self.stage = stage
        self.start_time = time.time()
        self.end_time: Optional[float] = None
        self.overall_result = ValidationResult()
        self.item_results: list[ValidationResult] = []
        self.summary_stats: dict[str, Any] = {}

    def add_item_result(self, result: ValidationResult) -> None:
        """Add validation result for an individual item."""
        self.item_results.append(result)
        if not result.is_valid:
            self.overall_result.is_valid = False
        self.overall_result.errors.extend(result.errors)
        self.overall_result.warnings.extend(result.warnings)

    def finalize(self) -> None:
        """Finalize the report and calculate summary statistics."""
        self.end_time = time.time()

        total_items = len(self.item_results)
        valid_items = sum(1 for r in self.item_results if r.is_valid)
        total_errors = sum(len(r.errors) for r in self.item_results)
        total_warnings = sum(len(r.warnings) for r in self.item_results)

        self.summary_stats = {
            "stage": self.stage.value,
            "duration_seconds": self.end_time - self.start_time,
            "total_items": total_items,
            "valid_items": valid_items,
            "invalid_items": total_items - valid_items,
            "total_errors": total_errors,
            "total_warnings": total_warnings,
            "success_rate": (valid_items / total_items) * 100 if total_items > 0 else 100,
        }

        self.overall_result.add_metadata("summary", self.summary_stats)

    def to_dict(self) -> dict[str, Any]:
        """Convert report to dictionary."""
        return {
            "stage": self.stage.value,
            "overall_valid": self.overall_result.is_valid,
            "summary": self.summary_stats,
            "errors": self.overall_result.errors,
            "warnings": self.overall_result.warnings,
        }


class ValidationMiddleware:
    """Middleware for integrating validation throughout the workflow."""

    def __init__(self, config: ValidationConfig):
        self.config = config
        self.reports: list[ValidationReport] = []
        self._validation_cache: dict[str, ValidationResult] = {}

    def validate_input_parameters(self, params: DesignParameters) -> ValidationReport:
        """Validate input design parameters."""
        stage = ValidationStage.INPUT
        report = ValidationReport(stage)

        if not self.config.is_enabled_for_stage(stage):
            logger.debug(f"Validation disabled for stage: {stage.value}")
            report.finalize()
            return report

        logger.info(f"Validating input parameters for stage: {stage.value}")

        # Validate parameter consistency
        param_result = ValidationUtils.validate_parameter_consistency(params)
        report.add_item_result(param_result)

        # Cross-validate Pydantic and Pandera constraints
        if self.config.validate_consistency:
            cross_val_result = ValidationUtils.cross_validate_pydantic_pandera()
            report.add_item_result(cross_val_result)

        report.finalize()
        self._handle_validation_result(stage, report)
        self.reports.append(report)

        return report

    def validate_transcripts(self, transcripts: list[TranscriptInfo]) -> ValidationReport:
        """Validate transcript data after retrieval."""
        stage = ValidationStage.TRANSCRIPT_RETRIEVAL
        report = ValidationReport(stage)

        if not self.config.is_enabled_for_stage(stage):
            report.finalize()
            return report

        logger.info(f"Validating {len(transcripts)} transcripts")

        for transcript in transcripts:
            result = ValidationResult()

            # Validate transcript has sequence
            if not transcript.sequence:
                result.add_error(f"Transcript {transcript.transcript_id} has no sequence")
            else:
                # Validate sequence composition
                if self.config.validate_sequences:
                    seq_result = ValidationUtils.validate_nucleotide_sequence(transcript.sequence, allow_ambiguous=True)
                    result.merge(seq_result)

                # Validate sequence length
                if len(transcript.sequence) < 100:
                    result.add_warning(
                        f"Transcript {transcript.transcript_id} is very short ({len(transcript.sequence)} nt)"
                    )

            # Validate required fields
            if not transcript.transcript_id:
                result.add_error("Transcript missing ID")

            report.add_item_result(result)

        report.finalize()
        self._handle_validation_result(stage, report)
        self.reports.append(report)

        return report

    def validate_design_results(self, design_result: DesignResult) -> ValidationReport:
        """Validate siRNA design results."""
        stage = ValidationStage.DESIGN
        report = ValidationReport(stage)

        if not self.config.is_enabled_for_stage(stage):
            report.finalize()
            return report

        logger.info(f"Validating design results with {len(design_result.candidates)} candidates")

        for candidate in design_result.candidates:
            result = ValidationResult()

            # Validate candidate consistency
            consistency_result = ValidationUtils.validate_candidate_consistency(candidate)
            result.merge(consistency_result)

            # Validate sequences
            if self.config.validate_sequences:
                guide_result = ValidationUtils.validate_nucleotide_sequence(
                    candidate.guide_sequence, allow_ambiguous=False
                )
                passenger_result = ValidationUtils.validate_nucleotide_sequence(
                    candidate.passenger_sequence, allow_ambiguous=False
                )
                length_result = ValidationUtils.validate_sirna_length(candidate.guide_sequence)

                result.merge(guide_result)
                result.merge(passenger_result)
                result.merge(length_result)

            # Validate biological constraints
            if self.config.validate_biology:
                bio_result = ValidationUtils.validate_biological_constraints(candidate)
                result.merge(bio_result)

            report.add_item_result(result)

        report.finalize()
        self._handle_validation_result(stage, report)
        self.reports.append(report)

        return report

    def validate_dataframe_output(self, df: pd.DataFrame, schema_type: str) -> ValidationReport:
        """Validate DataFrame output against pandera schemas."""
        stage = ValidationStage.OUTPUT
        report = ValidationReport(stage)

        if not self.config.is_enabled_for_stage(stage):
            report.finalize()
            return report

        logger.info(f"Validating {schema_type} DataFrame with {len(df)} rows")

        # Validate against schema
        schema_result = ValidationUtils.validate_dataframe_schema(df, schema_type)
        report.add_item_result(schema_result)

        report.finalize()
        self._handle_validation_result(stage, report)
        self.reports.append(report)

        return report

    def validate_transcript_id_consistency(
        self,
        transcripts: list[TranscriptInfo],
        candidates: list[SiRNACandidate],
        orf_data: Optional[pd.DataFrame] = None,
    ) -> ValidationReport:
        """Validate consistency of transcript IDs across datasets."""
        stage = ValidationStage.FILTERING  # Use filtering stage for consistency checks
        report = ValidationReport(stage)

        if not self.config.is_enabled_for_stage(stage) or not self.config.validate_consistency:
            report.finalize()
            return report

        logger.info("Validating transcript ID consistency across datasets")

        # Extract ID sets
        transcript_ids = {t.transcript_id for t in transcripts}
        candidate_ids = {c.transcript_id for c in candidates}
        orf_ids = set(orf_data["transcript_id"]) if orf_data is not None else set()

        # Validate consistency
        consistency_result = ValidationUtils.validate_transcript_ids_consistency(candidate_ids, orf_ids, transcript_ids)
        report.add_item_result(consistency_result)

        report.finalize()
        self._handle_validation_result(stage, report)
        self.reports.append(report)

        return report

    def save_validation_report(self, output_path: Path) -> None:
        """Save comprehensive validation report."""
        report_data = {
            "validation_config": {
                "default_level": self.config.default_level.value,
                "validate_sequences": self.config.validate_sequences,
                "validate_biology": self.config.validate_biology,
                "validate_consistency": self.config.validate_consistency,
            },
            "stage_reports": [report.to_dict() for report in self.reports],
            "summary": self._generate_summary(),
        }

        with output_path.open("w") as f:
            json.dump(report_data, f, indent=2, default=str)

        logger.info(f"Validation report saved to: {output_path}")

    def _generate_summary(self) -> dict[str, Any]:
        """Generate overall validation summary."""
        total_errors = sum(len(report.overall_result.errors) for report in self.reports)
        total_warnings = sum(len(report.overall_result.warnings) for report in self.reports)
        stages_validated = len(self.reports)
        stages_passed = sum(1 for report in self.reports if report.overall_result.is_valid)

        return {
            "total_stages": stages_validated,
            "stages_passed": stages_passed,
            "stages_failed": stages_validated - stages_passed,
            "total_errors": total_errors,
            "total_warnings": total_warnings,
            "overall_success": stages_passed == stages_validated,
        }

    def _handle_validation_result(self, stage: ValidationStage, report: ValidationReport) -> None:
        """Handle validation results based on configuration."""
        level = self.config.get_level_for_stage(stage)

        if level == ValidationLevel.DISABLED:
            return

        # Log results
        if report.overall_result.errors:
            if level == ValidationLevel.STRICT:
                logger.error(f"Validation failed for stage {stage.value}: {len(report.overall_result.errors)} errors")
                for error in report.overall_result.errors[:5]:  # Log first 5 errors
                    logger.error(f"  - {error}")
            else:
                logger.warning(f"Validation issues in stage {stage.value}: {len(report.overall_result.errors)} errors")

        if report.overall_result.warnings:
            logger.warning(
                f"Validation warnings for stage {stage.value}: {len(report.overall_result.warnings)} warnings"
            )

        # Handle failures
        if not report.overall_result.is_valid and self.config.should_fail_on_error(stage):
            error_summary = f"Validation failed for stage {stage.value}: {len(report.overall_result.errors)} errors"
            raise ValueError(error_summary)
