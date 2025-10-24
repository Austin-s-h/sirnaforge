"""Pydantic models for siRNA design data structures."""

from enum import Enum
from typing import Any, Callable, Optional, TypeVar, Union

import pandas as pd
import pandera as pa
from pandera.typing import DataFrame
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator

from sirnaforge.models.modifications import StrandMetadata, StrandRole
from sirnaforge.models.schemas import SiRNACandidateSchema
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import get_modification_summary

logger = get_logger(__name__)

# mypy-friendly typed alias for pydantic's untyped decorator factory
F = TypeVar("F", bound=Callable[..., Any])
FieldValidatorFactory = Callable[..., Callable[[F], F]]
field_validator_typed: FieldValidatorFactory = field_validator


class FilterCriteria(BaseModel):
    """Quality filters for siRNA candidate selection based on thermodynamic and empirical criteria."""

    # GC content filters (updated to match documentation: optimal 35-60%)
    gc_min: float = Field(
        default=35.0, ge=0, le=100, description="Minimum GC content % (balance stability/accessibility)"
    )
    gc_max: float = Field(default=60.0, ge=0, le=100, description="Maximum GC content % (prevent over-stabilization)")

    # Sequence composition filters
    max_poly_runs: int = Field(
        default=3, ge=1, description="Max consecutive identical nucleotides (avoid synthesis issues)"
    )

    # Secondary structure filters
    max_paired_fraction: float = Field(
        default=0.6, ge=0, le=1, description="Max secondary structure pairing (prevent rigid structures)"
    )

    # Thermodynamic asymmetry filters
    min_asymmetry_score: float = Field(
        default=0.65,
        ge=0.3,
        le=1,
        description=(
            "Minimum thermodynamic asymmetry score for guide strand selection into RISC. "
            "Higher values (0.65-0.85) promote correct 5' end instability for effective strand loading."
        ),
    )

    # Minimum Free Energy filters (optimal: -2 to -8 kcal/mol)
    mfe_min: Optional[float] = Field(
        default=-8.0, description="Minimum MFE threshold in kcal/mol (more negative = too stable)"
    )
    mfe_max: Optional[float] = Field(
        default=-2.0, description="Maximum MFE threshold in kcal/mol (less negative = too unstable)"
    )

    # Duplex stability filters (optimal: -15 to -25 kcal/mol)
    duplex_stability_min: Optional[float] = Field(
        default=-25.0, description="Minimum duplex ΔG threshold in kcal/mol (more negative = too stable)"
    )
    duplex_stability_max: Optional[float] = Field(
        default=-15.0, description="Maximum duplex ΔG threshold in kcal/mol (less negative = too unstable)"
    )

    # Melting temperature filters (optimal: 60-78°C for human cells)
    melting_temp_min: Optional[float] = Field(default=60.0, description="Minimum melting temperature in °C")
    melting_temp_max: Optional[float] = Field(default=78.0, description="Maximum melting temperature in °C")

    # End asymmetry filters (optimal: +2 to +6 kcal/mol)
    delta_dg_end_min: Optional[float] = Field(
        default=2.0, description="Minimum end asymmetry ΔΔG (dg_3p - dg_5p) in kcal/mol"
    )
    delta_dg_end_max: Optional[float] = Field(
        default=6.0, description="Maximum end asymmetry ΔΔG (dg_3p - dg_5p) in kcal/mol"
    )

    # Off-target filters
    max_off_target_count: Optional[int] = Field(
        default=3, ge=0, description="Maximum allowed off-target sites (goal: ≤3)"
    )

    @field_validator_typed("gc_max")
    @classmethod
    def gc_max_greater_than_min(cls, v: float, info: ValidationInfo) -> float:
        if "gc_min" in info.data and v < info.data["gc_min"]:
            raise ValueError("gc_max must be greater than or equal to gc_min")
        return v


class ScoringWeights(BaseModel):
    """Relative weights for composite siRNA scoring components."""

    asymmetry: float = Field(
        default=0.25, ge=0, le=1, description="Thermodynamic asymmetry weight (guide strand selection)"
    )
    gc_content: float = Field(
        default=0.20, ge=0, le=1, description="GC content optimization weight (stability balance)"
    )
    accessibility: float = Field(
        default=0.25, ge=0, le=1, description="Target accessibility weight (secondary structure)"
    )
    off_target: float = Field(default=0.20, ge=0, le=1, description="Off-target avoidance weight (specificity)")
    empirical: float = Field(
        default=0.10, ge=0, le=1, description="Empirical design rules weight (established patterns)"
    )

    @field_validator_typed("empirical")
    @classmethod
    def weights_sum_to_one(cls, v: float, info: ValidationInfo) -> float:
        total = sum(info.data.values()) + v
        if not (0.95 <= total <= 1.05):  # Allow small floating point errors
            raise ValueError(f"Scoring weights must sum to 1.0, got {total}")
        return v


class DesignParameters(BaseModel):
    """Complete configuration parameters for siRNA design workflow."""

    model_config = ConfigDict(extra="forbid")

    # Basic parameters
    sirna_length: int = Field(default=21, ge=19, le=23, description="siRNA duplex length in nucleotides")
    top_n: int = Field(default=50, ge=1, le=1000, description="Number of top-scoring candidates to return")

    # Filtering criteria
    filters: FilterCriteria = Field(default_factory=FilterCriteria, description="Quality control filters")

    # Scoring weights
    scoring: ScoringWeights = Field(default_factory=ScoringWeights, description="Component score weights")

    # Optional analysis parameters
    avoid_snps: bool = Field(default=True, description="Exclude regions with known SNPs")
    check_off_targets: bool = Field(default=True, description="Perform genome-wide off-target analysis")
    predict_structure: bool = Field(default=True, description="Calculate RNA secondary structures")

    # Chemical modification parameters
    apply_modifications: bool = Field(
        default=True, description="Automatically apply chemical modification patterns to designed siRNAs"
    )
    modification_pattern: str = Field(
        default="standard_2ome",
        description="Modification pattern to apply (standard_2ome, minimal_terminal, maximal_stability, none)",
    )
    default_overhang: str = Field(default="dTdT", description="Default overhang sequence (dTdT for DNA, UU for RNA)")

    # File paths (optional)
    # TODO: review snp incorporation feature
    snp_file: Optional[str] = Field(default=None, description="Path to SNP VCF file for avoidance")
    # Review genome index passing / FASTA selection
    genome_index: Optional[str] = Field(default=None, description="Path to BWA genome index for off-target search")


class SequenceType(str, Enum):
    """Categories of input sequence types for siRNA design."""

    TRANSCRIPT = "transcript"  # Full transcript sequence (mRNA)
    GENOMIC = "genomic"  # Genomic DNA sequence
    CDS = "cds"  # Protein-coding sequence only
    UTR = "utr"  # Untranslated region sequence


class SiRNACandidate(BaseModel):
    """Individual siRNA candidate with computed thermodynamic and efficacy properties."""

    model_config = ConfigDict(extra="forbid")

    # Identity
    id: str = Field(description="Unique siRNA candidate identifier")
    transcript_id: str = Field(description="Source transcript ID (e.g., ENST00000123456)")
    position: int = Field(ge=1, description="1-based start position in transcript")

    # Sequences
    guide_sequence: str = Field(min_length=19, max_length=23, description="Guide strand (antisense, loaded into RISC)")
    passenger_sequence: str = Field(
        min_length=19, max_length=23, description="Passenger strand (sense, typically degraded)"
    )

    # Basic properties
    gc_content: float = Field(ge=0, le=100, description="GC content % (optimal: 35-60%)")
    length: int = Field(ge=19, le=23, description="siRNA duplex length in nucleotides")

    # Thermodynamic properties
    asymmetry_score: float = Field(
        ge=0, le=1, description="Thermodynamic asymmetry score for RISC loading (optimal: ≥0.65)"
    )
    duplex_stability: Optional[float] = Field(default=None, description="Duplex formation ΔG in kcal/mol")

    # Secondary structure
    structure: Optional[str] = Field(default=None, description="RNA secondary structure (dot-bracket notation)")
    mfe: Optional[float] = Field(default=None, description="Minimum free energy in kcal/mol (optimal: -2 to -8)")
    paired_fraction: float = Field(default=0.0, ge=0, le=1, description="Fraction of paired bases (optimal: 0.4-0.6)")

    # Off-target analysis
    off_target_count: int = Field(default=0, ge=0, description="Number of potential off-target sites (goal: ≤3)")
    off_target_penalty: float = Field(default=0.0, ge=0, description="Off-target penalty score (lower is better)")

    # Transcript hit metrics (how many input transcripts this guide hits)
    transcript_hit_count: int = Field(
        default=1, ge=0, description="Number of input transcripts containing this guide sequence"
    )
    transcript_hit_fraction: float = Field(
        default=1.0, ge=0, le=1, description="Fraction of input transcripts targeted by this guide (1.0 = all)"
    )

    # Composite scoring
    component_scores: dict[str, float] = Field(default_factory=dict, description="Individual scoring component values")
    composite_score: float = Field(ge=0, le=100, description="Overall siRNA quality score (higher is better)")

    # Quality flags
    class FilterStatus(str, Enum):
        """Filter status codes for quality control."""

        # PASS is a domain status label, NOT a password. Bandit B105 false positive. # nosec B105
        PASS = "PASS"
        GC_OUT_OF_RANGE = "GC_OUT_OF_RANGE"
        POLY_RUNS = "POLY_RUNS"
        EXCESS_PAIRING = "EXCESS_PAIRING"
        LOW_ASYMMETRY = "LOW_ASYMMETRY"

    # Either True (passed) or one of the FilterStatus reasons (failed)
    passes_filters: Union[bool, FilterStatus] = Field(
        default=True, description="PASS if all filters passed, otherwise specific failure reason"
    )
    quality_issues: list[str] = Field(default_factory=list, description="List of detected quality concerns")

    # Optional chemical modification metadata
    guide_metadata: Optional[StrandMetadata] = Field(
        default=None,
        description="Optional StrandMetadata for guide strand with chemical modifications",
    )
    passenger_metadata: Optional[StrandMetadata] = Field(
        default=None,
        description="Optional StrandMetadata for passenger strand with chemical modifications",
    )

    @field_validator_typed("guide_sequence", "passenger_sequence")
    @classmethod
    def validate_nucleotide_sequence(cls, v: str) -> str:
        valid_bases = set("ATCGU")
        if not all(base.upper() in valid_bases for base in v):
            raise ValueError(f"Sequence contains invalid nucleotides: {v}")
        return v.upper()

    @field_validator_typed("passenger_sequence")
    @classmethod
    def sequences_same_length(cls, v: str, info: ValidationInfo) -> str:
        if "guide_sequence" in info.data and len(v) != len(info.data["guide_sequence"]):
            raise ValueError("Guide and passenger sequences must be the same length")
        return v

    def to_fasta(self, include_metadata: bool = False) -> str:
        """Return FASTA format representation of the guide sequence.

        Args:
            include_metadata: If True and guide_metadata is present, include it in the header

        Returns:
            FASTA-formatted string with candidate ID as header and guide sequence.
        """
        if include_metadata and self.guide_metadata:
            header = self.guide_metadata.to_fasta_header(target_gene=self.transcript_id, strand_role=StrandRole.GUIDE)
            # Extract just the header content after '>'
            header_content = header[1:] if header.startswith(">") else header
            return f">{header_content}\n{self.guide_sequence}\n"
        return f">{self.id}\n{self.guide_sequence}\n"


class DesignResult(BaseModel):
    """Complete results from siRNA design workflow with metadata and statistics."""

    model_config = ConfigDict(extra="forbid")

    # Input information
    input_file: str = Field(description="Path to input FASTA file processed")
    parameters: DesignParameters = Field(description="Design parameters used for this run")

    # Results
    candidates: list[SiRNACandidate] = Field(description="All generated siRNA candidates")
    top_candidates: list[SiRNACandidate] = Field(description="Top-scoring candidates (filtered and ranked)")

    # Summary statistics
    total_sequences: int = Field(ge=0, description="Number of input sequences processed")
    total_candidates: int = Field(ge=0, description="Total siRNA candidates generated")
    filtered_candidates: int = Field(ge=0, description="Candidates passing quality filters")

    # Processing metadata
    processing_time: float = Field(ge=0, description="Total processing time in seconds")
    tool_versions: dict[str, str] = Field(default_factory=dict, description="Software versions used in analysis")

    @pa.check_types
    def save_csv(self, filepath: str) -> DataFrame[SiRNACandidateSchema]:
        """Save siRNA candidates to CSV file with comprehensive validation.

        Exports all candidates to CSV format with full thermodynamic metrics.
        The DataFrame is validated against SiRNACandidateSchema before saving
        to ensure data integrity and proper column types.

        Args:
            filepath: Output CSV file path

        Returns:
            Validated DataFrame conforming to SiRNACandidateSchema

        Raises:
            pandera.errors.SchemaError: If data validation fails
        """
        df_data = []
        for candidate in self.candidates:
            cs = candidate.component_scores or {}

            # Get modification summary if modifications were applied
            mod_summary = get_modification_summary(candidate) if candidate.guide_metadata else {}

            row = {
                "id": candidate.id,
                "transcript_id": candidate.transcript_id,
                "position": candidate.position,
                "guide_sequence": candidate.guide_sequence,
                "passenger_sequence": candidate.passenger_sequence,
                "gc_content": candidate.gc_content,
                "asymmetry_score": candidate.asymmetry_score,
                # Thermodynamics and structure
                "structure": getattr(candidate, "structure", None),
                "mfe": getattr(candidate, "mfe", None),
                "paired_fraction": candidate.paired_fraction,
                "duplex_stability_dg": candidate.duplex_stability,
                "duplex_stability_score": cs.get("duplex_stability_score"),
                "dg_5p": cs.get("dg_5p"),
                "dg_3p": cs.get("dg_3p"),
                "delta_dg_end": cs.get("delta_dg_end"),
                "melting_temp_c": cs.get("melting_temp_c"),
                "off_target_count": candidate.off_target_count,
                "transcript_hit_count": candidate.transcript_hit_count,
                "transcript_hit_fraction": candidate.transcript_hit_fraction,
                "composite_score": candidate.composite_score,
                "passes_filters": (
                    candidate.passes_filters.value
                    if hasattr(candidate.passes_filters, "value")
                    else candidate.passes_filters
                ),
                # Chemical modifications
                "guide_overhang": mod_summary.get("guide_overhang", ""),
                "guide_modifications": mod_summary.get("guide_modifications", ""),
                "passenger_overhang": mod_summary.get("passenger_overhang", ""),
                "passenger_modifications": mod_summary.get("passenger_modifications", ""),
            }
            df_data.append(row)

        df = pd.DataFrame(df_data)

        # Validate DataFrame against schema - let failures bubble up
        logger.debug(f"Validating siRNA candidates DataFrame with {len(df)} rows")
        validated_df = SiRNACandidateSchema.validate(df)
        logger.info(f"siRNA candidates schema validation passed for {len(validated_df)} candidates")

        # Note: do not append design parameters as per-row columns to the candidates CSV.
        # Full design parameters are available in workflow metadata (`workflow_summary.json`).

        # Save validated DataFrame (with appended params if available)
        validated_df.to_csv(filepath, index=False)

        return validated_df

    def get_summary(self) -> dict[str, Any]:
        """Generate summary statistics for the design results.

        Returns:
            Dictionary containing key metrics including sequence counts,
            processing time, best score, and tool versions used.
        """
        return {
            "input_sequences": self.total_sequences,
            "total_candidates": self.total_candidates,
            "filtered_candidates": self.filtered_candidates,
            "top_candidates": len(self.top_candidates),
            "processing_time": f"{self.processing_time:.2f}s",
            "best_score": max([c.composite_score for c in self.top_candidates]) if self.top_candidates else 0,
            "tool_versions": self.tool_versions,
        }
