"""Pydantic models for siRNA design data structures."""

from enum import Enum
from typing import Any, Callable, Optional, TypeVar

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator

# mypy-friendly typed alias for pydantic's untyped decorator factory
F = TypeVar("F", bound=Callable[..., Any])
FieldValidatorFactory = Callable[..., Callable[[F], F]]
field_validator_typed: FieldValidatorFactory = field_validator


class FilterCriteria(BaseModel):
    """Criteria for filtering siRNA candidates."""

    gc_min: float = Field(default=30.0, ge=0, le=100, description="Minimum GC content percentage")
    gc_max: float = Field(default=52.0, ge=0, le=100, description="Maximum GC content percentage")
    max_poly_runs: int = Field(default=3, ge=1, description="Maximum consecutive identical nucleotides")
    max_paired_fraction: float = Field(default=0.6, ge=0, le=1, description="Maximum secondary structure pairing")
    min_asymmetry_score: float = Field(default=0.0, ge=0, le=1, description="Minimum thermodynamic asymmetry")

    @field_validator_typed("gc_max")
    @classmethod
    def gc_max_greater_than_min(cls, v: float, info: ValidationInfo) -> float:
        if "gc_min" in info.data and v < info.data["gc_min"]:
            raise ValueError("gc_max must be greater than or equal to gc_min")
        return v


class ScoringWeights(BaseModel):
    """Weights for composite scoring components."""

    asymmetry: float = Field(default=0.25, ge=0, le=1, description="Thermodynamic asymmetry weight")
    gc_content: float = Field(default=0.20, ge=0, le=1, description="GC content optimization weight")
    accessibility: float = Field(default=0.25, ge=0, le=1, description="Target accessibility weight")
    off_target: float = Field(default=0.20, ge=0, le=1, description="Off-target avoidance weight")
    empirical: float = Field(default=0.10, ge=0, le=1, description="Empirical rules weight")

    @field_validator_typed("empirical")
    @classmethod
    def weights_sum_to_one(cls, v: float, info: ValidationInfo) -> float:
        total = sum(info.data.values()) + v
        if not (0.95 <= total <= 1.05):  # Allow small floating point errors
            raise ValueError(f"Scoring weights must sum to 1.0, got {total}")
        return v


class DesignParameters(BaseModel):
    """Complete parameters for siRNA design."""

    model_config = ConfigDict(extra="forbid")

    # Basic parameters
    sirna_length: int = Field(default=21, ge=19, le=23, description="siRNA length in nucleotides")
    top_n: int = Field(default=50, ge=1, le=1000, description="Number of top candidates to return")

    # Filtering criteria
    filters: FilterCriteria = Field(default_factory=FilterCriteria)

    # Scoring weights
    scoring: ScoringWeights = Field(default_factory=ScoringWeights)

    # Optional analysis parameters
    avoid_snps: bool = Field(default=True, description="Avoid known SNP positions")
    check_off_targets: bool = Field(default=True, description="Perform off-target analysis")
    predict_structure: bool = Field(default=True, description="Predict secondary structures")

    # File paths (optional)
    # TODO: review snp incorporation feature
    snp_file: Optional[str] = Field(default=None, description="Path to SNP VCF file")
    # Review genome index passing / FASTA selection
    genome_index: Optional[str] = Field(default=None, description="Path to genome index for off-targets")


class SequenceType(str, Enum):
    """Types of input sequences."""

    TRANSCRIPT = "transcript"
    GENOMIC = "genomic"
    CDS = "cds"
    UTR = "utr"


class SiRNACandidate(BaseModel):
    """Individual siRNA candidate with all computed properties."""

    model_config = ConfigDict(extra="forbid")

    # Identity
    id: str = Field(description="Unique identifier for this candidate")
    transcript_id: str = Field(description="Source transcript identifier")
    position: int = Field(ge=1, description="1-based position in transcript")

    # Sequences
    guide_sequence: str = Field(min_length=19, max_length=23, description="Guide (antisense) sequence")
    passenger_sequence: str = Field(min_length=19, max_length=23, description="Passenger (sense) sequence")

    # Basic properties
    gc_content: float = Field(ge=0, le=100, description="GC content percentage")
    length: int = Field(ge=19, le=23, description="Sequence length")

    # Thermodynamic properties
    asymmetry_score: float = Field(ge=0, le=1, description="Thermodynamic asymmetry score")
    duplex_stability: Optional[float] = Field(default=None, description="Duplex stability (ΔG)")

    # Secondary structure
    structure: Optional[str] = Field(default=None, description="Predicted secondary structure")
    mfe: Optional[float] = Field(default=None, description="Minimum free energy")
    paired_fraction: float = Field(default=0.0, ge=0, le=1, description="Fraction of paired bases")

    # Off-target analysis
    off_target_count: int = Field(default=0, ge=0, description="Number of potential off-targets")
    off_target_penalty: float = Field(default=0.0, ge=0, description="Off-target penalty score")

    # Transcript hit metrics (how many input transcripts this guide hits)
    transcript_hit_count: int = Field(
        default=1, ge=0, description="Number of input transcripts that contain this guide sequence"
    )
    transcript_hit_fraction: float = Field(
        default=1.0, ge=0, le=1, description="Fraction of input transcripts hit by this guide sequence"
    )

    # Composite scoring
    component_scores: dict[str, float] = Field(default_factory=dict, description="Individual component scores")
    composite_score: float = Field(ge=0, le=100, description="Final composite score")

    # Quality flags
    passes_filters: bool = Field(default=True, description="Passes all quality filters")
    quality_issues: list[str] = Field(default_factory=list, description="List of quality concerns")

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

    def to_fasta(self) -> str:
        """Return FASTA format representation."""
        return f">{self.id}\n{self.guide_sequence}\n"


class DesignResult(BaseModel):
    """Complete results from siRNA design process."""

    model_config = ConfigDict(extra="forbid")

    # Input information
    input_file: str = Field(description="Path to input FASTA file")
    parameters: DesignParameters = Field(description="Parameters used for design")

    # Results
    candidates: list[SiRNACandidate] = Field(description="All siRNA candidates")
    top_candidates: list[SiRNACandidate] = Field(description="Top-scoring candidates")

    # Summary statistics
    total_sequences: int = Field(ge=0, description="Number of input sequences processed")
    total_candidates: int = Field(ge=0, description="Total candidates generated")
    filtered_candidates: int = Field(ge=0, description="Candidates passing filters")

    # Processing metadata
    processing_time: float = Field(ge=0, description="Processing time in seconds")
    tool_versions: dict[str, str] = Field(default_factory=dict, description="Tool versions used")

    def save_csv(self, filepath: str) -> None:
        """Save results to CSV file."""
        df_data = []
        for candidate in self.candidates:
            row = {
                "id": candidate.id,
                "transcript_id": candidate.transcript_id,
                "position": candidate.position,
                "guide_sequence": candidate.guide_sequence,
                "passenger_sequence": candidate.passenger_sequence,
                "gc_content": candidate.gc_content,
                "asymmetry_score": candidate.asymmetry_score,
                "paired_fraction": candidate.paired_fraction,
                "off_target_count": candidate.off_target_count,
                "transcript_hit_count": candidate.transcript_hit_count,
                "transcript_hit_fraction": candidate.transcript_hit_fraction,
                "composite_score": candidate.composite_score,
                "passes_filters": candidate.passes_filters,
            }
            df_data.append(row)

        df = pd.DataFrame(df_data)
        df.to_csv(filepath, index=False)

    def get_summary(self) -> dict[str, Any]:
        """Get summary statistics."""
        return {
            "input_sequences": self.total_sequences,
            "total_candidates": self.total_candidates,
            "filtered_candidates": self.filtered_candidates,
            "top_candidates": len(self.top_candidates),
            "processing_time": f"{self.processing_time:.2f}s",
            "best_score": max([c.composite_score for c in self.top_candidates]) if self.top_candidates else 0,
            "tool_versions": self.tool_versions,
        }
