"""Pandera schemas for siRNAforge data validation.

This module defines pandera schemas for validating the structure and content
of various table-like outputs from the siRNAforge pipeline.

Modern schemas using class-based approach with type annotations for improved
type safety, error reporting, and maintainability.

Use schemas: MySchema.validate(df) - validation errors provide detailed feedback.
"""

import re
from typing import Any, Callable, TypeVar, cast

import pandas as pd
import pandera.pandas as pa
from pandera.pandas import DataFrameModel, Field
from pandera.typing.pandas import Series

# Typed alias for pandera's dataframe_check decorator to satisfy mypy
F = TypeVar("F", bound=Callable[..., Any])
# Pandera's dataframe_check has a complex decorator signature; cast for mypy.
dataframe_check_typed = cast(Callable[[F], F], pa.dataframe_check)


# Custom validation functions for bioinformatics data
def valid_nucleotide_sequence(sequence: str) -> bool:
    """Validate nucleotide sequence contains only valid bases."""
    return bool(re.match(r"^[ATCGUN-]*$", sequence.upper())) if isinstance(sequence, str) else False


def valid_rna_sequence(sequence: str) -> bool:
    """Validate RNA sequence contains only valid RNA bases."""
    return bool(re.match(r"^[AUCG]+$", sequence.upper())) if isinstance(sequence, str) else False


def valid_strand(strand: str) -> bool:
    """Validate strand orientation."""
    return strand in ["+", "-"] if isinstance(strand, str) else False


def valid_codon(codon: str) -> bool:
    """Validate start/stop codon sequences."""
    valid_start = ["ATG"]
    valid_stop = ["TAA", "TAG", "TGA"]
    return codon.upper() in (valid_start + valid_stop) if isinstance(codon, str) else False


def sirna_length_range(sequence: str) -> bool:
    """Validate siRNA sequence length is in typical range."""
    return 19 <= len(sequence) <= 23 if isinstance(sequence, str) else False


# Schema configuration for better error reporting
class SchemaConfig:
    """Common configuration for all schemas."""

    coerce = True
    strict = True  # Ensure no unexpected columns
    ordered = False  # Allow columns in any order


class SiRNACandidateSchema(DataFrameModel):
    """Schema for siRNA candidate results (CSV output).

    Validates the structure and content of siRNA design results with
    comprehensive checks for biological validity and data integrity.
    """

    class Config(SchemaConfig):
        """Schema configuration with improved error reporting."""

        description = "siRNA candidate validation schema"
        title = "SiRNA Design Results"

    # Identity fields
    id: Series[str] = Field(description="Unique identifier for siRNA candidate")
    transcript_id: Series[str] = Field(description="Source transcript identifier")
    position: Series[int] = Field(ge=1, description="1-based position in transcript")

    # Sequence fields with validation
    guide_sequence: Series[str] = Field(description="Guide (antisense) sequence")
    passenger_sequence: Series[str] = Field(description="Passenger (sense) sequence")

    # Quantitative properties
    gc_content: Series[float] = Field(ge=0.0, le=100.0, description="GC content percentage")
    asymmetry_score: Series[float] = Field(ge=0.0, le=1.0, description="Thermodynamic asymmetry score")
    paired_fraction: Series[float] = Field(
        ge=0.0, le=1.0, description="Fraction of paired bases in secondary structure"
    )

    # Thermodynamic details (nullable if backend not available)
    structure: Series[Any] = Field(description="Predicted secondary structure (dot-bracket)", nullable=True)
    mfe: Series[float] = Field(description="Minimum free energy of guide structure", nullable=True)
    duplex_stability_dg: Series[float] = Field(description="Guide:passenger duplex stability ΔG", nullable=True)
    duplex_stability_score: Series[float] = Field(
        ge=0.0, le=1.0, description="Normalized duplex stability score [0-1]", nullable=True
    )
    dg_5p: Series[float] = Field(description="5' end duplex ΔG (positions 1-7)", nullable=True)
    dg_3p: Series[float] = Field(description="3' end duplex ΔG (positions 15-21)", nullable=True)
    delta_dg_end: Series[float] = Field(description="ΔΔG = dg_5p - dg_3p", nullable=True)
    melting_temp_c: Series[float] = Field(description="Estimated duplex melting temperature (°C)", nullable=True)

    # Off-target analysis results
    off_target_count: Series[int] = Field(ge=0, description="Number of potential off-targets")

    # Transcript hit metrics
    transcript_hit_count: Series[int] = Field(ge=0, description="Number of input transcripts containing this guide")
    transcript_hit_fraction: Series[float] = Field(
        ge=0.0, le=1.0, description="Fraction of input transcripts hit by this guide"
    )

    # Scoring results
    composite_score: Series[float] = Field(ge=0.0, le=100.0, description="Final composite score")

    # Quality control: allow legacy booleans or new status strings
    passes_filters: Series[Any] = Field(description="Filter status: True/False (legacy) or PASS/reason code (new)")

    @dataframe_check_typed
    def check_passes_filters_values(cls, df: pd.DataFrame) -> bool:
        """Ensure passes_filters contains allowed values: bool or allowed status strings."""
        allowed = {"PASS", "GC_OUT_OF_RANGE", "POLY_RUNS", "EXCESS_PAIRING", "LOW_ASYMMETRY"}
        series = df["passes_filters"]

        def _ok(v: Any) -> bool:
            if isinstance(v, bool):
                return True
            return isinstance(v, str) and v in allowed

        return bool(series.map(_ok).all())

    @dataframe_check_typed
    def check_sequence_lengths(cls, df: pd.DataFrame) -> bool:
        """Validate siRNA sequence lengths are in acceptable range."""
        guide_lengths = df["guide_sequence"].str.len()
        passenger_lengths = df["passenger_sequence"].str.len()
        return bool(guide_lengths.between(19, 23).all() and passenger_lengths.between(19, 23).all())

    @dataframe_check_typed
    def check_nucleotide_sequences(cls, df: pd.DataFrame) -> bool:
        """Validate sequences contain only valid nucleotide bases (DNA or RNA)."""
        guide_valid = df["guide_sequence"].str.match(r"^[ATCGU]+$").all()
        passenger_valid = df["passenger_sequence"].str.match(r"^[ATCGU]+$").all()
        return bool(guide_valid and passenger_valid)


class ORFValidationSchema(DataFrameModel):
    """Schema for ORF validation report (tab-delimited output).

    Validates open reading frame analysis results with proper handling
    of nullable fields and bioinformatics constraints.
    """

    class Config(SchemaConfig):
        """Schema configuration."""

        description = "ORF validation analysis schema"
        title = "ORF Analysis Results"
        strict = False  # Allow different dtypes for nullable fields

    # Basic sequence information
    transcript_id: Series[str] = Field(description="Transcript identifier")
    sequence_length: Series[int] = Field(ge=1, description="Sequence length in nucleotides")
    gc_content: Series[float] = Field(ge=0.0, le=100.0, description="Overall GC content percentage")

    # ORF detection results
    orfs_found: Series[int] = Field(ge=0, description="Number of ORFs detected")
    has_valid_orf: Series[bool] = Field(description="Whether transcript has valid ORF")

    # Longest ORF details (nullable if no ORF found) - allowing flexible types
    longest_orf_start: Series[Any] = Field(description="Start position of longest ORF (1-based)", nullable=True)
    longest_orf_end: Series[Any] = Field(description="End position of longest ORF (1-based)", nullable=True)
    longest_orf_length: Series[Any] = Field(description="Length of longest ORF in nucleotides", nullable=True)
    longest_orf_frame: Series[Any] = Field(description="Reading frame of longest ORF (0, 1, or 2)", nullable=True)

    # Codon information (nullable)
    start_codon: Series[Any] = Field(description="Start codon of longest ORF", nullable=True)
    stop_codon: Series[Any] = Field(description="Stop codon of longest ORF", nullable=True)

    # ORF-specific GC content
    orf_gc_content: Series[Any] = Field(description="GC content of longest ORF", nullable=True)

    # UTR/CDS characterization can be present in outputs but is not required by schema.
    # We intentionally omit these from the schema so tests with legacy columns still pass,
    # while Config.strict=False allows extra columns like utr5_length, utr3_length, etc.


class OffTargetHitsSchema(DataFrameModel):
    """Schema for off-target analysis results (TSV output).

    Validates off-target prediction results with relaxed constraints
    to accommodate various external tool outputs.
    """

    class Config(SchemaConfig):
        """Schema configuration with relaxed strictness for external tool outputs."""

        description = "Off-target analysis results schema"
        title = "Off-target Prediction Results"
        strict = False  # More lenient for external tool outputs

    # Query information
    qname: Series[str] = Field(description="Query sequence name/ID")

    # Target identification (nullable for no-hit cases)
    target_id: Series[Any] = Field(description="Target sequence identifier", nullable=True)
    species: Series[Any] = Field(description="Target species", nullable=True)

    # Genomic location (nullable)
    chromosome: Series[Any] = Field(description="Chromosome/contig", nullable=True)
    position: Series[Any] = Field(description="Genomic position", nullable=True)
    strand: Series[Any] = Field(description="Strand orientation", nullable=True)

    # Alignment metrics (nullable)
    mismatches: Series[Any] = Field(description="Number of mismatches", nullable=True)
    alignment_score: Series[Any] = Field(description="Alignment score", nullable=True)
    offtarget_score: Series[Any] = Field(description="Off-target penalty score", nullable=True)

    # Target sequence with alignment (nullable)
    target_sequence: Series[Any] = Field(description="Target sequence with alignment", nullable=True)
