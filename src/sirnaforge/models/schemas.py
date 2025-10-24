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
    """Validate nucleotide sequence contains only standard DNA/RNA bases.

    Args:
        sequence: Nucleotide sequence to validate

    Returns:
        True if sequence contains only A, T, C, G, U, N, or - characters
    """
    return bool(re.match(r"^[ATCGUN-]*$", sequence.upper())) if isinstance(sequence, str) else False


def valid_rna_sequence(sequence: str) -> bool:
    """Validate RNA sequence contains only standard RNA bases.

    Args:
        sequence: RNA sequence to validate

    Returns:
        True if sequence contains only A, U, C, G characters
    """
    return bool(re.match(r"^[AUCG]+$", sequence.upper())) if isinstance(sequence, str) else False


def valid_strand(strand: str) -> bool:
    """Validate genomic strand orientation notation.

    Args:
        strand: Strand indicator to validate

    Returns:
        True if strand is "+" or "-"
    """
    return strand in ["+", "-"] if isinstance(strand, str) else False


def valid_codon(codon: str) -> bool:
    """Validate start or stop codon sequences.

    Args:
        codon: Three-nucleotide codon to validate

    Returns:
        True if codon is a valid start (ATG) or stop (TAA, TAG, TGA) codon
    """
    valid_start = ["ATG"]
    valid_stop = ["TAA", "TAG", "TGA"]
    return codon.upper() in (valid_start + valid_stop) if isinstance(codon, str) else False


def sirna_length_range(sequence: str) -> bool:
    """Validate siRNA sequence length is in functional range.

    Args:
        sequence: siRNA sequence to validate

    Returns:
        True if sequence length is between 19-23 nucleotides (typical siRNA range)
    """
    return 19 <= len(sequence) <= 23 if isinstance(sequence, str) else False


# Schema configuration for better error reporting
class SchemaConfig:
    """Common configuration settings for all pandera schemas.

    Provides consistent validation behavior across all siRNAforge data schemas
    with type coercion, strict column checking, and flexible column ordering.
    """

    coerce = True
    strict = True  # Ensure no unexpected columns
    ordered = False  # Allow columns in any order


class SiRNACandidateSchema(DataFrameModel):
    """Validation schema for siRNA candidate results (CSV output).

    Ensures data integrity and biological validity of siRNA design results with
    comprehensive checks for sequence composition, thermodynamic parameters,
    and scoring metrics. Includes optimal value ranges for key metrics based
    on research-backed thermodynamic principles.

    Expected columns include sequences, thermodynamic scores (asymmetry, MFE,
    duplex stability), off-target counts, and composite quality scores.
    """

    class Config(SchemaConfig):
        """Schema configuration with improved error reporting."""

        description = "siRNA candidate validation schema"
        title = "SiRNA Design Results"
        # Allow DataFrames without modification columns (add as empty/null)
        add_missing_columns = True
        strict = False  # Don't reject DataFrames missing modification columns

    # Identity fields
    id: Series[str] = Field(description="Unique siRNA candidate identifier")
    transcript_id: Series[str] = Field(description="Source transcript ID (e.g., ENST00000123456)")
    position: Series[int] = Field(ge=1, description="1-based start position in transcript")

    # Sequence fields with validation
    guide_sequence: Series[str] = Field(description="Guide strand sequence (antisense, 19-23 nt)")
    passenger_sequence: Series[str] = Field(description="Passenger strand sequence (sense, 19-23 nt)")

    # Quantitative properties
    gc_content: Series[float] = Field(ge=0.0, le=100.0, description="GC content % (optimal: 35-60%)")
    asymmetry_score: Series[float] = Field(ge=0.0, le=1.0, description="Thermodynamic asymmetry score (optimal: ≥0.65)")
    paired_fraction: Series[float] = Field(
        ge=0.0, le=1.0, description="Fraction of paired bases in secondary structure (optimal: 0.4-0.8)"
    )

    # Thermodynamic details (nullable if backend not available)
    structure: Series[Any] = Field(description="RNA secondary structure in dot-bracket notation", nullable=True)
    mfe: Series[float] = Field(description="Minimum free energy in kcal/mol (optimal: -2 to -8)", nullable=True)
    duplex_stability_dg: Series[float] = Field(
        description="siRNA duplex ΔG in kcal/mol (optimal: -15 to -25)", nullable=True
    )
    duplex_stability_score: Series[float] = Field(
        ge=0.0, le=1.0, description="Normalized duplex stability score [0-1]", nullable=True
    )
    dg_5p: Series[float] = Field(description="5' end ΔG kcal/mol (positions 1-7)", nullable=True)
    dg_3p: Series[float] = Field(description="3' end ΔG kcal/mol (positions 15-21)", nullable=True)
    delta_dg_end: Series[float] = Field(
        description="End asymmetry ΔΔG = dg_3p - dg_5p (optimal: +2 to +6)", nullable=True
    )
    melting_temp_c: Series[float] = Field(description="Duplex melting temperature °C (optimal: 60-78°C)", nullable=True)

    # Off-target analysis results
    off_target_count: Series[int] = Field(ge=0, description="Number of potential off-target sites (goal: ≤3)")

    # Transcript hit metrics
    transcript_hit_count: Series[int] = Field(ge=0, description="Number of input transcripts containing this guide")
    transcript_hit_fraction: Series[float] = Field(
        ge=0.0, le=1.0, description="Fraction of input transcripts hit by this guide (1.0 = all transcripts)"
    )

    # Scoring results
    composite_score: Series[float] = Field(
        ge=0.0, le=100.0, description="Overall siRNA quality score (higher is better)"
    )

    # Quality control: allow legacy booleans or new status strings
    passes_filters: Series[Any] = Field(description="Filter result: PASS or failure reason (GC_OUT_OF_RANGE, etc.)")

    # Chemical modification columns (optional, nullable)
    # Using add_missing_columns to auto-add with null values
    guide_overhang: Series[str] = Field(
        description="Guide strand 3' overhang sequence (e.g., dTdT, UU)",
        nullable=True,
        coerce=True,
    )
    guide_modifications: Series[str] = Field(
        description="Guide strand modification summary",
        nullable=True,
        coerce=True,
    )
    passenger_overhang: Series[str] = Field(
        description="Passenger strand 3' overhang sequence",
        nullable=True,
        coerce=True,
    )
    passenger_modifications: Series[str] = Field(
        description="Passenger strand modification summary",
        nullable=True,
        coerce=True,
    )

    @dataframe_check_typed
    def check_passes_filters_values(cls, df: pd.DataFrame) -> bool:
        """Ensure passes_filters contains allowed filter status values."""
        allowed = {
            "PASS",
            "GC_OUT_OF_RANGE",
            "POLY_RUNS",
            "EXCESS_PAIRING",
            "LOW_ASYMMETRY",
        }
        series = df["passes_filters"]

        def _ok(v: Any) -> bool:
            if isinstance(v, bool):
                return True
            return isinstance(v, str) and v in allowed

        return bool(series.map(_ok).all())

    @dataframe_check_typed
    def check_sequence_lengths(cls, df: pd.DataFrame) -> bool:
        """Validate siRNA sequences are in functional range (19-23 nt)."""
        guide_lengths = df["guide_sequence"].str.len()
        passenger_lengths = df["passenger_sequence"].str.len()
        return bool(guide_lengths.between(19, 23).all() and passenger_lengths.between(19, 23).all())

    @dataframe_check_typed
    def check_nucleotide_sequences(cls, df: pd.DataFrame) -> bool:
        """Validate sequences contain only valid RNA/DNA bases."""
        guide_valid = df["guide_sequence"].str.match(r"^[ATCGU]+$").all()
        passenger_valid = df["passenger_sequence"].str.match(r"^[ATCGU]+$").all()
        return bool(guide_valid and passenger_valid)


class ORFValidationSchema(DataFrameModel):
    """Validation schema for open reading frame analysis results (tab-delimited output).

    Validates ORF detection and characterization results with proper handling
    of nullable fields for cases where no valid ORF is found. Includes metrics
    for transcript composition, ORF boundaries, codon usage, and GC content
    within coding regions.

    Used to validate outputs from ORF analysis tools and ensure data consistency
    for downstream siRNA target validation.
    """

    class Config(SchemaConfig):
        """Schema configuration."""

        description = "ORF validation analysis schema"
        title = "ORF Analysis Results"
        strict = False  # Allow different dtypes for nullable fields

    # Basic sequence information
    transcript_id: Series[str] = Field(description="Transcript identifier")
    sequence_length: Series[int] = Field(ge=1, description="Total transcript length in nucleotides")
    gc_content: Series[float] = Field(ge=0.0, le=100.0, description="Overall transcript GC content %")

    # ORF detection results
    orfs_found: Series[int] = Field(ge=0, description="Total number of open reading frames detected")
    has_valid_orf: Series[bool] = Field(description="True if transcript contains a valid protein-coding ORF")

    # Longest ORF details (nullable if no ORF found) - allowing flexible types
    longest_orf_start: Series[Any] = Field(description="Start position of longest ORF (1-based)", nullable=True)
    longest_orf_end: Series[Any] = Field(description="End position of longest ORF (1-based)", nullable=True)
    longest_orf_length: Series[Any] = Field(description="Longest ORF length in nucleotides", nullable=True)
    longest_orf_frame: Series[Any] = Field(description="Reading frame of longest ORF (0, 1, or 2)", nullable=True)

    # Codon information (nullable)
    start_codon: Series[Any] = Field(description="Start codon of longest ORF (usually ATG)", nullable=True)
    stop_codon: Series[Any] = Field(description="Stop codon of longest ORF (TAA, TAG, or TGA)", nullable=True)

    # ORF-specific GC content
    orf_gc_content: Series[Any] = Field(description="GC content % of the longest ORF region", nullable=True)

    # UTR/CDS characterization can be present in outputs but is not required by schema.
    # We intentionally omit these from the schema so tests with legacy columns still pass,
    # while Config.strict=False allows extra columns like utr5_length, utr3_length, etc.


class OffTargetHitsSchema(DataFrameModel):
    """Validation schema for off-target analysis results (TSV output).

    Validates off-target prediction results from external alignment tools with
    relaxed constraints to accommodate diverse tool outputs. Handles nullable
    fields for cases with no significant off-target matches and various
    alignment scoring systems.

    Supports results from BWA, BLAST, and other sequence similarity tools
    used for genome-wide off-target analysis.
    """

    class Config(SchemaConfig):
        """Schema configuration with relaxed strictness for external tool outputs."""

        description = "Off-target analysis results schema"
        title = "Off-target Prediction Results"
        strict = False  # More lenient for external tool outputs

    # Query information
    qname: Series[str] = Field(description="Query siRNA sequence identifier")

    # Target identification (nullable for no-hit cases)
    target_id: Series[Any] = Field(description="Off-target sequence/gene identifier", nullable=True)
    species: Series[Any] = Field(description="Target organism/species name", nullable=True)

    # Genomic location (nullable)
    chromosome: Series[Any] = Field(description="Chromosome or contig name", nullable=True)
    position: Series[Any] = Field(description="Genomic coordinate of potential off-target", nullable=True)
    strand: Series[Any] = Field(description="Strand orientation (+ or -)", nullable=True)

    # Alignment metrics (nullable)
    mismatches: Series[Any] = Field(description="Number of base mismatches with target", nullable=True)
    alignment_score: Series[Any] = Field(description="Sequence alignment score", nullable=True)
    offtarget_score: Series[Any] = Field(description="Off-target risk penalty score", nullable=True)

    # Target sequence with alignment (nullable)
    target_sequence: Series[Any] = Field(description="Aligned target sequence with mismatch notation", nullable=True)


class ModificationSummarySchema(DataFrameModel):
    """Validation schema for chemical modification summary data.

    Validates modification summary dictionaries returned by get_modification_summary()
    when converted to DataFrame format. Used for validating modification-enhanced
    siRNA candidate outputs.
    """

    class Config(SchemaConfig):
        """Schema configuration."""

        description = "Chemical modification summary schema"
        title = "Modification Summary"
        strict = False  # Allow extra fields

    guide_overhang: Series[str] = Field(description="Guide strand 3' overhang (e.g., dTdT, UU)")
    guide_modifications: Series[str] = Field(description="Guide modifications (e.g., 2OMe(11)+PS(2))")
    passenger_overhang: Series[str] = Field(description="Passenger strand 3' overhang")
    passenger_modifications: Series[str] = Field(description="Passenger modifications")

    @dataframe_check_typed
    def check_modification_format(cls, df: pd.DataFrame) -> bool:
        """Validate modification summary format."""
        # Valid formats: "none", "", "2OMe(11)", "2OMe(11)+PS(2)"
        for col in ["guide_modifications", "passenger_modifications"]:
            if col in df.columns:
                series = df[col]
                # Empty strings and "none" are valid
                valid_pattern = r"^(none|)$|^([A-Z0-9]+\(\d+\)(\+[A-Z0-9]+\(\d+\))*)$"
                if not series.str.match(valid_pattern, na=False).all():
                    return False
        return True
