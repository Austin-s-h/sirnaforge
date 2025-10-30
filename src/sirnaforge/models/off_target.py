"""Pydantic models for off-target analysis data structures.

This module provides validated data models for:
- BWA alignment results (both genome and miRNA)
- Aggregated off-target summaries
- Analysis metadata and statistics

Using Pydantic ensures type safety, automatic validation, and clean
serialization to JSON/TSV formats.
"""

import re
from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Optional

from pydantic import BaseModel, Field, field_validator


class AlignmentStrand(str, Enum):
    """Genomic strand orientation."""

    FORWARD = "+"
    REVERSE = "-"


class AnalysisMode(str, Enum):
    """BWA alignment analysis mode."""

    MIRNA_SEED = "mirna_seed"
    TRANSCRIPTOME = "transcriptome"


# Base classes for shared functionality
class BaseAlignmentHit(BaseModel, ABC):
    """Base class for alignment hits with common fields and validators.

    This abstract base class contains all shared fields and validation logic
    for both off-target and miRNA alignment hits.
    """

    model_config = {"frozen": False, "validate_assignment": True}

    # Query information (common to all hits)
    qname: str = Field(description="Query siRNA sequence identifier")
    qseq: str = Field(description="Query sequence (siRNA candidate)")

    # Alignment details (common to all hits)
    coord: int = Field(ge=0, description="Alignment coordinate (0-based)")
    strand: AlignmentStrand = Field(description="Strand orientation")
    cigar: str = Field(description="CIGAR string describing alignment")
    mapq: int = Field(ge=0, le=255, description="Mapping quality (0-255)")

    # Scoring (common to all hits)
    as_score: Optional[int] = Field(default=None, description="Alignment score (AS tag)")
    nm: int = Field(ge=0, description="Edit distance / number of mismatches")
    seed_mismatches: int = Field(ge=0, description="Mismatches in seed region (positions 2-8)")
    offtarget_score: float = Field(ge=0.0, description="Off-target penalty score")

    @field_validator("qseq")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        """Ensure sequence contains only valid nucleotide characters."""
        if not v:
            raise ValueError("Query sequence cannot be empty")
        valid_chars = set("ATCGUN-")
        if not set(v.upper()).issubset(valid_chars):
            raise ValueError(f"Invalid nucleotide characters in sequence: {v}")
        return v.upper()

    @field_validator("cigar")
    @classmethod
    def validate_cigar(cls, v: str) -> str:
        """Basic CIGAR string validation."""
        if not v or v == "*":
            return v
        # CIGAR should be numbers followed by operation letters
        if not re.match(r"^(\d+[MIDNSHPX=])+$", v):
            raise ValueError(f"Invalid CIGAR string format: {v}")
        return v

    @abstractmethod
    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for TSV/JSON serialization."""
        pass

    def to_tsv_row(self) -> str:
        """Convert to TSV row format."""
        d = self.to_dict()
        return "\t".join(str(d[k]) for k in d)

    @classmethod
    @abstractmethod
    def tsv_header(cls) -> str:
        """Get TSV header line."""
        pass


class OffTargetHit(BaseAlignmentHit):
    """Single off-target alignment hit from BWA analysis.

    Represents one potential off-target binding site identified by
    sequence alignment against a reference genome or transcriptome.
    """

    # Target information (specific to genome/transcriptome hits)
    rname: str = Field(description="Reference sequence identifier (chromosome, transcript, etc.)")

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for TSV/JSON serialization."""
        return {
            "qname": self.qname,
            "qseq": self.qseq,
            "rname": self.rname,
            "coord": self.coord,
            "strand": self.strand.value,
            "cigar": self.cigar,
            "mapq": self.mapq,
            "as_score": self.as_score if self.as_score is not None else "NA",
            "nm": self.nm,
            "seed_mismatches": self.seed_mismatches,
            "offtarget_score": self.offtarget_score,
        }

    @classmethod
    def tsv_header(cls) -> str:
        """Get TSV header line."""
        return "qname\tqseq\trname\tcoord\tstrand\tcigar\tmapq\tas_score\tnm\tseed_mismatches\tofftarget_score"


class MiRNAHit(BaseAlignmentHit):
    """Single miRNA seed match hit from BWA analysis.

    Represents a potential miRNA-like seed match identified by
    alignment against miRNA databases.
    """

    # miRNA-specific information
    species: str = Field(description="Species of miRNA database")
    database: str = Field(description="miRNA database name (mirgenedb, mirbase, etc.)")
    mirna_id: str = Field(description="miRNA identifier (e.g., hsa-miR-21-5p)")

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for TSV/JSON serialization."""
        return {
            "qname": self.qname,
            "qseq": self.qseq,
            "species": self.species,
            "database": self.database,
            "mirna_id": self.mirna_id,
            "coord": self.coord,
            "strand": self.strand.value,
            "cigar": self.cigar,
            "mapq": self.mapq,
            "as_score": self.as_score if self.as_score is not None else "NA",
            "nm": self.nm,
            "seed_mismatches": self.seed_mismatches,
            "offtarget_score": self.offtarget_score,
        }

    @classmethod
    def tsv_header(cls) -> str:
        """Get TSV header line."""
        return (
            "qname\tqseq\tspecies\tdatabase\tmirna_id\tcoord\tstrand\tcigar\tmapq\tas_score\t"
            "nm\tseed_mismatches\tofftarget_score"
        )


class BaseSummary(BaseModel):
    """Base class for analysis summary statistics with common metadata fields."""

    model_config = {"frozen": False}

    candidate_id: str = Field(description="Candidate identifier")
    total_sequences: int = Field(ge=0, description="Number of sequences analyzed")
    total_hits: int = Field(ge=0, description="Total hits found")

    # Common metadata
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Analysis timestamp")
    status: str = Field(default="completed", description="Analysis status")


class AnalysisSummary(BaseSummary):
    """Summary statistics for a single candidate's off-target analysis."""

    species: str = Field(description="Target species analyzed")
    mode: AnalysisMode = Field(description="Analysis mode used")

    # Alignment statistics
    mean_mapq: Optional[float] = Field(default=None, ge=0.0, le=255.0, description="Mean mapping quality")
    mean_mismatches: Optional[float] = Field(default=None, ge=0.0, description="Mean number of mismatches")
    mean_seed_mismatches: Optional[float] = Field(default=None, ge=0.0, description="Mean seed region mismatches")


class MiRNASummary(BaseSummary):
    """Summary statistics for miRNA seed match analysis."""

    mirna_database: str = Field(description="miRNA database used")
    species_analyzed: list[str] = Field(description="List of species analyzed")
    hits_per_species: dict[str, int] = Field(default_factory=dict, description="Hit counts by species")

    # Analysis parameters
    parameters: dict[str, Any] = Field(
        default_factory=lambda: {"seed_start": 2, "seed_end": 8, "mode": "mirna_seed"},
        description="Analysis parameters",
    )


class BaseAggregatedSummary(BaseModel):
    """Base class for aggregated analysis summaries with common fields."""

    model_config = {"frozen": False}

    species_analyzed: list[str] = Field(description="List of species analyzed")
    analysis_files_processed: int = Field(ge=0, description="Number of analysis files processed")

    # Common file paths
    combined_tsv: Optional[Path] = Field(default=None, description="Path to combined TSV file")
    combined_json: Optional[Path] = Field(default=None, description="Path to combined JSON file")
    summary_file: Optional[Path] = Field(default=None, description="Path to summary file")

    # Common metadata
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat(), description="Aggregation timestamp")


class AggregatedOffTargetSummary(BaseAggregatedSummary):
    """Summary of aggregated off-target results across multiple candidates and genomes."""

    total_results: int = Field(ge=0, description="Total off-target hits aggregated")


class AggregatedMiRNASummary(BaseAggregatedSummary):
    """Summary of aggregated miRNA results across multiple candidates."""

    total_mirna_hits: int = Field(ge=0, description="Total miRNA seed matches")
    mirna_database: str = Field(description="miRNA database used")

    hits_per_species: dict[str, int] = Field(default_factory=dict, description="Hit counts by species")
    hits_per_candidate: dict[str, int] = Field(default_factory=dict, description="Hit counts by candidate")

    total_candidates: int = Field(ge=0, description="Total candidates analyzed")
