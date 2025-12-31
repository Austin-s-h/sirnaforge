"""Pydantic models for transcript annotation data structures."""

from typing import Optional

from pydantic import BaseModel, ConfigDict, Field

from sirnaforge.config.reference_policy import ReferenceChoice


class Interval(BaseModel):
    """Genomic interval with start, end, and optional strand information."""

    seq_region_name: str = Field(description="Chromosome or sequence region name")
    start: int = Field(ge=0, description="Start position (0-based or 1-based depending on source)")
    end: int = Field(ge=0, description="End position (inclusive)")
    strand: Optional[int] = Field(
        default=None, description="Strand: 1 for forward, -1 for reverse, None for unstranded"
    )

    model_config = ConfigDict(frozen=True)

    def __hash__(self) -> int:
        """Make Interval hashable for use in sets/dicts."""
        return hash((self.seq_region_name, self.start, self.end, self.strand))


class TranscriptAnnotation(BaseModel):
    """Comprehensive transcript annotation from genomic databases.

    Contains transcript metadata, genomic coordinates, exon/CDS structure,
    and source provenance for reproducibility.
    """

    # Core identifiers
    transcript_id: str = Field(description="Transcript stable identifier (e.g., ENST00000269305)")
    gene_id: str = Field(description="Parent gene stable identifier (e.g., ENSG00000141510)")
    symbol: Optional[str] = Field(default=None, description="Gene symbol (e.g., TP53)")
    biotype: Optional[str] = Field(default=None, description="Transcript biotype (e.g., protein_coding)")

    # Genomic coordinates
    seq_region_name: str = Field(description="Chromosome or contig name (e.g., '17', 'chr17')")
    start: int = Field(ge=1, description="Genomic start position (1-based)")
    end: int = Field(ge=1, description="Genomic end position (inclusive)")
    strand: int = Field(description="Strand: 1 for forward/sense, -1 for reverse/antisense")

    # Structural features
    exons: list[Interval] = Field(default_factory=list, description="List of exon intervals")
    cds: list[Interval] = Field(
        default_factory=list, description="List of coding sequence (CDS) intervals, empty if non-coding"
    )

    # Source metadata
    provider: str = Field(description="Annotation provider (e.g., 'ensembl_rest', 'vep')")
    endpoint: Optional[str] = Field(default=None, description="API endpoint or base URL used for retrieval")
    reference_choice: Optional[str] = Field(
        default=None, description="Reference assembly/release info (e.g., 'GRCh38.p13', 'ensembl_release_110')"
    )

    model_config = ConfigDict(use_enum_values=True)


class TranscriptAnnotationBundle(BaseModel):
    """Collection of transcript annotations with resolution tracking.

    Bundles multiple transcript annotations from a single query,
    tracks which IDs were successfully resolved, and maintains
    reference provenance.
    """

    transcripts: dict[str, TranscriptAnnotation] = Field(
        default_factory=dict, description="Map of transcript_id -> annotation"
    )
    unresolved: list[str] = Field(default_factory=list, description="IDs that could not be resolved")
    reference_choice: ReferenceChoice = Field(
        description="Reference selection metadata (explicit, default, or disabled)"
    )

    model_config = ConfigDict(arbitrary_types_allowed=True, use_enum_values=True)

    @property
    def resolved_count(self) -> int:
        """Number of successfully resolved transcripts."""
        return len(self.transcripts)

    @property
    def unresolved_count(self) -> int:
        """Number of IDs that could not be resolved."""
        return len(self.unresolved)

    @property
    def total_requested(self) -> int:
        """Total number of IDs requested (resolved + unresolved)."""
        return self.resolved_count + self.unresolved_count
