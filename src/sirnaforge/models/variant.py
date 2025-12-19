"""Pydantic models for genomic variant data structures."""

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, ConfigDict, Field


class VariantMode(str, Enum):
    """Mode for how variants should be handled in siRNA design."""

    TARGET = "target"  # Design siRNAs specifically targeting the variant allele
    AVOID = "avoid"  # Avoid designing siRNAs that overlap variant positions
    BOTH = "both"  # Generate candidates for both reference and alternate alleles


class VariantSource(str, Enum):
    """Trusted sources for variant data, ordered by priority."""

    CLINVAR = "clinvar"  # Priority 1: Clinical variant database
    ENSEMBL = "ensembl"  # Priority 2: Ensembl Variation database
    DBSNP = "dbsnp"  # Priority 3: dbSNP reference SNP database
    LOCAL_VCF = "local-vcf"  # Local VCF file


class ClinVarSignificance(str, Enum):
    """ClinVar clinical significance classifications."""

    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain significance"
    LIKELY_BENIGN = "Likely benign"
    BENIGN = "Benign"
    CONFLICTING = "Conflicting interpretations of pathogenicity"
    OTHER = "Other"


class VariantRecord(BaseModel):
    """Complete variant record with annotations from multiple sources."""

    model_config = ConfigDict(extra="forbid")

    # Core variant identity
    id: Optional[str] = Field(default=None, description="Variant identifier (e.g., rs12345, ClinVar accession)")
    chr: str = Field(description="Chromosome (normalized to 'chrN' format for GRCh38)")
    pos: int = Field(ge=1, description="1-based genomic position")
    ref: str = Field(min_length=1, description="Reference allele")
    alt: str = Field(min_length=1, description="Alternate allele")

    # Assembly
    assembly: str = Field(default="GRCh38", description="Reference genome assembly (only GRCh38 supported)")

    # Source tracking
    sources: list[VariantSource] = Field(
        default_factory=list, description="List of sources where this variant was found"
    )

    # ClinVar annotations
    clinvar_significance: Optional[ClinVarSignificance] = Field(
        default=None, description="ClinVar clinical significance classification"
    )

    # Population frequency
    af: Optional[float] = Field(
        default=None, ge=0.0, le=1.0, description="Global allele frequency (from gnomAD/Ensembl/dbSNP)"
    )

    # Population-specific allele frequencies (optional, for geographic targeting)
    population_afs: dict[str, float] = Field(
        default_factory=dict,
        description="Population-specific allele frequencies (e.g., {'AFR': 0.15, 'EUR': 0.02, 'EAS': 0.08})",
    )

    # Functional annotations
    annotations: dict[str, Any] = Field(
        default_factory=dict, description="Additional functional annotations (gene, transcript consequences, etc.)"
    )

    # Data provenance
    provenance: dict[str, Any] = Field(
        default_factory=dict, description="Source-specific metadata (query timestamp, API version, etc.)"
    )

    def to_vcf_style(self) -> str:
        """Return variant in VCF-style coordinate format: chr:pos:ref:alt."""
        return f"{self.chr}:{self.pos}:{self.ref}:{self.alt}"

    def get_max_population_af(self) -> Optional[float]:
        """Get the maximum allele frequency across all populations.

        Returns:
            Maximum population-specific AF, or None if no population data available
        """
        if not self.population_afs:
            return None
        return max(self.population_afs.values())

    def get_effective_af_for_mode(self, mode: "VariantMode") -> Optional[float]:  # noqa: F821
        """Get the effective allele frequency based on variant mode.

        For 'avoid' mode: Use max population AF if available (to avoid SNPs
        prevalent in any geographic group), otherwise use global AF.

        For 'target' or 'both' modes: Use global AF (targets most common alleles).

        Args:
            mode: Variant mode (target/avoid/both)

        Returns:
            Effective allele frequency for filtering, or None if no AF data
        """
        if mode == VariantMode.AVOID:
            # In avoid mode, use max population AF to avoid SNPs prevalent
            # in any geographic group (e.g., >10% in one population even if
            # global AF is only 1%)
            max_pop_af = self.get_max_population_af()
            if max_pop_af is not None:
                return max_pop_af
        # For target/both modes, or if no population data, use global AF
        return self.af

    def get_primary_source(self) -> Optional[VariantSource]:
        """Get the highest priority source for this variant."""
        priority_order = [VariantSource.CLINVAR, VariantSource.ENSEMBL, VariantSource.DBSNP, VariantSource.LOCAL_VCF]
        for source in priority_order:
            if source in self.sources:
                return source
        return self.sources[0] if self.sources else None


class VariantQueryType(str, Enum):
    """Types of variant query identifiers."""

    RSID = "rsid"  # rsID format: rs12345
    COORDINATE = "coordinate"  # VCF-style: chr:pos:ref:alt
    HGVS = "hgvs"  # HGVS notation: NM_000546.6:c.215C>G


class VariantQuery(BaseModel):
    """Parsed variant query with normalized components."""

    model_config = ConfigDict(extra="forbid")

    raw_input: str = Field(description="Original user-provided query string")
    query_type: VariantQueryType = Field(description="Detected query type")

    # Parsed components (populated based on query type)
    rsid: Optional[str] = Field(default=None, description="rsID if query_type is RSID")
    chr: Optional[str] = Field(default=None, description="Chromosome if coordinate or HGVS")
    pos: Optional[int] = Field(default=None, ge=1, description="Position if coordinate")
    ref: Optional[str] = Field(default=None, description="Reference allele if coordinate")
    alt: Optional[str] = Field(default=None, description="Alternate allele if coordinate")
    hgvs: Optional[str] = Field(default=None, description="HGVS string if query_type is HGVS")

    # Assembly
    assembly: str = Field(default="GRCh38", description="Target assembly (only GRCh38 supported)")
