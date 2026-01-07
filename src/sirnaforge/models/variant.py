"""Pydantic models for genomic variant data structures."""

from enum import Enum
from typing import Any

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


class EnsemblMapping(BaseModel):
    """Ensembl variation mapping information.

    Represents genomic mapping data for a variant from the Ensembl Variation API.
    Contains coordinate and allele information for a specific genomic location.
    """

    model_config = ConfigDict(extra="allow")  # Allow additional fields from API

    location: str = Field(description="Genomic location in chr:start-end format")
    allele_string: str = Field(description="Allele string (e.g., 'C/A' or 'C/A/G/T')")
    assembly_name: str = Field(description="Genome assembly name (e.g., 'GRCh38')")
    seq_region_name: str = Field(description="Sequence region (chromosome) name")
    strand: int = Field(description="Strand orientation (1 = forward, -1 = reverse)")
    start: int = Field(description="Start position (1-based)")
    end: int = Field(description="End position (1-based)")
    coord_system: str = Field(description="Coordinate system (e.g., 'chromosome')")


class EnsemblPopulationFrequency(BaseModel):
    """Population frequency data from Ensembl.

    Contains allele frequency information for a specific population from
    sources like 1000 Genomes, gnomAD, ExAC, etc.
    """

    model_config = ConfigDict(extra="allow")

    population: str = Field(description="Population identifier (e.g., 'AFR', 'EUR', '1000GENOMES:phase_3:AFR')")
    frequency: float | None = Field(description="Allele frequency in this population", ge=0.0, le=1.0)
    allele_count: int | None = Field(default=None, description="Count of alleles observed", ge=0)
    allele_number: int | None = Field(default=None, description="Total number of alleles", ge=0)
    allele: str | None = Field(default=None, description="Allele for which frequency is reported")


class EnsemblVariationResponse(BaseModel):
    """Response model for Ensembl Variation API.

    Represents the complete response from the Ensembl REST API variation endpoint.
    Contains variant metadata, mappings, frequencies, and clinical information.
    """

    model_config = ConfigDict(extra="allow")  # Allow additional fields

    name: str = Field(description="Variant name (rsID)")
    var_class: str = Field(description="Variant class (SNP, insertion, deletion, etc.)")
    source: str = Field(description="Data source description")
    most_severe_consequence: str | None = Field(description="Most severe functional consequence")
    MAF: float | None = Field(description="Minor allele frequency", ge=0.0, le=1.0)
    minor_allele: str | None = Field(description="Minor allele")
    ambiguity: str = Field(description="IUPAC ambiguity code")
    mappings: list[EnsemblMapping] = Field(description="Genomic mappings")
    clinical_significance: list[str] | None = Field(description="Clinical significance terms")
    synonyms: list[str] | None = Field(description="Alternative identifiers")
    evidence: list[str] | None = Field(description="Supporting evidence types")
    populations: list[EnsemblPopulationFrequency | None] = Field(
        default_factory=list, description="Population frequency data"
    )


class ClinVarVariationResponse(BaseModel):
    """Response model for ClinVar variation summary.

    Represents the response from NCBI ClinVar E-utilities esummary endpoint.
    Contains clinical significance, variation details, and associated conditions.
    """

    model_config = ConfigDict(extra="allow")

    uid: str = Field(description="ClinVar variation ID")
    obj_type: str | None = Field(description="Object type (e.g., 'single nucleotide variant')")
    accession: str | None = Field(description="Accession number (e.g., 'VCV000376655')")
    title: str | None = Field(description="Variation title with HGVS notation")
    germline_classification: dict | None = Field(description="Germline classification data")
    clinical_impact_classification: dict | None = Field(description="Clinical impact classification")
    variation_set: list[dict] | None = Field(description="Variation set information")
    genes: list[dict] | None = Field(description="Associated genes")
    molecular_consequence_list: list[str] | None = Field(description="Molecular consequences")


class VariantRecord(BaseModel):
    """Complete variant record with annotations from multiple sources."""

    model_config = ConfigDict(extra="forbid")

    # Core variant identity
    id: str | None = Field(default=None, description="Variant identifier (e.g., rs12345, ClinVar accession)")
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
    clinvar_significance: ClinVarSignificance | None = Field(
        default=None, description="ClinVar clinical significance classification"
    )

    # Population frequency
    af: float | None = Field(
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

    def get_max_population_af(self) -> float | None:
        """Get the maximum allele frequency across all populations.

        Returns:
            Maximum population-specific AF, or None if no population data available
        """
        if not self.population_afs:
            return None
        return max(self.population_afs.values())

    def get_effective_af_for_mode(self, mode: "VariantMode") -> float | None:  # noqa: F821
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

    def get_primary_source(self) -> VariantSource | None:
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
    rsid: str | None = Field(default=None, description="rsID if query_type is RSID")
    chr: str | None = Field(default=None, description="Chromosome if coordinate or HGVS")
    pos: int | None = Field(default=None, ge=1, description="Position if coordinate")
    ref: str | None = Field(default=None, description="Reference allele if coordinate")
    alt: str | None = Field(default=None, description="Alternate allele if coordinate")
    hgvs: str | None = Field(default=None, description="HGVS string if query_type is HGVS")

    # Assembly
    assembly: str = Field(default="GRCh38", description="Target assembly (only GRCh38 supported)")
