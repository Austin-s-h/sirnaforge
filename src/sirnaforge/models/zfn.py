"""Pydantic models for Zinc Finger Nuclease (ZFN) design and off-target search."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Literal, cast

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo

from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.typed_decorators import field_validator_typed

logger = get_logger(__name__)


class Strand(str, Enum):
    """Genome strand labels for reported sites."""

    PLUS = "+"
    MINUS = "-"


class DimerMode(str, Enum):
    """What dimerization outcomes to consider when enumerating cut sites."""

    HETERODIMER_ONLY = "heterodimer_only"
    INCLUDE_HOMODIMERS = "include_homodimers"


class ZFNAlgorithm(str, Enum):
    """Off-target ranking model family for ZFN sites."""

    HOMOLOGY = "homology"
    CONSERVED_G = "conserved_g"
    ZFN_V2 = "zfn_v2"


class MatchOrientation(str, Enum):
    """Orientation/ordering of half-sites at a genomic locus."""

    LR = "L...R"
    RL = "R...L"
    LL = "L...L"
    RR = "R...R"


class IUPACMode(str, Enum):
    """How degenerate bases are handled in the query."""

    NONE = "none"
    ALLOW_IUPAC = "allow_iupac"
    EXPAND_IUPAC = "expand_iupac"


class ZFNHalfSiteConstraints(BaseModel):
    """Constraints on one ZFN half-site."""

    model_config = ConfigDict(extra="forbid")

    min_len: int = Field(default=9, ge=6, description="Minimum half-site length in bp")
    max_len: int = Field(default=18, ge=6, description="Maximum half-site length in bp")
    max_mismatches: int = Field(default=2, ge=0, le=6, description="Max mismatches allowed per half-site search")
    seed_len_from_fokI: int | None = Field(
        default=6,
        ge=1,
        description="If set, enforce stricter mismatch handling nearest FokI",
    )
    seed_max_mismatches: int | None = Field(
        default=1,
        ge=0,
        description="Max mismatches allowed in the seed region nearest FokI",
    )
    iupac_mode: IUPACMode = Field(
        default=IUPACMode.ALLOW_IUPAC,
        description="Interpretation of IUPAC ambiguity codes in query half-sites",
    )
    window_stride: int = Field(
        default=1,
        ge=1,
        le=50,
        description="Sliding-window stride in bp for exhaustive genomic scanning (1 = fully exhaustive)",
    )

    @field_validator_typed("max_len")
    @classmethod
    def max_len_ge_min_len(cls, v: int, info: ValidationInfo) -> int:
        """Validate max_len >= min_len."""
        if "min_len" in info.data and v < info.data["min_len"]:
            raise ValueError("max_len must be >= min_len")
        return v


class ZFNSpacerConstraints(BaseModel):
    """Spacer constraints between Left and Right half-sites at a cleavage site."""

    model_config = ConfigDict(extra="forbid")

    allowed_spacer_lengths: list[int] = Field(
        default_factory=lambda: [5, 6, 7],
        min_length=1,
        description="Allowed spacer lengths between half-sites",
    )
    require_opposite_strands: bool = Field(
        default=True,
        description="Require half-sites on opposite strands in canonical ZFN configuration",
    )

    @field_validator_typed("allowed_spacer_lengths")
    @classmethod
    def normalize_spacer_lengths(cls, v: list[int]) -> list[int]:
        """Normalize spacer lengths to a sorted, unique positive list."""
        parsed_values: set[int] = set()
        for length in v:
            parsed = int(length)
            if parsed > 0:
                parsed_values.add(parsed)
        cleaned = sorted(parsed_values)
        if not cleaned:
            raise ValueError("allowed_spacer_lengths must contain at least one positive integer")
        return cleaned


class ZFNOffTargetFilterCriteria(BaseModel):
    """Criteria to fail a candidate ZFN due to predicted off-target burden."""

    model_config = ConfigDict(extra="forbid")

    max_total_sites: int | None = Field(default=None, ge=0, description="Max total predicted cut sites allowed")
    max_exonic_sites: int | None = Field(default=None, ge=0, description="Max predicted sites in exons")
    max_promoter_sites: int | None = Field(default=None, ge=0, description="Max predicted sites in promoters")
    min_site_score_to_count: float = Field(
        default=70.0,
        ge=0,
        le=100,
        description="Only count predicted sites with score >= this threshold",
    )


class ZFNScoringWeights(BaseModel):
    """Weights for composite ZFN candidate scoring."""

    model_config = ConfigDict(extra="forbid")

    on_target_quality: float = Field(default=0.35, ge=0, le=1, description="On-target site preference weight")
    off_target_specificity: float = Field(default=0.45, ge=0, le=1, description="Off-target avoidance weight")
    manufacturability: float = Field(default=0.20, ge=0, le=1, description="Constructability / assembly weight")

    @field_validator_typed("manufacturability")
    @classmethod
    def weights_sum_to_one(cls, v: float, info: ValidationInfo) -> float:
        """Validate scoring weights sum to ~1.0."""
        total = sum(info.data.values()) + v
        if not (0.95 <= total <= 1.05):
            raise ValueError(f"ZFN scoring weights must sum to 1.0, got {total}")
        return v


class ZFNMutationType(str, Enum):
    """Allowed mutation categories for ZFN sub-finger constraints."""

    SUBSTITUTION = "substitution"
    TRANSITION = "transition"
    TRANSVERSION = "transversion"
    INSERTION = "insertion"
    DELETION = "deletion"


class ZFNSubfingerMutationConstraint(BaseModel):
    """Mutation allowance definition for one ZFN sub-finger."""

    subfinger_index: int = Field(ge=1, description="1-based sub-finger index")
    max_mutations: int = Field(ge=0, description="Maximum allowed mutations for this sub-finger")
    mutation_types: list[ZFNMutationType] = Field(
        min_length=1,
        description="Allowed mutation types for this sub-finger",
    )


class ZFNDefaultSubfingerMutationConstraint(BaseModel):
    """Default mutation allowance applied to each sub-finger."""

    max_mutations: int = Field(ge=0, description="Maximum allowed mutations for each sub-finger")
    mutation_types: list[ZFNMutationType] = Field(
        min_length=1,
        description="Allowed mutation types for each sub-finger",
    )


class ZFNOverallMutationConstraint(BaseModel):
    """Global mutation allowance applied across all sub-fingers."""

    max_mutations: int = Field(ge=0, description="Maximum allowed mutations across all sub-fingers")
    mutation_types: list[ZFNMutationType] = Field(
        min_length=1,
        description="Allowed mutation types for the global budget",
    )


class ZFNMutationConstraints(BaseModel):
    """Composite container for all ZFN sub-finger and global mutation budgets."""

    model_config = ConfigDict(extra="forbid")

    subfinger_mutations: list[ZFNSubfingerMutationConstraint] = Field(
        default_factory=lambda: cast(list[ZFNSubfingerMutationConstraint], []),
        description="Per-sub-finger mutation allowances",
    )
    default_subfinger_mutation: ZFNDefaultSubfingerMutationConstraint | None = Field(
        default=None,
        description="Default mutation allowance applied to each sub-finger",
    )
    overall_mutations: list[ZFNOverallMutationConstraint] = Field(
        default_factory=lambda: cast(list[ZFNOverallMutationConstraint], []),
        description="Global mutation budgets applied across all sub-fingers",
    )


class ZFNShardingConfig(BaseModel):
    """Optional chromosome/chunk sharding controls for scalable ZFN search."""

    model_config = ConfigDict(extra="forbid")

    enabled: bool = Field(default=True, description="Enable chromosome/chunk sharding")
    chunk_size_bp: int = Field(default=20_000_000, ge=1, description="Nominal chunk size in base pairs")
    overlap_bp: int = Field(
        default=50,
        ge=0,
        description="Chunk overlap in base pairs (auto-raised to safe minimum when needed)",
    )
    chromosomes: list[str] = Field(
        default_factory=lambda: cast(list[str], []),
        description=(
            "Optional chromosome filter tokens; empty means all contigs. "
            "Supports exact labels (chr3/3), ranges (1-5), groups (autosomes, sex), "
            "and glob patterns (chrUn_*)."
        ),
    )
    max_workers: int = Field(default=1, ge=1, le=128, description="Parallel shard workers for in-process search")

    @field_validator_typed("chromosomes")
    @classmethod
    def validate_chromosomes(cls, v: list[str]) -> list[str]:
        """Normalize chromosome names while preserving input order."""
        cleaned: list[str] = []
        for raw in v:
            token = raw.strip()
            if token and token not in cleaned:
                cleaned.append(token)
        return cleaned


class ZFNDesignParameters(BaseModel):
    """Top-level configuration for ZFN pair evaluation and off-target search."""

    model_config = ConfigDict(extra="forbid")

    search_space_reference: str | None = Field(
        default="ensembl_human_hg38_primary",
        description=(
            "Cache source key for retrieval via GenomeManager "
            "(e.g., ensembl_human_hg38_primary); ignored when search_space_fasta is set"
        ),
    )
    search_space_fasta: str | None = Field(
        default=None,
        description="Explicit FASTA path/URL for genomic search space (retrieved and cached when needed)",
    )
    search_space_index: str | None = Field(
        default=None,
        description="Optional pre-built index prefix for search space (reserved for indexed backends)",
    )
    left_half_site: str = Field(description="Left ZFN half-site sequence")
    right_half_site: str = Field(description="Right ZFN half-site sequence")
    half_site_constraints: ZFNHalfSiteConstraints = Field(default_factory=ZFNHalfSiteConstraints)
    spacer_constraints: ZFNSpacerConstraints = Field(default_factory=ZFNSpacerConstraints)
    dimer_mode: DimerMode = Field(default=DimerMode.HETERODIMER_ONLY)
    algorithm: ZFNAlgorithm = Field(default=ZFNAlgorithm.ZFN_V2)
    top_n_sites: int = Field(default=5000, ge=1, description="Maximum number of retained ranked sites")
    report_n_sites: int = Field(default=200, ge=1, description="Number of sites to include in reports")
    off_target_filters: ZFNOffTargetFilterCriteria = Field(default_factory=ZFNOffTargetFilterCriteria)
    scoring: ZFNScoringWeights = Field(default_factory=ZFNScoringWeights)
    mutation_constraints: ZFNMutationConstraints | None = Field(
        default=None,
        description="Optional per-sub-finger and global mutation budgets for manufacturability scoring",
    )
    sharding: ZFNShardingConfig = Field(default_factory=ZFNShardingConfig)

    @field_validator_typed("left_half_site", "right_half_site")
    @classmethod
    def validate_dna_bases_and_iupac(cls, v: str) -> str:
        """Validate half-sites as DNA/IUPAC strings."""
        allowed = set("ACGTNRYWSKMBDHV")
        vv = v.upper().replace(" ", "")
        if not vv:
            raise ValueError("Half-site sequence cannot be empty")
        if any(base not in allowed for base in vv):
            raise ValueError(f"Invalid base(s) in half-site: {v}")
        return vv

    @field_validator_typed("right_half_site")
    @classmethod
    def validate_half_site_lengths(cls, v: str, info: ValidationInfo) -> str:
        """Validate half-site lengths against configured range."""
        constraints = info.data.get("half_site_constraints")
        if constraints is None:
            return v

        left = info.data.get("left_half_site", "")
        for label, seq in (("left_half_site", left), ("right_half_site", v)):
            if seq and not (constraints.min_len <= len(seq) <= constraints.max_len):
                raise ValueError(
                    f"{label} length {len(seq)} out of range [{constraints.min_len}, {constraints.max_len}]"
                )
        return v

    @field_validator_typed("search_space_reference")
    @classmethod
    def validate_search_space_source(cls, v: str | None) -> str | None:
        """Normalize search-space reference keys."""
        if v is None:
            return None
        normalized = v.strip()
        if not normalized:
            raise ValueError("search_space_reference cannot be empty when provided")
        return normalized

    @field_validator_typed("search_space_fasta")
    @classmethod
    def validate_search_space_fasta(cls, v: str | None) -> str | None:
        """Normalize explicit FASTA path/URL strings."""
        if v is None:
            return None
        normalized = v.strip()
        if not normalized:
            raise ValueError("search_space_fasta cannot be empty when provided")
        return normalized

    def canonical_search_contract(self) -> ZFNSearchContract:
        """Return canonical, JSON-serializable ZFN search contract for reports/workflows."""
        return ZFNSearchContract(
            left_half_site=self.left_half_site,
            right_half_site=self.right_half_site,
            allowed_spacer_lengths=list(self.spacer_constraints.allowed_spacer_lengths),
            max_mismatches_per_half_site=self.half_site_constraints.max_mismatches,
            seed_len_from_foki=self.half_site_constraints.seed_len_from_fokI,
            seed_max_mismatches=self.half_site_constraints.seed_max_mismatches,
            dimer_mode=self.dimer_mode,
            algorithm=self.algorithm,
            require_opposite_strands=self.spacer_constraints.require_opposite_strands,
            orientation_convention="L...R genomic ordering",
        )


class ZFNSearchContract(BaseModel):
    """Canonical internal contract for ZFN search/ranking behavior."""

    model_config = ConfigDict(extra="forbid")

    left_half_site: str = Field(description="Canonical uppercase left half-site sequence")
    right_half_site: str = Field(description="Canonical uppercase right half-site sequence")
    allowed_spacer_lengths: list[int] = Field(description="Sorted unique spacer lengths in bp")
    max_mismatches_per_half_site: int = Field(ge=0, description="Mismatch budget per half-site")
    seed_len_from_foki: int | None = Field(default=None, ge=1)
    seed_max_mismatches: int | None = Field(default=None, ge=0)
    dimer_mode: DimerMode = Field(description="Pairing mode for hit assembly")
    algorithm: ZFNAlgorithm = Field(description="Algorithm family used for scoring")
    require_opposite_strands: bool = Field(description="Whether same-strand pairs are disallowed")
    orientation_convention: str = Field(description="Human-readable orientation convention")


class GenomicAnnotationConfig(BaseModel):
    """Optional annotation configuration for off-target site classification."""

    model_config = ConfigDict(extra="forbid")

    annotation_path: str | None = Field(default=None, description="Local GTF/GFF path for exon/intron annotation")
    annotation_reference: str | None = Field(
        default=None,
        description=(
            "Cache key/URI for annotation retrieval. Current implementation resolves local paths and cached downloads."
        ),
    )
    cache_dir: str | None = Field(default=None, description="Optional cache directory for resolved annotation assets")
    promoter_upstream_bp: int = Field(default=2000, ge=0, description="Promoter window upstream of TSS")
    promoter_downstream_bp: int = Field(default=200, ge=0, description="Promoter window downstream of TSS")

    @field_validator_typed("annotation_path", "annotation_reference", "cache_dir")
    @classmethod
    def validate_optional_paths(cls, v: str | None) -> str | None:
        """Normalize optional string paths/keys."""
        if v is None:
            return None
        normalized = v.strip()
        if not normalized:
            raise ValueError("Path/reference values cannot be empty strings")
        return normalized

    def resolved_annotation_path(self) -> Path | None:
        """Return resolved local annotation path when available."""
        if self.annotation_path:
            return Path(self.annotation_path)
        if self.annotation_reference:
            candidate = Path(self.annotation_reference)
            if candidate.exists():
                return candidate
        return None


class ZFNOffTargetSite(BaseModel):
    """One predicted cleavage site for a ZFN pair at a genomic locus."""

    model_config = ConfigDict(extra="forbid")

    site_id: str = Field(description="Unique predicted site identifier")
    chrom: str = Field(description="Chromosome/contig name")
    start_1based: int = Field(ge=1, description="1-based start coordinate of full site")
    end_1based: int = Field(ge=1, description="1-based end coordinate of full site")
    strand: Strand = Field(description="Reported strand for this site")
    orientation: MatchOrientation = Field(description="Half-site orientation at locus")
    spacer_len: int = Field(ge=0, description="Spacer length between half-sites")
    sequence: str = Field(description="Full site sequence: half-site + spacer + half-site")
    left_mismatches: int = Field(ge=0)
    right_mismatches: int = Field(ge=0)
    left_seed_mismatches: int = Field(default=0, ge=0)
    right_seed_mismatches: int = Field(default=0, ge=0)
    left_mismatch_positions: list[int] = Field(default_factory=list)
    right_mismatch_positions: list[int] = Field(default_factory=list)
    total_mismatches: int = Field(ge=0)
    score: float = Field(ge=0, le=100, description="Predicted off-target similarity/cleavage likelihood score")
    score_components: dict[str, float] = Field(default_factory=dict)
    dimer_compatible: bool = Field(default=True)
    region: Literal["exon", "promoter", "intron", "intergenic", "unknown"] = Field(default="unknown")
    nearest_gene: str | None = Field(default=None)
    left_aligned: str = Field(description="Alignment string for left half-site")
    right_aligned: str = Field(description="Alignment string for right half-site")


class ZFNCandidate(BaseModel):
    """A ZFN pair candidate with summary metrics and composite score."""

    model_config = ConfigDict(extra="forbid")

    id: str = Field(description="Unique candidate identifier")
    left_half_site: str = Field(description="Selected left half-site")
    right_half_site: str = Field(description="Selected right half-site")
    allowed_spacers: list[int] = Field(description="Spacer lengths considered")
    predicted_sites_total: int = Field(default=0, ge=0)
    predicted_sites_exonic: int = Field(default=0, ge=0)
    predicted_sites_promoter: int = Field(default=0, ge=0)
    worst_site_score: float | None = Field(default=None, ge=0, le=100)
    best_offtarget_score: float | None = Field(default=None, ge=0, le=100)
    passes_offtarget_filters: bool | str = Field(default=True, description="PASS if within thresholds, else reason")
    component_scores: dict[str, float] = Field(default_factory=dict)
    composite_score: float = Field(ge=0, le=100, description="Overall design score")
    top_offtargets: list[ZFNOffTargetSite] = Field(default_factory=lambda: cast(list[ZFNOffTargetSite], []))


class ZFNDesignResult(BaseModel):
    """Complete results for ZFN pair evaluation and off-target search."""

    model_config = ConfigDict(extra="forbid")

    parameters: ZFNDesignParameters = Field(description="Parameters used")
    annotation: GenomicAnnotationConfig | None = Field(default=None)
    candidates: list[ZFNCandidate] = Field(default_factory=lambda: cast(list[ZFNCandidate], []))
    off_target_sites: list[ZFNOffTargetSite] = Field(default_factory=lambda: cast(list[ZFNOffTargetSite], []))
    processing_time_s: float = Field(ge=0)
    tool_versions: dict[str, str] = Field(default_factory=dict)

    def save_offtargets_csv(self, filepath: str) -> pd.DataFrame:
        """Save off-target site table to CSV."""
        rows: list[dict[str, Any]] = []
        for site in self.off_target_sites:
            rows.append(
                {
                    "site_id": site.site_id,
                    "chrom": site.chrom,
                    "start_1based": site.start_1based,
                    "end_1based": site.end_1based,
                    "strand": site.strand.value,
                    "orientation": site.orientation.value,
                    "spacer_len": site.spacer_len,
                    "sequence": site.sequence,
                    "left_mismatches": site.left_mismatches,
                    "right_mismatches": site.right_mismatches,
                    "left_seed_mismatches": site.left_seed_mismatches,
                    "right_seed_mismatches": site.right_seed_mismatches,
                    "left_mismatch_positions": ",".join(str(pos) for pos in site.left_mismatch_positions),
                    "right_mismatch_positions": ",".join(str(pos) for pos in site.right_mismatch_positions),
                    "total_mismatches": site.total_mismatches,
                    "score": site.score,
                    "dimer_compatible": site.dimer_compatible,
                    "score_components": site.score_components,
                    "region": site.region,
                    "nearest_gene": site.nearest_gene,
                    "left_aligned": site.left_aligned,
                    "right_aligned": site.right_aligned,
                }
            )
        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)
        return df

    def get_summary(self) -> dict[str, Any]:
        """Generate a compact run summary."""
        return {
            "candidates": len(self.candidates),
            "off_target_sites": len(self.off_target_sites),
            "processing_time_s": round(self.processing_time_s, 3),
            "algorithm": self.parameters.algorithm.value,
            "search_space_reference": self.parameters.search_space_reference,
        }
