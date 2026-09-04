"""Pydantic models for siRNA design data structures."""

import json
from enum import Enum
from typing import Any

import pandas as pd
from pandera.typing import DataFrame
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo

from sirnaforge.models.modifications import StrandMetadata, StrandRole
from sirnaforge.models.schemas import SiRNACandidateSchema
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import get_modification_summary
from sirnaforge.utils.typed_decorators import check_types_typed, field_validator_typed, model_validator_typed

logger = get_logger(__name__)

# Sequence-length bounds for observed/input siRNA-like sequences the off-target
# engine analyzes (SiRNACandidate). The classic siRNA guide is ~19-23 nt, but
# the engine also handles longer species (e.g. Dicer 3' read-through isoforms,
# extended guides). We accept up to ENGINE_MAX_GUIDE_LEN and only *warn* above
# the RECOMMENDED biological max. Note: DesignParameters.sirna_length (what the
# tool *designs*) stays capped at 23 — this only relaxes what it will *analyze*.
MIN_GUIDE_LEN = 19
RECOMMENDED_MAX_GUIDE_LEN = 23
ENGINE_MAX_GUIDE_LEN = 40

# Default RISC-loading asymmetry gate, shared with
# ThermodynamicCalculator.is_thermodynamically_favorable so the two cannot drift.
DEFAULT_MIN_ASYMMETRY_SCORE = 0.65

# Attainable range of SiRNADesigner._calculate_empirical_score. The simplified
# Reynolds rule adjusts a 0.5 base score by +/-0.1 per criterion, so it can never
# reach 1.0; min_empirical_score is bounded by these so an unsatisfiable
# threshold fails at construction instead of rejecting every candidate.
EMPIRICAL_SCORE_MIN = 0.4
EMPIRICAL_SCORE_MAX = 0.7
DEFAULT_MIN_EMPIRICAL_SCORE = 0.5

# Canonical composite-score term names, shared by ScoringWeights and the scorer.
COMPOSITE_TERM_NAMES = (
    "asymmetry",
    "gc_content",
    "accessibility",
    "empirical",
    "off_target",
    "isoform_coverage",
    "conservation",
)


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
        default=DEFAULT_MIN_ASYMMETRY_SCORE,
        ge=0.3,
        le=1,
        description=(
            "Minimum thermodynamic asymmetry score for guide strand selection into RISC. "
            "Applied to SiRNACandidate.asymmetry_score. "
            "Higher values (0.65-0.85) promote correct 5' end instability for effective strand loading."
        ),
    )

    # Empirical (simplified Reynolds) design-rule filter
    min_empirical_score: float = Field(
        default=DEFAULT_MIN_EMPIRICAL_SCORE,
        ge=EMPIRICAL_SCORE_MIN,
        le=EMPIRICAL_SCORE_MAX,
        description=(
            "Minimum empirical design-rule score. Applied to the 'empirical' component score, "
            f"whose attainable range is {EMPIRICAL_SCORE_MIN}-{EMPIRICAL_SCORE_MAX}; the default rejects only "
            "candidates penalised at guide position 19 with no G/C at position 1."
        ),
    )

    @field_validator_typed("gc_max")
    @classmethod
    def gc_max_greater_than_min(cls, v: float, info: ValidationInfo) -> float:
        """Validate that gc_max is greater than or equal to gc_min."""
        if "gc_min" in info.data and v < info.data["gc_min"]:
            raise ValueError("gc_max must be greater than or equal to gc_min")
        return v


class OffTargetFilterCriteria(BaseModel):
    """Filtering criteria for off-target analysis results.

    Controls which siRNA candidates fail due to excessive off-target potential.
    """

    # Genuine off-target threshold (on-target, ortholog and repeat-mediated hits excluded)
    max_off_target_count: int | None = Field(
        default=3,
        ge=0,
        description="Maximum genuine off-target sites (goal: ≤3). Excludes on-target, ortholog and repeat hits.",
    )

    # Transcriptome off-target thresholds
    # The three transcriptome thresholds below count GENUINE off-target hits only: on-target,
    # ortholog and repeat-mediated hits are classified out before they reach these counters.
    max_transcriptome_hits_0mm: int | None = Field(
        default=1,
        ge=0,
        description="Maximum perfect-match genuine off-target transcriptome hits (excludes on-target isoforms)",
    )
    max_transcriptome_hits_1mm: int | None = Field(
        default=10, ge=0, description="Maximum 1-mismatch genuine off-target hits (typical: 5-10, None = no limit)"
    )
    max_transcriptome_hits_2mm: int | None = Field(
        default=50, ge=0, description="Maximum 2-mismatch genuine off-target hits (typical: 20-50, None = no limit)"
    )
    # Read against transcriptome_hits_seed_0mm. This is the only *targeted* gate on a partial
    # (clipped or gapped) hit whose seed paired perfectly: such a hit has a guide-level nm > 2, so
    # it lands in no mismatch stratum and max_transcriptome_hits_{0,1,2}mm never see it. The only
    # other thing that counts it is the blunt max_off_target_count ceiling (default 3).
    # Unlike the three mismatch thresholds this one counts ALL screened species, matching the
    # reported transcriptome_hits_seed_0mm column exactly. Defaults to None (off) because a safe
    # ceiling has not been calibrated against truth data; set it to opt in.
    max_transcriptome_seed_perfect: int | None = Field(
        default=None, ge=0, description="Maximum transcriptome hits with perfect seed (positions 2-8, None = no limit)"
    )

    # miRNA off-target thresholds
    max_mirna_perfect_seed: int | None = Field(
        default=0, ge=0, description="Maximum perfect miRNA seed matches (typical: 3-5, None = no limit)"
    )
    max_mirna_1mm_seed: int | None = Field(
        default=10, ge=0, description="Maximum 1-mismatch miRNA seed hits (typical: 10-20, None = no limit)"
    )
    fail_on_high_risk_mirna: bool = Field(
        default=True, description="Fail if high-risk miRNA hits detected (perfect seed + offtarget_score < 5.0)"
    )

    # Combined off-target threshold
    max_total_offtarget_hits: int | None = Field(
        default=None, ge=0, description="Maximum total off-target hits (transcriptome + miRNA, None = no limit)"
    )


class ScoringWeights(BaseModel):
    """Relative weights for composite siRNA scoring components."""

    asymmetry: float = Field(
        default=0.12, ge=0, le=1, description="Thermodynamic asymmetry weight (guide strand selection)"
    )
    gc_content: float = Field(
        default=0.10, ge=0, le=1, description="GC content optimization weight (stability balance)"
    )
    accessibility: float = Field(
        default=0.13, ge=0, le=1, description="Target accessibility weight (secondary structure)"
    )
    empirical: float = Field(
        default=0.15, ge=0, le=1, description="Empirical design rules weight (established patterns)"
    )
    off_target: float = Field(
        default=0.25,
        ge=0,
        le=1,
        description="Post-screen genuine off-target specificity weight (on-target, ortholog and repeat excluded)",
    )
    isoform_coverage: float = Field(
        default=0.15, ge=0, le=1, description="Protein-coding isoform coverage weight (targeting completeness)"
    )
    conservation: float = Field(
        default=0.10, ge=0, le=1, description="Cross-species ortholog conservation weight (specificity check)"
    )

    @model_validator_typed(mode="after")
    def weights_sum_to_one(self) -> "ScoringWeights":
        """Validate that scoring weights sum to approximately 1.0."""
        total = sum(getattr(self, term) for term in COMPOSITE_TERM_NAMES)
        if not (0.95 <= total <= 1.05):
            raise ValueError(f"Scoring weights must sum to 1.0, got {total:.3f}")
        return self


class DesignMode(str, Enum):
    """Design mode for siRNA/miRNA-biogenesis-aware workflows."""

    SIRNA = "sirna"  # Standard siRNA design mode
    MIRNA = "mirna"  # miRNA-biogenesis-aware design mode
    ZFN = "zfn"  # Zinc-finger nuclease pair evaluation with exhaustive off-target search


class MiRNADesignConfig(BaseModel):
    """Configuration preset for miRNA-biogenesis-aware siRNA design.

    This config encapsulates thresholds, defaults, and scoring weights
    optimized for miRNA-like processing (Drosha/Dicer recognition,
    Argonaute loading preferences, seed-based off-target analysis).
    """

    model_config = ConfigDict(extra="forbid")

    # Conservative thermodynamic thresholds for miRNA mode
    gc_min: float = Field(default=30.0, ge=0, le=100, description="Minimum GC content % for miRNA mode")
    gc_max: float = Field(default=52.0, ge=0, le=100, description="Maximum GC content % for miRNA mode")
    asymmetry_min: float = Field(default=0.65, ge=0, le=1, description="Minimum asymmetry score for Argonaute loading")
    max_homopolymer: int = Field(default=3, ge=1, description="Maximum homopolymer run length")

    # Canonical duplex format defaults
    overhang: str = Field(default="UU", description="Default overhang for miRNA mode (UU for RNA)")
    modifications: str = Field(
        default="standard_2ome", description="Default chemical modification pattern for miRNA mode"
    )

    # Off-target preset
    off_target_preset: str = Field(
        default="MIRNA_SEED_7_8", description="Off-target analysis preset (seed-based matching)"
    )

    # Scoring weights for miRNA-specific features
    scoring_weights: dict[str, float] = Field(
        default_factory=lambda: {
            "ago_start_bonus": 0.1,  # Bonus for A/U at guide position 1
            "pos1_mismatch_bonus": 0.05,  # Bonus for G:U wobble or mismatch at position 1
            "seed_clean_bonus": 0.15,  # Bonus for clean seed region (positions 2-8)
            "supp_13_16_bonus": 0.1,  # Bonus for 3' supplementary pairing potential
            "five_p_end_destabilization_bonus": 0.1,  # Bonus for destabilized 5' guide end
        },
        description="Scoring weight modifiers for miRNA-specific features",
    )

    # Pri-miRNA hairpin validation (enabled only when hairpin context is provided)
    enable_pri_hairpin_validation: bool = Field(
        default=False,
        description="Enable pri-miRNA hairpin structure validation (requires hairpin input)",
    )


class DesignParameters(BaseModel):
    """Complete configuration parameters for siRNA design workflow."""

    model_config = ConfigDict(extra="forbid")

    # Design mode selection
    design_mode: DesignMode = Field(
        default=DesignMode.SIRNA,
        description="Design mode: sirna (default), mirna (miRNA-biogenesis-aware), or zfn",
    )

    # Basic parameters
    sirna_length: int = Field(default=21, ge=19, le=23, description="siRNA duplex length in nucleotides")
    top_n: int = Field(
        default=500,
        ge=1,
        description=(
            "Number of top-ranked candidates reported in top_candidates. Screening is applied to "
            "every distinct candidate sequence regardless of this value, so it no longer gates "
            "off-target analysis, and it does not change how many candidates are enumerated."
        ),
    )

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
    snp_file: str | None = Field(default=None, description="Path to SNP VCF file for avoidance")
    # Review genome index passing / FASTA selection
    genome_index: str | None = Field(default=None, description="Path to BWA genome index for off-target search")


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

    # Sequences (accept up to the engine max; 19-23 nt is the recommended range)
    guide_sequence: str = Field(
        min_length=MIN_GUIDE_LEN,
        max_length=ENGINE_MAX_GUIDE_LEN,
        description="Guide strand (antisense, loaded into RISC). 19-23 nt recommended; up to 40 nt analyzed.",
    )
    passenger_sequence: str = Field(
        min_length=MIN_GUIDE_LEN,
        max_length=ENGINE_MAX_GUIDE_LEN,
        description="Passenger strand (sense, typically degraded). 19-23 nt recommended; up to 40 nt analyzed.",
    )

    # Basic properties
    gc_content: float = Field(ge=0, le=100, description="GC content % (optimal: 35-60%)")
    length: int = Field(
        ge=MIN_GUIDE_LEN, le=ENGINE_MAX_GUIDE_LEN, description="siRNA length in nucleotides (19-23 recommended)"
    )

    # Thermodynamic properties
    asymmetry_score: float = Field(
        ge=0, le=1, description="Thermodynamic asymmetry score for RISC loading (optimal: ≥0.65)"
    )
    duplex_stability: float | None = Field(default=None, description="Duplex formation ΔG in kcal/mol")

    # Secondary structure
    structure: str | None = Field(default=None, description="RNA secondary structure (dot-bracket notation)")
    mfe: float | None = Field(default=None, description="Minimum free energy in kcal/mol (optimal: -2 to -8)")
    paired_fraction: float = Field(default=0.0, ge=0, le=1, description="Fraction of paired bases (optimal: 0.4-0.6)")

    # Off-target analysis
    off_target_screened: bool = Field(
        default=False,
        description=(
            "True once off-target screening has run for this candidate. Distinguishes "
            "'screened, no hits' from 'never screened' -- both leave the hit counts below at 0."
        ),
    )
    off_target_count: int = Field(
        default=0,
        ge=0,
        description="Number of genuine off-target sites (on-target, ortholog and repeat hits excluded, goal: ≤3)",
    )
    # Reporting only -- nothing scores or filters on this field, and its direction depends on
    # which stage wrote it last, so do not compare values across candidates screened differently:
    #   * design time (SiRNADesigner._calculate_off_target_score) writes an internal-repeat 7-mer
    #     penalty, where HIGHER is worse;
    #   * after screening, _integrate_offtarget_results overwrites it with the MAXIMUM
    #     ``offtarget_score`` over the candidate's hits, where higher is SAFER -- 0.0 is reserved
    #     for a full-length exact match (the highest-risk hit there is). Taking the max therefore
    #     reports the candidate's *least* worrying hit, and since ``nm`` became a guide-level
    #     distance the values got wider (a clipped minus-strand partial hit reports ~76-98 where
    #     it used to report 0.0).
    off_target_penalty: float = Field(
        default=0.0,
        ge=0,
        description=(
            "Reporting only, direction depends on provenance: design-time internal-repeat penalty "
            "(higher = worse), overwritten post-screen by max offtarget_score (higher = safer, "
            "0.0 = perfect match). Use off_target_count / the hit strata to judge risk."
        ),
    )

    # Detailed transcriptome off-target metrics. _total counts every genuine off-target hit at
    # any mismatch count (so it always agrees with off_target_count); _0mm/_1mm/_2mm are
    # stratified nm<=2 SUBSETS of _total, not addends -- nm>=3 hits land in _total only.
    transcriptome_hits_total: int = Field(
        default=0, ge=0, description="Total genuine off-target transcriptome hits (any mismatch count)"
    )
    transcriptome_hits_0mm: int = Field(
        default=0, ge=0, description="Perfect-match subset of transcriptome_hits_total (0 mismatches)"
    )
    transcriptome_hits_1mm: int = Field(default=0, ge=0, description="1-mismatch subset of transcriptome_hits_total")
    transcriptome_hits_2mm: int = Field(default=0, ge=0, description="2-mismatch subset of transcriptome_hits_total")
    # The only counter that sees a PARTIAL hit whose seed paired perfectly. Because nm is a
    # guide-level distance, a clipped or gapped hit (e.g. 6S15M/NM:i:0 on the minus strand, where
    # the clip lands on guide positions 16-21 and leaves the seed intact) carries nm=6: it is
    # counted here and in transcriptome_hits_total / off_target_count, but falls in NO mismatch
    # stratum, so max_transcriptome_hits_{0,1,2}mm cannot gate it. Only
    # max_transcriptome_seed_perfect (default None) and max_off_target_count (default 3) do.
    # Unlike the _0mm/_1mm/_2mm counters this one is not split by species.
    transcriptome_hits_seed_0mm: int = Field(
        default=0, ge=0, description="Transcriptome hits with perfect seed match (positions 2-8), all species"
    )
    on_target_confirmed: bool = Field(
        default=False,
        description="Whether any hit was recognised as the query gene (see on_target_hits for the count)",
    )

    # Hit classification metrics (four-way classifier: on-target, ortholog, repeat, off-target)
    on_target_hits: int = Field(
        default=0, ge=0, description="Hits classified as on-target (query gene in query species)"
    )
    ortholog_hits: int = Field(default=0, ge=0, description="Hits classified as ortholog (same gene, other species)")
    repeat_hits: int = Field(default=0, ge=0, description="Hits classified as repeat element")
    ortholog_species: str = Field(
        default="", description="Comma-separated canonical species with at least one ortholog hit"
    )

    # Repeat detection (design-time k-mer frequency check)
    repeat_flagged: bool = Field(
        default=False, description="True if guide exceeds repeat transcript-fraction threshold at design time"
    )
    repeat_transcript_fraction: float = Field(
        default=0.0, ge=0, le=1, description="Fraction of reference transcripts containing this guide"
    )

    # miRNA off-target metrics
    mirna_hits_total: int = Field(default=0, ge=0, description="Total miRNA seed match hits")
    mirna_hits_0mm_seed: int = Field(default=0, ge=0, description="Perfect miRNA seed matches (positions 2-8)")
    mirna_hits_1mm_seed: int = Field(default=0, ge=0, description="miRNA seed matches with 1 mismatch in seed")
    mirna_hits_high_risk: int = Field(
        default=0, ge=0, description="High-risk miRNA hits (perfect seed + low offtarget_score)"
    )

    # miRNA-specific fields (populated when design_mode == "mirna")
    guide_pos1_base: str | None = Field(
        default=None, description="Nucleotide at guide position 1 (for Argonaute selection scoring)"
    )
    pos1_pairing_state: str | None = Field(
        default=None, description="Pairing state at position 1: perfect, wobble, or mismatch"
    )
    seed_class: str | None = Field(default=None, description="Seed match class: 6mer, 7mer-m8, 7mer-a1, or 8mer")
    supp_13_16_score: float | None = Field(
        default=None, ge=0, le=1, description="3' supplementary pairing score (positions 13-16)"
    )
    seed_7mer_hits: int | None = Field(
        default=None, ge=0, description="Number of 7mer seed matches in off-target analysis"
    )
    seed_8mer_hits: int | None = Field(
        default=None, ge=0, description="Number of 8mer seed matches in off-target analysis"
    )
    seed_hits_weighted: float | None = Field(
        default=None, ge=0, description="Weighted seed hits by 3' UTR abundance (if expression data provided)"
    )
    off_target_seed_risk_class: str | None = Field(
        default=None, description="Off-target risk classification: low, medium, high"
    )

    # Transcript hit metrics (how many input transcripts this guide hits)
    transcript_hit_count: int = Field(
        default=1, ge=0, description="Number of input transcripts containing this guide sequence"
    )
    transcript_hit_fraction: float = Field(
        default=1.0, ge=0, le=1, description="Fraction of input transcripts targeted by this guide (1.0 = all)"
    )

    # Post-screen sub-scores (isoform coverage and conservation)
    isoform_coverage: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description="Protein-coding isoform coverage sub-score (hit/total, inactive if no protein-coding isoforms)",
    )
    conservation_score: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description="Cross-species conservation sub-score (ortholog species hit / requested, inactive in single-species)",
    )

    # Composite scoring
    component_scores: dict[str, float] = Field(default_factory=dict, description="Individual scoring component values")
    composite_score: float = Field(ge=0, le=100, description="Overall siRNA quality score (higher is better)")
    score_asymmetry: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of asymmetry term to composite score"
    )
    score_gc_content: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of GC content term to composite score"
    )
    score_accessibility: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of accessibility term to composite score"
    )
    score_empirical: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of empirical term to composite score"
    )
    score_off_target: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of off-target term to composite score"
    )
    score_isoform_coverage: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of isoform coverage term to composite score"
    )
    score_conservation: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of conservation term to composite score"
    )
    scored_after_screening: bool = Field(
        default=False,
        description="True if composite score includes post-screen terms (off-target, isoform, conservation)",
    )
    weight_set_version: str = Field(
        default="", description="Scoring weight set version that produced composite_score (empty = not yet scored)"
    )

    # Quality flags
    class FilterStatus(str, Enum):
        """Filter status codes for quality control."""

        # PASS is a domain status label, NOT a password. Bandit B105 false positive. # nosec B105
        PASS = "PASS"
        GC_OUT_OF_RANGE = "GC_OUT_OF_RANGE"
        POLY_RUNS = "POLY_RUNS"
        EXCESS_PAIRING = "EXCESS_PAIRING"
        LOW_ASYMMETRY = "LOW_ASYMMETRY"
        LOW_EMPIRICAL_SCORE = "LOW_EMPIRICAL_SCORE"
        DIRTY_CONTROL = "DIRTY_CONTROL"
        REPEAT_ELEMENT = "REPEAT_ELEMENT"
        EXCESS_OFF_TARGETS = "EXCESS_OFF_TARGETS"
        TRANSCRIPTOME_PERFECT_MATCH = "TRANSCRIPTOME_PERFECT_MATCH"
        TRANSCRIPTOME_1MM = "TRANSCRIPTOME_1MM"
        TRANSCRIPTOME_2MM = "TRANSCRIPTOME_2MM"
        TRANSCRIPTOME_SEED_PERFECT = "TRANSCRIPTOME_SEED_PERFECT"
        MIRNA_PERFECT_SEED = "MIRNA_PERFECT_SEED"
        HIGH_RISK_MIRNA = "HIGH_RISK_MIRNA"
        TOTAL_OFFTARGETS = "TOTAL_OFFTARGETS"

    # Either True (passed) or one of the FilterStatus reasons (failed)
    passes_filters: bool | FilterStatus = Field(
        default=True, description="PASS if all filters passed, otherwise specific failure reason"
    )
    quality_issues: list[str] = Field(default_factory=list, description="List of detected quality concerns")

    # Variant-specific fields (populated when variant targeting/avoidance is enabled)
    overlapped_variants: list[dict[str, Any]] = Field(
        default_factory=list,
        description="Variants that overlap this siRNA candidate position (serialized VariantRecord)",
    )
    allele_specific: bool = Field(
        default=False, description="Whether this candidate is specific to a particular allele (ref or alt)"
    )
    targeted_alleles: list[str] = Field(
        default_factory=list,
        description="Which alleles this candidate targets (e.g., ['ref'], ['alt'], or ['ref', 'alt'])",
    )
    variant_mode: str | None = Field(
        default=None, description="Variant handling mode used: 'target', 'avoid', or 'both'"
    )

    # Optional chemical modification metadata
    guide_metadata: StrandMetadata | None = Field(
        default=None,
        description="Optional StrandMetadata for guide strand with chemical modifications",
    )
    passenger_metadata: StrandMetadata | None = Field(
        default=None,
        description="Optional StrandMetadata for passenger strand with chemical modifications",
    )

    @field_validator_typed("guide_sequence", "passenger_sequence")
    @classmethod
    def validate_nucleotide_sequence(cls, v: str) -> str:
        """Validate that sequence contains only valid nucleotides."""
        valid_bases = set("ATCGU")
        if not all(base.upper() in valid_bases for base in v):
            raise ValueError(f"Sequence contains invalid nucleotides: {v}")
        return v.upper()

    @field_validator_typed("guide_sequence", "passenger_sequence")
    @classmethod
    def warn_beyond_recommended_length(cls, v: str) -> str:
        """Warn (do not fail) when a sequence exceeds the recommended biological max.

        19-23 nt is the recommended siRNA range; longer sequences (e.g. Dicer 3'
        read-through isoforms) are still analyzed up to ENGINE_MAX_GUIDE_LEN.
        """
        if len(v) > RECOMMENDED_MAX_GUIDE_LEN:
            logger.warning(
                "Sequence length %d nt exceeds the recommended max (%d nt); analyzing anyway (engine max %d nt).",
                len(v),
                RECOMMENDED_MAX_GUIDE_LEN,
                ENGINE_MAX_GUIDE_LEN,
            )
        return v

    @field_validator_typed("passenger_sequence")
    @classmethod
    def sequences_same_length(cls, v: str, info: ValidationInfo) -> str:
        """Validate that passenger sequence is same length as guide sequence."""
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


def build_candidate_row(candidate: SiRNACandidate) -> dict[str, Any]:
    """Map one SiRNACandidate to its canonical output-row dict.

    The single source of truth for candidate CSV columns, shared by DesignResult.save_csv (the
    `sirnaforge design` path) and SiRNAWorkflow.step6_generate_reports (the `sirnaforge workflow`
    path) so the two writers cannot drift on which columns they emit. Optional attributes use a
    tolerant getattr since the workflow path feeds candidates from several producers.
    """
    cs = candidate.component_scores or {}
    mod_summary = get_modification_summary(candidate) if candidate.guide_metadata else {}
    pass_state = candidate.passes_filters
    passes_filters = pass_state.value if hasattr(pass_state, "value") else pass_state

    def _maybe_attr(name: str, default: Any = None) -> Any:
        return getattr(candidate, name, default)

    return {
        "id": candidate.id,
        "transcript_id": candidate.transcript_id,
        "position": candidate.position,
        "guide_sequence": candidate.guide_sequence,
        "passenger_sequence": candidate.passenger_sequence,
        "gc_content": candidate.gc_content,
        "asymmetry_score": candidate.asymmetry_score,
        # Thermodynamics and structure
        "structure": _maybe_attr("structure"),
        "mfe": _maybe_attr("mfe"),
        "paired_fraction": candidate.paired_fraction,
        "duplex_stability_dg": candidate.duplex_stability,
        "duplex_stability_score": cs.get("duplex_stability_score"),
        "dg_5p": cs.get("dg_5p"),
        "dg_3p": cs.get("dg_3p"),
        "delta_dg_end": cs.get("delta_dg_end"),
        "melting_temp_c": cs.get("melting_temp_c"),
        "off_target_screened": candidate.off_target_screened,
        "off_target_count": candidate.off_target_count,
        "off_target_penalty": candidate.off_target_penalty,
        # Hit classification metrics (issue #80)
        "on_target_hits": candidate.on_target_hits,
        "ortholog_hits": candidate.ortholog_hits,
        "repeat_hits": candidate.repeat_hits,
        "ortholog_species": candidate.ortholog_species,
        "repeat_flagged": candidate.repeat_flagged,
        "repeat_transcript_fraction": candidate.repeat_transcript_fraction,
        # Legacy transcriptome/miRNA hit metrics
        "transcriptome_hits_total": _maybe_attr("transcriptome_hits_total", 0),
        "transcriptome_hits_0mm": _maybe_attr("transcriptome_hits_0mm", 0),
        "transcriptome_hits_1mm": _maybe_attr("transcriptome_hits_1mm", 0),
        "transcriptome_hits_2mm": _maybe_attr("transcriptome_hits_2mm", 0),
        "transcriptome_hits_seed_0mm": _maybe_attr("transcriptome_hits_seed_0mm", 0),
        "on_target_confirmed": _maybe_attr("on_target_confirmed", False),
        "mirna_hits_total": _maybe_attr("mirna_hits_total", 0),
        "mirna_hits_0mm_seed": _maybe_attr("mirna_hits_0mm_seed", 0),
        "mirna_hits_1mm_seed": _maybe_attr("mirna_hits_1mm_seed", 0),
        "mirna_hits_high_risk": _maybe_attr("mirna_hits_high_risk", 0),
        # miRNA-specific columns (nullable)
        "guide_pos1_base": _maybe_attr("guide_pos1_base"),
        "pos1_pairing_state": _maybe_attr("pos1_pairing_state"),
        "seed_class": _maybe_attr("seed_class"),
        "supp_13_16_score": _maybe_attr("supp_13_16_score"),
        "seed_7mer_hits": _maybe_attr("seed_7mer_hits"),
        "seed_8mer_hits": _maybe_attr("seed_8mer_hits"),
        "seed_hits_weighted": _maybe_attr("seed_hits_weighted"),
        "off_target_seed_risk_class": _maybe_attr("off_target_seed_risk_class"),
        # Transcript hit metrics
        "transcript_hit_count": candidate.transcript_hit_count,
        "transcript_hit_fraction": candidate.transcript_hit_fraction,
        # Post-screen sub-scores
        "isoform_coverage": candidate.isoform_coverage,
        "conservation_score": candidate.conservation_score,
        # Composite scoring
        "composite_score": candidate.composite_score,
        "score_asymmetry": candidate.score_asymmetry,
        "score_gc_content": candidate.score_gc_content,
        "score_accessibility": candidate.score_accessibility,
        "score_empirical": candidate.score_empirical,
        "score_off_target": candidate.score_off_target,
        "score_isoform_coverage": candidate.score_isoform_coverage,
        "score_conservation": candidate.score_conservation,
        "scored_after_screening": candidate.scored_after_screening,
        "weight_set_version": candidate.weight_set_version,
        "passes_filters": passes_filters,
        # Chemical modifications
        "guide_overhang": mod_summary.get("guide_overhang", ""),
        "guide_modifications": mod_summary.get("guide_modifications", ""),
        "passenger_overhang": mod_summary.get("passenger_overhang", ""),
        "passenger_modifications": mod_summary.get("passenger_modifications", ""),
        # Variant-aware annotations
        "variant_mode": _maybe_attr("variant_mode"),
        "allele_specific": _maybe_attr("allele_specific", False),
        "targeted_alleles": json.dumps(_maybe_attr("targeted_alleles", [])),
        "overlapped_variants": json.dumps(_maybe_attr("overlapped_variants", [])),
    }


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
    rejected_candidates: list[SiRNACandidate] = Field(
        default_factory=list,
        description=("Candidates discarded during initial filtering (used for dirty controls and auditing)"),
    )

    @check_types_typed
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
        df_data = [build_candidate_row(candidate) for candidate in self.candidates]
        df = pd.DataFrame(df_data)

        # Convert nullable integer columns to pandas Int64 dtype for proper None handling
        nullable_int_cols = ["seed_7mer_hits", "seed_8mer_hits"]
        for col in nullable_int_cols:
            if col in df.columns:
                df[col] = df[col].astype("Int64")

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
