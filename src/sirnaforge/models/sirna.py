"""Pydantic models for siRNA design data structures."""

import json
from enum import Enum
from typing import Any, ClassVar

import pandas as pd
from pandera.typing import DataFrame
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo

from sirnaforge.core.repeat_detection import DEFAULT_REPEAT_TRANSCRIPT_FRACTION
from sirnaforge.models.modifications import StrandMetadata, StrandRole
from sirnaforge.models.policy import DECLARED_FILTER_IDS, FilterAction, FilterEvaluation
from sirnaforge.models.schemas import SiRNACandidateSchema
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import get_modification_summary
from sirnaforge.utils.typed_decorators import check_types_typed, field_validator_typed, model_validator_typed

logger = get_logger(__name__)

#: What a filter that did not run writes, rather than a blank a reader mistakes for a verdict.
_NOT_EVALUATED = FilterEvaluation.NOT_EVALUATED.value

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
# The max is 0.6, not 0.7: issue #96 deleted the G/C-at-guide-position-1 clause
# (it contradicted the A/U-at-position-1 biogenesis rule), leaving {0.4, 0.5, 0.6}.
EMPIRICAL_SCORE_MIN = 0.4
EMPIRICAL_SCORE_MAX = 0.6
# Default 0.4 == EMPIRICAL_SCORE_MIN, i.e. the gate is inert by default.
# It was 0.5, which after the position-1 clause was deleted rejected every candidate carrying C at
# guide position 19 -- a 12.9% cut of the passing pool on a reference TP53 run (264 -> 230). Measured
# against 2,816 siRNAs with knockdown data that rule is BACKWARDS: C at guide position 19 associates
# with MORE knockdown (mean inhibition 0.735 vs 0.668, p = 2.4e-13), and the A/U preference the
# rubric rewards there associates with LESS (0.650 vs 0.719, p = 3.7e-17).
# The cause is a strand error. Reynolds' criteria are numbered on the SENSE strand, so they land on
# the guide 5' end, not the guide 3' end: A/U count at guide positions 1-5 tracks efficacy at
# rho +0.378 (monotonic, 0.455 -> 0.849), while at positions 17-21 it is null (rho -0.017, p = 0.37).
# Raise this only once the rubric is scored at the end of the guide that carries the signal.
DEFAULT_MIN_EMPIRICAL_SCORE = 0.4

# Every scored term name, in reporting order. This is a *union* over the named weight vectors
# below, used for column ordering and for iterating contributions -- it is NOT a weight vector
# and nothing validates against it. Each vector validates against its own TERM_NAMES, because
# the three vectors score different terms (3, 4 and 7 of them) and a single global tuple used to
# make a missing term look like a licence to renormalise -- which is no longer permitted anywhere.
# `pos1_mismatch` is absent: issue #102 removed it from the only vector that held it, so no vector
# scores it and `score_pos1_mismatch` -- a *contribution* column -- is now always null. The pairing
# state itself is still reported, on `guide_pos1_base` and `pos1_pairing_state`.
COMPOSITE_TERM_NAMES = (
    "off_target",
    "target_accessibility",
    "asymmetry",
    "gc_content",
    "ago_start",
    "supp_13_16",
)

# Hand-authored weight vectors must already sum to 1.0. The tolerance covers float
# representation of two-decimal literals only -- it is not room to be approximately normalised.
WEIGHT_SUM_TOLERANCE = 1e-9

# RNAplfold parameters for target-site accessibility. W=150/L=100 sits near the plateau of the
# benchmark correlation (rho +0.249 at W=40 rising to +0.269 at W=240, n=2,779 siRNAs with measured
# knockdown) without paying W=240's cost. Configurable rather than constant because they move a
# site's accessibility percentile substantially.
DEFAULT_PLFOLD_WINDOW = 150
DEFAULT_PLFOLD_MAX_BP_SPAN = 100

# Log-probability floor for normalising P(seed-end 8-mer unpaired) into a [0, 1] feature.
# Derived from the benchmark distribution of log10 P at W=150/L=100 (p1 = -4.93, p50 = -1.72,
# p99 = -0.18), so -5.0 captures ~99% of observed sites. Deliberately a fixed floor, not a
# per-transcript one: self-calibrating normalisation would make scores incomparable between
# targets and between runs.
DEFAULT_ACCESSIBILITY_LOG_FLOOR = -5.0


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
    max_repeat_transcript_fraction: float = Field(
        default=DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
        gt=0,
        le=1,
        description="Ceiling on the fraction of reference transcripts containing the guide; fails REPEAT_ELEMENT",
    )
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

    # Empirical (simplified Reynolds) design-rule filter. Since issue #96 the empirical score is
    # gate-only: it is computed and reported, and this threshold is the only thing that reads it.
    min_empirical_score: float = Field(
        default=DEFAULT_MIN_EMPIRICAL_SCORE,
        ge=EMPIRICAL_SCORE_MIN,
        le=EMPIRICAL_SCORE_MAX,
        description=(
            "Minimum empirical design-rule score. Applied to the 'empirical' component score, "
            f"whose attainable range is {EMPIRICAL_SCORE_MIN}-{EMPIRICAL_SCORE_MAX}. Defaults to "
            f"{EMPIRICAL_SCORE_MIN}, which rejects nothing: the rubric's positional rules are applied to "
            "the guide 3' end, where measured knockdown data shows no signal, so gating on them is not "
            "justified. Gate only -- 'empirical' is not a scored term."
        ),
    )

    # Protein-coding isoform coverage gate. Defaults to None (off) so default behaviour is
    # unchanged: coverage is computed and reported on every candidate either way, and no ceiling
    # has been calibrated against truth data. Only readable post-screening, since the numerator
    # comes from the guide -> source-transcript map built during screening.
    min_isoform_coverage: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description=(
            "Minimum protein-coding isoform coverage fraction gating LOW_ISOFORM_COVERAGE "
            "(None = no gate, the default). Applied post-screening to SiRNACandidate.isoform_coverage; "
            "a candidate whose coverage could not be computed is never failed by it."
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
    # 15 is calibrated against one internal 94-design reference set: it is the lowest cap at which
    # the gate enriches for expert-chosen guides rather than depleting them. At 3 the gate was
    # depleted (p = 0.89) and at 10 it carried no information (p = 0.51). Single-target calibration.
    max_off_target_count: int | None = Field(
        default=15,
        ge=0,
        description="Maximum genuine off-target sites (default 15). Excludes on-target, ortholog and repeat hits.",
    )

    # Transcriptome off-target thresholds
    # The three transcriptome thresholds below count GENUINE off-target hits only: on-target,
    # ortholog and repeat-mediated hits are classified out before they reach these counters.
    max_transcriptome_hits_0mm: int | None = Field(
        default=1,
        ge=0,
        description=(
            "Maximum perfect-match genuine off-target transcriptome hits (excludes on-target "
            "isoforms). HUMAN-STRATIFIED: read against hits whose species is human or unlabelled, so it does not match the identically named all-species column"
        ),
    )
    max_transcriptome_hits_1mm: int | None = Field(
        default=10,
        ge=0,
        description=(
            "Maximum 1-mismatch genuine off-target hits (typical: 5-10, None = no limit). HUMAN-STRATIFIED: read against hits whose species is human or unlabelled, so it does not match the identically named all-species column"
        ),
    )
    max_transcriptome_hits_2mm: int | None = Field(
        default=50,
        ge=0,
        description=(
            "Maximum 2-mismatch genuine off-target hits (typical: 20-50, None = no limit). HUMAN-STRATIFIED: read against hits whose species is human or unlabelled, so it does not match the identically named all-species column"
        ),
    )
    # Read against transcriptome_hits_seed_0mm. This is the only *targeted* gate on a partial
    # (clipped or gapped) hit whose seed paired perfectly: such a hit has a guide-level nm > 2, so
    # it lands in no mismatch stratum and max_transcriptome_hits_{0,1,2}mm never see it. The only
    # other thing that counts it is the blunt max_off_target_count ceiling (default 15).
    # Unlike the three mismatch thresholds this one counts ALL screened species, matching the
    # reported transcriptome_hits_seed_0mm column exactly. Defaults to None (off) because a safe
    # ceiling has not been calibrated against truth data; set it to opt in.
    max_transcriptome_seed_perfect: int | None = Field(
        default=None, ge=0, description="Maximum transcriptome hits with perfect seed (positions 2-8, None = no limit)"
    )

    # miRNA off-target thresholds
    max_mirna_perfect_seed: int | None = Field(
        default=0,
        ge=0,
        description=(
            "Maximum perfect miRNA seed matches (typical: 3-5, None = no limit). HUMAN-STRATIFIED: read against hits labelled human ONLY -- unlike the transcriptome gates, an unlabelled hit is not counted -- so it does not match the identically named all-species column"
        ),
    )
    # None, not 10: the policy resolves this filter to `off`, and a threshold beside an off action
    # reads as a limit that applies. The manifest published `action: off` and `threshold: 10` side by
    # side, which is two different answers to "is 11 rejected?".
    max_mirna_1mm_seed: int | None = Field(
        default=None,
        ge=0,
        description=(
            "Maximum 1-mismatch miRNA seed hits (typical: 10-20, None = no limit). NOT ENFORCED in "
            "0.7.1: no gate reads it, so the run policy resolves this filter to off rather than to passing"
        ),
    )
    fail_on_high_risk_mirna: bool = Field(
        default=True,
        description=(
            "Fail if high-risk miRNA hits detected (perfect seed + offtarget_score < 5.0). HUMAN-STRATIFIED: read against hits labelled human ONLY -- unlike the transcriptome gates, an unlabelled hit is not counted -- so it does not match the identically named all-species column"
        ),
    )

    # Combined off-target threshold
    max_total_offtarget_hits: int | None = Field(
        default=None,
        ge=0,
        description=(
            "Maximum total off-target hits (transcriptome + miRNA, None = no limit). HUMAN-STRATIFIED: read against human-or-unlabelled transcriptome hits plus human-labelled-only miRNA hits, so it does not match the identically named all-species column"
        ),
    )


class WeightVector(BaseModel):
    """A named, hand-authored weight vector over one explicit term set.

    Nothing in siRNAforge scales, renormalises or divides these numbers at runtime. A vector that
    does not already sum to 1.0, or that carries no name, is a construction-time error -- which is
    what makes the weight recorded in the manifest the weight that actually applied.

    Subclasses declare ``VECTOR_NAME`` and ``TERM_NAMES``; the validator below reads both, so a
    subclass cannot forget either and still be usable.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    VECTOR_NAME: ClassVar[str] = ""
    TERM_NAMES: ClassVar[tuple[str, ...]] = ()

    @model_validator_typed(mode="after")
    def named_and_normalised(self) -> "WeightVector":
        """Reject an unnamed, empty or mis-summed vector at construction."""
        cls = type(self)
        if not cls.VECTOR_NAME:
            raise ValueError(f"{cls.__name__} declares no VECTOR_NAME; every weight vector must be named")
        if not cls.TERM_NAMES:
            raise ValueError(f"{cls.__name__} declares no TERM_NAMES; a vector must name the terms it scores")
        total = sum(getattr(self, term) for term in cls.TERM_NAMES)
        if abs(total - 1.0) > WEIGHT_SUM_TOLERANCE:
            raise ValueError(
                f"Weight vector '{cls.VECTOR_NAME}' must sum to exactly 1.0, got {total:.6f}. "
                "Weights are never renormalised at runtime, so a vector that does not sum to 1.0 "
                "would silently rescale every score it produced."
            )
        return self

    @property
    def name(self) -> str:
        """Vector name recorded on every row and in the run manifest."""
        return type(self).VECTOR_NAME

    @property
    def terms(self) -> tuple[str, ...]:
        """Exactly the terms this vector scores; the scorer requires all of them."""
        return type(self).TERM_NAMES

    def as_mapping(self) -> dict[str, float]:
        """Term name -> weight, in declaration order."""
        return {term: float(getattr(self, term)) for term in type(self).TERM_NAMES}


class DesignWeights(WeightVector):
    """``design_v4``: the design-stage vector, over the three terms computable before screening.

    ``design_score`` is NOT comparable with ``composite_score``: it is a different vector over a
    different term set. It is **not** systematically the larger of the two, and issue #102 deleted
    that claim -- the two vectors weight their shared terms differently, so which is higher depends
    on the candidate. All features at 0.5 with ``off_target = 1.0`` gives ``design_v4`` 50.0 against
    ``postscreen_sirna_v4``'s 62.5. Compare neither to the other; compare each within its stage.
    """

    VECTOR_NAME: ClassVar[str] = "design_v4"
    TERM_NAMES: ClassVar[tuple[str, ...]] = ("target_accessibility", "asymmetry", "gc_content")

    # asymmetry takes the top slot rather than target_accessibility. On 900 benchmark siRNAs with
    # measured knockdown the two are statistically indistinguishable (Spearman rho +0.273 vs +0.267),
    # so this is not a large evidential gap -- but asymmetry additionally correlates rho +0.53 with
    # A/U content at guide positions 1-5 (rho +0.438 on that 900-siRNA folding subsample; +0.378 on
    # the full 2,816, where #102 measured guide position 1 alone higher still at +0.415, so A/U(1-5)
    # is not the strongest single feature). Accessibility explaining ~7% of rank variance did not
    # justify 0.40.
    target_accessibility: float = Field(
        default=0.35,
        ge=0,
        le=1,
        description=(
            "Target-site accessibility weight: log-scaled P(the 8 nt of the mRNA target site "
            "pairing guide positions 1-8 are unpaired), from RNAplfold on the transcript"
        ),
    )
    asymmetry: float = Field(
        default=0.40, ge=0, le=1, description="Thermodynamic asymmetry weight (guide strand selection)"
    )
    gc_content: float = Field(
        default=0.25, ge=0, le=1, description="GC content optimization weight (stability balance)"
    )


class PostScreenSiRNAWeights(WeightVector):
    """``postscreen_sirna_v4``: ``design_v4``'s terms plus ``off_target``, one extra term.

    ``off_target`` holds its nominal 0.25 while the scored budget shrank from six terms to four, so
    its share of the scored budget rises from 0.25/0.60 to 0.25/1.00.

    Issue #102 corrected the justification that used to sit here. It read "the term measured 2.24x
    its nominal share of composite variance on the reference run"; that 2.24x is from an internal
    run under weight set 2.0.0 and has never been re-derived, and the frozen public baseline
    measures the opposite -- **0.39x nominal, the least influential of the four scored terms** on a
    human-only screen, its contribution compressed against the 25.0 ceiling. Screening one more
    species doubles the share to 0.19, so that number describes the screening scope, not the term.
    Whether 0.25 is the right weight is therefore open, and the honest statement is that it is a
    declared prior awaiting a full-reference run.
    """

    VECTOR_NAME: ClassVar[str] = "postscreen_sirna_v4"
    TERM_NAMES: ClassVar[tuple[str, ...]] = ("off_target", "target_accessibility", "asymmetry", "gc_content")

    off_target: float = Field(
        default=0.25,
        ge=0,
        le=1,
        description="Post-screen genuine off-target specificity weight (on-target, ortholog and repeat excluded)",
    )
    target_accessibility: float = Field(
        default=0.30, ge=0, le=1, description="Target-site accessibility weight (RNAplfold opening probability)"
    )
    asymmetry: float = Field(
        default=0.25, ge=0, le=1, description="Thermodynamic asymmetry weight (guide strand selection)"
    )
    gc_content: float = Field(
        default=0.20, ge=0, le=1, description="GC content optimization weight (stability balance)"
    )


class PostScreenMiRNAWeights(WeightVector):
    """``postscreen_mirna_v4``: the miRNA-biogenesis-aware post-screen vector, 6 declared terms.

    Hand-authored, not derived from ``postscreen_sirna_v4``: deriving it by scaling would be the
    1.25 divisor in a new costume. The biogenesis terms replaced an undeclared bonus that was folded
    into the score and then divided out of it, so ``--design-mode mirna`` used to move every weight
    by a factor absent from the manifest.

    Issue #102 removed ``pos1_mismatch``, which held 0.05 here while being **exactly constant at
    0.0** -- one distinct value over all 13,415 scored candidates of the public baseline, 0.000% of
    variance. It is constant by construction, not by accident: the passenger is the exact reverse
    complement of the guide, so guide position 1 always forms a Watson-Crick pair. A term that ranks
    nothing must not consume weight. ``biogenesis_features`` still computes it and the miRNA design
    path still writes it to ``component_scores``, and the pairing state stays on the row as
    ``guide_pos1_base`` and ``pos1_pairing_state``, so a run stays auditable and the term can return
    if the designer ever builds deliberately mismatched passengers. Its *contribution* column
    ``score_pos1_mismatch`` is now always null, because no vector scores it; see
    ``models/scoring_profile.py``.

    The 0.05 it released was **not** reassigned by judgement. The four shared terms are now exactly
    ``0.80 x postscreen_sirna_v4`` -- 0.20 / 0.24 / 0.20 / 0.16 -- which is the proportional-scaling
    rule this docstring already declared, now reachable at 2 dp without rounding and therefore
    without the off_target/asymmetry tie-break the rounding used to force (the source vector ties
    those two at 0.25 apiece; this one ties them at 0.20). ``ago_start`` and ``supp_13_16`` keep
    their declared 0.10 each, so the biogenesis budget is 0.20 and the shared budget 0.80.

    Declared expert priors, reviewed and accepted, and **not fitted to any dataset**; a run scored
    under different numbers is not comparable, so bump SCORING_WEIGHT_SET_VERSION if they change.
    """

    VECTOR_NAME: ClassVar[str] = "postscreen_mirna_v4"
    TERM_NAMES: ClassVar[tuple[str, ...]] = (
        "off_target",
        "target_accessibility",
        "asymmetry",
        "gc_content",
        "ago_start",
        "supp_13_16",
    )

    off_target: float = Field(default=0.20, ge=0, le=1, description="Post-screen genuine off-target specificity weight")
    target_accessibility: float = Field(
        default=0.24, ge=0, le=1, description="Target-site accessibility weight (RNAplfold opening probability)"
    )
    asymmetry: float = Field(default=0.20, ge=0, le=1, description="Thermodynamic asymmetry weight")
    gc_content: float = Field(default=0.16, ge=0, le=1, description="GC content optimization weight")
    ago_start: float = Field(
        default=0.10, ge=0, le=1, description="Argonaute loading preference weight (A/U at guide position 1)"
    )
    supp_13_16: float = Field(
        default=0.10, ge=0, le=1, description="3' supplementary pairing weight (guide positions 13-16)"
    )


class ScoringWeights(BaseModel):
    """The named weight vectors this run may score with, one per (stage, design mode).

    A vector is chosen, never combined: ``vector_for`` returns exactly one, its name is stamped on
    the candidate and written to the manifest, and the scorer requires every term it declares. There
    is deliberately no flat weight attribute here -- a single flat vector is what let the scorer
    renormalise over "whichever terms happened to be populated".
    """

    model_config = ConfigDict(extra="forbid")

    design: DesignWeights = Field(
        default_factory=DesignWeights, description="Design-stage vector (design_v4), scores design_score"
    )
    postscreen_sirna: PostScreenSiRNAWeights = Field(
        default_factory=PostScreenSiRNAWeights,
        description="Post-screen siRNA vector (postscreen_sirna_v4), scores composite_score",
    )
    postscreen_mirna: PostScreenMiRNAWeights = Field(
        default_factory=PostScreenMiRNAWeights,
        description="Post-screen miRNA-biogenesis-aware vector (postscreen_mirna_v4), scores composite_score",
    )

    def vector_for(self, *, post_screen: bool, design_mode: "DesignMode | None" = None) -> WeightVector:
        """Return the one vector that applies, by stage and design mode.

        The design stage uses a single vector in every mode: the miRNA biogenesis terms need
        screening-independent inputs but are only meaningful against the post-screen term set, and
        giving miRNA mode its own design vector would make two design_scores incomparable.
        """
        if not post_screen:
            return self.design
        if design_mode == DesignMode.MIRNA:
            return self.postscreen_mirna
        return self.postscreen_sirna

    def all_vectors(self) -> tuple[WeightVector, ...]:
        """Every declared vector, for manifest recording and whole-config validation."""
        return (self.design, self.postscreen_sirna, self.postscreen_mirna)

    def as_manifest(self) -> dict[str, dict[str, float]]:
        """Vector name -> its weights, so a row's ``weight_vector`` resolves to the numbers used."""
        return {vector.name: vector.as_mapping() for vector in self.all_vectors()}


class TargetAccessibilityConfig(BaseModel):
    """RNAplfold settings for the target-site accessibility scoring term.

    All three fields change the numeric scale of `target_accessibility`, so a run that
    alters them is not comparable with a default run. They are recorded in the run
    manifest via DesignParameters.
    """

    model_config = ConfigDict(extra="forbid")

    window_size: int = Field(
        default=DEFAULT_PLFOLD_WINDOW,
        ge=20,
        le=1000,
        description=(
            "RNAplfold averaging window W in nucleotides (default 150). Larger is mildly better "
            "on the benchmark (rho +0.249 at W=40 to +0.269 at W=240) and W=150-240 is the plateau."
        ),
    )
    max_bp_span: int = Field(
        default=DEFAULT_PLFOLD_MAX_BP_SPAN,
        ge=10,
        le=1000,
        description="RNAplfold maximum base-pair span L in nucleotides (default 100); must not exceed window_size",
    )
    log_floor: float = Field(
        default=DEFAULT_ACCESSIBILITY_LOG_FLOOR,
        lt=0,
        description=(
            "log10 probability treated as zero accessibility when normalising the term into [0, 1] "
            "(default -5.0, which captures ~99% of benchmark sites). A fixed floor on purpose: a "
            "per-transcript floor would make composite scores incomparable between targets."
        ),
    )

    @model_validator_typed(mode="after")
    def span_within_window(self) -> "TargetAccessibilityConfig":
        """RNAplfold silently clamps L to W, so reject the combination instead of misreporting it."""
        if self.max_bp_span > self.window_size:
            raise ValueError(f"max_bp_span ({self.max_bp_span}) must not exceed window_size ({self.window_size})")
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

    # The miRNA-specific scoring weights live in PostScreenMiRNAWeights, not here. This config
    # carries thresholds and format defaults only. `seed_clean_bonus` and
    # `five_p_end_destabilization_bonus` used to be declared here and read nowhere in src/ --
    # 0.25 of declared bonus weight that did nothing -- and were deleted in issue #96.

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
    # None means report every ranked candidate. Kept as an optional int for backwards compatibility
    # with callers that pass an explicit ceiling; every consumer slices with ``[: top_n]``, and
    # ``list[:None]`` is a no-op, so None needs no special-casing downstream.
    top_n: int | None = Field(
        default=None,
        ge=1,
        description=(
            "Number of top-ranked candidates reported in top_candidates (None = no limit, the "
            "default). Screening is applied to every distinct candidate sequence regardless of this "
            "value, so it does not gate off-target analysis or change how many candidates are "
            "enumerated -- it only truncates what is reported."
        ),
    )

    # Filtering criteria
    filters: FilterCriteria = Field(default_factory=FilterCriteria, description="Quality control filters")
    # Read by SiRNAWorkflow when gating on off-target screening results. Without this field the
    # getattr in _check_offtarget_filters could never resolve, so the criteria were unreachable.
    offtarget_filters: OffTargetFilterCriteria = Field(
        default_factory=OffTargetFilterCriteria, description="Off-target rejection thresholds"
    )

    # Scoring weights
    scoring: ScoringWeights = Field(default_factory=ScoringWeights, description="Component score weights")

    # Target-site accessibility (RNAplfold) settings for the target_accessibility term
    target_accessibility: TargetAccessibilityConfig = Field(
        default_factory=TargetAccessibilityConfig,
        description="RNAplfold window and normalisation settings for the target_accessibility term",
    )

    # Optional analysis parameters
    # Defaults to False because nothing reads it. No gate, scorer or enumerator consults this flag,
    # so a True default put "avoid_snps: true" in every manifest ever written for a run that did no
    # such thing -- a claim about the run that was never true. It stays in the model as the switch a
    # variant-aware enumerator would read; until then, off is the only honest default.
    avoid_snps: bool = Field(
        default=False,
        description="Exclude regions with known SNPs. NOT IMPLEMENTED in 0.7.1: no code reads this flag",
    )
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

    # Guide self-structure. Reported, and the EXCESS_PAIRING gate input; not a scoring term.
    # mfe == 0.0 with an all-dots structure is the open chain -- the physical floor of the MFE,
    # not a failed fold and not a sentinel.
    structure: str | None = Field(default=None, description="Guide-strand secondary structure (dot-bracket notation)")
    mfe: float | None = Field(
        default=None, description="Guide-strand minimum free energy in kcal/mol (0.0 = open chain, the floor)"
    )
    paired_fraction: float = Field(
        default=0.0,
        ge=0,
        le=1,
        description=(
            "Fraction of guide bases paired in its own MFE structure. Quantised to 2k/length by the "
            "dot-bracket. Reported and gates EXCESS_PAIRING; no longer a composite scoring term."
        ),
    )

    # Target-site accessibility (RNAplfold on the transcript, all windows anchored on the site's
    # 3' end, which is the end the guide seed pairs). None means it could not be computed, and
    # since no weight is ever redistributed the candidate then carries no score at all -- the
    # value is never substituted with a default.
    target_accessibility_p: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description="P(the 8 nt of the target site pairing guide positions 1-8 are unpaired). Scored term input.",
    )
    target_accessibility_p_17mer: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description="P(the 3'-most 17 nt of the target site are unpaired). Reported only, never scored.",
    )
    target_accessibility_p_site: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description="P(the entire target site is unpaired). Reported only, never scored.",
    )

    # Off-target analysis
    off_target_screened: bool = Field(
        default=False,
        description=(
            "True once this candidate's screen produced usable evidence: it reached the aligner "
            "and its query species was aligned. Distinguishes 'screened, no hits' from 'never "
            "screened' -- both leave the hit counts below at 0. It does NOT promise every "
            "requested species was aligned: an alignment for some OTHER species can fail while "
            "this stays True, so the counts below are always a lower bound, never a guaranteed "
            "total (offtarget_summary.filtering_stats.unscreened_species names the shortfall). "
            "False means the screen yielded no usable evidence for this candidate at all -- never "
            "run, never submitted, or its query species never aligned -- in which case the counts "
            "are unknown rather than zero, and any hits that were found are still reported."
        ),
    )
    off_target_count: int = Field(
        default=0,
        ge=0,
        description=(
            "Sites counted as liabilities: on-target, ortholog and repeat hits excluded; hits "
            "whose class could not be decided INCLUDED (see undetermined_hits), so a missing "
            "reference cannot loosen the screen. The enforced ceiling is "
            "OffTargetFilterCriteria.max_off_target_count (default 15)"
        ),
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
            "SUPERSEDED by score_off_target; kept for continuity, do not gate or rank on it. "
            "Reporting only, direction depends on provenance: design-time internal-repeat penalty "
            "(higher = worse), overwritten post-screen by max offtarget_score (higher = safer, "
            "0.0 = perfect match). Measured on a 40,079-candidate run it barely tracks the quantity "
            "a reader assumes it means: Pearson +0.20 against off_target_count, with 43% of rows "
            "pinned at its 132 ceiling. Use off_target_count / the hit strata to judge risk."
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
    # max_transcriptome_seed_perfect (default None) and max_off_target_count (default 15) do.
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
    undetermined_hits: int = Field(
        default=0,
        ge=0,
        description=(
            "Hits whose class could not be decided because the hit species has no transcript "
            "index. Included in off_target_count and reported here so the unqualified share is visible"
        ),
    )
    ortholog_species: str = Field(
        default="", description="Comma-separated canonical species with at least one ortholog hit"
    )
    # #103's join key. The aligner is handed one FASTA record per DISTINCT guide sequence, so a hit
    # row's qname is the representative's id, not this candidate's: joining hits on `id` silently
    # attributes one guide's evidence to one of its median-6 (max-34 on the frozen baseline)
    # candidates. Joining on the guide sequence instead works only while the spellings are
    # byte-identical, which U-vs-T guides are not.
    screen_query_id: str | None = Field(
        default=None,
        description=(
            "Query id this candidate was screened under (qname on its hit rows); None when the candidate "
            "was not submitted to the aligner's input FASTA. A candidate that was submitted but whose "
            "alignment never ran still carries the key — off_target_screened=False is the signal there."
        ),
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

    # Reported, non-scoring evidence. Both left the composite in issue #96 -- they are still
    # computed on every candidate and still written to every row, but no weight reads them.
    # isoform_coverage additionally feeds the optional FilterCriteria.min_isoform_coverage gate.
    isoform_coverage: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description=(
            "Protein-coding isoform coverage (hit/total, None if no protein-coding isoforms). "
            "Reported and the optional LOW_ISOFORM_COVERAGE gate input; not a scoring term."
        ),
    )
    conservation_score: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description=(
            "Cross-species conservation fraction (ortholog species hit / requested, None in "
            "single-species runs). Reported only; not a scoring term."
        ),
    )

    # Composite scoring. Two scores on two declared vectors, deliberately not one field:
    #   design_score    -- design_v4, 3 terms, available before screening
    #   composite_score -- postscreen_{sirna,mirna}_v4, available only after screening
    # They are NOT comparable: different term sets, and the two vectors also weight their *shared*
    # terms differently, so neither score is systematically the larger. Issue #102 deleted the claim
    # that design_score is the more optimistic number -- 0.5 features with off_target = 1.0 give
    # design_v4 50.0 against postscreen_sirna_v4's 62.5.
    component_scores: dict[str, float] = Field(default_factory=dict, description="Individual scoring component values")
    design_score: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description=(
            "Design-stage score on the design_v4 vector (target_accessibility, asymmetry, "
            "gc_content). None when a term could not be computed. Not comparable with composite_score."
        ),
    )
    composite_score: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description=(
            "Post-screen siRNA quality score, higher is better. None until off-target screening has "
            "produced usable evidence for this candidate -- it is not computable before that."
        ),
    )
    score_off_target: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of off-target term to composite score"
    )
    score_target_accessibility: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of target-site accessibility term to the score"
    )
    score_asymmetry: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of asymmetry term to the score"
    )
    score_gc_content: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of GC content term to the score"
    )
    score_ago_start: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of the Argonaute-start term (miRNA mode only)"
    )
    score_pos1_mismatch: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of the position-1 pairing term (miRNA mode only)"
    )
    score_supp_13_16: float | None = Field(
        default=None, ge=0, le=100, description="Contribution of the 3' supplementary pairing term (miRNA mode only)"
    )
    scored_after_screening: bool = Field(
        default=False,
        description="True if composite_score was computed post-screening (the only stage that computes it)",
    )
    weight_set_version: str = Field(
        default="", description="Scoring weight set version that produced the score (empty = not yet scored)"
    )
    weight_vector: str = Field(
        default="",
        description=(
            "Name of the hand-authored weight vector that produced the score (design_v4, "
            "postscreen_sirna_v4 or postscreen_mirna_v4), so a row traces to the exact weights used"
        ),
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
        LOW_ISOFORM_COVERAGE = "LOW_ISOFORM_COVERAGE"
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

    # One verdict per declared filter, so a gate's outcome is not lost when another gate rejects the
    # same candidate. `passes_filters` holds a SINGLE label and each gate overwrites it, so a label
    # count was never a rejection count: on a 40,079-candidate run 14,266 candidates carrying an
    # off-target label had already failed a design gate, and the design gates' own counts read 3x
    # lower than the number of candidates they actually rejected.
    #
    # These also carry the value each gate compared. That is not redundancy: six gates read
    # human-stratified counters that exist only as locals at the gate site, and `max_poly_runs` reads
    # a length that is not a candidate field at all, so for those filters the row cannot otherwise be
    # used to check the verdict. `evidence_exported` on the descriptor says which ones those are.
    filter_verdicts: dict[str, str] = Field(
        default_factory=dict,
        description="filter_id -> FilterEvaluation value (pass/fail/unknown/not_evaluated)",
    )
    filter_observed: dict[str, float | None] = Field(
        default_factory=dict,
        description="filter_id -> the value the gate compared against its threshold",
    )

    # Either True (passed) or one of the FilterStatus reasons (failed). Reflects FAIL-action filters
    # only: a warn-action filter records its verdict above and leaves this alone, which is what makes
    # "record it, do not reject" representable.
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

    def record_filter_verdict(
        self,
        filter_id: str,
        *,
        observed: float | None,
        passed: bool,
        action: FilterAction,
        status: "SiRNACandidate.FilterStatus | None" = None,
    ) -> None:
        """Record one gate's own outcome, and reject only if that gate's action says to.

        This is the split that makes ``FilterAction.WARN`` mean anything. Before it, a gate had
        exactly one way to express a failure -- overwrite ``passes_filters`` -- so "record it, do not
        reject" was unrepresentable and ``WARN`` sat in the vocabulary unapplied.

        ``passes_filters`` keeps first-failure-wins: it is a single label and the first FAIL-action
        gate to reject a candidate owns it. That is now only a display choice, because every gate's
        verdict survives in ``filter_verdicts`` regardless of which label won.

        Args:
            filter_id: The declared filter this verdict belongs to.
            observed: The value compared against the threshold. ``None`` records ``unknown``: the gate
                was in force but its evidence was unavailable, which is not the same as passing.
            passed: Whether the comparison succeeded.
            action: The resolved action for this filter. Only ``FAIL`` may reject.
            status: The label to put on ``passes_filters`` when this gate rejects. Required for a
                FAIL-action gate; a warn-action gate never needs one.
        """
        self.filter_observed[filter_id] = observed
        if observed is None:
            self.filter_verdicts[filter_id] = FilterEvaluation.UNKNOWN.value
            return

        self.filter_verdicts[filter_id] = FilterEvaluation.PASS.value if passed else FilterEvaluation.FAIL.value
        if passed or action is not FilterAction.FAIL:
            return
        if status is not None and self.passes_filters is True:
            self.passes_filters = status

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


def ranking_score(candidate: SiRNACandidate) -> float:
    """The number a candidate is currently ranked by.

    ``composite_score`` when screening produced it, otherwise the design-stage ``design_score``.
    The two are not comparable, so callers that may hold a mixture must keep them apart (see
    SiRNAWorkflow._apply_post_screen_ranking); this helper only answers "which number does this
    candidate have". An unscored candidate sorts last at 0.0 rather than raising: rejected
    candidates and pre-designed guides legitimately have neither.
    """
    if candidate.composite_score is not None:
        return candidate.composite_score
    return candidate.design_score if candidate.design_score is not None else 0.0


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
    verdicts = candidate.filter_verdicts or {}
    observed = candidate.filter_observed or {}

    def _maybe_attr(name: str, default: Any = None) -> Any:
        return getattr(candidate, name, default)

    return {
        "id": candidate.id,
        # The join key for row-level off-target evidence: equals qname in the hit table.
        "screen_query_id": _maybe_attr("screen_query_id"),
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
        # Target-site accessibility: the first is scored, the other two are reported only
        "target_accessibility_p": _maybe_attr("target_accessibility_p"),
        "target_accessibility_p_17mer": _maybe_attr("target_accessibility_p_17mer"),
        "target_accessibility_p_site": _maybe_attr("target_accessibility_p_site"),
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
        "undetermined_hits": _maybe_attr("undetermined_hits", 0),
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
        # Scoring. design_score and composite_score are different vectors over different term
        # sets, so both are emitted and neither is filled in from the other.
        "design_score": candidate.design_score,
        "composite_score": candidate.composite_score,
        "score_off_target": candidate.score_off_target,
        "score_target_accessibility": candidate.score_target_accessibility,
        "score_asymmetry": candidate.score_asymmetry,
        "score_gc_content": candidate.score_gc_content,
        "score_ago_start": _maybe_attr("score_ago_start"),
        # score_pos1_mismatch is not exported: #102 removed pos1_mismatch from the only vector that
        # held it, so the contribution is null for every candidate of every run. An always-null column
        # is not neutral -- a reader who finds it in the header tries to use it. The field stays on the
        # model for whichever vector scores the term next; re-add the column with that vector.
        "score_supp_13_16": _maybe_attr("score_supp_13_16"),
        "empirical_score": cs.get("empirical"),
        "scored_after_screening": candidate.scored_after_screening,
        "weight_set_version": candidate.weight_set_version,
        "weight_vector": _maybe_attr("weight_vector", ""),
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
        # Per-filter verdicts, appended last so every column above keeps its position and meaning.
        # One pair per declared filter, always the same columns in the same order: a filter that did
        # not run writes `not_evaluated` rather than a blank, because a blank cell is read as a verdict.
        **{f"{filter_id}_verdict": verdicts.get(filter_id, _NOT_EVALUATED) for filter_id in DECLARED_FILTER_IDS},
        **{f"{filter_id}_observed": observed.get(filter_id) for filter_id in DECLARED_FILTER_IDS},
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
            "best_score": max([ranking_score(c) for c in self.top_candidates]) if self.top_candidates else 0,
            "tool_versions": self.tool_versions,
        }
