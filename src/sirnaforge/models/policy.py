"""Shared run-policy vocabulary: run modes, filter verdicts, scope and target intent.

Data contracts only. This module deliberately imports nothing from ``workflow`` and resolves,
defaults and applies nothing: the resolver lives in ``config/run_policy.py`` (#99), the evidence
producers in ``core``/Nextflow (#100), and intent evaluation in the classification path (#101).
Its only job is to give those packages one spelling for each concept.

Owned by #99; extended in place by #100 (verdict records) and #101 (intent fields, scope). A branch
needing a field this module lacks adds it here rather than building a parallel object elsewhere.

Two things are deliberately expressed as data rather than as field names. A filter's *direction*
lived only in the ``max_*``/``min_*`` prefix and in the body of ``_check_offtarget_filters``, so
:class:`FilterComparator` names it; a filter's *scope* lived only in which counter a gate happened
to read, so :class:`FilterScope` names the species set, mismatch ceiling and hit classes it counts.
Together they let a consumer regenerate a gate from the model instead of restating it.
"""

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, model_validator

DEFAULT_TARGET_SPECIES = "human"


class RunMode(str, Enum):
    """How much evidence a run is claiming to have.

    Attributes:
        DESIGN_ONLY: No screening was requested; only design-time evidence exists.
        EXPLORATORY: Screening ran, but incomplete evidence may be retained and labelled.
        QUALIFIED: Required evidence must be complete for a candidate to qualify.
    """

    DESIGN_ONLY = "design_only"
    EXPLORATORY = "exploratory"
    QUALIFIED = "qualified"


class FilterAction(str, Enum):
    """What a filter does with its own verdict.

    Attributes:
        OFF: The filter is not evaluated at all.
        WARN: The verdict is recorded and reported, but does not reject a candidate.
        FAIL: The verdict rejects a candidate.
    """

    OFF = "off"
    WARN = "warn"
    FAIL = "fail"


#: Every declared filter, in registry order. Lives here rather than beside ``FILTER_SPECS`` because
#: the candidate row needs it to emit one fixed verdict column per filter, and ``models`` cannot
#: import ``config.run_policy`` -- that module imports this package. ``FILTER_SPECS`` is the single
#: source of truth for what a filter *is*; this is only the ordered set of ids, and a unit test pins
#: the two against each other so a filter cannot be added in one place and forgotten in the other.
DECLARED_FILTER_IDS: tuple[str, ...] = (
    "gc_content_min",
    "gc_content_max",
    "max_poly_runs",
    "max_paired_fraction",
    "min_asymmetry_score",
    "min_empirical_score",
    "min_isoform_coverage",
    "max_off_target_count",
    "max_transcriptome_hits_0mm",
    "max_transcriptome_hits_1mm",
    "max_transcriptome_hits_2mm",
    "max_transcriptome_seed_perfect",
    "max_mirna_perfect_seed",
    "max_mirna_1mm_seed",
    "fail_on_high_risk_mirna",
    "max_total_offtarget_hits",
)


class FilterEvaluation(str, Enum):
    """The verdict itself, which is distinct from the action taken on it.

    Four outcomes, because three different things were previously spelled as a blank or a
    borrowed ``PASS``: the filter did not run, the filter ran but its evidence was missing, and
    the filter ran and decided. ``NOT_EVALUATED`` exists because #100 emits fixed columns -- an
    off filter must write a word, not an empty cell that a reader mistakes for a verdict.

    Attributes:
        PASS: Evidence was available and the candidate met the threshold.
        FAIL: Evidence was available and the candidate did not meet the threshold.
        UNKNOWN: The filter was in force but could not be decided: its evidence is unavailable.
        NOT_EVALUATED: The filter was not applied at all, so it makes no claim either way.
    """

    PASS = "pass"
    FAIL = "fail"
    UNKNOWN = "unknown"
    NOT_EVALUATED = "not_evaluated"


class UnknownEvidenceAction(str, Enum):
    """What an ``UNKNOWN`` evaluation does to a candidate.

    A deliberately narrower vocabulary than :class:`FilterAction`: ``OFF`` means "never
    evaluated", so it cannot be the response to a filter that *was* evaluated and came back
    undecided. The two spellings match ``FilterAction``'s so the serialised words agree.

    Attributes:
        WARN: The undecided candidate is retained and labelled.
        FAIL: The undecided candidate is rejected.
    """

    WARN = "warn"
    FAIL = "fail"


class FilterStage(str, Enum):
    """When a filter's input exists, which decides whether it can be evaluated at all.

    A design-only run holds no screening evidence, so every post-screen filter in it is
    ``NOT_EVALUATED`` -- not passed. Keeping the stage as data is what lets the resolver say that
    without a per-filter conditional.

    Attributes:
        DESIGN: The input is computed while candidates are enumerated.
        POST_SCREEN: The input comes from off-target screening.
    """

    DESIGN = "design"
    POST_SCREEN = "post_screen"


class SettingSource(str, Enum):
    """Where a resolved setting's value came from, recorded per setting in the manifest.

    Precedence is built-in profile < config file < explicit caller value; ``RUN_MODE_RULE`` and
    ``DESIGN_MODE_PRESET`` are the two derivations, and both are recorded so a value nobody typed
    can still be traced to the rule that produced it.

    Attributes:
        BUILTIN_PROFILE: The versioned built-in profile's baseline.
        DESIGN_MODE_PRESET: The design mode's preset, applied only where the caller said nothing.
        CONFIG_FILE: A configuration file.
        EXPLICIT: The caller stated it (a CLI option or an API keyword argument).
        RUN_MODE_RULE: Derived from the resolved run mode.
    """

    BUILTIN_PROFILE = "builtin_profile"
    DESIGN_MODE_PRESET = "design_mode_preset"
    CONFIG_FILE = "config_file"
    EXPLICIT = "explicit"
    RUN_MODE_RULE = "run_mode_rule"


class SettingProvenance(BaseModel):
    """One setting, its resolved value, and which authority supplied it.

    Attributes:
        key: Stable setting name (a filter id, or a named policy setting).
        value: The value that applies, JSON-representable.
        source: Which authority supplied it.
        detail: Why, when the source alone does not say.
    """

    key: str = Field(min_length=1, description="Stable setting name")
    value: bool | int | float | str | None = Field(description="Resolved value; None means no value applies")
    source: SettingSource = Field(description="Authority that supplied the value")
    detail: str | None = Field(default=None, description="Why this source won, when that is not obvious")

    model_config = ConfigDict(frozen=True, extra="forbid")


class ProfileIdentity(BaseModel):
    """Which built-in profile a run resolved against, and a hash of its contents.

    The name and version identify the profile; the hash is what makes two runs comparable, because
    a profile edited in place would otherwise report the same identity as its predecessor.

    Attributes:
        name: Profile name.
        version: Profile version.
        content_hash: ``sha256:`` digest over the profile's canonical serialisation.
        experimental: Whether the profile's numbers are calibrated. Every 0.7.1 profile is not.
        description: What the profile is for.
    """

    name: str = Field(min_length=1, description="Profile name")
    version: str = Field(min_length=1, description="Profile version")
    content_hash: str = Field(min_length=1, description="sha256 digest of the profile's canonical serialisation")
    experimental: bool = Field(description="False only for a profile whose numbers are calibrated")
    description: str = Field(min_length=1, description="What the profile is for")

    model_config = ConfigDict(frozen=True, extra="forbid")


class FilterComparator(str, Enum):
    """The passing predicate a filter applies, read as ``observed <comparator> threshold``.

    Named as the *pass* direction rather than the fail direction so a descriptor reads the same
    way in Python and in a generic client-side evaluator. ``max_*`` fields are :attr:`LE`,
    ``min_*`` fields are :attr:`GE`; a boolean "fail if any" flag is :attr:`LE` against 0. Ordering
    comparators only -- over a numeric threshold an equality gate is a pair of bounds, so adding
    ``eq``/``ne`` would only widen the table a client evaluator has to restate.

    Attributes:
        LE: Passes when the observed value is at most the threshold (a ceiling).
        LT: Passes when the observed value is below the threshold.
        GE: Passes when the observed value is at least the threshold (a floor).
        GT: Passes when the observed value is above the threshold.
    """

    LE = "le"
    LT = "lt"
    GE = "ge"
    GT = "gt"

    def passes(self, observed: float, threshold: float) -> bool:
        """Apply this comparator, so its meaning has exactly one implementation in Python.

        A client-side evaluator has to restate the table in its own language; nothing else does.
        """
        if self is FilterComparator.LE:
            return observed <= threshold
        if self is FilterComparator.LT:
            return observed < threshold
        if self is FilterComparator.GE:
            return observed >= threshold
        return observed > threshold


class FilterScope(BaseModel):
    """Which alignments a threshold counts, as data instead of as a counter's name.

    Every axis is independently unrestricted-by-default, and an unrestricted axis means "counts
    everything the screen saw on this axis" -- a well-defined scope, not an undeclared one. What
    the *default off-target scope should be* is a policy decision (D3: human, nm <= 2) and belongs
    to the resolver, not to this type.

    Attributes:
        species: Canonical species names counted. Empty means every screened species.
        max_mismatches: Highest guide-level mismatch count in scope. ``None`` means no ceiling.
        hit_classes: ``HitClass`` values counted (``core/hit_classification.py``), held as strings
            so this contract stays a leaf module and #101 can add members without editing it.
            Empty means every class, including on-target isoform and repeat alignments.
    """

    species: frozenset[str] = Field(default=frozenset(), description="Species counted; empty means all screened")
    max_mismatches: int | None = Field(default=None, ge=0, description="Mismatch ceiling in scope; None means no cap")
    hit_classes: frozenset[str] = Field(
        default=frozenset(), description="HitClass values counted; empty means every class"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def no_blank_members(self) -> "FilterScope":
        """Reject a blank species or hit class, which would silently narrow a scope to nothing real."""
        for axis, values in (("species", self.species), ("hit_classes", self.hit_classes)):
            if any(not value.strip() for value in values):
                raise ValueError(f"{axis} contains a blank entry; omit the axis instead of naming an empty value")
        return self


class FilterDescriptor(BaseModel):
    """One filter's identity and the comparison it makes, independent of any candidate.

    This is the unit #103 emits so a report can re-apply a gate without reimplementing it:
    ``{column, comparator, threshold, action, scope}`` generated from the model rather than
    restated in JavaScript. ``threshold`` may be ``None`` -- an undeclared threshold is why a
    filter can be in force and still not be evaluable.

    Attributes:
        filter_id: Stable machine identity, also the column stem for exported per-filter columns.
        column: The candidate column or feature the comparison reads.
        comparator: The passing direction.
        threshold: The value compared against, or ``None`` when none is declared.
        scope: Which alignments the threshold counts.
        action: What the filter does with its own verdict.
        stage: When the compared input exists. Defaults to :attr:`FilterStage.DESIGN`, which is the
            0.7.1 behaviour of every gate -- always evaluable -- so the default asserts nothing new;
            a filter whose input comes from screening must say so, and the resolver always does.
    """

    filter_id: str = Field(min_length=1, description="Stable machine identity of the filter")
    column: str = Field(min_length=1, description="Candidate column or feature read")
    comparator: FilterComparator = Field(description="Passing direction of the comparison")
    threshold: int | float | None = Field(description="Value compared against; None means undeclared")
    scope: FilterScope = Field(description="Alignments the threshold counts; state it even when unrestricted")
    action: FilterAction = Field(description="What the filter does with its verdict")
    stage: FilterStage = Field(default=FilterStage.DESIGN, description="When the compared input exists")

    model_config = ConfigDict(frozen=True, extra="forbid")


class FilterVerdict(BaseModel):
    """What one filter observed on one candidate, and what it concluded.

    The descriptor and the outcome travel together so a reported verdict can never be read
    against a threshold other than the one it was decided under, and a decided verdict must be the
    one its own descriptor produces. That is what lets a client re-apply a gate and get the
    pipeline's answer back: a descriptor that cannot reproduce its verdict is rejected here rather
    than discovered as a disagreement in a report. A gate that is not one comparison against one
    threshold does not get to publish a descriptor claiming it is.

    Whether an ``UNKNOWN`` rejects a candidate is deliberately not answered here: that needs
    :attr:`EvidenceRequirements.unknown_evidence_action`, which the caller holds.

    Attributes:
        descriptor: The filter as configured for this candidate.
        observed: The value the gate actually compared -- not a rounded display value, which could
            reproduce the wrong verdict. ``None`` when nothing was observed; present alongside
            ``NOT_EVALUATED`` is legitimate, being a measured value whose gate did not run.
        evaluation: The outcome.
    """

    descriptor: FilterDescriptor = Field(description="The filter as configured when this verdict was reached")
    observed: int | float | None = Field(default=None, description="Value read; None means nothing was observed")
    evaluation: FilterEvaluation = Field(description="Outcome of the comparison")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def verdict_is_supported(self) -> "FilterVerdict":
        """An off filter must not decide, and a decided verdict must be the one its descriptor produces."""
        if self.descriptor.action is FilterAction.OFF and self.evaluation is not FilterEvaluation.NOT_EVALUATED:
            raise ValueError(f"filter {self.descriptor.filter_id} is off, so its evaluation must be NOT_EVALUATED")
        if self.evaluation not in (FilterEvaluation.PASS, FilterEvaluation.FAIL):
            return self
        if self.observed is None or self.descriptor.threshold is None:
            raise ValueError(
                f"filter {self.descriptor.filter_id} reports {self.evaluation.value} without both an observed "
                "value and a threshold; use UNKNOWN or NOT_EVALUATED instead"
            )
        passed = self.descriptor.comparator.passes(self.observed, self.descriptor.threshold)
        if passed is not (self.evaluation is FilterEvaluation.PASS):
            raise ValueError(
                f"filter {self.descriptor.filter_id} reports {self.evaluation.value}, but "
                f"{self.observed} {self.descriptor.comparator.value} {self.descriptor.threshold} is {passed}"
            )
        return self


class ScreeningChannel(str, Enum):
    """An independent source of off-target liability evidence.

    Attributes:
        TRANSCRIPTOME: Full-length/near-full-length transcript alignment.
        MIRNA_SEED: Resemblance to a known miRNA seed.
        TRANSCRIPT_SEED: Complementary seed sites in transcript sequence.
    """

    TRANSCRIPTOME = "transcriptome"
    MIRNA_SEED = "mirna_seed"
    TRANSCRIPT_SEED = "transcript_seed"


class TargetSelectivity(str, Enum):
    """Whether every same-gene transcript is wanted, or only the declared ones.

    A same-gene hit is not automatically desirable, so this cannot be inferred from the target
    set: the two intents differ on exactly the transcripts neither set names.

    Attributes:
        PAN_ISOFORM: Knockdown of any transcript of the target gene is intended.
        ISOFORM_SELECTIVE: Only the required transcripts are intended; other same-gene
            transcripts are liabilities unless separately declared.
    """

    PAN_ISOFORM = "pan_isoform"
    ISOFORM_SELECTIVE = "isoform_selective"


class OrthologyEvidenceSource(str, Enum):
    """How a cross-species correspondence was established.

    The published hit table spells the same two tiers as ``gene_id``/``symbol_heuristic``
    (:class:`sirnaforge.core.hit_classification.OrthologEvidence`, resolved by
    :mod:`sirnaforge.data.orthology`); nothing converts between the two vocabularies yet.

    Attributes:
        EXPLICIT_MAPPING: A declared orthology mapping or database assertion.
        SYMBOL_EQUALITY: Gene symbols matched. A heuristic, and never a claim of validated
            orthology on its own.
    """

    EXPLICIT_MAPPING = "explicit_mapping"
    SYMBOL_EQUALITY = "symbol_equality"


class OrthologyAssertion(BaseModel):
    """One recorded cross-species correspondence, with how it was established.

    Attributes:
        species: Canonical species the correspondence is asserted in.
        gene_id: Version-stripped gene ID in that species, when known.
        gene_symbol: Gene symbol in that species, when known.
        source: How the correspondence was established.
        confidence: Reported confidence, when the source provides one.
    """

    species: str = Field(min_length=1, description="Canonical species the correspondence is asserted in")
    gene_id: str | None = Field(default=None, description="Version-stripped gene ID, when known")
    gene_symbol: str | None = Field(default=None, description="Gene symbol, when known")
    source: OrthologyEvidenceSource = Field(description="How the correspondence was established")
    confidence: float | None = Field(default=None, ge=0.0, le=1.0, description="Reported confidence, when available")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def identifies_a_gene(self) -> "OrthologyAssertion":
        """An assertion naming neither a gene ID nor a symbol asserts nothing."""
        if self.gene_id is None and self.gene_symbol is None:
            raise ValueError(f"orthology assertion for {self.species} names neither a gene ID nor a symbol")
        return self

    @property
    def is_validated_orthology(self) -> bool:
        """True only for an explicit mapping; symbol equality is a heuristic and says so here."""
        return self.source is OrthologyEvidenceSource.EXPLICIT_MAPPING


class TargetIntent(BaseModel):
    """What the user intends to knock down, and where they intend to look for liabilities.

    The two species sets are deliberately separate. A species may be a target, a screening
    subject, both, or neither: screening mouse does not make mouse knockdown a requirement, and
    requiring human knockdown does not by itself request a human screen. Callers that need one
    set to imply the other must say so explicitly.

    The three transcript-set fields are separate because coverage silently shrinking to the
    isoforms that survived design filtering is the defect they exist to prevent:
    ``annotation_universe`` is what the annotation source knew, ``coverage_denominator`` is what
    coverage is reported against and is frozen before design filtering, and
    ``enumeration_inputs`` is what guides were actually enumerated on.

    ``selectivity`` has no default: pan-isoform and isoform-selective runs disagree on exactly
    the transcripts neither the required nor the excluded set names, and inferring one is how
    "all retrieved same-gene transcripts are intended" became an unexamined assumption.

    Attributes:
        selectivity: Whether same-gene transcripts outside the required set are intended.
        target_species: Canonical species names whose transcripts are meant to be knocked down.
        offtarget_screen_species: Canonical species names to screen for off-target liabilities.
        required_transcript_ids: Transcripts a guide must cover, versionless.
        excluded_transcript_ids: Transcripts a guide must not cover, versionless. A hit on one of
            these is on-target by gene taxonomy and still unacceptable.
        annotation_universe: Every transcript the annotation source offered, versionless.
        coverage_denominator: The frozen denominator coverage is reported against. Empty means no
            denominator was frozen, so coverage is unknown rather than complete.
        enumeration_inputs: Transcripts actually submitted to guide enumeration.
        annotation_provenance: Identity of the annotation source behind the sets above.
        orthology_evidence: Recorded cross-species correspondences supporting non-query members
            of ``target_species``. Absence is unknown, not a denial.
    """

    selectivity: TargetSelectivity = Field(description="Pan-isoform or isoform-selective; the caller must declare it")
    target_species: frozenset[str] = Field(
        default=frozenset({DEFAULT_TARGET_SPECIES}),
        description="Species whose transcripts are the intended knockdown targets",
    )
    offtarget_screen_species: frozenset[str] = Field(
        default=frozenset(),
        description="Species to screen for off-target liabilities; independent of the target set",
    )
    required_transcript_ids: frozenset[str] = Field(
        default=frozenset(), description="Transcript IDs a guide is required to cover"
    )
    excluded_transcript_ids: frozenset[str] = Field(
        default=frozenset(), description="Transcript IDs a guide must not cover"
    )
    annotation_universe: frozenset[str] = Field(
        default=frozenset(), description="Every transcript the annotation source offered"
    )
    coverage_denominator: frozenset[str] = Field(
        default=frozenset(), description="Frozen coverage denominator; empty means coverage is unknown"
    )
    enumeration_inputs: frozenset[str] = Field(
        default=frozenset(), description="Transcripts actually submitted to guide enumeration"
    )
    annotation_provenance: str | None = Field(
        default=None, description="Identity of the annotation source behind the transcript sets"
    )
    orthology_evidence: tuple[OrthologyAssertion, ...] = Field(
        default=(), description="Cross-species correspondences supporting non-query target species"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def enumeration_was_filtered(self) -> bool | None:
        """Whether enumeration saw less than the coverage denominator, or ``None`` if unknowable.

        ``None`` when either set is undeclared, because an undeclared denominator cannot prove
        that enumeration was complete.
        """
        if not self.coverage_denominator or not self.enumeration_inputs:
            return None
        return bool(self.coverage_denominator - self.enumeration_inputs)


class Requiredness(str, Enum):
    """Whether one unit of evidence must be complete for a candidate to qualify.

    Attributes:
        REQUIRED: Incompleteness blocks qualification.
        EXPLORATORY: Incompleteness is reportable but not disqualifying.
    """

    REQUIRED = "required"
    EXPLORATORY = "exploratory"


class ChannelRequirement(BaseModel):
    """Requiredness for one channel in one species.

    Per channel *and* species because the two axes do not factorise: a transcriptome screen may be
    required in human and exploratory in mouse while the miRNA channel is exploratory in both.
    Enumerating pairs is what makes #101's guarantee hold -- adding or removing an exploratory
    pair cannot change any required pair, because required pairs are named individually.

    Attributes:
        channel: Liability channel.
        species: Canonical species name.
        requiredness: Whether this pair's completion is required.
    """

    channel: ScreeningChannel = Field(description="Liability channel")
    species: str = Field(min_length=1, description="Canonical species name")
    requiredness: Requiredness = Field(description="Required or exploratory for this pair")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def key(self) -> tuple[str, str]:
        """The channel/species pair this requirement speaks for."""
        return (self.channel.value, self.species)


class EvidenceRequirements(BaseModel):
    """Which screening evidence a run requires, and how an undecided filter is treated.

    The single authority on requiredness. Plan and evidence entries deliberately carry no
    ``required`` flag of their own: a producer that could restate requiredness could contradict
    this object, and the contradiction would be invisible.

    ``unknown_evidence_action`` has no default on purpose: choosing what an ``UNKNOWN``
    evaluation does to a candidate is a policy decision owned by the resolver, not a fallback
    this contract is entitled to pick.

    Attributes:
        channel_requirements: One entry per channel x species pair the run takes a position on.
        unknown_evidence_action: What to do with a filter that evaluated to ``UNKNOWN``.
    """

    channel_requirements: tuple[ChannelRequirement, ...] = Field(
        default=(), description="Requiredness per channel x species pair"
    )
    unknown_evidence_action: UnknownEvidenceAction = Field(
        description="Action applied when a filter's evaluation is UNKNOWN; the caller must decide"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def one_position_per_pair(self) -> "EvidenceRequirements":
        """Two entries for one channel/species pair would reintroduce the ambiguity this type removes."""
        keys = [entry.key for entry in self.channel_requirements]
        if len(keys) != len(set(keys)):
            duplicates = sorted({key for key in keys if keys.count(key) > 1})
            raise ValueError(f"duplicate channel/species requirements: {duplicates}")
        return self

    def requiredness_of(self, channel: ScreeningChannel, species: str) -> Requiredness | None:
        """Requiredness declared for one pair, or ``None`` when the run declared none.

        ``None`` is not a synonym for exploratory: what an undeclared pair means is the caller's
        decision, and it is a different decision in each run mode.
        """
        for entry in self.channel_requirements:
            if entry.channel is channel and entry.species == species:
                return entry.requiredness
        return None

    @property
    def required_pairs(self) -> frozenset[tuple[str, str]]:
        """Channel/species pairs whose completion is required."""
        return frozenset(
            entry.key for entry in self.channel_requirements if entry.requiredness is Requiredness.REQUIRED
        )
