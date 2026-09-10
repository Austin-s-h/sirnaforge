"""One resolution point for run mode, filter actions and versioned numerical defaults.

``resolve_run_policy`` is **pure**: it downloads nothing, folds nothing and aligns nothing. It
resolves precedence once -- built-in versioned profile < design-mode preset < config file <
explicit CLI/API value -- validates the result by constructing
:class:`sirnaforge.models.sirna.DesignParameters`, and returns an immutable
:class:`ResolvedRunPolicy` that the CLI, the Python API and the off-target-only path all receive.

Two properties are the point of the module:

* **A value that nobody typed is traceable to the rule that produced it.** Every setting carries a
  :class:`~sirnaforge.models.policy.SettingProvenance` record, so "resolved to 52 because the miRNA
  preset applied" and "resolved to 52 because you asked for 52" are different facts in the manifest.
* **Provenance, not value equality, decides whether a preset applies.** ``cli.py`` used to compare a
  supplied value against the siRNA default, so an explicit ``--gc-max 60`` in miRNA mode was
  rewritten to 52. A field-set test cannot make that mistake: an omitted option and an option typed
  with the default value are different inputs here.

What this module deliberately does **not** do: resolve the screening *reference* (that is the
transcriptome resolver, tracked separately), decide candidate eligibility from evidence (#100), or
promote any number. Every threshold below is the one 0.7.1 already applied; the ``legacy`` profile
carries them verbatim so a resolved run reproduces an unresolved one.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, get_args

import tomllib
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.models.policy import (
    DEFAULT_TARGET_SPECIES,
    ChannelRequirement,
    EvidenceRequirements,
    FilterAction,
    FilterComparator,
    FilterDescriptor,
    FilterScope,
    FilterStage,
    ProfileIdentity,
    Requiredness,
    RunMode,
    ScreeningChannel,
    SettingProvenance,
    SettingSource,
    UnknownEvidenceAction,
)
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    FilterCriteria,
    MiRNADesignConfig,
    OffTargetFilterCriteria,
    ScoringWeights,
    TargetAccessibilityConfig,
)

LEGACY_PROFILE_NAME = "legacy"
DEFAULT_PROFILE_NAME = LEGACY_PROFILE_NAME
POLICY_SCHEMA_VERSION = "0.7.1"


class RunPolicyError(ValueError):
    """A run cannot proceed as configured, raised before any expensive work happens.

    A ``ValueError`` so one CLI boundary catches it alongside Pydantic's own errors.
    """


class EntryPoint(str, Enum):
    """Which public surface asked for a policy, which decides only the *default* run mode.

    Attributes:
        DESIGN_COMMAND: ``sirnaforge design`` / a direct designer call. Design-only by default.
        SCREENING_WORKFLOW: ``sirnaforge workflow`` / ``run_sirna_workflow``. Qualified by default.
        OFFTARGET_ONLY: ``sirnaforge offtarget`` / ``run_offtarget_only_workflow``. Qualified.
    """

    DESIGN_COMMAND = "design_command"
    SCREENING_WORKFLOW = "screening_workflow"
    OFFTARGET_ONLY = "offtarget_only"


ENTRY_POINT_DEFAULT_RUN_MODE: Mapping[EntryPoint, RunMode] = {
    EntryPoint.DESIGN_COMMAND: RunMode.DESIGN_ONLY,
    EntryPoint.SCREENING_WORKFLOW: RunMode.QUALIFIED,
    EntryPoint.OFFTARGET_ONLY: RunMode.QUALIFIED,
}


@dataclass(frozen=True)
class SettingSpec:
    """Where one public setting name lands on the validated models.

    ``key`` is the public name used by the CLI, the config file, the provenance record and the
    manifest; ``model`` and ``field`` are where it is validated. The indirection exists because
    three public names (``plfold_window`` and friends) do not match their model fields, and a
    manifest keyed on model fields would not be greppable against a CLI option.
    """

    key: str
    model: str
    field: str


# Every setting the resolver accepts. A key absent from here is an undeclared term and is rejected,
# so a typo in a config file cannot read as an omitted option.
SETTING_SPECS: tuple[SettingSpec, ...] = (
    SettingSpec("gc_min", "filters", "gc_min"),
    SettingSpec("gc_max", "filters", "gc_max"),
    SettingSpec("max_poly_runs", "filters", "max_poly_runs"),
    SettingSpec("max_paired_fraction", "filters", "max_paired_fraction"),
    SettingSpec("min_asymmetry_score", "filters", "min_asymmetry_score"),
    SettingSpec("min_empirical_score", "filters", "min_empirical_score"),
    SettingSpec("min_isoform_coverage", "filters", "min_isoform_coverage"),
    SettingSpec("max_off_target_count", "offtarget_filters", "max_off_target_count"),
    SettingSpec("max_transcriptome_hits_0mm", "offtarget_filters", "max_transcriptome_hits_0mm"),
    SettingSpec("max_transcriptome_hits_1mm", "offtarget_filters", "max_transcriptome_hits_1mm"),
    SettingSpec("max_transcriptome_hits_2mm", "offtarget_filters", "max_transcriptome_hits_2mm"),
    SettingSpec("max_transcriptome_seed_perfect", "offtarget_filters", "max_transcriptome_seed_perfect"),
    SettingSpec("max_mirna_perfect_seed", "offtarget_filters", "max_mirna_perfect_seed"),
    SettingSpec("max_mirna_1mm_seed", "offtarget_filters", "max_mirna_1mm_seed"),
    SettingSpec("fail_on_high_risk_mirna", "offtarget_filters", "fail_on_high_risk_mirna"),
    SettingSpec("max_total_offtarget_hits", "offtarget_filters", "max_total_offtarget_hits"),
    SettingSpec("plfold_window", "target_accessibility", "window_size"),
    SettingSpec("plfold_max_bp_span", "target_accessibility", "max_bp_span"),
    SettingSpec("accessibility_log_floor", "target_accessibility", "log_floor"),
    SettingSpec("sirna_length", "design", "sirna_length"),
    SettingSpec("top_n", "design", "top_n"),
    SettingSpec("check_off_targets", "design", "check_off_targets"),
    SettingSpec("predict_structure", "design", "predict_structure"),
    SettingSpec("avoid_snps", "design", "avoid_snps"),
    SettingSpec("apply_modifications", "design", "apply_modifications"),
    SettingSpec("modification_pattern", "design", "modification_pattern"),
    SettingSpec("default_overhang", "design", "default_overhang"),
)

SETTING_BY_KEY: Mapping[str, SettingSpec] = {spec.key: spec for spec in SETTING_SPECS}

_MODEL_FOR_TARGET: Mapping[str, type[BaseModel]] = {
    "filters": FilterCriteria,
    "offtarget_filters": OffTargetFilterCriteria,
    "target_accessibility": TargetAccessibilityConfig,
    "design": DesignParameters,
}

# The one place the ``legacy`` profile knowingly differs from a model field default. The CLI has
# always sent 30.0, so 30.0 is what "today's numbers verbatim" means for a resolved run; the
# FilterCriteria field default of 35.0 applies only to code that constructs the model directly.
# Which of the two is biologically right is a filter-default decision and is not settled here.
LEGACY_DEFAULT_EXCEPTIONS: Mapping[str, tuple[Any, str]] = {
    "gc_min": (30.0, "delivered CLI default; FilterCriteria's field default of 35.0 is not what any run applied"),
}


def _owner(parameters: DesignParameters, spec: SettingSpec) -> BaseModel:
    """The model instance one setting lives on; the pseudo-model ``design`` is the root object."""
    return parameters if spec.model == "design" else getattr(parameters, spec.model)


def _switch_off_value(spec: SettingSpec) -> Any:
    """The value that removes one setting's threshold, or raise if the field has no such state.

    A gate is switched off by clearing the number it compares against, which is a state the model
    already has for every nullable threshold (``_check_offtarget_filters`` skips a ``None``) and for
    the one boolean flag. The six design-stage thresholds are plain floats with no "absent" value, so
    the resolver refuses rather than fabricating an inert number and reporting it as a threshold.
    """
    field = _MODEL_FOR_TARGET[spec.model].model_fields[spec.field]
    if type(None) in get_args(field.annotation):
        return None
    if field.annotation is bool:
        return False
    raise RunPolicyError(
        f"{spec.key} has no 'no threshold' state in 0.7.1, so the filter reading it cannot be "
        "switched off; widen the threshold instead"
    )


def _model_default(spec: SettingSpec) -> Any:
    """The declared default of the field a setting lands on, so a profile cannot restate a number."""
    return _MODEL_FOR_TARGET[spec.model].model_fields[spec.field].get_default(call_default_factory=False)


@dataclass(frozen=True)
class RunPolicyProfile:
    """A named, versioned baseline for every declared setting.

    Baselines are *derived* from the model field defaults rather than retyped, with the exceptions
    above declared explicitly, so a profile and the models it validates against cannot drift.
    """

    name: str
    version: str
    description: str
    experimental: bool
    baseline: Mapping[str, Any]
    exceptions: Mapping[str, str]

    @staticmethod
    def legacy() -> RunPolicyProfile:
        """0.7.1's shipped numbers, exactly as an unresolved run applied them."""
        baseline = {spec.key: _model_default(spec) for spec in SETTING_SPECS}
        exceptions: dict[str, str] = {}
        for key, (value, reason) in LEGACY_DEFAULT_EXCEPTIONS.items():
            baseline[key] = value
            exceptions[key] = reason
        return RunPolicyProfile(
            name=LEGACY_PROFILE_NAME,
            version=POLICY_SCHEMA_VERSION,
            description=(
                "0.7.1's shipped thresholds, carried verbatim. Not calibrated: most numbers are "
                "expert priors, and the two that were measured were measured on one target."
            ),
            experimental=True,
            baseline=baseline,
            exceptions=exceptions,
        )

    def identity(self) -> ProfileIdentity:
        """Name, version and a content hash over the baseline this profile actually applied."""
        payload = json.dumps(
            {"name": self.name, "version": self.version, "baseline": self.baseline, "exceptions": self.exceptions},
            sort_keys=True,
            default=str,
        )
        digest = hashlib.sha256(payload.encode()).hexdigest()
        return ProfileIdentity(
            name=self.name,
            version=self.version,
            content_hash=f"sha256:{digest}",
            experimental=self.experimental,
            description=self.description,
        )


BUILTIN_PROFILES: Mapping[str, RunPolicyProfile] = {LEGACY_PROFILE_NAME: RunPolicyProfile.legacy()}


def _mirna_preset() -> Mapping[str, Any]:
    """miRNA-mode setting overrides, read off ``MiRNADesignConfig`` so the two cannot diverge.

    Only settings the preset genuinely moves are listed. ``modification_pattern`` is deliberately
    absent: ``MiRNADesignConfig.modifications`` equals the siRNA default, so the sentinel that used
    to "apply" it never changed anything.
    """
    preset = MiRNADesignConfig()
    return {
        "gc_min": preset.gc_min,
        "gc_max": preset.gc_max,
        "min_asymmetry_score": preset.asymmetry_min,
        "max_poly_runs": preset.max_homopolymer,
        "default_overhang": preset.overhang,
        "modification_pattern": preset.modifications,
    }


class ResolvedFilter(BaseModel):
    """One filter as this run will apply it, plus what a reader needs to interpret it.

    The descriptor is the machine-readable half (#103 re-applies it client-side). The rest is the
    half that keeps a report honest: ``definition`` says what the number counts, ``transform``
    names any conversion between the model field and the comparison, and ``evidence_exported``
    says whether the column the comparison reads is actually in the candidate CSV today.
    """

    descriptor: FilterDescriptor = Field(description="The gate as configured for this run")
    setting_key: str = Field(min_length=1, description="Public setting name carrying the threshold")
    definition: str = Field(min_length=1, description="What the compared number counts, in words")
    transform: str | None = Field(default=None, description="Conversion applied between field and comparison")
    evidence_exported: bool = Field(
        description="Whether the column the comparison reads is present in the candidate CSV today"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def filter_id(self) -> str:
        """Stable identity of the gate."""
        return self.descriptor.filter_id

    @property
    def is_evaluated(self) -> bool:
        """False when the gate is off or has no threshold: either way it makes no claim."""
        return self.descriptor.action is not FilterAction.OFF and self.descriptor.threshold is not None


@dataclass(frozen=True)
class _FilterSpec:
    """Static description of one gate: everything about it that is not the resolved threshold."""

    filter_id: str
    setting_key: str
    column: str
    comparator: FilterComparator
    stage: FilterStage
    definition: str
    default_action: FilterAction
    scope_species: tuple[str, ...] = ()
    scope_max_mismatches: int | None = None
    scope_hit_classes: tuple[str, ...] = ()
    transform: str | None = None
    evidence_exported: bool = True


# Liability classes as the pipeline counts them: an undecidable hit is counted as a liability, so a
# missing reference cannot loosen a screen (``core.hit_annotation.LIABILITY_CLASSES``). Held as
# strings because ``FilterScope`` is a leaf contract.
_LIABILITY_CLASSES: tuple[str, ...] = ("off_target", "undetermined")

# Why several off-target gates carry an explicit ``{human}`` scope: the counters they read are
# built from hits whose species is human *or* unlabelled, so the gate is not the all-species number
# the equally-named candidate column reports. Exporting the human-stratified counters is #101's
# filter-scope work; naming the scope here is what makes the discrepancy legible instead of silent.
_HUMAN_STRATIFIED = ("human",)
_HUMAN_STRATIFIED_NOTE = (
    "counts hits whose species is human or unlabelled, not every screened species; the identically "
    "named candidate column is the all-species number, so the two disagree on a multi-species run"
)

FILTER_SPECS: tuple[_FilterSpec, ...] = (
    _FilterSpec(
        filter_id="gc_content_min",
        setting_key="gc_min",
        column="gc_content",
        comparator=FilterComparator.GE,
        stage=FilterStage.DESIGN,
        definition="Guide GC content percentage floor; fails GC_OUT_OF_RANGE.",
        default_action=FilterAction.FAIL,
    ),
    _FilterSpec(
        filter_id="gc_content_max",
        setting_key="gc_max",
        column="gc_content",
        comparator=FilterComparator.LE,
        stage=FilterStage.DESIGN,
        definition="Guide GC content percentage ceiling; fails GC_OUT_OF_RANGE.",
        default_action=FilterAction.FAIL,
    ),
    _FilterSpec(
        filter_id="max_poly_runs",
        setting_key="max_poly_runs",
        column="max_poly_run_length",
        comparator=FilterComparator.LE,
        stage=FilterStage.DESIGN,
        definition="Longest run of one nucleotide in the guide; fails POLY_RUNS.",
        default_action=FilterAction.FAIL,
        evidence_exported=False,
    ),
    _FilterSpec(
        filter_id="max_paired_fraction",
        setting_key="max_paired_fraction",
        column="paired_fraction",
        comparator=FilterComparator.LE,
        stage=FilterStage.DESIGN,
        definition="Fraction of guide bases paired in its own MFE structure; fails EXCESS_PAIRING.",
        default_action=FilterAction.FAIL,
        transform="quantised to 2k/length by the dot-bracket, so few values are attainable on a 21mer",
    ),
    _FilterSpec(
        filter_id="min_asymmetry_score",
        setting_key="min_asymmetry_score",
        column="asymmetry_score",
        comparator=FilterComparator.GE,
        stage=FilterStage.DESIGN,
        definition="Thermodynamic asymmetry floor for RISC loading; fails LOW_ASYMMETRY.",
        default_action=FilterAction.FAIL,
    ),
    _FilterSpec(
        filter_id="min_empirical_score",
        setting_key="min_empirical_score",
        column="empirical_score",
        comparator=FilterComparator.GE,
        stage=FilterStage.DESIGN,
        definition="Simplified Reynolds design-rule floor; fails LOW_EMPIRICAL_SCORE.",
        default_action=FilterAction.FAIL,
        transform="attainable range is {0.4, 0.5, 0.6}; the default of 0.4 rejects nothing",
    ),
    _FilterSpec(
        filter_id="min_isoform_coverage",
        setting_key="min_isoform_coverage",
        column="isoform_coverage",
        comparator=FilterComparator.GE,
        stage=FilterStage.POST_SCREEN,
        definition="Protein-coding isoform coverage floor; fails LOW_ISOFORM_COVERAGE.",
        default_action=FilterAction.FAIL,
    ),
    _FilterSpec(
        filter_id="max_off_target_count",
        setting_key="max_off_target_count",
        column="off_target_count",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=(
            "Genuine off-target sites over every screened species and every mismatch count, with "
            "on-target, ortholog and repeat hits classified out and undecidable hits counted in; "
            "fails EXCESS_OFF_TARGETS."
        ),
        default_action=FilterAction.FAIL,
        scope_hit_classes=_LIABILITY_CLASSES,
    ),
    _FilterSpec(
        filter_id="max_transcriptome_hits_0mm",
        setting_key="max_transcriptome_hits_0mm",
        column="transcriptome_hits_0mm_human",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=f"Perfect-match genuine off-target transcriptome hits: {_HUMAN_STRATIFIED_NOTE}.",
        default_action=FilterAction.FAIL,
        scope_species=_HUMAN_STRATIFIED,
        scope_max_mismatches=0,
        scope_hit_classes=_LIABILITY_CLASSES,
        evidence_exported=False,
    ),
    _FilterSpec(
        filter_id="max_transcriptome_hits_1mm",
        setting_key="max_transcriptome_hits_1mm",
        column="transcriptome_hits_1mm_human",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=f"1-mismatch genuine off-target hits: {_HUMAN_STRATIFIED_NOTE}.",
        default_action=FilterAction.FAIL,
        scope_species=_HUMAN_STRATIFIED,
        scope_max_mismatches=1,
        scope_hit_classes=_LIABILITY_CLASSES,
        evidence_exported=False,
    ),
    _FilterSpec(
        filter_id="max_transcriptome_hits_2mm",
        setting_key="max_transcriptome_hits_2mm",
        column="transcriptome_hits_2mm_human",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=f"2-mismatch genuine off-target hits: {_HUMAN_STRATIFIED_NOTE}.",
        default_action=FilterAction.FAIL,
        scope_species=_HUMAN_STRATIFIED,
        scope_max_mismatches=2,
        scope_hit_classes=_LIABILITY_CLASSES,
        evidence_exported=False,
    ),
    _FilterSpec(
        filter_id="max_transcriptome_seed_perfect",
        setting_key="max_transcriptome_seed_perfect",
        column="transcriptome_hits_seed_0mm",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=(
            "Transcriptome hits whose seed (positions 2-8) paired perfectly, over every screened "
            "species. The only targeted gate on a clipped or gapped hit with an intact seed, which "
            "carries a guide-level nm > 2 and so lands in no mismatch stratum; fails "
            "TRANSCRIPTOME_SEED_PERFECT."
        ),
        default_action=FilterAction.FAIL,
    ),
    _FilterSpec(
        filter_id="max_mirna_perfect_seed",
        setting_key="max_mirna_perfect_seed",
        column="mirna_hits_0mm_seed_human",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=f"Perfect miRNA seed matches: {_HUMAN_STRATIFIED_NOTE}; fails MIRNA_PERFECT_SEED.",
        default_action=FilterAction.FAIL,
        scope_species=_HUMAN_STRATIFIED,
        evidence_exported=False,
    ),
    _FilterSpec(
        filter_id="max_mirna_1mm_seed",
        setting_key="max_mirna_1mm_seed",
        column="mirna_hits_1mm_seed",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=(
            "1-mismatch miRNA seed hits. Declared with a threshold of 10 but read by no gate in "
            "0.7.1, so it resolves to off: the number is reported and nothing acts on it."
        ),
        default_action=FilterAction.OFF,
    ),
    _FilterSpec(
        filter_id="fail_on_high_risk_mirna",
        setting_key="fail_on_high_risk_mirna",
        column="mirna_hits_high_risk_human",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=(
            f"High-risk miRNA hits (perfect seed and offtarget_score < 5.0): {_HUMAN_STRATIFIED_NOTE}; "
            "fails HIGH_RISK_MIRNA."
        ),
        default_action=FilterAction.FAIL,
        scope_species=_HUMAN_STRATIFIED,
        transform="a boolean flag, expressed as a ceiling of 0 so it carries a comparator like every other gate",
        evidence_exported=False,
    ),
    _FilterSpec(
        filter_id="max_total_offtarget_hits",
        setting_key="max_total_offtarget_hits",
        column="total_offtarget_hits_human",
        comparator=FilterComparator.LE,
        stage=FilterStage.POST_SCREEN,
        definition=(f"Transcriptome plus miRNA hits combined: {_HUMAN_STRATIFIED_NOTE}; fails TOTAL_OFFTARGETS."),
        default_action=FilterAction.FAIL,
        evidence_exported=False,
        scope_species=_HUMAN_STRATIFIED,
    ),
)

FILTER_SPEC_BY_ID: Mapping[str, _FilterSpec] = {spec.filter_id: spec for spec in FILTER_SPECS}


class ResolvedRunPolicy(BaseModel):
    """Everything one run resolved to, and where each part came from.

    Immutable. ``design_parameters`` is the validated owner of every number -- this object adds the
    run mode, the per-filter actions, the evidence requirements and the provenance trail, and does
    not hold a second copy of any threshold.
    """

    schema_version: str = Field(default=POLICY_SCHEMA_VERSION, description="Policy payload version")
    entry_point: EntryPoint = Field(description="Public surface that resolved this policy")
    run_mode: RunMode = Field(description="How much evidence this run claims to have")
    design_mode: DesignMode = Field(description="siRNA/miRNA/ZFN design mode; orthogonal to run_mode")
    profile: ProfileIdentity = Field(description="Built-in profile this run resolved against")
    design_parameters: DesignParameters = Field(description="Validated parameters; the only holder of thresholds")
    filters: tuple[ResolvedFilter, ...] = Field(description="Every declared gate, in reporting order")
    evidence_requirements: EvidenceRequirements = Field(description="Which screening evidence this run requires")
    requested: tuple[SettingProvenance, ...] = Field(description="What the caller and config file actually stated")
    resolved: tuple[SettingProvenance, ...] = Field(description="Every declared setting and the authority that won")

    model_config = ConfigDict(frozen=True, extra="forbid")

    def descriptor(self, filter_id: str) -> FilterDescriptor:
        """The descriptor for one gate, so a consumer never rebuilds a threshold by hand."""
        for resolved in self.filters:
            if resolved.filter_id == filter_id:
                return resolved.descriptor
        raise KeyError(f"no filter named {filter_id!r}")

    @property
    def evaluated_filters(self) -> tuple[ResolvedFilter, ...]:
        """Gates that will actually decide something. The rest are ``NOT_EVALUATED``, not passed."""
        return tuple(resolved for resolved in self.filters if resolved.is_evaluated)

    def source_of(self, key: str) -> SettingSource:
        """Which authority supplied one setting's value."""
        for record in self.resolved:
            if record.key == key:
                return record.source
        raise KeyError(f"no setting named {key!r}")

    def value_of(self, key: str) -> Any:
        """The value one setting resolved to, read from the provenance trail."""
        for record in self.resolved:
            if record.key == key:
                return record.value
        raise KeyError(f"no setting named {key!r}")

    def as_manifest(self) -> dict[str, Any]:
        """The manifest block: requested and resolved settings, actions, terms and profile identity."""
        return {
            "schema_version": self.schema_version,
            "entry_point": self.entry_point.value,
            "run_mode": self.run_mode.value,
            "design_mode": self.design_mode.value,
            "profile": self.profile.model_dump(mode="json"),
            "requested_settings": [record.model_dump(mode="json") for record in self.requested],
            "resolved_settings": [record.model_dump(mode="json") for record in self.resolved],
            "filters": [
                {
                    "filter_id": resolved.filter_id,
                    "setting_key": resolved.setting_key,
                    "column": resolved.descriptor.column,
                    "comparator": resolved.descriptor.comparator.value,
                    "threshold": resolved.descriptor.threshold,
                    "action": resolved.descriptor.action.value,
                    "stage": resolved.descriptor.stage.value,
                    "scope": resolved.descriptor.scope.model_dump(mode="json"),
                    "definition": resolved.definition,
                    "transform": resolved.transform,
                    "evidence_exported": resolved.evidence_exported,
                    "evaluated": resolved.is_evaluated,
                }
                for resolved in self.filters
            ],
            "evidence_requirements": self.evidence_requirements.model_dump(mode="json"),
        }


@dataclass(frozen=True)
class _Layer:
    """One precedence layer: the values it supplies and the authority it speaks with."""

    source: SettingSource
    values: Mapping[str, Any]
    detail: str | None = None


def load_policy_config(path: Path | str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    """Read a policy config file, returning ``(settings, run_directives)``.

    JSON and TOML only, chosen by suffix. ``settings`` are declared setting keys; ``run_directives``
    are the file's ``run_mode`` / ``design_mode`` / ``profile`` / ``filter_actions`` keys. Anything
    else is an undeclared term and raises here rather than being ignored.
    """
    config_path = Path(path)
    if not config_path.is_file():
        raise RunPolicyError(f"policy config file not found: {config_path}")
    text = config_path.read_text()
    try:
        payload = tomllib.loads(text) if config_path.suffix.lower() == ".toml" else json.loads(text)
    except (json.JSONDecodeError, tomllib.TOMLDecodeError) as exc:
        raise RunPolicyError(f"policy config file {config_path} is not valid {config_path.suffix.lstrip('.')}: {exc}")
    if not isinstance(payload, dict):
        raise RunPolicyError(f"policy config file {config_path} must contain a mapping at the top level")

    directives = {
        key: payload.pop(key) for key in ("run_mode", "design_mode", "profile", "filter_actions") if key in payload
    }
    settings = payload.pop("settings", {})
    if payload:
        raise RunPolicyError(
            f"policy config file {config_path} declares unknown top-level keys: {sorted(payload)}. "
            "Settings belong under 'settings'."
        )
    if not isinstance(settings, dict):
        raise RunPolicyError(f"policy config file {config_path}: 'settings' must be a mapping")
    _reject_undeclared(settings, f"policy config file {config_path}")
    return settings, directives


def _reject_undeclared(values: Mapping[str, Any], origin: str) -> None:
    """An unknown setting name must fail loudly: silently ignoring it reads as an applied override."""
    unknown = sorted(set(values) - set(SETTING_BY_KEY))
    if unknown:
        raise RunPolicyError(
            f"{origin} declares settings this build does not have: {unknown}. Declared settings: "
            f"{sorted(SETTING_BY_KEY)}"
        )


def _coerce_action(value: Any, filter_id: str) -> FilterAction:
    """Parse a per-filter action, naming the filter in the error so a typo is findable."""
    try:
        return FilterAction(str(value).lower())
    except ValueError as exc:
        raise RunPolicyError(
            f"filter {filter_id!r} was given action {value!r}; choose one of "
            f"{[action.value for action in FilterAction]}"
        ) from exc


def _resolve_actions(
    requested: Mapping[str, Any] | None,
) -> dict[str, FilterAction]:
    """Per-filter actions the caller stated, validated against the declared filter set."""
    if not requested:
        return {}
    unknown = sorted(set(requested) - set(FILTER_SPEC_BY_ID))
    if unknown:
        raise RunPolicyError(
            f"filter actions name filters this build does not have: {unknown}. Declared filters: "
            f"{sorted(FILTER_SPEC_BY_ID)}"
        )
    return {filter_id: _coerce_action(value, filter_id) for filter_id, value in requested.items()}


def _parse_run_mode(value: Any) -> RunMode:
    try:
        return value if isinstance(value, RunMode) else RunMode(str(value).lower())
    except ValueError as exc:
        raise RunPolicyError(f"invalid run mode {value!r}; choose one of {[mode.value for mode in RunMode]}") from exc


def parse_design_mode(value: Any) -> DesignMode:
    """Parse a design-mode string. Design mode is orthogonal to run mode and implies no evidence."""
    try:
        return value if isinstance(value, DesignMode) else DesignMode(str(value).lower())
    except ValueError as exc:
        raise RunPolicyError(
            f"invalid design mode {value!r}; choose one of {[mode.value for mode in DesignMode]}"
        ) from exc


def _resolve_run_mode(
    *,
    entry_point: EntryPoint,
    requested_run_mode: RunMode | str | None,
    stated: Mapping[str, Any],
    legacy_skip_screening: bool | None,
    require_screening_completeness: bool | None,
) -> tuple[RunMode, list[SettingProvenance], dict[str, Any]]:
    """Resolve the run mode and the screening switch together, since each constrains the other.

    Returns the mode, the provenance records explaining it, and any derived setting values (today
    only ``check_off_targets``) that the mode itself implies.
    """
    explicit_mode = _parse_run_mode(requested_run_mode) if requested_run_mode is not None else None
    stated_screening = stated.get("check_off_targets")
    records: list[SettingProvenance] = []
    derived: dict[str, Any] = {}

    if legacy_skip_screening and stated_screening is True:
        raise RunPolicyError(
            "skip-screening was requested while check_off_targets=True was stated explicitly; "
            "these are the same switch and must agree"
        )
    screening_off = legacy_skip_screening is True or stated_screening is False

    if explicit_mode is RunMode.QUALIFIED and screening_off:
        raise RunPolicyError(
            "run_mode=qualified requires off-target screening, but screening was explicitly "
            "disabled. Use run_mode=design_only, or leave screening enabled"
        )
    if explicit_mode is RunMode.DESIGN_ONLY and stated_screening is True:
        raise RunPolicyError(
            "run_mode=design_only cannot be combined with check_off_targets=True: a design-only run "
            "holds no screening evidence by definition"
        )
    if explicit_mode is RunMode.QUALIFIED and require_screening_completeness is False:
        raise RunPolicyError(
            "run_mode=qualified cannot waive required screening completeness: that is what "
            "qualified means. Use run_mode=exploratory to keep incomplete evidence and say so"
        )

    if explicit_mode is not None:
        mode = explicit_mode
        records.append(SettingProvenance(key="run_mode", value=mode.value, source=SettingSource.EXPLICIT, detail=None))
    elif legacy_skip_screening is True:
        mode = RunMode.DESIGN_ONLY
        records.append(
            SettingProvenance(
                key="run_mode",
                value=mode.value,
                source=SettingSource.RUN_MODE_RULE,
                detail="legacy skip-screening flag maps to design-only",
            )
        )
    elif stated_screening is False:
        mode = RunMode.DESIGN_ONLY
        records.append(
            SettingProvenance(
                key="run_mode",
                value=mode.value,
                source=SettingSource.RUN_MODE_RULE,
                detail="check_off_targets=False was stated, so no screening evidence exists",
            )
        )
    elif require_screening_completeness is False:
        mode = RunMode.EXPLORATORY
        records.append(
            SettingProvenance(
                key="run_mode",
                value=mode.value,
                source=SettingSource.RUN_MODE_RULE,
                detail="screening completeness was waived, which is exploratory rather than qualified",
            )
        )
    else:
        mode = ENTRY_POINT_DEFAULT_RUN_MODE[entry_point]
        records.append(
            SettingProvenance(
                key="run_mode",
                value=mode.value,
                source=SettingSource.BUILTIN_PROFILE,
                detail=f"default run mode for the {entry_point.value} entry point",
            )
        )

    if mode is RunMode.DESIGN_ONLY and "check_off_targets" not in stated and legacy_skip_screening is None:
        derived["check_off_targets"] = False
    if screening_off:
        derived["check_off_targets"] = False
    return mode, records, derived


def _evidence_requirements(
    *, run_mode: RunMode, query_species: str, screen_species: Sequence[str]
) -> EvidenceRequirements:
    """Declare which channel/species pairs must complete, and what an undecided filter costs.

    Qualified requires the transcriptome channel in the *query* species -- the one whose alignment
    decides on-target membership and post-screen scoring. Every other pair is exploratory: a
    secondary species that fails to align is reportable, not disqualifying. Nothing here enforces
    the requirement; eligibility against it is #100's.
    """
    if run_mode is RunMode.DESIGN_ONLY:
        return EvidenceRequirements(unknown_evidence_action=UnknownEvidenceAction.WARN)

    requiredness = Requiredness.REQUIRED if run_mode is RunMode.QUALIFIED else Requiredness.EXPLORATORY
    pairs: list[ChannelRequirement] = [
        ChannelRequirement(channel=ScreeningChannel.TRANSCRIPTOME, species=query_species, requiredness=requiredness)
    ]
    for species in dict.fromkeys(normalize_species_name(name) for name in screen_species):
        if species == query_species:
            continue
        pairs.append(
            ChannelRequirement(
                channel=ScreeningChannel.TRANSCRIPTOME, species=species, requiredness=Requiredness.EXPLORATORY
            )
        )
    pairs.append(
        ChannelRequirement(
            channel=ScreeningChannel.MIRNA_SEED, species=query_species, requiredness=Requiredness.EXPLORATORY
        )
    )
    return EvidenceRequirements(
        channel_requirements=tuple(pairs),
        unknown_evidence_action=(
            UnknownEvidenceAction.FAIL if run_mode is RunMode.QUALIFIED else UnknownEvidenceAction.WARN
        ),
    )


PASSTHROUGH_FIELDS: frozenset[str] = frozenset({"snp_file", "genome_index", "scoring"})


def _build_design_parameters(
    values: Mapping[str, Any], design_mode: DesignMode, passthrough: Mapping[str, Any]
) -> DesignParameters:
    """Construct the validated parameter tree, translating Pydantic errors into policy errors.

    Constructing rather than copying is deliberate: it is what keeps every bound and every
    cross-field rule (``gc_max >= gc_min``, ``max_bp_span <= window_size``) in force on a resolved
    run, and it happens before any caller has created a directory or downloaded a reference.

    ``passthrough`` carries the DesignParameters fields that are inputs rather than policy -- a SNP
    file path, a prebuilt index, an explicit weight set. They are not resolved and not overridable
    by a profile; they are forwarded so a resolved run does not silently drop them.
    """
    unknown = sorted(set(passthrough) - PASSTHROUGH_FIELDS)
    if unknown:
        raise RunPolicyError(f"these are not pass-through parameters: {unknown}")
    grouped: dict[str, dict[str, Any]] = {
        "filters": {},
        "offtarget_filters": {},
        "target_accessibility": {},
        "design": {},
    }
    for key, value in values.items():
        spec = SETTING_BY_KEY[key]
        grouped[spec.model][spec.field] = value
    try:
        return DesignParameters(
            design_mode=design_mode,
            filters=FilterCriteria(**grouped["filters"]),
            offtarget_filters=OffTargetFilterCriteria(**grouped["offtarget_filters"]),
            target_accessibility=TargetAccessibilityConfig(**grouped["target_accessibility"]),
            **{key: value for key, value in passthrough.items() if value is not None},
            **grouped["design"],
        )
    except ValidationError as exc:
        raise RunPolicyError(format_validation_error(exc)) from exc


def format_validation_error(exc: ValidationError) -> str:
    """Render a Pydantic error as one line per problem, without a traceback.

    Cross-field validators (``max_bp_span`` versus ``window_size``) previously escaped the CLI as an
    unhandled traceback because they raise at construction, not at option-parsing time.
    """
    lines: list[str] = []
    for error in exc.errors():
        location = ".".join(str(part) for part in error.get("loc", ())) or exc.title
        message = str(error.get("msg", "")).removeprefix("Value error, ")
        lines.append(f"{location}: {message}")
    return "; ".join(lines) or str(exc)


def _check_finite_weights(scoring: ScoringWeights) -> None:
    """A non-finite weight must fail before expensive work, not produce a NaN score at the end."""
    for vector in scoring.all_vectors():
        for term, weight in vector.as_mapping().items():
            if not math.isfinite(weight):
                raise RunPolicyError(f"weight vector '{vector.name}' has a non-finite weight for term '{term}'")


def _resolve_filters(
    *, parameters: DesignParameters, run_mode: RunMode, actions: Mapping[str, FilterAction]
) -> tuple[ResolvedFilter, ...]:
    """Attach each gate's resolved threshold and action to its static description.

    Three things turn a gate off, and they are different facts: the caller said ``off``; the gate
    has no declared threshold, so it cannot be evaluated; or the run holds no evidence of the kind
    the gate reads. All three produce ``FilterAction.OFF``, which is what makes a disabled filter
    ``NOT_EVALUATED`` rather than passed.
    """
    resolved: list[ResolvedFilter] = []
    for spec in FILTER_SPECS:
        setting = SETTING_BY_KEY[spec.setting_key]
        threshold = getattr(_owner(parameters, setting), setting.field)
        if isinstance(threshold, bool):
            # A boolean flag is a ceiling of zero: enabled means "at most zero hits", disabled is off.
            action = spec.default_action if threshold else FilterAction.OFF
            threshold = 0
        else:
            action = spec.default_action
        if threshold is None:
            action = FilterAction.OFF
        if spec.stage is FilterStage.POST_SCREEN and run_mode is RunMode.DESIGN_ONLY:
            action = FilterAction.OFF
        override = actions.get(spec.filter_id)
        if override is not None:
            if override is not FilterAction.OFF and action is FilterAction.OFF and threshold is None:
                raise RunPolicyError(
                    f"filter {spec.filter_id!r} was set to {override.value} but has no threshold; "
                    f"set {spec.setting_key} as well, or leave the filter off"
                )
            action = override
        resolved.append(
            ResolvedFilter(
                descriptor=FilterDescriptor(
                    filter_id=spec.filter_id,
                    column=spec.column,
                    comparator=spec.comparator,
                    threshold=threshold,
                    scope=FilterScope(
                        species=frozenset(spec.scope_species),
                        max_mismatches=spec.scope_max_mismatches,
                        hit_classes=frozenset(spec.scope_hit_classes),
                    ),
                    action=action,
                    stage=spec.stage,
                ),
                setting_key=spec.setting_key,
                definition=spec.definition,
                transform=spec.transform,
                evidence_exported=spec.evidence_exported,
            )
        )
    return tuple(resolved)


def resolve_run_policy(
    *,
    entry_point: EntryPoint,
    design_mode: DesignMode | str | None = None,
    run_mode: RunMode | str | None = None,
    profile_name: str = DEFAULT_PROFILE_NAME,
    config_file: Path | str | None = None,
    stated: Mapping[str, Any] | None = None,
    filter_actions: Mapping[str, Any] | None = None,
    passthrough: Mapping[str, Any] | None = None,
    legacy_skip_screening: bool | None = None,
    require_screening_completeness: bool | None = None,
    query_species: str | None = None,
    screen_species: Sequence[str] = (),
) -> ResolvedRunPolicy:
    """Resolve one run's policy. Pure: no downloads, no folding, no alignment.

    Args:
        entry_point: Which public surface is asking. Sets only the default run mode.
        design_mode: ``sirna``/``mirna``/``zfn``. Orthogonal to ``run_mode``: choosing miRNA design
            says nothing about how complete the screening evidence has to be.
        run_mode: Explicit run mode. Omitted, it follows the entry point and the screening switch.
        profile_name: Built-in versioned profile supplying the baseline.
        config_file: Optional JSON/TOML file, which beats the profile and loses to explicit values.
        stated: Settings the caller *actually stated* -- a field set, not a value comparison. Keys
            must be declared in :data:`SETTING_SPECS`; a value of ``None`` is a stated ``None``.
        filter_actions: ``filter_id -> off|warn|fail``.
        passthrough: DesignParameters fields that are inputs rather than policy (``snp_file``,
            ``genome_index``, an explicit ``scoring`` weight set). Forwarded, never defaulted.
        legacy_skip_screening: The legacy skip flag. True maps visibly to design-only.
        require_screening_completeness: Set False to keep incomplete evidence; incompatible with an
            explicit qualified mode, and on its own it resolves to exploratory.
        query_species: Species the target transcripts belong to; the pair that qualified mode
            requires. Defaults to human, as the rest of the pipeline does.
        screen_species: Species being screened against, recorded as exploratory pairs.

    Returns:
        The immutable resolved policy.

    Raises:
        RunPolicyError: For an undeclared setting, an invalid bound, a non-finite weight or an
            incompatible mode/settings combination -- all before any expensive work.
    """
    stated = dict(stated or {})
    _reject_undeclared(stated, "the caller")
    resolved_design_mode = parse_design_mode(design_mode) if design_mode is not None else DesignMode.SIRNA

    profile = BUILTIN_PROFILES.get(profile_name)
    if profile is None:
        raise RunPolicyError(f"unknown profile {profile_name!r}; available profiles: {sorted(BUILTIN_PROFILES)}")

    config_settings: Mapping[str, Any] = {}
    config_actions: Mapping[str, Any] = {}
    if config_file is not None:
        config_settings, directives = load_policy_config(config_file)
        config_actions = directives.get("filter_actions") or {}
        if run_mode is None and "run_mode" in directives:
            run_mode = directives["run_mode"]
        if design_mode is None and "design_mode" in directives:
            resolved_design_mode = parse_design_mode(directives["design_mode"])
        if "profile" in directives and directives["profile"] != profile.name:
            raise RunPolicyError(
                f"policy config file selects profile {directives['profile']!r}, which must be passed as "
                "profile_name so the caller and the file cannot disagree about the baseline"
            )

    resolved_run_mode, mode_records, mode_derived = _resolve_run_mode(
        entry_point=entry_point,
        requested_run_mode=run_mode,
        stated=stated,
        legacy_skip_screening=legacy_skip_screening,
        require_screening_completeness=require_screening_completeness,
    )

    preset = _mirna_preset() if resolved_design_mode is DesignMode.MIRNA else {}
    layers = (
        _Layer(SettingSource.BUILTIN_PROFILE, profile.baseline, f"{profile.name} profile {profile.version}"),
        _Layer(SettingSource.DESIGN_MODE_PRESET, preset, "miRNA design-mode preset" if preset else None),
        _Layer(SettingSource.CONFIG_FILE, config_settings, str(config_file) if config_file else None),
        _Layer(SettingSource.EXPLICIT, stated, "stated by the caller"),
        _Layer(SettingSource.RUN_MODE_RULE, mode_derived, f"implied by run_mode={resolved_run_mode.value}"),
    )

    values: dict[str, Any] = {}
    winners: dict[str, _Layer] = {}
    for layer in layers:
        for key, value in layer.values.items():
            values[key] = value
            winners[key] = layer

    # A gate the caller switched off must stop comparing, not merely record an intent: clearing its
    # threshold is what the existing gate code already reads as "no gate".
    actions = _resolve_actions({**config_actions, **(filter_actions or {})})
    for filter_id, action in actions.items():
        if action is not FilterAction.OFF:
            continue
        spec = SETTING_BY_KEY[FILTER_SPEC_BY_ID[filter_id].setting_key]
        values[spec.key] = _switch_off_value(spec)
        winners[spec.key] = _Layer(SettingSource.EXPLICIT, {}, f"cleared because filter {filter_id} was switched off")

    parameters = _build_design_parameters(values, resolved_design_mode, passthrough or {})
    _check_finite_weights(parameters.scoring)

    filters = _resolve_filters(parameters=parameters, run_mode=resolved_run_mode, actions=actions)

    requested: list[SettingProvenance] = [
        SettingProvenance(key=key, value=value, source=SettingSource.CONFIG_FILE, detail=str(config_file))
        for key, value in config_settings.items()
    ]
    requested += [
        SettingProvenance(key=key, value=value, source=SettingSource.EXPLICIT, detail=None)
        for key, value in stated.items()
    ]
    requested += [
        SettingProvenance(key=f"action.{filter_id}", value=action.value, source=SettingSource.EXPLICIT, detail=None)
        for filter_id, action in actions.items()
    ]

    resolved_records: list[SettingProvenance] = list(mode_records)
    resolved_records.append(
        SettingProvenance(
            key="design_mode",
            value=resolved_design_mode.value,
            source=SettingSource.EXPLICIT if design_mode is not None else SettingSource.BUILTIN_PROFILE,
            detail=None if design_mode is not None else "no design mode stated",
        )
    )
    for spec in SETTING_SPECS:
        layer = winners[spec.key]
        detail = layer.detail
        if layer.source is SettingSource.BUILTIN_PROFILE and spec.key in profile.exceptions:
            detail = f"{detail}; declared exception to the model field default ({profile.exceptions[spec.key]})"
        resolved_records.append(
            SettingProvenance(key=spec.key, value=values[spec.key], source=layer.source, detail=detail)
        )
    for resolved_filter in filters:
        resolved_records.append(
            SettingProvenance(
                key=f"action.{resolved_filter.filter_id}",
                value=resolved_filter.descriptor.action.value,
                source=(
                    SettingSource.EXPLICIT
                    if resolved_filter.filter_id in actions
                    else SettingSource.RUN_MODE_RULE
                    if resolved_filter.descriptor.stage is FilterStage.POST_SCREEN
                    and resolved_run_mode is RunMode.DESIGN_ONLY
                    else SettingSource.BUILTIN_PROFILE
                ),
                detail=(
                    "post-screen gate in a design-only run: not evaluated"
                    if resolved_filter.descriptor.stage is FilterStage.POST_SCREEN
                    and resolved_run_mode is RunMode.DESIGN_ONLY
                    else None
                ),
            )
        )

    return ResolvedRunPolicy(
        entry_point=entry_point,
        run_mode=resolved_run_mode,
        design_mode=resolved_design_mode,
        profile=profile.identity(),
        design_parameters=parameters,
        filters=filters,
        evidence_requirements=_evidence_requirements(
            run_mode=resolved_run_mode,
            query_species=normalize_species_name(query_species) if query_species else DEFAULT_TARGET_SPECIES,
            screen_species=screen_species,
        ),
        requested=tuple(requested),
        resolved=tuple(resolved_records),
    )


def describe_parameters(
    parameters: DesignParameters,
    *,
    entry_point: EntryPoint,
    run_mode: RunMode | str | None = None,
    query_species: str | None = None,
    screen_species: Sequence[str] = (),
) -> ResolvedRunPolicy:
    """Wrap an already-built ``DesignParameters`` in a resolved policy without changing a number.

    The adapter for callers that construct ``DesignParameters`` themselves -- direct
    ``WorkflowConfig`` users and older keyword-argument APIs. Every setting is reported as stated by
    the caller, because handing over a fully specified parameter object *is* stating it; nothing is
    defaulted, so an existing run resolves to exactly the numbers it already had.
    """
    stated: dict[str, Any] = {}
    for spec in SETTING_SPECS:
        stated[spec.key] = getattr(_owner(parameters, spec), spec.field)
    return resolve_run_policy(
        entry_point=entry_point,
        design_mode=parameters.design_mode,
        run_mode=run_mode,
        stated=stated,
        passthrough={
            "snp_file": parameters.snp_file,
            "genome_index": parameters.genome_index,
            "scoring": parameters.scoring,
        },
        query_species=query_species,
        screen_species=screen_species,
    )


def default_for(key: str) -> Any:
    """The ``legacy`` profile's value for one setting, for CLI help and documentation tables.

    Help text and parameter tables read this instead of restating a number, which is what stopped
    the documented off-target cap (3) and the documented design weights (0.40/0.35) from drifting
    away from the code (15, and 0.35/0.40) unnoticed.
    """
    if key not in SETTING_BY_KEY:
        raise KeyError(f"no setting named {key!r}")
    return BUILTIN_PROFILES[DEFAULT_PROFILE_NAME].baseline[key]


def mirna_preset_default_for(key: str) -> Any:
    """The miRNA design-mode preset's value for one setting, or the profile's when it moves nothing.

    CLI help for an option the preset moves must name both numbers, or a miRNA user reads a default
    that never applies to them.
    """
    if key not in SETTING_BY_KEY:
        raise KeyError(f"no setting named {key!r}")
    return _mirna_preset().get(key, default_for(key))


def declared_filter_ids() -> tuple[str, ...]:
    """Every gate this build declares, in reporting order."""
    return tuple(spec.filter_id for spec in FILTER_SPECS)


def switchable_filter_ids() -> tuple[str, ...]:
    """Gates that can be switched off, which is those whose threshold has an absent state.

    The other six read a plain float with no "no threshold" value, so switching them off would need
    the gate application itself to read the action -- which it does not do in 0.7.1.
    """
    switchable: list[str] = []
    for spec in FILTER_SPECS:
        try:
            _switch_off_value(SETTING_BY_KEY[spec.setting_key])
        except RunPolicyError:
            continue
        switchable.append(spec.filter_id)
    return tuple(switchable)


__all__ = [
    "BUILTIN_PROFILES",
    "DEFAULT_PROFILE_NAME",
    "EntryPoint",
    "FILTER_SPECS",
    "LEGACY_PROFILE_NAME",
    "POLICY_SCHEMA_VERSION",
    "ResolvedFilter",
    "ResolvedRunPolicy",
    "RunPolicyError",
    "RunPolicyProfile",
    "SETTING_SPECS",
    "SettingSpec",
    "declared_filter_ids",
    "switchable_filter_ids",
    "describe_parameters",
    "default_for",
    "format_validation_error",
    "load_policy_config",
    "mirna_preset_default_for",
    "parse_design_mode",
    "resolve_run_policy",
]
