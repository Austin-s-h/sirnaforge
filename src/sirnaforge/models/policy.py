"""Shared run-policy vocabulary: run modes, filter actions, evaluations and target intent.

Data contracts only. This module deliberately imports nothing from ``workflow`` and resolves,
defaults and applies nothing: the resolver lives in ``config/run_policy.py`` (#99), the evidence
producers in ``core``/Nextflow (#100), and intent evaluation in the classification path (#101).
Its only job is to give those packages one spelling for each concept.
"""

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field

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


class FilterEvaluation(str, Enum):
    """The verdict itself, which is distinct from the action taken on it.

    A filter with ``FilterAction.OFF`` produces no evaluation at all rather than ``PASS``:
    "not evaluated" and "evaluated and passed" are different claims. ``UNKNOWN`` is the third
    outcome, for a filter whose evidence is missing or censored; it must never be reported as
    ``PASS``.

    Attributes:
        PASS: Evidence was available and the candidate met the threshold.
        FAIL: Evidence was available and the candidate did not meet the threshold.
        UNKNOWN: The filter could not be decided because its evidence is unavailable.
    """

    PASS = "pass"
    FAIL = "fail"
    UNKNOWN = "unknown"


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


class TargetIntent(BaseModel):
    """What the user intends to knock down, and where they intend to look for liabilities.

    The two species sets are deliberately separate. A species may be a target, a screening
    subject, both, or neither: screening mouse does not make mouse knockdown a requirement, and
    requiring human knockdown does not by itself request a human screen. Callers that need one
    set to imply the other must say so explicitly.

    Nothing here is resolved or defaulted beyond ``target_species``, whose documented default is
    human alone. Empty sets mean "not declared", and the consumer decides what that implies.

    Attributes:
        target_species: Canonical species names whose transcripts are meant to be knocked down.
        offtarget_screen_species: Canonical species names to screen for off-target liabilities.
        required_transcript_ids: Transcripts a guide must cover, versionless.
        excluded_transcript_ids: Transcripts a guide must not cover, versionless. A hit on one of
            these is on-target by gene taxonomy and still unacceptable.
    """

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

    model_config = ConfigDict(frozen=True)


class EvidenceRequirements(BaseModel):
    """Which screening evidence a run requires, and how an undecided filter is treated.

    ``unknown_evidence_action`` has no default on purpose: choosing what an ``UNKNOWN``
    evaluation does to a candidate is a policy decision owned by the resolver, not a fallback
    this contract is entitled to pick.

    Attributes:
        required_channels: Channels whose evidence must be complete.
        exploratory_channels: Channels whose evidence is descriptive; incompleteness is reportable
            but not disqualifying.
        required_species: Species whose evidence must be complete.
        exploratory_species: Species whose evidence is descriptive.
        unknown_evidence_action: What to do with a filter that evaluated to ``UNKNOWN``.
    """

    required_channels: frozenset[ScreeningChannel] = Field(
        default=frozenset(), description="Channels whose evidence must be complete"
    )
    exploratory_channels: frozenset[ScreeningChannel] = Field(
        default=frozenset(), description="Channels whose evidence is descriptive only"
    )
    required_species: frozenset[str] = Field(default=frozenset(), description="Species whose evidence must be complete")
    exploratory_species: frozenset[str] = Field(
        default=frozenset(), description="Species whose evidence is descriptive only"
    )
    unknown_evidence_action: FilterAction = Field(
        description="Action applied when a filter's evaluation is UNKNOWN; the caller must decide"
    )

    model_config = ConfigDict(frozen=True)
