"""Screening plan and screening evidence contracts.

Types only, per #100: no producers, no consumers, no reconciliation. The plan says what a run
intended to screen; the evidence says what each channel × species × reference actually did. They
are separate objects because "we never asked" and "we asked and it failed" must not collapse into
the same shape, and because an absent output file is not evidence of completion.

The counts model exists so an unknown can never be published as a zero: ``hits_retained`` is
``None`` until something observed it, and a retained count that was capped or truncated says so.
"""

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field

from sirnaforge.models.policy import ScreeningChannel

EVIDENCE_SCHEMA_VERSION = "1"
"""Version of the plan/evidence payload shape. Bumping it is a visible, breaking change."""


class EvidenceStatus(str, Enum):
    """Outcome of one planned screening unit.

    Attributes:
        COMPLETE: The configured search ran to completion. This is a statement about the search,
            not a claim that every biologically possible hit was found.
        FAILED: The search was requested and did not complete (aligner error, missing reference).
        NOT_REQUESTED: The search was never asked for, so there is nothing to reconcile.
        CENSORED: The search ran but its results were withheld or truncated to the point that the
            observed counts are a lower bound rather than a total.
    """

    COMPLETE = "complete"
    FAILED = "failed"
    NOT_REQUESTED = "not_requested"
    CENSORED = "censored"


class ObservedCounts(BaseModel):
    """What a screening unit observed, with the flags that stop an unknown becoming a zero.

    Attributes:
        hits_retained: Hits kept after the channel's own filtering. ``None`` means unobserved --
            never write 0 for "we do not know".
        is_lower_bound: True when more hits may exist than were retained, for any reason.
        hit_cap: The per-query retention cap in force, if any.
        truncated: True when the cap or another limit actually discarded hits.
    """

    hits_retained: int | None = Field(default=None, ge=0, description="Retained hit count; None means unobserved")
    is_lower_bound: bool = Field(default=False, description="True when the retained count may understate the total")
    hit_cap: int | None = Field(default=None, ge=0, description="Per-query retention cap in force, if any")
    truncated: bool = Field(default=False, description="True when hits were actually discarded by a cap or limit")

    model_config = ConfigDict(frozen=True)


class ScreeningPlanEntry(BaseModel):
    """One channel × species × reference the run intended to screen.

    ``search_settings`` is a mapping, so the model itself is not hashable; ``key`` is the hashable
    identity to index or join on.

    Attributes:
        channel: Which liability channel this entry covers.
        species: Canonical species name, normalized by the caller before construction.
        reference_id: Identity of the reference searched (index name, database release). ``None``
            when the channel has no reference identity to report.
        required: True when this entry's completion is required rather than exploratory.
        guide_set_digest: Digest of the submitted guide set, so evidence can be tied to the guides
            it was actually produced from.
        search_settings: The search settings relevant to reproducing this entry.
    """

    channel: ScreeningChannel = Field(description="Liability channel")
    species: str = Field(description="Canonical species name")
    reference_id: str | None = Field(default=None, description="Reference identity searched")
    required: bool = Field(default=False, description="True when completion is required, not exploratory")
    guide_set_digest: str | None = Field(default=None, description="Digest of the submitted guide set")
    search_settings: dict[str, str | int | float | bool | None] = Field(
        default_factory=dict, description="Search settings relevant to this entry"
    )

    model_config = ConfigDict(frozen=True)

    @property
    def key(self) -> tuple[str, str, str | None]:
        """Identity shared with the matching evidence entry: channel, species, reference."""
        return (self.channel.value, self.species, self.reference_id)


class ScreeningPlan(BaseModel):
    """What a run intended to screen, recorded before any reference can be dropped.

    Attributes:
        schema_version: Payload shape version.
        entries: The planned units, in declaration order.
    """

    schema_version: str = Field(default=EVIDENCE_SCHEMA_VERSION, description="Payload shape version")
    entries: tuple[ScreeningPlanEntry, ...] = Field(default=(), description="Planned screening units")

    model_config = ConfigDict(frozen=True)


class ScreeningEvidenceEntry(BaseModel):
    """What one planned screening unit actually did.

    Attributes:
        channel: Liability channel, matching the plan entry.
        species: Canonical species name, matching the plan entry.
        reference_id: Reference identity, matching the plan entry.
        required: True when this entry's completion is required rather than exploratory.
        status: Completion outcome.
        counts: Observed counts and their lower-bound/cap flags.
        submitted_guide_digest: Digest of the guide set actually submitted, for comparison against
            the plan's digest.
        submitted_guides: Number of guides submitted, or ``None`` when unobserved.
        processed_guides: Number of guides the channel reported processing, or ``None``.
        detail: Human-readable reason, required for FAILED and CENSORED to be actionable.
    """

    channel: ScreeningChannel = Field(description="Liability channel")
    species: str = Field(description="Canonical species name")
    reference_id: str | None = Field(default=None, description="Reference identity searched")
    required: bool = Field(default=False, description="True when completion is required, not exploratory")
    status: EvidenceStatus = Field(description="Completion outcome")
    counts: ObservedCounts = Field(default_factory=ObservedCounts, description="Observed counts and their flags")
    submitted_guide_digest: str | None = Field(default=None, description="Digest of the submitted guide set")
    submitted_guides: int | None = Field(default=None, ge=0, description="Guides submitted; None means unobserved")
    processed_guides: int | None = Field(default=None, ge=0, description="Guides processed; None means unobserved")
    detail: str | None = Field(default=None, description="Reason for a failed, censored or absent result")

    model_config = ConfigDict(frozen=True)

    @property
    def key(self) -> tuple[str, str, str | None]:
        """Identity shared with the matching plan entry: channel, species, reference."""
        return (self.channel.value, self.species, self.reference_id)


class ScreeningEvidence(BaseModel):
    """What a run actually screened, to be reconciled against a plan by its consumer.

    Attributes:
        schema_version: Payload shape version.
        entries: The observed units, one per plan entry the producer could speak to.
    """

    schema_version: str = Field(default=EVIDENCE_SCHEMA_VERSION, description="Payload shape version")
    entries: tuple[ScreeningEvidenceEntry, ...] = Field(default=(), description="Observed screening units")

    model_config = ConfigDict(frozen=True)
