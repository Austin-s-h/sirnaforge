"""Screening plan and screening evidence contracts.

Types only, per #100: no producers, no consumers, no reconciliation. The plan says what a run
intended to screen; the evidence says what each channel × species × reference actually did. They
are separate objects because "we never asked" and "we asked and it failed" must not collapse into
the same shape, and because an absent output file is not evidence of completion.

Owned by #100; extended in place by #101 (per-count cap semantics, the count matrix). Requiredness
is deliberately absent from every entry here: :class:`~sirnaforge.models.policy.EvidenceRequirements`
is its single authority, and an entry that could restate it could contradict it invisibly.

The counts model exists so an unknown can never be published as a zero: every count is ``None``
until something observed it, and a count that was capped or truncated says so on its own rather
than through one flag shared across three quantities.
"""

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, model_validator

from sirnaforge.models.policy import FilterScope, ScreeningChannel

EVIDENCE_SCHEMA_VERSION = "2"
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


class ObservedCount(BaseModel):
    """One observed quantity with its own cap and truncation semantics.

    Attributes:
        value: The count observed. ``None`` means unobserved -- never write 0 for "we do not know".
        is_lower_bound: True when more may exist than were counted, for any reason.
        cap: The retention or reporting cap in force for this quantity, if any.
        truncated: True when the cap or another limit actually discarded members.
    """

    value: int | None = Field(default=None, ge=0, description="Observed count; None means unobserved")
    is_lower_bound: bool = Field(default=False, description="True when the count may understate the total")
    cap: int | None = Field(default=None, ge=0, description="Cap in force for this quantity, if any")
    truncated: bool = Field(default=False, description="True when a cap or limit actually discarded members")

    model_config = ConfigDict(frozen=True, extra="forbid")


class ObservedCounts(BaseModel):
    """What a screening unit observed, at the three aggregation units that are reported separately.

    Sites, distinct transcripts and distinct genes have their own caps and their own truncation
    behaviour: a per-query alignment cap truncates sites, while a gene-level cap over unresolved
    transcripts truncates genes. One shared flag cannot describe both, and collapsing them is how
    a gene-level limit silently drops the transcripts it could not resolve.

    Attributes:
        sites: Retained alignments or seed sites -- the hit-level count for alignment channels.
        distinct_transcripts: Distinct transcripts carrying at least one retained site.
        distinct_genes: Distinct genes carrying at least one retained site.
        unresolved_gene_sites: Retained sites whose transcript could not be mapped to a gene, so a
            gene-level count is a lower bound rather than a total.
    """

    sites: ObservedCount = Field(default_factory=ObservedCount, description="Retained alignments or seed sites")
    distinct_transcripts: ObservedCount = Field(
        default_factory=ObservedCount, description="Distinct transcripts with a retained site"
    )
    distinct_genes: ObservedCount = Field(
        default_factory=ObservedCount, description="Distinct genes with a retained site"
    )
    unresolved_gene_sites: ObservedCount = Field(
        default_factory=ObservedCount, description="Retained sites with no resolvable gene"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")


class HitCountCell(BaseModel):
    """One cell of a species × mismatch-class × hit-class count.

    Attributes:
        species: Canonical species name.
        mismatches: Guide-level mismatch class this cell counts.
        hit_class: ``HitClass`` value (``core/hit_classification.py``), held as a string so this
            contract stays independent of the classifier's membership.
        sites: Alignments or seed sites in this cell.
    """

    species: str = Field(min_length=1, description="Canonical species name")
    mismatches: int = Field(ge=0, description="Guide-level mismatch class")
    hit_class: str = Field(min_length=1, description="HitClass value this cell counts")
    sites: int = Field(ge=0, description="Alignments or seed sites in this cell")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def key(self) -> tuple[str, int, str]:
        """The cell's coordinates."""
        return (self.species, self.mismatches, self.hit_class)


class HitCountMatrix(BaseModel):
    """A complete count over a declared scope, so an absent cell is a real zero.

    ``scope`` is what makes the matrix readable: inside it an absent cell means zero, and outside
    it the matrix says nothing at all. A consumer that embeds hit rows for part of the scope and
    counts for the rest needs exactly this distinction -- a guide with liabilities beyond the
    embedded rows must not render like a clean one.

    Attributes:
        scope: The species, mismatch ceiling and hit classes this matrix is complete over.
        cells: The non-zero cells, in any order.
    """

    scope: FilterScope = Field(description="Scope the matrix is complete over; absent cells inside it are zero")
    cells: tuple[HitCountCell, ...] = Field(default=(), description="Non-zero cells")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def cells_are_unique(self) -> "HitCountMatrix":
        """Two cells at one coordinate would make the total depend on iteration order."""
        keys = [cell.key for cell in self.cells]
        if len(keys) != len(set(keys)):
            duplicates = sorted({key for key in keys if keys.count(key) > 1})
            raise ValueError(f"duplicate hit-count cells: {duplicates}")
        return self


class ScreeningPlanEntry(BaseModel):
    """One channel × species × reference × guide-set the run intended to screen.

    ``search_settings`` is a mapping, so the model itself is not hashable; ``key`` is the hashable
    identity to index or join on. The guide-set digest is part of that identity because two
    screens of the same reference with different guide sets are different evidence, and joining
    them would attribute one screen's counts to the other's guides.

    Attributes:
        channel: Which liability channel this entry covers.
        species: Canonical species name, normalized by the caller before construction.
        reference_id: Identity of the reference searched (index name, database release). ``None``
            when the channel has no reference identity to report.
        guide_set_digest: Digest of the guide set this entry plans to submit. No default: it is
            half the join identity, and a ``None`` identity collides with every other ``None``.
        search_settings: The search settings relevant to reproducing this entry.
    """

    channel: ScreeningChannel = Field(description="Liability channel")
    species: str = Field(min_length=1, description="Canonical species name")
    reference_id: str | None = Field(default=None, description="Reference identity searched")
    guide_set_digest: str = Field(min_length=1, description="Digest of the guide set planned for submission")
    search_settings: dict[str, str | int | float | bool | None] = Field(
        default_factory=dict, description="Search settings relevant to this entry"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def key(self) -> tuple[str, str, str | None, str]:
        """Identity shared with the matching evidence entry: channel, species, reference, guide set."""
        return (self.channel.value, self.species, self.reference_id, self.guide_set_digest)


class ScreeningPlan(BaseModel):
    """What a run intended to screen, recorded before any reference can be dropped.

    Attributes:
        schema_version: Payload shape version.
        entries: The planned units, in declaration order.
    """

    schema_version: str = Field(default=EVIDENCE_SCHEMA_VERSION, description="Payload shape version")
    entries: tuple[ScreeningPlanEntry, ...] = Field(default=(), description="Planned screening units")

    model_config = ConfigDict(frozen=True, extra="forbid")


class ScreeningEvidenceEntry(BaseModel):
    """What one planned screening unit actually did.

    Attributes:
        channel: Liability channel, matching the plan entry.
        species: Canonical species name, matching the plan entry.
        reference_id: Reference identity, matching the plan entry.
        guide_set_digest: Digest of the guide set this entry answers for, matching the plan entry.
        status: Completion outcome.
        counts: Observed counts at each reported aggregation unit.
        submitted_guide_digest: Digest of the guide set actually submitted. Differing from
            ``guide_set_digest`` is a real defect and is representable so it can be detected.
        submitted_guides: Number of guides submitted, or ``None`` when unobserved.
        processed_guides: Number of guides the channel reported processing, or ``None``.
        detail: Human-readable reason, required for FAILED and CENSORED to be actionable.
    """

    channel: ScreeningChannel = Field(description="Liability channel")
    species: str = Field(min_length=1, description="Canonical species name")
    reference_id: str | None = Field(default=None, description="Reference identity searched")
    guide_set_digest: str = Field(min_length=1, description="Digest of the guide set this entry answers for")
    status: EvidenceStatus = Field(description="Completion outcome")
    counts: ObservedCounts = Field(default_factory=ObservedCounts, description="Observed counts and their flags")
    submitted_guide_digest: str | None = Field(default=None, description="Digest of the guide set actually submitted")
    submitted_guides: int | None = Field(default=None, ge=0, description="Guides submitted; None means unobserved")
    processed_guides: int | None = Field(default=None, ge=0, description="Guides processed; None means unobserved")
    detail: str | None = Field(default=None, description="Reason for a failed, censored or absent result")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @property
    def key(self) -> tuple[str, str, str | None, str]:
        """Identity shared with the matching plan entry: channel, species, reference, guide set."""
        return (self.channel.value, self.species, self.reference_id, self.guide_set_digest)


class ScreeningEvidence(BaseModel):
    """What a run actually screened, to be reconciled against a plan by its consumer.

    Attributes:
        schema_version: Payload shape version.
        entries: The observed units, one per plan entry the producer could speak to.
    """

    schema_version: str = Field(default=EVIDENCE_SCHEMA_VERSION, description="Payload shape version")
    entries: tuple[ScreeningEvidenceEntry, ...] = Field(default=(), description="Observed screening units")

    model_config = ConfigDict(frozen=True, extra="forbid")
