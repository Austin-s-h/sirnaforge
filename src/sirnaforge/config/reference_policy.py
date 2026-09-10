"""Reference resolution for the screening path: one resolver, one typed reference.

Screening vocabulary is ``transcriptome`` throughout (#99). ``genome`` survives only as the ZFN
:class:`ReferenceKind`, where it means genomic DNA: siRNA and miRNA act on mRNA, so the
transcriptome is the correct reference for both on- and off-target, and ``kind`` is keyed to
modality and nothing else.

Two rules this module exists to enforce:

- **A reference's species comes from the resolved reference, never from the parameter it arrived
  on.** :class:`SpeciesAuthority` records which authority supplied it, in the order declared >
  bundled-source registry > the reference's own headers > unresolved.
- **One door.** An explicit index override and a resolved default are both
  :class:`ReferenceRequest` values that differ in ``state``/``reason``/``form`` -- their
  provenance -- and resolve to the same :class:`ScreeningReference` shape, so neither can skip the
  metadata the classifier needs.

Imports of ``sirnaforge.models`` and ``sirnaforge.data`` are deliberately lazy: both packages'
``__init__`` reach back into this module, so a module-level import would be circular.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sirnaforge.models.evidence import ScreeningPlan, ScreeningPlanEntry
    from sirnaforge.models.policy import FilterScope
    from sirnaforge.models.sirna import DesignMode

DEFAULT_TRANSCRIPTOME_SOURCES: tuple[str, ...] = (
    "ensembl_human_cdna",
    "ensembl_mouse_cdna",
    "ensembl_rat_cdna",
    "ensembl_macaque_cdna",
)
DEFAULT_TRANSCRIPTOME_SOURCE = DEFAULT_TRANSCRIPTOME_SOURCES[0]

DEFAULT_MIRNA_SOURCE = "mirgenedb"
DEFAULT_MIRNA_CANONICAL_SPECIES: tuple[str, ...] = (
    "chicken",
    "pig",
    "rat",
    "mouse",
    "human",
    "rhesus",
    "macaque",
)


class ReferenceState(str, Enum):
    """Describe how a reference input was selected."""

    EXPLICIT = "explicit"
    DEFAULT = "default"
    DISABLED = "disabled"


@dataclass(frozen=True)
class ReferenceChoice:
    """Normalized representation of a resolved reference input."""

    value: str | None
    state: ReferenceState
    reason: str

    @property
    def enabled(self) -> bool:
        """Return True when a usable reference has been selected."""
        return self.value is not None and self.state is not ReferenceState.DISABLED

    @staticmethod
    def explicit(value: str, reason: str = "user-provided") -> ReferenceChoice:
        """Create an explicit user-selected reference choice."""
        return ReferenceChoice(value=value, state=ReferenceState.EXPLICIT, reason=reason)

    @staticmethod
    def default(value: str, reason: str = "auto-selected") -> ReferenceChoice:
        """Create a default-sourced reference choice."""
        return ReferenceChoice(value=value, state=ReferenceState.DEFAULT, reason=reason)

    @staticmethod
    def disabled(reason: str) -> ReferenceChoice:
        """Create a disabled reference choice with context."""
        return ReferenceChoice(value=None, state=ReferenceState.DISABLED, reason=reason)

    def to_metadata(self) -> dict[str, str | None | bool]:
        """Return a serializable snapshot for logs/JSON summaries."""
        return {
            "value": self.value,
            "state": self.state.value,
            "reason": self.reason,
            "enabled": self.enabled,
        }


@dataclass(frozen=True)
class ReferenceSelection:
    """Container describing zero or more resolved references."""

    choices: tuple[ReferenceChoice, ...] = field(default_factory=tuple)
    disabled_reason: str | None = None

    @property
    def enabled(self) -> bool:
        """Return True when at least one reference is configured."""
        return bool(self.choices)

    @staticmethod
    def disabled(reason: str) -> ReferenceSelection:
        """Create a disabled selection with a descriptive reason."""
        return ReferenceSelection(choices=(), disabled_reason=reason)

    def to_metadata(self) -> dict[str, object]:
        """Render selection metadata for logging."""
        return {
            "enabled": self.enabled,
            "disabled_reason": self.disabled_reason,
            "choices": [choice.to_metadata() for choice in self.choices],
        }


@dataclass(frozen=True)
class WorkflowInputSpec:
    """Raw workflow inputs prior to policy resolution."""

    input_fasta: str | None = None
    transcriptome_argument: str | None = None
    default_transcriptomes: Sequence[str] = field(default_factory=lambda: DEFAULT_TRANSCRIPTOME_SOURCES)
    design_only: bool = False
    # Stays False by default: resolving the defaults for a caller who supplied their own FASTA
    # downloads and indexes multi-gigabyte cDNA references they never requested. Callers that
    # want that must opt in, and naming a reference via transcriptome_argument always wins.
    allow_transcriptome_for_input_fasta: bool = False


class ReferencePolicyResolver:
    """Resolve workflow defaults while preserving intent metadata."""

    def __init__(self, spec: WorkflowInputSpec):
        """Create a resolver for a specific workflow input specification."""
        self.spec = spec

    def resolve_transcriptomes(self) -> ReferenceSelection:
        """Return one or more transcriptome references."""
        if self.spec.design_only:
            return ReferenceSelection.disabled("design-only mode requested")

        transcription_arg = (self.spec.transcriptome_argument or "").strip()
        if transcription_arg:
            choice = ReferenceChoice.explicit(transcription_arg, reason="explicit transcriptome override")
            return ReferenceSelection(choices=(choice,))

        if self.spec.input_fasta and not self.spec.allow_transcriptome_for_input_fasta:
            return ReferenceSelection.disabled(
                "input FASTA without an explicit transcriptome reference: design-only mode"
            )

        defaults = tuple(ref.strip() for ref in self.spec.default_transcriptomes if ref and ref.strip())
        if defaults:
            choices = tuple(ReferenceChoice.default(ref, reason="auto transcriptome default") for ref in defaults)
            return ReferenceSelection(choices=choices)

        return ReferenceSelection.disabled("no transcriptome defaults configured")


class ReferenceKind(str, Enum):
    """What a screening reference contains, keyed to modality and nothing else.

    Attributes:
        TRANSCRIPTOME: mRNA/cDNA sequence. The reference for siRNA and miRNA screening, on- and
            off-target alike.
        GENOME: Genomic DNA. ZFN only -- ZFNs bind DNA. No genomic screening mode exists for siRNA.
    """

    TRANSCRIPTOME = "transcriptome"
    GENOME = "genome"


class ReferenceForm(str, Enum):
    """How a requested reference arrives, which decides how it is materialised.

    Not a second door: both forms resolve to one :class:`ScreeningReference` and both build the
    species label and the transcript index. The form only says what work materialising it takes.

    Attributes:
        REFERENCE: A bundled source name, local path or URL to fetch, cache and index.
        PREBUILT_INDEX: A ``species:index_prefix`` the caller has already built.
    """

    REFERENCE = "reference"
    PREBUILT_INDEX = "prebuilt_index"


class SpeciesAuthority(str, Enum):
    """Which authority supplied a reference's species label.

    Ordered strongest first. Header inference is a *fallback*, not the authority: a caller who
    names the species means it, and only silence hands the decision to the file.

    Attributes:
        DECLARED: The caller stated it (``species:`` prefix on a reference or index entry).
        SOURCE_REGISTRY: A bundled source's own registry entry stated it.
        REFERENCE_HEADERS: Inferred from the reference's own Ensembl headers.
        UNRESOLVED: Nothing said. The label is :data:`UNRESOLVED_SPECIES`, which matches no species,
            so hits against it cannot be qualified.
    """

    DECLARED = "declared"
    SOURCE_REGISTRY = "source_registry"
    REFERENCE_HEADERS = "reference_headers"
    UNRESOLVED = "unresolved"


#: Search settings carried on a plan entry, matching ``ScreeningPlanEntry.search_settings``.
SearchSettings = dict[str, str | int | float | bool | None]

#: Label for a reference whose species nothing could establish. Deliberately not "transcriptome":
#: that is the *kind*, and using it as a species is the defect this module exists to close.
UNRESOLVED_SPECIES = "unknown"


class ReferenceKindError(ValueError):
    """A reference's kind does not match the modality's, raised before any expensive work."""


@dataclass(frozen=True)
class ReferenceRequest:
    """One screening reference a run asked for, before anything is fetched.

    Attributes:
        value: Bundled source name, path, URL or index prefix.
        form: How it is materialised.
        state: Provenance -- explicit caller override or resolved default.
        reason: Why this request exists, carried into the resolved reference.
        declared_species: Species the caller named for it, if any.
    """

    value: str
    form: ReferenceForm
    state: ReferenceState
    reason: str
    declared_species: str | None = None

    @property
    def is_explicit(self) -> bool:
        """Whether the caller asked for this reference by name."""
        return self.state is ReferenceState.EXPLICIT


@dataclass(frozen=True)
class ReferenceRejection:
    """A requested reference that could not be used, and why.

    Kept beside the resolved references rather than dropped: a species with no usable reference is
    a completeness fact about the run, and erasing it is how "unscreened" came to read as "clean".

    Attributes:
        species: Species the request was for, or :data:`UNRESOLVED_SPECIES` when even that is unknown.
        identity: The reference that failed.
        reason: What went wrong, in terms a user can act on.
    """

    species: str
    identity: str
    reason: str


@dataclass(frozen=True)
class ScreeningReference:
    """One resolved screening reference: ``{species, kind, identity, index}`` plus its provenance.

    Attributes:
        species: Canonical species, taken from the reference itself (see :class:`SpeciesAuthority`).
        kind: Transcriptome or genome, keyed to the run's modality.
        identity: What was resolved -- the source name, path, URL or index prefix as requested.
        index: What the aligner is handed: an index prefix, or the FASTA to index when
            ``needs_index_build`` is set.
        species_authority: Which authority supplied ``species``.
        form: How the reference was materialised.
        state: Provenance -- explicit override or resolved default.
        reason: Why this reference is in the run.
        fasta: The sequence file the classifier's transcript index was built from, when there is one.
        needs_index_build: True when no index exists yet, so the pipeline must build one. It decides
            which pipeline parameter carries this reference: a FASTA named as an index prefix aligns
            nothing, which the pipeline then reports as a completed screen.
    """

    species: str
    kind: ReferenceKind
    identity: str
    index: str
    species_authority: SpeciesAuthority
    form: ReferenceForm
    state: ReferenceState
    reason: str
    fasta: str | None = None
    needs_index_build: bool = False

    @property
    def index_entry(self) -> str:
        """The ``species:path`` token the pipeline parameter is built from."""
        return f"{self.species}:{self.index}"

    def to_metadata(self) -> dict[str, str | None]:
        """Serializable snapshot for the published reference summary."""
        return {
            "species": self.species,
            "kind": self.kind.value,
            "identity": self.identity,
            "index": self.index,
            "species_authority": self.species_authority.value,
            "form": self.form.value,
            "state": self.state.value,
            "reason": self.reason,
            "fasta": self.fasta,
            "needs_index_build": str(self.needs_index_build),
        }

    def plan_entry(self, *, guide_set_digest: str, search_settings: SearchSettings | None = None) -> ScreeningPlanEntry:
        """The :class:`~sirnaforge.models.evidence.ScreeningPlanEntry` this reference plans to screen.

        The channel is always ``TRANSCRIPTOME``: only a transcriptome reference ever resolves, since a
        genome-kind run refuses every screening reference before it is materialised.
        """
        from sirnaforge.models.evidence import ScreeningPlanEntry  # noqa: PLC0415
        from sirnaforge.models.policy import ScreeningChannel  # noqa: PLC0415

        return ScreeningPlanEntry(
            channel=ScreeningChannel.TRANSCRIPTOME,
            species=self.species,
            reference_id=self.identity,
            guide_set_digest=guide_set_digest,
            search_settings=dict(search_settings or {}),
        )


@dataclass(frozen=True)
class ScreeningReferenceSet:
    """Every reference one run resolved, the requests it could not use, and the scope that covers.

    Attributes:
        kind: The modality's reference kind. Every member shares it, by construction.
        references: Resolved references, in resolution order.
        rejections: Requests that produced no usable reference.
        requested_species: Species the run asked to screen, before resolution.
    """

    kind: ReferenceKind
    references: tuple[ScreeningReference, ...] = ()
    rejections: tuple[ReferenceRejection, ...] = ()
    requested_species: tuple[str, ...] = ()

    @property
    def species(self) -> tuple[str, ...]:
        """Species with a usable reference, deduplicated, in resolution order."""
        return tuple(dict.fromkeys(reference.species for reference in self.references))

    @property
    def index_parameter(self) -> str:
        """The ``species:index_prefix,...`` value for references whose index already exists."""
        return ",".join(
            dict.fromkeys(reference.index_entry for reference in self.references if not reference.needs_index_build)
        )

    @property
    def fasta_parameter(self) -> str:
        """The ``species:path,...`` value for references the pipeline still has to index."""
        return ",".join(
            dict.fromkeys(reference.index_entry for reference in self.references if reference.needs_index_build)
        )

    @property
    def enabled(self) -> bool:
        """Whether anything resolved."""
        return bool(self.references)

    def scope(self) -> FilterScope:
        """The species this screen resolved references for, as an explicit set rather than a counter's name.

        Planned, not covered: it is fixed when the references resolve, and a species whose alignment
        later published nothing is reported as a shortfall rather than removed from here.
        """
        from sirnaforge.models.policy import FilterScope  # noqa: PLC0415

        return FilterScope(species=frozenset(self.species))

    def plan(self, *, guide_set_digest: str, search_settings: SearchSettings | None = None) -> ScreeningPlan:
        """The screening plan these references intend to execute, one entry per reference."""
        from sirnaforge.models.evidence import ScreeningPlan  # noqa: PLC0415

        return ScreeningPlan(
            entries=tuple(
                reference.plan_entry(guide_set_digest=guide_set_digest, search_settings=search_settings)
                for reference in self.references
            )
        )

    def to_metadata(self) -> dict[str, object]:
        """Serializable snapshot: what resolved, what did not, and over which species."""
        return {
            "kind": self.kind.value,
            "requested_species": list(self.requested_species),
            "resolved_species": list(self.species),
            "references": [reference.to_metadata() for reference in self.references],
            "rejections": [
                {"species": rejection.species, "identity": rejection.identity, "reason": rejection.reason}
                for rejection in self.rejections
            ],
        }


def screening_kind_for_design_mode(design_mode: DesignMode | str) -> ReferenceKind:
    """The reference kind a design mode screens against.

    ``DesignMode`` is read by value rather than imported: ``sirnaforge.models`` imports this module,
    so importing it back at module scope would be circular.
    """
    mode = str(getattr(design_mode, "value", design_mode)).strip().lower()
    return ReferenceKind.GENOME if mode == "zfn" else ReferenceKind.TRANSCRIPTOME


def split_declared_species(value: str) -> tuple[str | None, str]:
    """Split an optional ``species:`` prefix off a reference, returning ``(species, remainder)``.

    The prefix is only read as a species when it is a name the species registry recognises, so a
    URL scheme (``https:``) or a Windows drive letter is never mistaken for one.
    """
    from sirnaforge.data.species_registry import CANONICAL_SPECIES_ALIAS_MAP  # noqa: PLC0415

    entry = value.strip()
    head, separator, remainder = entry.partition(":")
    if not separator or not remainder.strip():
        return None, entry
    if head.strip().lower() not in CANONICAL_SPECIES_ALIAS_MAP:
        return None, entry
    return head.strip(), remainder.strip()


def parse_index_entries(raw: str | None, *, option: str, reason: str) -> tuple[ReferenceRequest, ...]:
    """Parse ``species:index_prefix`` entries into explicit prebuilt-index requests.

    Args:
        raw: Comma-separated entries, or None.
        option: Option name quoted in the error, so a malformed entry is findable.
        reason: Why these requests exist, carried into the resolved reference.

    Raises:
        ValueError: An entry does not name a species, or names one the registry does not recognise.
            An index prefix carries no headers to fall back on, so a declared species here is the
            only authority for that reference's label -- an unrecognised one would be accepted as
            the species (``transcriptome:`` and typos alike), and orthology would silently resolve
            against nothing. ``--species`` refuses an unknown name for the same reason.
    """
    from sirnaforge.data.species_registry import CANONICAL_SPECIES_ALIAS_MAP  # noqa: PLC0415

    if not raw:
        return ()

    requests: list[ReferenceRequest] = []
    for token in raw.split(","):
        entry = token.strip()
        if not entry:
            continue
        species, separator, prefix = entry.partition(":")
        if not separator or not species.strip() or not prefix.strip():
            raise ValueError(f"{option} entries must be in species:/index_prefix form; got {entry!r}")
        if species.strip().lower() not in CANONICAL_SPECIES_ALIAS_MAP:
            supported = ", ".join(sorted(set(CANONICAL_SPECIES_ALIAS_MAP.values())))
            raise ValueError(
                f"{option} entry {entry!r} names an unsupported species {species.strip()!r}. "
                f"Supported canonical species: {supported}"
            )
        requests.append(
            ReferenceRequest(
                value=prefix.strip(),
                form=ReferenceForm.PREBUILT_INDEX,
                state=ReferenceState.EXPLICIT,
                reason=reason,
                declared_species=species.strip(),
            )
        )
    return tuple(requests)


def kind_of_bundled_source(identity: str) -> ReferenceKind | None:
    """The kind of a bundled Ensembl source name, or None when the name is not a bundled source.

    This is what makes a modality mismatch catchable before any download: a genomic assembly name
    handed to a transcriptome screen is recognisable from the source table alone.
    """
    from sirnaforge.data.ensembl_references import build_genome_sources, build_transcriptome_sources  # noqa: PLC0415

    name = identity.strip()
    if name in build_transcriptome_sources():
        return ReferenceKind.TRANSCRIPTOME
    if name in build_genome_sources():
        return ReferenceKind.GENOME
    return None


def resolve_reference_species(
    *,
    declared: str | None,
    source_species: str | None,
    header_species: str | None,
) -> tuple[str, SpeciesAuthority, str | None]:
    """Resolve a reference's species and say which authority supplied it.

    Returns ``(species, authority, conflict)``. ``conflict`` is a message when a weaker authority
    disagrees with the winner -- reported rather than silently discarded, because a declared species
    that contradicts the file is usually a mistake in one of the two.
    """
    from sirnaforge.data.species_registry import normalize_species_name  # noqa: PLC0415

    candidates = (
        (declared, SpeciesAuthority.DECLARED),
        (source_species, SpeciesAuthority.SOURCE_REGISTRY),
        (header_species, SpeciesAuthority.REFERENCE_HEADERS),
    )
    winner: tuple[str, SpeciesAuthority] | None = None
    for value, authority in candidates:
        label = (value or "").strip()
        if not label:
            continue
        normalized = normalize_species_name(label)
        if winner is None:
            winner = (normalized, authority)
        elif normalized != winner[0]:
            return (
                *winner,
                f"{authority.value} says {normalized!r} but {winner[1].value} says {winner[0]!r}; using {winner[0]!r}",
            )
    if winner is None:
        return UNRESOLVED_SPECIES, SpeciesAuthority.UNRESOLVED, None
    return (*winner, None)


def build_screening_requests(
    *,
    kind: ReferenceKind,
    selection: ReferenceSelection,
    index_requests: Sequence[ReferenceRequest] = (),
) -> tuple[tuple[ReferenceRequest, ...], tuple[ReferenceRejection, ...]]:
    """Turn a selection and any explicit index requests into one ordered request list.

    The single place a screening reference can enter a run. Explicit index requests come first: they
    are an override, and the species they declare decide the screen.

    Raises:
        ReferenceKindError: An explicitly requested reference's kind does not match ``kind``. A
            *default* that does not match is disabled instead and returned as a rejection, because a
            default nobody asked for must not fail a run.
    """
    requests = list(index_requests)
    for choice in selection.choices:
        if not choice.value:
            continue
        requests.append(
            ReferenceRequest(
                value=choice.value,
                form=ReferenceForm.REFERENCE,
                state=choice.state,
                reason=choice.reason,
            )
        )

    kept: list[ReferenceRequest] = []
    rejections: list[ReferenceRejection] = []
    for request in requests:
        declared, value = (
            (request.declared_species, request.value)
            if request.form is ReferenceForm.PREBUILT_INDEX
            else split_declared_species(request.value)
        )
        request = ReferenceRequest(  # noqa: PLW2901 - normalized in place; the split is part of parsing
            value=value,
            form=request.form,
            state=request.state,
            reason=request.reason,
            declared_species=declared,
        )
        source_kind = kind_of_bundled_source(request.value)
        if kind is ReferenceKind.GENOME:
            # ZFN screens genomic DNA. A transcriptome reference reaching it is a mismatch whether or
            # not it names a bundled source, so the check cannot rest on the source table alone -- and
            # a prebuilt index arrives on --transcriptome-indices, so it is a transcriptome index by
            # construction and must be refused here too, not only in its REFERENCE form.
            mismatch = source_kind is not ReferenceKind.GENOME
        else:
            mismatch = source_kind is not None and source_kind is not kind
        if not mismatch:
            kept.append(request)
            continue
        detail = (
            f"reference {request.value!r} is a {(source_kind or ReferenceKind.TRANSCRIPTOME).value} reference, "
            f"but this run screens against a {kind.value}"
        )
        if request.is_explicit:
            raise ReferenceKindError(
                f"{detail}. siRNA and miRNA screening always uses a transcriptome and ZFN always uses a "
                "genome, so this reference cannot be screened as requested."
            )
        rejections.append(
            ReferenceRejection(
                species=UNRESOLVED_SPECIES,
                identity=request.value,
                reason=f"{detail}, so the default was not resolved",
            )
        )
    return tuple(kept), tuple(rejections)


def render_reference_selection_label(selection: ReferenceSelection) -> str:
    """Render a stable, human-readable label for CLI/config summaries."""
    if selection.enabled:
        rendered_choices = [f"{choice.value} ({choice.state.value})" for choice in selection.choices if choice.value]
        return ", ".join(rendered_choices)

    reason = selection.disabled_reason or "not available"
    return f"disabled ({reason})"


__all__ = [
    "DEFAULT_TRANSCRIPTOME_SOURCE",
    "DEFAULT_TRANSCRIPTOME_SOURCES",
    "DEFAULT_MIRNA_SOURCE",
    "DEFAULT_MIRNA_CANONICAL_SPECIES",
    "UNRESOLVED_SPECIES",
    "ReferenceChoice",
    "ReferenceForm",
    "ReferenceKind",
    "ReferenceKindError",
    "ReferenceRejection",
    "ReferenceRequest",
    "ReferenceSelection",
    "ReferencePolicyResolver",
    "ReferenceState",
    "ScreeningReference",
    "ScreeningReferenceSet",
    "SpeciesAuthority",
    "WorkflowInputSpec",
    "build_screening_requests",
    "kind_of_bundled_source",
    "parse_index_entries",
    "render_reference_selection_label",
    "resolve_reference_species",
    "screening_kind_for_design_mode",
    "split_declared_species",
]
