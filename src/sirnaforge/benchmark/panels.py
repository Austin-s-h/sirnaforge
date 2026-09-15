"""Benchmark panel registry: architectures, compatibility verdicts, honestly-labelled fixtures.

Issue #109. A "panel" is a named siRNA/miRNA efficacy dataset the benchmark surface can prepare an
artifact from. This module is the one place a panel is described: its citation, its declared duplex
architecture, the columns a source table must carry, and whether this repository actually vendors
its bytes. ``sirnaforge.benchmark.prepare`` reads a descriptor and a raw source row and calls
:func:`derive_observation` here for every field the panel/architecture determines; it does not
re-derive a slice or a compatibility verdict itself, so the rule lives in exactly one place.

Two verified facts this module works around rather than papering over (full accounting in
``tests/data/benchmarks/README.md``, which is the authority for the vendored bytes):

1. The issue's cited source, ``docs/prd_benchmark_artifacts_and_variable_length.md``, does not exist
   in this repository -- not in the working tree, not in any branch's history. The issue body is the
   entire specification; nothing below claims to summarise a document that was read.
2. Four of the panels #109 named are now vendored, and every descriptor below is re-derived from
   those bytes rather than from the issue's prose. ``tests/data/benchmarks/oligogym/records.csv``
   (4,113 rows, 22 columns, written by ``scripts/prepare_oligo_benchmarks.py``) carries the measured
   guide/passenger sequences and labels for ``ichihara`` (2,850 rows over two OligoGym datasets),
   ``martinelli`` (907) and ``shmushkovich`` (356), each row naming its own DOI in ``source_url``.
   What is *not* vendored is the upstream OligoGym extract those rows were derived from
   (``tests/data/external/oligogym/*.csv.gz``, named in every row's ``source_file``, absent from this
   tree) and any primary-source PDF; so a DOI, a per-dataset row count and the sequences themselves
   are the whole of what this repository can show, and the citations below say exactly that.

Every ``oligogym/records.csv`` row is ``context_type="synthetic_neutral_flanks"``: the transcript it
was written against is ``70 x A + reverse-complement(guide) + 70 x A``, so its ``site_start``/
``site_end`` (1-based 71..) are positions in a *fabricated* context. That is a real coordinate in a
real file and a useless one for accessibility, which is why it is recorded under its own target
identity value, :attr:`TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL`, rather than being pushed into
``panel_local`` (which claims a coordinate in the panel's own measured target) or ``confirmed``
(a native transcript, which #110 owns and which ``tests/data/benchmarks/oligogym_native_design/``
-- not read here -- is the mapped counterpart of).

Architectures (:class:`PanelArchitecture`) are the three duplex geometries #109's fixed-length
design surface can consume. ``ASYMMETRIC`` is excluded from every compatibility verdict by
construction, not merely under-supported: #109's scope excludes variable-length/asymmetric design
entirely (#110 owns giving it a real, native path), and silently truncating or padding a 15/20
duplex into a 19-23 nt paired core would misrepresent a measured sequence, which the parent issue
also forbids. :func:`derive_observation` therefore never slices an asymmetric record -- it returns
the verbatim guide with ``compatibility_status="incompatible"``. All 356 vendored Shmushkovich rows
take that path.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, model_validator

#: The predeclared development/held-out split rule this module implements. See
#: :func:`predeclared_split`.
PREDECLARED_SPLIT_RULE_ID = "sha256_accession_parity_v1"


class PanelArchitecture(str, Enum):
    """Duplex geometry a panel's guide/passenger pair was measured as.

    Attributes:
        PAIRED_CORE_WITH_OVERHANG: A fixed-length complementary core (this repo: 19 nt) plus a
            measured 3' overhang on the guide -- the canonical Tuschl-rule 21-mer duplex Huesken and
            Ichihara both used.
        FULLY_COMPLEMENTARY: Guide and passenger pair over their full, equal length; no overhang.
        ASYMMETRIC: Guide and passenger differ in length (e.g. Shmushkovich's 15/20 duplex). Never
            compatible with #109's fixed-length interface; see #110.
    """

    PAIRED_CORE_WITH_OVERHANG = "paired_core_with_overhang"
    FULLY_COMPLEMENTARY = "fully_complementary"
    ASYMMETRIC = "asymmetric"


class DuplexPairingStatus(str, Enum):
    """Per-record pairing geometry, as stated by the panel. A superset of :class:`PanelArchitecture`.

    ``UNSTATED`` is the fourth value the artifact schema declares, for a panel that measured a guide
    without ever recording how -- or whether -- it pairs a passenger. No descriptor in
    :data:`PANEL_REGISTRY` emits it, because all six declare an architecture; it is declared here
    because the vocabulary is shared with ``sirnaforge.benchmark.artifact``'s observation schema, and
    a future panel legitimately needs it.
    """

    PAIRED_CORE_WITH_OVERHANG = "paired_core_with_overhang"
    FULLY_COMPLEMENTARY = "fully_complementary"
    ASYMMETRIC = "asymmetric"
    UNSTATED = "unstated"


class CompatibilityStatus(str, Enum):
    """Whether one observation can enter #109's fixed-length paired-core design path."""

    COMPATIBLE = "compatible"
    INCOMPATIBLE = "incompatible"


class TargetIdentityStatus(str, Enum):
    """How well this repository knows *where* an observation's target site is.

    Extends the two-value vocabulary ``sirnaforge.benchmark.artifact.TargetIdentityStatus`` declared
    (``unavailable | panel_local``) with a third value, because the bytes vendored since #109 was
    settled fit neither: ``tests/data/benchmarks/oligogym/records.csv`` carries a real, exact
    ``site_start``/``site_end`` -- and that coordinate is a position in a synthetic transcript the
    adapter fabricated (``70 x A + revcomp(guide) + 70 x A``), not in anything the panel measured.
    Calling that ``panel_local`` would overclaim (it reads as "the panel's own target coordinate")
    and calling it ``unavailable`` would throw away a coordinate the artifact can honestly reproduce,
    so #109 records it as what it is.

    There is deliberately no ``confirmed`` member. A native transcript identity is #110's to
    establish (``tests/data/benchmarks/oligogym_native_design/``), so #109 cannot emit one even by
    accident -- the value does not exist in this vocabulary to be emitted.

    Attributes:
        UNAVAILABLE: The panel states no target coordinate at all, in any context.
        PANEL_LOCAL: A coordinate in the panel's own measured target sequence.
        SYNTHETIC_CONTEXT_LOCAL: A coordinate in a synthetic context this repository generated. Exact
            and reproducible; evidence about design enumeration only, never about native
            accessibility.
    """

    UNAVAILABLE = "unavailable"
    PANEL_LOCAL = "panel_local"
    SYNTHETIC_CONTEXT_LOCAL = "synthetic_context_local"


def predeclared_split(accession: str) -> str:
    """Development/held-out split for one accession, by SHA-256 parity.

    Reuses -- rather than reimplements from scratch -- the exact rule
    ``scripts/validate_scoring_profiles.py::split_of`` pinned for the Huesken panel (issues #97,
    #102): ``sha256(accession)[0] % 2 == 0 -> development``. That script is a standalone entry point,
    not an importable package member, so the rule is duplicated here under its own id,
    :data:`PREDECLARED_SPLIT_RULE_ID`, rather than imported.
    ``tests/unit/test_benchmark_panels.py`` pins this function against the exact accession lists
    published in ``tests/unit/data/README.md`` so the duplication cannot silently drift from its
    named authority.

    Deterministic and independent of every efficacy value and every feature: computable from the
    accession string alone, and fixed before any number derived from it was computed.
    """
    digest = hashlib.sha256(accession.encode()).digest()
    return "development" if digest[0] % 2 == 0 else "held_out"


class PanelColumnMapping(BaseModel):
    """Which columns of a panel's raw source table hold each measured field.

    ``passenger_column``, ``accession_column`` and ``measured_value_column`` are all optional: a
    panel may measure a guide only (no passenger), record no numeric endpoint, or -- for a panel with
    no predeclared split -- carry no accession a split could be computed from.

    For the three vendored OligoGym-derived panels the mapping onto
    ``tests/data/benchmarks/oligogym/records.csv``'s 22 columns is, field by field:

    ==========================  =================================================================
    observation field           ``records.csv`` column
    ==========================  =================================================================
    ``full_guide_sequence``     ``guide_sequence`` (the antisense strand; ``guide_source`` says
                                which OligoGym component it came from)
    ``passenger_sequence``      ``passenger_sequence``
    ``measured_value``          ``efficacy_higher_is_better`` for Ichihara/Martinelli, where it is
                                identical to ``label_raw``/``label_processed`` in all 3,757 rows;
                                ``label_processed`` for Shmushkovich, whose
                                ``efficacy_higher_is_better`` is a derived ``100 - label`` flip
                                (verified on all 356 rows) that no descriptor here reads
    ``measured_endpoint``       not a column: ``label_semantics``/``label_direction`` are constant
                                per dataset, and the descriptor names the endpoint they state
    ``target_identity_status``  ``context_type`` (constant ``synthetic_neutral_flanks``), checked
                                per row against the descriptor rather than trusted
    ``target_start_1based``     ``site_start`` (1-based, 71 on every row)
    ``target_end_1based``       ``site_end``
    row selection               ``dataset`` -- see :class:`PanelRowSelector`; one file holds all four
    ==========================  =================================================================

    Some columns are deliberately read by nothing. ``benchmark_id``/``transcript_id`` are ids the
    adapter minted (``bench_<dataset>_row<n>``); ``input_fasta`` names the synthetic FASTA rather than
    a transcript; ``source_row`` is 1-based where ``BenchmarkObservation.source_row_index`` is 0-based
    within the *selected* rows, so equating them would silently mislabel provenance; and
    ``source_file``/``source_url``/``guide_source``/``label_semantics`` are provenance the descriptor
    states once for the panel rather than re-reading per row.
    ``target`` (a gene symbol, empty on all 1,263 Martinelli/Shmushkovich rows and some Ichihara ones)
    is *not* mapped to ``accession_column``: it is neither an accession nor unique per site, and
    ``predeclared_split`` is pinned to accessions (see :func:`predeclared_split`), so splitting on it
    would be a new, unaudited rule wearing an audited rule's id.
    """

    guide_column: str = Field(min_length=1)
    passenger_column: str | None = None
    accession_column: str | None = None
    measured_value_column: str | None = None
    context_type_column: str | None = None
    target_site_start_column: str | None = None
    target_site_end_column: str | None = None

    model_config = ConfigDict(frozen=True, extra="forbid")


class AssayLabelSource(BaseModel):
    """Where a record's ``assay_label`` comes from: one source-table column, or a fixed constant.

    Exactly one of the two must be set. A per-row column wins when present, because it is more
    specific than a panel-wide constant; the constant exists because every panel in
    :data:`PANEL_REGISTRY` currently runs a single assay that does not vary row to row.
    """

    column: str | None = None
    constant: str | None = None

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def _exactly_one_source(self) -> AssayLabelSource:
        if (self.column is None) == (self.constant is None):
            raise ValueError("AssayLabelSource requires exactly one of column or constant")
        return self

    def resolve(self, row: Mapping[str, str]) -> str:
        """Assay label for one raw row. The artifact schema requires it non-empty."""
        label = row[self.column] if self.column is not None else self.constant
        if not label:
            raise ValueError("assay_label resolved empty; the artifact schema requires it non-empty")
        return label


class PanelRowSelector(BaseModel):
    """Which rows of a shared source table belong to one panel.

    Needed because the vendored bytes arrived as a single table:
    ``tests/data/benchmarks/oligogym/records.csv`` holds all four OligoGym datasets, distinguished
    only by its ``dataset`` column. Without a selector, ``benchmark prepare --panel ichihara`` would
    ingest Martinelli's and Shmushkovich's rows too and stamp Ichihara's declared architecture on
    them -- the exact relabelling #109 exists to prevent -- so the descriptor states its own rows and
    the reader honours it.

    ``values`` is a closed set, never a prefix or a pattern: a new OligoGym dataset appearing in a
    regenerated ``records.csv`` must be declared here deliberately, not swept into whichever panel's
    name it happens to start with.
    """

    column: str = Field(min_length=1)
    values: tuple[str, ...] = Field(min_length=1)

    model_config = ConfigDict(frozen=True, extra="forbid")

    def matches(self, row: Mapping[str, str]) -> bool:
        """Whether one raw row belongs to this panel.

        Raises:
            KeyError: The selector column is absent -- a table that cannot say which panel a row
                belongs to must fail loudly, not silently contribute every row or none.
        """
        return row[self.column] in self.values


class PanelDescriptor(BaseModel):
    """Everything the benchmark surface needs to know about one named panel.

    Declarative only: no file I/O, no pandas. ``sirnaforge.benchmark.prepare`` reads a descriptor and
    a raw ``Mapping[str, str]`` row (however it obtained the row) and passes both to
    :func:`derive_observation` for the biology; this model states the contract, it does not apply it.

    Attributes:
        panel_id: Registry key, ``[a-z0-9_]+`` -- also the artifact directory's panel component.
        display_name: Human-readable name for CLI/manifest output.
        citation: Primary citation, or an explicit statement that none has been independently
            verified in this repository (see the module docstring's second verified fact).
        redistribution: Third-party redistribution the bytes actually came from, when not the
            primary source. ``None`` when obtained from the primary source, or when no bytes are
            vendored at all (``data_present=False``).
        architecture: Declared duplex geometry; see :class:`PanelArchitecture`.
        declared_paired_length: The fixed-length paired core this panel's architecture implies.
            ``None`` only for :attr:`PanelArchitecture.ASYMMETRIC`, which has none.
        columns: Source-table column mapping; see :class:`PanelColumnMapping`.
        assay_label: Where ``assay_label`` comes from; see :class:`AssayLabelSource`.
        measured_endpoint: Name of the numeric quantity ``columns.measured_value_column`` holds, e.g.
            ``"inhibition_fraction"``. Required non-empty even when no column is configured, because
            a number with no endpoint is uninterpretable and this documents what the panel measures
            (or would measure) regardless of whether it is vendored.
        split_rule_id: Predeclared split rule id (currently only :data:`PREDECLARED_SPLIT_RULE_ID`),
            or ``None`` when this panel declares no split.
        data_present: Whether this repository vendors real bytes for this panel. True for
            ``huesken_subset`` and for the four OligoGym-derived ids vendored by ``f4beab7`` (three
            ingestible panels plus the ``oligogym`` aggregate they share a table with); see the module
            docstring for what "present" does and does not cover.
        vendored_csv: Repo-relative path of those bytes, or ``None`` when ``data_present`` is False.
            Kept here rather than only in ``prepare.py`` so one file answers "does this repository
            ship this panel, and from where" -- and so ``content_hash`` changes when the answer does.
        row_selector: Which rows of ``vendored_csv`` (or of a ``--panel-csv``) are this panel's; see
            :class:`PanelRowSelector`. ``None`` means every row.
        target_identity_status: How well the source table locates the target site; see
            :class:`TargetIdentityStatus`. Anything but ``UNAVAILABLE`` requires the coordinate
            columns that back it.
        aggregate_of: Panel ids this descriptor is only the shared redistribution of. Set for
            ``oligogym`` alone, whose one table spans three duplex geometries, so no single
            ``architecture`` could describe it; :func:`derive_observation` refuses an aggregate rather
            than deriving a row under a geometry that is wrong for 405 of its 4,113 rows.
    """

    panel_id: str = Field(pattern=r"^[a-z0-9_]+$")
    display_name: str = Field(min_length=1)
    citation: str = Field(min_length=1)
    redistribution: str | None = None
    architecture: PanelArchitecture
    declared_paired_length: int | None = Field(default=None, ge=19, le=23)
    columns: PanelColumnMapping
    assay_label: AssayLabelSource
    measured_endpoint: str = Field(min_length=1)
    split_rule_id: str | None = None
    data_present: bool
    vendored_csv: str | None = None
    row_selector: PanelRowSelector | None = None
    target_identity_status: TargetIdentityStatus = TargetIdentityStatus.UNAVAILABLE
    aggregate_of: tuple[str, ...] | None = None

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def _paired_length_matches_architecture(self) -> PanelDescriptor:
        if self.architecture is PanelArchitecture.ASYMMETRIC:
            if self.declared_paired_length is not None:
                raise ValueError("an asymmetric panel declares no fixed paired length")
        elif self.declared_paired_length is None:
            raise ValueError(f"a {self.architecture.value} panel must declare a paired length")
        return self

    @model_validator(mode="after")
    def _data_present_iff_a_vendored_path_is_named(self) -> PanelDescriptor:
        """``data_present`` is a claim about bytes, so it must name them (or claim nothing).

        The two halves drifting apart is the failure this catches: a ``True`` with no path is a claim
        no reader can check, and a path with ``False`` is bytes the manifest would then disown.
        """
        if self.data_present and self.vendored_csv is None:
            raise ValueError(f"panel {self.panel_id!r} claims data_present but names no vendored_csv")
        if not self.data_present and self.vendored_csv is not None:
            raise ValueError(f"panel {self.panel_id!r} names vendored_csv {self.vendored_csv!r} but claims no data")
        return self

    @model_validator(mode="after")
    def _target_identity_has_the_columns_it_claims(self) -> PanelDescriptor:
        """A located target site must say which columns locate it, and in which context.

        ``SYNTHETIC_CONTEXT_LOCAL`` additionally requires ``context_type_column``: the whole point of
        that value is that the context is fabricated, and :func:`derive_observation` verifies that
        claim per row instead of trusting the descriptor, which it can only do if the row says so.
        """
        located = self.target_identity_status is not TargetIdentityStatus.UNAVAILABLE
        has_coordinates = (
            self.columns.target_site_start_column is not None and self.columns.target_site_end_column is not None
        )
        if located and not has_coordinates:
            raise ValueError(
                f"panel {self.panel_id!r} declares target_identity_status "
                f"{self.target_identity_status.value!r} but maps no site start/end columns"
            )
        if not located and has_coordinates:
            raise ValueError(
                f"panel {self.panel_id!r} maps site start/end columns but declares target_identity_status "
                "'unavailable'; a coordinate this repository can read is a coordinate it must label"
            )
        if self.target_identity_status is TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL and (
            self.columns.context_type_column is None
        ):
            raise ValueError(
                f"panel {self.panel_id!r} declares a synthetic-context site but maps no "
                "context_type_column, so no row could be checked against that claim"
            )
        return self

    def selects_row(self, row: Mapping[str, str]) -> bool:
        """Whether one raw row is this panel's. Always True when the panel declares no selector.

        The reader (``sirnaforge.benchmark.prepare``) must call this before deriving a row, because a
        vendored table may be shared between panels (:class:`PanelRowSelector`).
        """
        return True if self.row_selector is None else self.row_selector.matches(row)

    def content_hash(self) -> str:
        """``sha256:`` digest over the descriptor, mirroring ``ProfileIdentity.content_hash``.

        Covers everything this descriptor decides -- citation, redistribution, architecture, paired
        length, column mapping, assay label source, split rule and ``data_present`` -- so an edit to
        any one of them (including flipping ``data_present`` the day real bytes are vendored) changes
        the hash a manifest records. The same guarantee ``config/run_policy.py::RunPolicyProfile
        .identity`` gives a resolved policy, and for the same reason: two artifacts sharing this hash
        must have been built from the same declared contract, not merely the same panel name.
        """
        payload = json.dumps(self.model_dump(mode="json"), sort_keys=True, default=str)
        digest = hashlib.sha256(payload.encode()).hexdigest()
        return f"sha256:{digest}"


class DerivedObservation(BaseModel):
    """One panel record's architecture-derived fields.

    The part of ``observations.csv`` (see ``sirnaforge.benchmark.artifact``) that a panel descriptor
    determines from a raw row, independent of source-file provenance and observation identity, which
    its caller already has (the source file path/hash, the row index, the panel id).

    ``target_identity_status``/``target_context_type``/``target_start_1based``/``target_end_1based``
    are here rather than in ``prepare.py`` because the columns that answer them are a fact about the
    panel's table, which is this module's subject. ``target_context_type`` is the row's own words for
    the context those coordinates index (``synthetic_neutral_flanks`` for every vendored OligoGym
    row); it is carried so a reader never has to infer the context from the status name.
    """

    architecture: PanelArchitecture
    assay_label: str
    measured_endpoint: str
    measured_value: float | None
    full_guide_sequence: str
    guide_length: int
    paired_guide_sequence: str
    paired_slice_start_1based: int
    paired_length: int
    passenger_sequence: str | None
    guide_3p_overhang: str | None
    passenger_3p_overhang: str | None
    duplex_pairing_status: DuplexPairingStatus
    target_identity_status: TargetIdentityStatus
    target_context_type: str | None
    target_start_1based: int | None
    target_end_1based: int | None
    compatibility_status: CompatibilityStatus
    compatibility_reason: str
    split: str | None

    model_config = ConfigDict(frozen=True, extra="forbid")


_ALLOWED_BASES = frozenset("ACGTUN")

#: The one ``context_type`` value ``tests/data/benchmarks/oligogym/records.csv`` uses, on all 4,113
#: rows: the transcript each site was written into is ``70 x A + revcomp(guide) + 70 x A``. Spelled
#: once here so the descriptor's declaration and the per-row check cannot disagree.
_SYNTHETIC_NEUTRAL_FLANKS = "synthetic_neutral_flanks"


def _clean_sequence(raw: str, *, field: str) -> str:
    """Upper-case and validate one sequence field. Never pads, trims or reorders it."""
    sequence = raw.strip().upper()
    if not sequence:
        raise ValueError(f"{field} is empty")
    invalid = sorted(set(sequence) - _ALLOWED_BASES)
    if invalid:
        raise ValueError(f"{field} has non-nucleotide characters: {invalid}")
    return sequence


def _check_compatibility(
    architecture: PanelArchitecture,
    guide_length: int,
    passenger_length: int | None,
    requested_paired_length: int,
) -> tuple[CompatibilityStatus, str]:
    """The one compatibility rule every panel and every fixture goes through.

    ``ASYMMETRIC`` is incompatible unconditionally -- #109 excludes asymmetric design entirely, so no
    ``requested_paired_length`` can make a 15/20 duplex compatible; #110 is the issue that gives it a
    real path. ``FULLY_COMPLEMENTARY`` requires the measured guide (and, when present, the passenger)
    to already equal the requested length, because inventing a slice would fabricate an overhang that
    was never measured. ``PAIRED_CORE_WITH_OVERHANG`` requires only that the guide be at least as
    long as the requested core, since the remainder is the (possibly zero-length) measured overhang.
    """
    if architecture is PanelArchitecture.ASYMMETRIC:
        duplex = (
            f"guide {guide_length} nt"
            if passenger_length is None
            else f"guide {guide_length} nt / passenger {passenger_length} nt"
        )
        return (
            CompatibilityStatus.INCOMPATIBLE,
            f"asymmetric duplex ({duplex}) is out of #109's fixed-length paired-core scope; "
            "#110 owns native asymmetric handling",
        )
    if architecture is PanelArchitecture.FULLY_COMPLEMENTARY:
        if passenger_length is not None and passenger_length != guide_length:
            return (
                CompatibilityStatus.INCOMPATIBLE,
                f"declared fully complementary but guide ({guide_length} nt) and passenger "
                f"({passenger_length} nt) differ in length",
            )
        if guide_length != requested_paired_length:
            return (
                CompatibilityStatus.INCOMPATIBLE,
                f"fully complementary duplex measured at {guide_length} nt; a "
                f"{requested_paired_length} nt paired core would invent an overhang that was never "
                "measured",
            )
        return CompatibilityStatus.COMPATIBLE, ""
    # PAIRED_CORE_WITH_OVERHANG
    if guide_length < requested_paired_length:
        return (
            CompatibilityStatus.INCOMPATIBLE,
            f"guide is {guide_length} nt, shorter than the requested {requested_paired_length} nt paired core",
        )
    return CompatibilityStatus.COMPATIBLE, ""


def _derive_target_identity(
    descriptor: PanelDescriptor, row: Mapping[str, str]
) -> tuple[TargetIdentityStatus, str | None, int | None, int | None]:
    """Target-site locality for one row: status, the context's own name, and its 1-based span.

    Verifies the descriptor's context claim against the row rather than asserting it: a
    ``SYNTHETIC_CONTEXT_LOCAL`` panel whose row says its context is anything other than
    ``synthetic_neutral_flanks`` raises, because the one thing #109 must never do is let a coordinate
    from a native context be recorded under a synthetic label -- or the reverse, which is how a
    fabricated position at 1-based 71 would end up read as evidence about a real transcript.
    """
    if descriptor.target_identity_status is TargetIdentityStatus.UNAVAILABLE:
        return TargetIdentityStatus.UNAVAILABLE, None, None, None

    context_type: str | None = None
    if descriptor.columns.context_type_column is not None:
        context_type = row[descriptor.columns.context_type_column].strip()
        if (
            descriptor.target_identity_status is TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL
            and context_type != _SYNTHETIC_NEUTRAL_FLANKS
        ):
            raise ValueError(
                f"panel {descriptor.panel_id!r} declares a {_SYNTHETIC_NEUTRAL_FLANKS} context but this "
                f"row states context_type {context_type!r}; refusing to record its site coordinates "
                "under the wrong context"
            )

    # Both columns are guaranteed present by PanelDescriptor's validator, so a KeyError here is a
    # malformed table (the caller reports the row), not an under-declared descriptor.
    start = int(row[str(descriptor.columns.target_site_start_column)])
    end = int(row[str(descriptor.columns.target_site_end_column)])
    if start < 1 or end < start:
        raise ValueError(f"target site span {start}..{end} is not a 1-based, non-empty interval")
    return descriptor.target_identity_status, context_type, start, end


def derive_observation(
    descriptor: PanelDescriptor,
    row: Mapping[str, str],
    *,
    requested_paired_length: int,
) -> DerivedObservation:
    """Architecture-derive one raw row into everything ``benchmark.prepare`` needs from a panel.

    ``requested_paired_length`` is the artifact directory's length tag (#109: equal to
    ``descriptor.declared_paired_length`` by default, but a caller may request a shorter core from a
    longer measured guide). Slicing is always 5'-anchored at position 1: every architecture #109
    supports measures its paired core from the guide's 5' end, so "the audited slice" never needs a
    second offset convention -- :attr:`DerivedObservation.paired_slice_start_1based` is always ``1``.

    An incompatible record is returned, never raised: :attr:`DerivedObservation.compatibility_status`
    carries the verdict so the caller can still write the observation row -- #109 requires every
    measured observation to appear in ``observations.csv``, incompatible or not. Its
    ``paired_guide_sequence`` equals ``full_guide_sequence`` verbatim in that case, never truncated or
    padded to ``requested_paired_length``: a compatibility failure is a statement about the duplex's
    geometry, not a licence to misrepresent the sequence that was measured.

    Raises:
        ValueError: The descriptor is an aggregate (:attr:`PanelDescriptor.aggregate_of`), the row's
            ``context_type`` contradicts the descriptor's, a sequence field is empty or non-nucleotide,
            or a declared split has no accession column.
    """
    if descriptor.aggregate_of is not None:
        raise ValueError(
            f"panel {descriptor.panel_id!r} is an aggregate of {', '.join(descriptor.aggregate_of)} and has no "
            "single duplex architecture of its own; prepare one of those panels instead"
        )
    full_guide_sequence = _clean_sequence(row[descriptor.columns.guide_column], field="full_guide_sequence")
    guide_length = len(full_guide_sequence)

    passenger_sequence: str | None = None
    if descriptor.columns.passenger_column is not None:
        raw_passenger = row.get(descriptor.columns.passenger_column)
        if raw_passenger:
            passenger_sequence = _clean_sequence(raw_passenger, field="passenger_sequence")
    passenger_length = None if passenger_sequence is None else len(passenger_sequence)

    status, reason = _check_compatibility(
        descriptor.architecture, guide_length, passenger_length, requested_paired_length
    )

    guide_3p_overhang: str | None
    passenger_3p_overhang: str | None
    if status is CompatibilityStatus.COMPATIBLE:
        paired_guide_sequence = full_guide_sequence[:requested_paired_length]
        guide_3p_overhang = full_guide_sequence[requested_paired_length:]
        passenger_3p_overhang = None if passenger_sequence is None else passenger_sequence[requested_paired_length:]
        paired_length = requested_paired_length
    else:
        paired_guide_sequence = full_guide_sequence
        guide_3p_overhang = None
        passenger_3p_overhang = None
        paired_length = guide_length

    measured_value: float | None = None
    if descriptor.columns.measured_value_column is not None:
        raw_value = row.get(descriptor.columns.measured_value_column)
        if raw_value not in (None, ""):
            measured_value = float(raw_value)

    split: str | None = None
    if descriptor.split_rule_id is not None:
        if descriptor.split_rule_id != PREDECLARED_SPLIT_RULE_ID:
            raise ValueError(f"unknown split rule id: {descriptor.split_rule_id}")
        if descriptor.columns.accession_column is None:
            raise ValueError("split_rule_id declared without an accession_column to split on")
        split = predeclared_split(row[descriptor.columns.accession_column])

    target_identity_status, target_context_type, target_start, target_end = _derive_target_identity(descriptor, row)

    return DerivedObservation(
        architecture=descriptor.architecture,
        assay_label=descriptor.assay_label.resolve(row),
        measured_endpoint=descriptor.measured_endpoint,
        measured_value=measured_value,
        full_guide_sequence=full_guide_sequence,
        guide_length=guide_length,
        paired_guide_sequence=paired_guide_sequence,
        paired_slice_start_1based=1,
        paired_length=paired_length,
        passenger_sequence=passenger_sequence,
        guide_3p_overhang=guide_3p_overhang,
        passenger_3p_overhang=passenger_3p_overhang,
        duplex_pairing_status=DuplexPairingStatus(descriptor.architecture.value),
        target_identity_status=target_identity_status,
        target_context_type=target_context_type,
        target_start_1based=target_start,
        target_end_1based=target_end,
        compatibility_status=status,
        compatibility_reason=reason,
        split=split,
    )


# --------------------------------------------------------------------------------------------------
# Registry
#
# Every citation below is something this repository can be made to show. Huesken's is reused verbatim
# from tests/unit/data/README.md and Ichihara's from docs/models_and_scoring.md, both vetted here
# before #109. For Martinelli and Shmushkovich the whole of what this repository holds is a DOI --
# recorded per row in tests/data/benchmarks/oligogym/records.csv's source_url column -- so the DOI is
# what is recorded, with no author list, title or journal inferred from it. That replaces the
# placeholder strings #109 shipped when no bytes were vendored; it is not a claim that anyone here
# has read the papers.
#
# Row counts, the guide convention and the synthetic-flank rule are read out of
# tests/data/benchmarks/oligogym/manifest.json and that directory's README.md, and the duplex
# geometries below are measured off records.csv itself (each descriptor states the count it rests
# on); tests/unit/test_benchmark_panels.py re-derives every one of those numbers from the vendored
# bytes so a regenerated records.csv cannot leave a stale claim standing here.
# --------------------------------------------------------------------------------------------------

#: The one vendored table the three OligoGym-derived panels share; see :class:`PanelRowSelector`.
_OLIGOGYM_RECORDS_CSV = "tests/data/benchmarks/oligogym/records.csv"

#: What every OligoGym-derived descriptor's ``redistribution`` has to say, before its per-panel part.
#: The synthetic flanks are the load-bearing disclosure: the rows are runnable design inputs, and
#: their coordinates are positions in a context this repository fabricated, not in a transcript.
_OLIGOGYM_REDISTRIBUTION_PREFIX = (
    "Third-party redistribution: OligoGym, adapted by scripts/prepare_oligo_benchmarks.py into "
    f"{_OLIGOGYM_RECORDS_CSV} (4,113 rows; see tests/data/benchmarks/README.md). Guides are "
    "OligoGym's second fasta component / HELM RNA2 antisense strand. Each row's design context is "
    "synthetic -- 70 A + reverse-complement(guide) + 70 A, site at 1-based 71 -- so labels and "
    "sequences are as measured while the coordinates are not native; the Ensembl-mapped counterpart "
    "is tests/data/benchmarks/oligogym_native_design/ (#110), which this module does not read. The "
    "upstream OligoGym extract each row names in source_file "
    "(tests/data/external/oligogym/*.csv.gz) is NOT vendored here, and no primary-source full text "
    "has been read in this repository. "
)

_HUESKEN_CITATION = (
    "Huesken D, Lange J, Mickanin C, Weiler J, Asselbergs F, Warner J, Meloon B, Engel S, Rosenberg "
    'A, Cohen D, Labow M, Reinhardt M, Natt F, Hall J. "Design of a genome-wide siRNA library using '
    'an artificial neural network." Nat Biotechnol. 2005 Aug;23(8):995-1001. doi:10.1038/nbt1118. '
    "PMID: 16025102."
)
_HUESKEN_REDISTRIBUTION = (
    "Third-party redistribution: https://github.com/apkrfi/unMod-siRNA-Pred, file "
    "raw_data/All_Dataset.csv. The primary paper is not open access, so the efficacy values were not "
    "verified against it -- see tests/unit/data/README.md."
)

HUESKEN_SUBSET = PanelDescriptor(
    panel_id="huesken_subset",
    display_name="Huesken 2005 subset (180-row vendored redistribution)",
    citation=_HUESKEN_CITATION,
    redistribution=_HUESKEN_REDISTRIBUTION,
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=PanelColumnMapping(
        guide_column="guide_sequence",
        accession_column="accession",
        measured_value_column="efficacy",
    ),
    assay_label=AssayLabelSource(constant="huesken_2005_knockdown_inhibition"),
    measured_endpoint="inhibition_fraction",
    split_rule_id=PREDECLARED_SPLIT_RULE_ID,
    data_present=True,
    vendored_csv="tests/unit/data/sirna_efficacy_subset.csv",
)

HUESKEN_FULL = PanelDescriptor(
    panel_id="huesken_full",
    display_name="Huesken 2005 full panel (2,816 rows / 41 accessions; untracked)",
    citation=_HUESKEN_CITATION,
    redistribution=(
        _HUESKEN_REDISTRIBUTION
        + " The full panel (work/sirna_bench.csv) is untracked in this repository; only the 180-row "
        "huesken_subset is vendored. scripts/validate_scoring_profiles.py reads the full panel under "
        "these column names for its own calibration record, which this descriptor's column mapping "
        "matches so the two agree on what the panel is."
    ),
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=PanelColumnMapping(
        guide_column="siRNA_seq",
        accession_column="B",
        measured_value_column="efficacy",
    ),
    assay_label=AssayLabelSource(constant="huesken_2005_knockdown_inhibition"),
    measured_endpoint="inhibition_fraction",
    split_rule_id=PREDECLARED_SPLIT_RULE_ID,
    data_present=False,
)


def _oligogym_columns(measured_value_column: str) -> PanelColumnMapping:
    """The records.csv mapping every OligoGym-derived panel shares, spelled once.

    Parameterised by the one column the panels disagree about: Ichihara and Martinelli read
    ``efficacy_higher_is_better``, which is byte-identical to their ``label_raw``/``label_processed``;
    Shmushkovich reads ``label_processed``, because for it that same column is a derived direction
    flip. See :class:`PanelColumnMapping`'s table for the field-by-field mapping.

    ``accession_column`` is deliberately absent: records.csv carries no accession, so none of these
    panels declares a split rule -- a split invented here would not be the predeclared one.
    """
    return PanelColumnMapping(
        guide_column="guide_sequence",
        passenger_column="passenger_sequence",
        measured_value_column=measured_value_column,
        context_type_column="context_type",
        target_site_start_column="site_start",
        target_site_end_column="site_end",
    )


ICHIHARA = PanelDescriptor(
    panel_id="ichihara",
    display_name="Ichihara 2007 thermodynamic duplex panel (2,850 OligoGym-derived rows)",
    citation=(
        "Ichihara M, Murakumo Y, Masuda A, Matsuura T, Asai N, Jijiwa M, Ishida M, Shinmi J, "
        "Yatsuya H, Qiao S, Takahashi M, Ohno K. Thermodynamic instability of siRNA duplex is a "
        "prerequisite for dependable prediction of siRNA activities. Nucleic Acids Research. "
        "2007;35(18):e123. Corroborated by the vendored bytes: every row of both datasets records "
        "source_url https://doi.org/10.1093/nar/gkm699, which resolves to this article."
    ),
    redistribution=(
        _OLIGOGYM_REDISTRIBUTION_PREFIX + "This panel is OligoGym's two Ichihara datasets, "
        "ichihara_2007_1 (2,431 rows) and ichihara_2007_2 (419 rows), which carry the same DOI."
    ),
    # Measured off records.csv, not assumed: all 2,850 rows are a 21 nt guide with a 19 nt passenger
    # equal to reverse-complement(guide[:19]) -- a 19 nt paired core plus a measured 2 nt 3' guide
    # overhang, the same geometry as the vendored Huesken subset.
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=_oligogym_columns("efficacy_higher_is_better"),
    # A constant, not the label_semantics column: that column holds one value across all 2,850 rows,
    # so a per-row read would add provenance noise without adding information.
    assay_label=AssayLabelSource(constant="ichihara_2007_percent_inhibition"),
    # Percent (the vendored rows span -27.8..134.1), not Huesken's 0-1 inhibition_fraction. Named so
    # the difference is legible in observations.csv: the two panels' values must never be pooled.
    measured_endpoint="percent_inhibition_relative_to_control_higher_is_better",
    split_rule_id=None,
    data_present=True,
    vendored_csv=_OLIGOGYM_RECORDS_CSV,
    row_selector=PanelRowSelector(column="dataset", values=("ichihara_2007_1", "ichihara_2007_2")),
    target_identity_status=TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL,
)

MARTINELLI = PanelDescriptor(
    panel_id="martinelli",
    display_name="Martinelli 21-mer duplex panel (907 OligoGym-derived rows)",
    citation=(
        "The only citation this repository can read out of its own bytes is the DOI every row of "
        "OligoGym's martinelli_2023_1 dataset carries in source_url: "
        "https://doi.org/10.1016/j.ygeno.2024.110815. No author list, title or journal is recorded "
        "here, because none has been read here; note also that the dataset name says 2023 while the "
        "DOI string says 2024, and this repository holds nothing that resolves which is the "
        "publication year."
    ),
    redistribution=(_OLIGOGYM_REDISTRIBUTION_PREFIX + "This panel is OligoGym's martinelli_2023_1 dataset (907 rows)."),
    # Re-declared against the vendored bytes, which contradict the fully-complementary 21-mer #109
    # assumed from the issue text alone: of the 907 rows, 858 have passenger[:19] ==
    # reverse-complement(guide[:19]) with both strands 21 nt -- a 19 nt paired core with a measured
    # 2 nt 3' overhang on each strand -- 18 are blunt-complementary over the full 21 nt, and 31 match
    # neither exactly (6 with 1-3 mismatches inside that 19 nt window, 25 in another register).
    # Declaring FULLY_COMPLEMENTARY/21 would have recorded guide_3p_overhang="" ("measured and blunt")
    # for 889 rows whose bytes say otherwise, which is the misrepresentation #109 forbids; the
    # paired-core declaration is the geometry 94.6% of the rows show and understates the rest rather
    # than inventing pairing for them.
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=_oligogym_columns("efficacy_higher_is_better"),
    assay_label=AssayLabelSource(constant="martinelli_2023_percent_knockdown"),
    measured_endpoint="percent_knockdown_of_target_mrna_higher_is_better",
    split_rule_id=None,
    data_present=True,
    vendored_csv=_OLIGOGYM_RECORDS_CSV,
    row_selector=PanelRowSelector(column="dataset", values=("martinelli_2023_1",)),
    target_identity_status=TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL,
)

SHMUSHKOVICH = PanelDescriptor(
    panel_id="shmushkovich",
    display_name="Shmushkovich asymmetric 15/20 hsiRNA panel (356 OligoGym-derived rows, all incompatible)",
    citation=(
        "The only citation this repository can read out of its own bytes is the DOI every row of "
        "OligoGym's shmushkovich_2018_1 dataset carries in source_url: "
        "https://doi.org/10.1093/nar/gky745. No author list, title or journal is recorded here, "
        "because none has been read here."
    ),
    redistribution=(
        _OLIGOGYM_REDISTRIBUTION_PREFIX + "This panel is OligoGym's shmushkovich_2018_1 dataset (356 rows). "
        "Every one of them is incompatible with #109's fixed-length interface by construction; they "
        "are vendored, and ingested, to be recorded and refused, never to be sliced into a "
        "19-23 nt core."
    ),
    # Confirmed against the vendored bytes: all 356 rows are a 20 nt guide with a 15 nt passenger
    # equal to reverse-complement(guide[:15]) -- the 15/20 asymmetric hsiRNA duplex #109 excludes,
    # with a 5 nt single-stranded guide tail. #110 owns giving it a real path.
    architecture=PanelArchitecture.ASYMMETRIC,
    declared_paired_length=None,
    # label_processed, the panel's own direction (percentage of target mRNA remaining, lower is
    # better). records.csv also carries efficacy_higher_is_better, which for this dataset alone is a
    # derived 100 - label_processed flip (exact on all 356 rows); reading it here would record a
    # transformation as a measurement, so this descriptor deliberately does not.
    columns=_oligogym_columns("label_processed"),
    assay_label=AssayLabelSource(constant="shmushkovich_2018_percent_target_remaining"),
    measured_endpoint="percent_target_mrna_remaining_lower_is_better",
    split_rule_id=None,
    data_present=True,
    vendored_csv=_OLIGOGYM_RECORDS_CSV,
    row_selector=PanelRowSelector(column="dataset", values=("shmushkovich_2018_1",)),
    target_identity_status=TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL,
)

OLIGOGYM = PanelDescriptor(
    panel_id="oligogym",
    display_name="OligoGym aggregate table (4,113 rows; prepare a member panel instead)",
    citation=(
        "OligoGym is the aggregator the three panels above were derived from, and this repository "
        "holds no citation for OligoGym itself -- no DOI, no URL, no upstream extract "
        "(tests/data/external/oligogym/ is absent). What it holds is the adapter that wrote the "
        "vendored table, scripts/prepare_oligo_benchmarks.py, and the four member DOIs its rows "
        "carry: https://doi.org/10.1093/nar/gkm699 (Ichihara, both datasets), "
        "https://doi.org/10.1016/j.ygeno.2024.110815 (Martinelli) and "
        "https://doi.org/10.1093/nar/gky745 (Shmushkovich)."
    ),
    redistribution=(
        _OLIGOGYM_REDISTRIBUTION_PREFIX + "This descriptor is the whole table rather than one "
        "dataset, so it is an aggregate: derive_observation refuses it and names its members."
    ),
    # Recorded because the schema requires an architecture, and never read: derive_observation
    # refuses an aggregate descriptor before it looks at one. It is the majority geometry (3,708 of
    # the 4,113 rows), not a description of the table -- the other 405 rows are Martinelli's 49
    # non-core-19 rows and all 356 asymmetric Shmushkovich rows, which is exactly why a single
    # architecture cannot stand for this file and why the members are the ingestible unit.
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=_oligogym_columns("efficacy_higher_is_better"),
    assay_label=AssayLabelSource(constant="oligogym_aggregate_not_ingested"),
    measured_endpoint="mixed_per_member_endpoint_see_member_panels",
    split_rule_id=None,
    data_present=True,
    vendored_csv=_OLIGOGYM_RECORDS_CSV,
    row_selector=None,
    target_identity_status=TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL,
    aggregate_of=("ichihara", "martinelli", "shmushkovich"),
)

#: Every declared panel, keyed by ``panel_id``. Six entries: five vendored (``huesken_subset`` plus
#: the four OligoGym-derived ones, of which ``oligogym`` is an aggregate of the other three) and one,
#: ``huesken_full``, whose table is untracked here.
PANEL_REGISTRY: Mapping[str, PanelDescriptor] = {
    d.panel_id: d for d in (HUESKEN_SUBSET, HUESKEN_FULL, ICHIHARA, MARTINELLI, SHMUSHKOVICH, OLIGOGYM)
}


def describe_panel(panel_id: str) -> PanelDescriptor:
    """Registry lookup that names every valid id on failure, rather than a bare ``KeyError``."""
    try:
        return PANEL_REGISTRY[panel_id]
    except KeyError as exc:
        valid = ", ".join(sorted(PANEL_REGISTRY))
        raise ValueError(f"unknown benchmark panel_id {panel_id!r}; valid ids: {valid}") from exc
