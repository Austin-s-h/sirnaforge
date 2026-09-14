"""Benchmark panel registry: architectures, compatibility verdicts, honestly-labelled fixtures.

Issue #109. A "panel" is a named siRNA/miRNA efficacy dataset the benchmark surface can prepare an
artifact from. This module is the one place a panel is described: its citation, its declared duplex
architecture, the columns a source table must carry, and whether this repository actually vendors
its bytes. ``sirnaforge.benchmark.prepare`` reads a descriptor and a raw source row and calls
:func:`derive_observation` here for every field the panel/architecture determines; it does not
re-derive a slice or a compatibility verdict itself, so the rule lives in exactly one place.

Two verified facts this module works around rather than papering over (full accounting in
``tests/unit/data/benchmark/README.md``):

1. The issue's cited source, ``docs/prd_benchmark_artifacts_and_variable_length.md``, does not exist
   in this repository -- not in the working tree, not in any branch's history. The issue body is the
   entire specification; nothing below claims to summarise a document that was read.
2. Of the five panels the issue names -- Huesken, Ichihara, Martinelli, Shmushkovich, OligoGym --
   only a 180-row Huesken redistribution (``huesken_subset``, ``tests/unit/data/sirna_efficacy_subset
   .csv``) ships real bytes in this repository. The other five descriptors below are declared so the
   registry, a future CLI and a manifest can name them and record ``data_present=False`` honestly;
   for those, ``columns`` is a *contract* a future ``--panel-csv`` would need to satisfy, not a
   description of a file anyone here has read. Nothing in this module reads or writes real Ichihara,
   Martinelli, Shmushkovich or OligoGym bytes, and none are vendored under ``tests/unit/data/``.

Architectures (:class:`PanelArchitecture`) are the three duplex geometries #109's fixed-length
design surface can consume. ``ASYMMETRIC`` is excluded from every compatibility verdict by
construction, not merely under-supported: #109's scope excludes variable-length/asymmetric design
entirely (#110 owns giving it a real, native path), and silently truncating or padding a 15/20
duplex into a 19-23 nt paired core would misrepresent a measured sequence, which the parent issue
also forbids. :func:`derive_observation` therefore never slices an asymmetric record -- it returns
the verbatim guide with ``compatibility_status="incompatible"``.
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
    """

    guide_column: str = Field(min_length=1)
    passenger_column: str | None = None
    accession_column: str | None = None
    measured_value_column: str | None = None

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
        data_present: Whether this repository vendors real bytes for this panel. ``False`` for every
            descriptor except ``huesken_subset``; see the module docstring.
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

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def _paired_length_matches_architecture(self) -> PanelDescriptor:
        if self.architecture is PanelArchitecture.ASYMMETRIC:
            if self.declared_paired_length is not None:
                raise ValueError("an asymmetric panel declares no fixed paired length")
        elif self.declared_paired_length is None:
            raise ValueError(f"a {self.architecture.value} panel must declare a paired length")
        return self

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
    compatibility_status: CompatibilityStatus
    compatibility_reason: str
    split: str | None

    model_config = ConfigDict(frozen=True, extra="forbid")


_ALLOWED_BASES = frozenset("ACGTUN")


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
    """
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
        compatibility_status=status,
        compatibility_reason=reason,
        split=split,
    )


# --------------------------------------------------------------------------------------------------
# Registry
#
# Citations below are reused verbatim from tests/unit/data/README.md (Huesken) and
# docs/models_and_scoring.md (Ichihara), both already vetted in this repository. Martinelli,
# Shmushkovich and OligoGym have no citation vetted anywhere in this repository: the strings below
# name them as the issue does and say so explicitly, rather than inventing a DOI or PMID this
# repository cannot back. See the module docstring's second verified fact and
# tests/unit/data/benchmark/README.md.
# --------------------------------------------------------------------------------------------------

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

ICHIHARA = PanelDescriptor(
    panel_id="ichihara",
    display_name="Ichihara 2007 thermodynamic duplex panel",
    citation=(
        "Ichihara M, Murakumo Y, Masuda A, Matsuura T, Asai N, Jijiwa M, Ishida M, Shinmi J, "
        "Yatsuya H, Qiao S, Takahashi M, Ohno K. Thermodynamic instability of siRNA duplex is a "
        "prerequisite for dependable prediction of siRNA activities. Nucleic Acids Research. "
        "2007;35(18):e123."
    ),
    redistribution=None,
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=PanelColumnMapping(
        guide_column="guide_sequence",
        accession_column="accession",
        measured_value_column="efficacy",
    ),
    assay_label=AssayLabelSource(constant="ichihara_2007_placeholder_assay"),
    measured_endpoint="inhibition_fraction",
    split_rule_id=None,
    data_present=False,
)

MARTINELLI = PanelDescriptor(
    panel_id="martinelli",
    display_name="Martinelli fully complementary duplex panel",
    citation=(
        "Martinelli et al. -- named in issue #109/#110 for its fully complementary 21-nt duplex "
        "records. No primary citation has been independently verified in this repository: no bytes "
        "are vendored, and this repository's history holds no prior reference to this panel."
    ),
    redistribution=None,
    architecture=PanelArchitecture.FULLY_COMPLEMENTARY,
    declared_paired_length=21,
    columns=PanelColumnMapping(
        guide_column="guide_sequence",
        passenger_column="passenger_sequence",
        accession_column="accession",
        measured_value_column="efficacy",
    ),
    assay_label=AssayLabelSource(constant="martinelli_placeholder_assay"),
    measured_endpoint="efficacy_unspecified_placeholder",
    split_rule_id=None,
    data_present=False,
)

SHMUSHKOVICH = PanelDescriptor(
    panel_id="shmushkovich",
    display_name="Shmushkovich asymmetric 15/20 duplex panel",
    citation=(
        "Shmushkovich et al. -- named in issue #109/#110 for its asymmetric 15/20 duplex "
        "architecture. No primary citation has been independently verified in this repository: no "
        "bytes are vendored, and this repository's history holds no prior reference to this panel."
    ),
    redistribution=None,
    architecture=PanelArchitecture.ASYMMETRIC,
    declared_paired_length=None,
    columns=PanelColumnMapping(
        guide_column="guide_sequence",
        passenger_column="passenger_sequence",
        accession_column="accession",
        measured_value_column="efficacy",
    ),
    assay_label=AssayLabelSource(constant="shmushkovich_placeholder_assay"),
    measured_endpoint="efficacy_unspecified_placeholder",
    split_rule_id=None,
    data_present=False,
)

OLIGOGYM = PanelDescriptor(
    panel_id="oligogym",
    display_name="OligoGym benchmark aggregator",
    citation=(
        "OligoGym -- named in issue #109/#110 as an aggregator panel; the issue text does not name a "
        "single architecture for it. No primary citation, source or bytes have been located in this "
        "repository, and no prior reference to it exists in this repository's history. The "
        "paired-core architecture and length declared below are a placeholder pending #110, not a "
        "description of a file anyone here has read."
    ),
    redistribution=None,
    architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
    declared_paired_length=19,
    columns=PanelColumnMapping(
        guide_column="guide_sequence",
        accession_column="accession",
        measured_value_column="efficacy",
    ),
    assay_label=AssayLabelSource(constant="oligogym_placeholder_assay"),
    measured_endpoint="efficacy_unspecified_placeholder",
    split_rule_id=None,
    data_present=False,
)

#: Every declared panel, keyed by ``panel_id``. Six entries: one real (``huesken_subset``) and five
#: declared-but-not-vendored, per the module docstring's second verified fact.
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
