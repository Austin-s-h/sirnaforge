"""Benchmark artifact/manifest schema and lossless CSV/JSON round trip (#109, bm-artifact slice).

An artifact is one directory per (compatible panel x paired length):
``<panel_id>__len<paired_length>/`` holding ``manifest.json``, ``observations.csv`` (written by
``benchmark prepare``) and, once ``benchmark design`` has run, ``design_inputs.fasta``,
``candidates_all.csv`` and ``accounting.csv``. This module owns the three row/manifest schemas
(:class:`BenchmarkObservation`, :class:`BenchmarkAccountingRow`, :class:`BenchmarkArtifactManifest`),
the two schema-version constants, the directory-name helper, the fixed inner filenames, and the
read/write helpers for the two CSVs plus the manifest. It does not read a panel, run a design or
build a manifest -- those are ``prepare``/``design``, in modules this file does not own.

Two properties make two independent ``prepare`` runs on the same inputs byte-identical, which is
the regression assertion the artifact buys: rows are always written sorted by their join key, and
no column anywhere carries a timestamp or an absolute path (``manifest.json``'s ``created_utc`` is
the one deliberate exception, and it lives outside both CSVs).

**Why the two overhang columns need a sentinel and nothing else does.** RFC 4180 does not
distinguish an empty field from a missing one -- ``csv`` hands back ``""`` for both, the same way
Postgres's ``COPY`` needed its own backslash-N marker for exactly this gap. Every optional column in
this schema uses that ambiguity harmlessly, because an empty cell and "not stated" are the same
fact everywhere except ``guide_3p_overhang``/``passenger_3p_overhang``: there, an empty string
(measured, blunt) and ``None`` (the panel states no overhang) are two different, both-real facts
that both have to round-trip, so ``_OVERHANG_UNSTATED`` marks the ``None`` case and a plain empty
cell is left to mean blunt.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from sirnaforge.config.run_policy import POLICY_SCHEMA_VERSION
from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation, SettingSource

#: The two CSV column contracts (``observations.csv``, ``accounting.csv``). Bump on any column
#: addition, removal, reorder or type change -- a reader keys its parsing off this, not off
#: sniffing the header.
BENCHMARK_ARTIFACT_SCHEMA_VERSION = "1.0"

#: ``manifest.json``'s own shape. Separate from the above because #110 is expected to add fields
#: to one without the other.
BENCHMARK_MANIFEST_SCHEMA_VERSION = "1.0"

#: Fixed inner filenames. ``manifest.json`` matches the name ``workflow.py::_build_fair_manifest``
#: already writes, so a reader has one habit across both artifact kinds.
MANIFEST_FILENAME = "manifest.json"
OBSERVATIONS_FILENAME = "observations.csv"
DESIGN_INPUTS_FASTA_FILENAME = "design_inputs.fasta"
CANDIDATES_ALL_FILENAME = "candidates_all.csv"
ACCOUNTING_FILENAME = "accounting.csv"

#: The design length range ``sirnaforge design`` already accepts; an artifact outside it could
#: never be consumed by the fixed-length path it exists to feed.
PAIRED_LENGTH_BOUNDS: tuple[int, int] = (19, 23)

#: Directory-name convention: ``benchmark/panels.py`` panel ids are lowercase-and-underscore by
#: construction, so a name outside this pattern is a caller error, not a panel this build knows.
_PANEL_ID_PATTERN_MESSAGE = "must match ^[a-z0-9_]+$ (lowercase letters, digits, underscore)"


def _panel_id_is_valid(panel_id: str) -> bool:
    """Whether ``panel_id`` matches the declared ``[a-z0-9_]+`` directory-name convention."""
    return bool(panel_id) and all(char in "abcdefghijklmnopqrstuvwxyz0123456789_" for char in panel_id)


class BenchmarkArtifactError(ValueError):
    """A benchmark artifact file cannot be trusted as written, raised before any row is silently defaulted.

    Kept a ``ValueError`` so it can share a catch clause with Pydantic's own :class:`ValidationError`
    at whichever CLI boundary later slices add, the way :class:`~sirnaforge.config.run_policy.RunPolicyError`
    does for policy resolution. Raised instead of filling in a placeholder for a missing declared
    column: a silent default for one panel's overhang or citation column would be indistinguishable
    from a real value in every report read off the artifact afterwards.
    """


def artifact_dir_name(panel_id: str, paired_length: int) -> str:
    """The fixed ``<panel_id>__len<paired_length>`` directory name for one artifact.

    One directory per (compatible panel x paired length) is the on-disk unit #109 defines; this is
    the single place that spells it, so ``prepare`` and ``design`` cannot drift on the convention.

    Args:
        panel_id: Registry panel id. Must match ``[a-z0-9_]+``.
        paired_length: The design length the artifact is paired to, 19-23.

    Returns:
        The directory name, e.g. ``huesken_full__len21``.

    Raises:
        BenchmarkArtifactError: ``panel_id`` or ``paired_length`` is not a value this build can
            build a directory for.
    """
    if not _panel_id_is_valid(panel_id):
        raise BenchmarkArtifactError(f"panel_id {panel_id!r} {_PANEL_ID_PATTERN_MESSAGE}")
    low, high = PAIRED_LENGTH_BOUNDS
    if not low <= paired_length <= high:
        raise BenchmarkArtifactError(f"paired_length {paired_length} must be between {low} and {high}")
    return f"{panel_id}__len{paired_length}"


# --------------------------------------------------------------------------------------------------
# Shared vocabulary. Declared once so BenchmarkObservation and BenchmarkArtifactManifest cannot each
# spell the same closed set of strings slightly differently.
# --------------------------------------------------------------------------------------------------

ArchitectureId = Literal["paired_core_with_overhang", "fully_complementary", "asymmetric"]
DuplexPairingStatus = Literal["fully_complementary", "paired_core_with_overhang", "asymmetric", "unstated"]
CompatibilityStatus = Literal["compatible", "incompatible"]
#: Deliberately has no ``confirmed`` member: #110 owns native transcript mapping, and #109 must be
#: structurally unable to claim it. ``synthetic_context_local`` is the honest third state for a site
#: whose coordinates are real and reproducible but locate it inside a FABRICATED context -- the
#: OligoGym design inputs are 70 A + revcomp(guide) + 70 A with the site at 1-based 71. ``panel_local``
#: would read as a coordinate in the panel's own measured target, which overclaims; ``unavailable``
#: would discard a coordinate the artifact can honestly reproduce, which underclaims.
TargetIdentityStatus = Literal["unavailable", "panel_local", "synthetic_context_local"]
DesignContextSource = Literal["panel_transcript", "measured_target_site"]
SplitLabel = Literal["development", "held_out"]
GuideMatch = Literal["exact", "paired_core_exact", "none"]

#: The sentinel that marks "the panel states no overhang" in ``guide_3p_overhang``/
#: ``passenger_3p_overhang``. Never a legal overhang value (those are empty or nucleotide/
#: modification codes), so it cannot collide with real data; a plain empty CSV cell is left free to
#: mean "measured, blunt".
_OVERHANG_UNSTATED = "<not_stated>"


class BenchmarkObservation(BaseModel):
    """One measured row from a benchmark panel, as it exists before any design candidate is built.

    Written by ``benchmark prepare`` into ``observations.csv``, one row per source row, compatible
    or not: an observation this build cannot pair to a design length is still recorded, with
    ``compatibility_status="incompatible"`` and a reason, rather than dropped. That is what keeps
    "the original observation when no design candidate is produced" true regardless of what design
    does later -- #109's acceptance criterion for that phrase is met by this row existing at all,
    not by anything ``accounting.csv`` adds on top of it.

    Frozen and ``extra="forbid"``: an artifact is read back into exactly the object that wrote it,
    or the reader raises, never guesses.
    """

    observation_id: str = Field(min_length=1, description="`{panel_id}:{source_row_index:06d}`")
    panel_id: str = Field(min_length=1, description="Registry panel id")
    source_row_index: int = Field(ge=0, description="0-based index in the source table")
    architecture: ArchitectureId = Field(description="Registry architecture id")
    assay_label: str = Field(min_length=1, description="Per-row label, or the descriptor's constant")
    measured_endpoint: str = Field(min_length=1, description="What measured_value counts, e.g. inhibition_fraction")
    measured_value: float | None = Field(description="None when the panel records no numeric endpoint")
    full_guide_sequence: str = Field(min_length=1, description="As measured, upper-cased only; never padded/trimmed")
    guide_length: int = Field(ge=1, description="len(full_guide_sequence)")
    paired_guide_sequence: str = Field(min_length=1, description="A slice of full_guide_sequence")
    paired_slice_start_1based: int = Field(ge=1, description="1-based start of the slice, for auditability")
    paired_length: int = Field(description="len(paired_guide_sequence); equals the artifact directory's len tag")
    passenger_sequence: str | None = Field(description="None = the panel measured only a guide")
    guide_3p_overhang: str | None = Field(description='"" = measured and blunt; None = panel states no overhang')
    passenger_3p_overhang: str | None = Field(description="Same convention as guide_3p_overhang")
    duplex_pairing_status: DuplexPairingStatus = Field(description="How the panel describes strand pairing")
    source_citation: str = Field(min_length=1, description="The primary citation")
    source_redistribution: str | None = Field(description="None only when obtained from the primary source")
    source_file: str = Field(min_length=1, description="Path of the source table as given")
    source_sha256: str = Field(min_length=1, description="SHA-256 of the source table's bytes")
    target_transcript_id: str | None = Field(description="#110 placeholder; None in #109")
    target_identity_status: TargetIdentityStatus = Field(
        description="unavailable, panel_local or synthetic_context_local; never confirmed -- #110 owns native"
    )
    target_start_1based: int | None = Field(description="Populated only for panel_local")
    target_end_1based: int | None = Field(description="Populated only for panel_local")
    target_strand: str | None = Field(description="Populated only for panel_local")
    design_context_id: str = Field(min_length=1, description="Record id in design_inputs.fasta")
    design_context_source: DesignContextSource = Field(description="panel_transcript or measured_target_site")
    split: SplitLabel | None = Field(description="Only when the descriptor declares a split rule")
    compatibility_status: CompatibilityStatus = Field(description="Whether this row entered the artifact's design")
    compatibility_reason: str = Field(description='"" when compatible')

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def _sequences_are_upper_case(self) -> BenchmarkObservation:
        """`full_guide_sequence`/`paired_guide_sequence` are stated "as measured, upper-cased only".

        A lower-case letter here would mean some caller normalised or transcribed the sequence
        rather than copying it, which is exactly the silent relabelling #109 rules out.
        """
        for name in ("full_guide_sequence", "paired_guide_sequence"):
            value = getattr(self, name)
            if value != value.upper():
                raise ValueError(f"{name} must be upper-cased only: {value!r}")
        return self

    @model_validator(mode="after")
    def _compatibility_reason_matches_status(self) -> BenchmarkObservation:
        """A compatible row states no reason; an incompatible one must state one.

        The pairing is what makes ``compatibility_reason`` machine-checkable rather than a free-text
        field a reader has to trust: "" claims nothing went wrong, so nothing went wrong is the only
        value it may hold.
        """
        if self.compatibility_status == "compatible" and self.compatibility_reason != "":
            raise ValueError("compatibility_reason must be '' when compatibility_status is 'compatible'")
        if self.compatibility_status == "incompatible" and self.compatibility_reason == "":
            raise ValueError("compatibility_reason must be non-empty when compatibility_status is 'incompatible'")
        return self


class BenchmarkAccountingRow(BaseModel):
    """One observation<->candidate join, carrying both filter verdict sets for that pairing.

    Deliberately narrow: candidate columns stay owned by
    :func:`sirnaforge.models.sirna.build_candidate_row` in ``candidates_all.csv``, and this row
    joins to it on ``candidate_id`` rather than re-emitting them, so the candidate column contract
    cannot drift by having two writers. An observation with no candidate still gets a row here
    (``entered_design=False``, ``candidate_id=None``), which is the other half of never dropping an
    observation -- ``observations.csv`` never omits the row, this file never omits its join.

    ``default_filter_status``/``benchmark_filter_status`` are re-derived, not two independently run
    screens: see :mod:`sirnaforge.benchmark.design`'s module docstring for why re-deriving the
    default verdict from the benchmark run's own recorded evidence is exact rather than approximate
    (widening is monotone, so a narrowing is refused), which is the property that makes emitting both
    sets here meaningful.
    """

    observation_id: str = Field(min_length=1, description="Join key into observations.csv")
    panel_id: str = Field(min_length=1, description="Registry panel id")
    paired_length: int = Field(description="Artifact directory's len tag")
    candidate_id: str | None = Field(description="None when no candidate was designed for this observation")
    designed_guide_sequence: str | None = Field(description="None when no candidate was designed")
    guide_match: GuideMatch = Field(description="RNA-normalised comparison against the measured guide")
    entered_design: bool = Field(description="Whether a candidate was built for this observation at all")
    default_filter_status: FilterEvaluation = Field(description="Aggregate verdict under the default policy")
    default_filter_reasons: str = Field(description='";"-joined failing filter_ids, registry order, "" if none')
    benchmark_filter_status: FilterEvaluation = Field(description="Aggregate verdict under the benchmark policy")
    benchmark_filter_reasons: str = Field(description='";"-joined failing filter_ids, registry order, "" if none')

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def _no_candidate_implies_not_evaluated(self) -> BenchmarkAccountingRow:
        """A gate cannot be blamed for a candidate that was never built (#109 filter-accounting rule).

        ``entered_design=False`` must carry ``candidate_id=None`` and both statuses
        ``not_evaluated`` with empty reasons -- the third fact (``entered_design``) is otherwise
        just a redundant restatement of ``candidate_id is None`` instead of the independent check it
        is meant to be.
        """
        if not self.entered_design:
            if self.candidate_id is not None or self.designed_guide_sequence is not None:
                raise ValueError("entered_design=False must carry candidate_id=None and designed_guide_sequence=None")
            for status, reasons in (
                (self.default_filter_status, self.default_filter_reasons),
                (self.benchmark_filter_status, self.benchmark_filter_reasons),
            ):
                if status is not FilterEvaluation.NOT_EVALUATED or reasons != "":
                    raise ValueError("entered_design=False must carry not_evaluated statuses and empty reason columns")
        return self


# --------------------------------------------------------------------------------------------------
# manifest.json
# --------------------------------------------------------------------------------------------------


class ManifestPanelBlock(BaseModel):
    """Which panel this artifact was prepared from, and whether its bytes are actually present."""

    panel_id: str = Field(min_length=1)
    display_name: str = Field(min_length=1)
    architecture: ArchitectureId
    citation: str = Field(min_length=1)
    redistribution: str | None = Field(description="None only when obtained from the primary source")
    data_present: bool = Field(description="False for a named-but-not-vendored panel (Ichihara/Martinelli/...)")
    descriptor_hash: str = Field(min_length=1, description="sha256:... over the registry descriptor")

    model_config = ConfigDict(frozen=True, extra="forbid")


class ManifestInputEntry(BaseModel):
    """One checksummed input the artifact was built from, keyed by a role unique within ``inputs``.

    ``role`` is a key, not a label: ``design_artifact`` *replaces* the entry carrying a role it is
    re-recording rather than appending a second one, so re-running ``benchmark design`` over an
    artifact leaves one entry per input instead of growing the provenance list on every run.
    """

    role: str = Field(min_length=1, description="Unique within `inputs`; e.g. panel_csv, prepared_observations_csv")
    path: str = Field(min_length=1)
    sha256: str = Field(min_length=1)
    size_bytes: int = Field(ge=0)

    model_config = ConfigDict(frozen=True, extra="forbid")


class ManifestOutputEntry(BaseModel):
    """One artifact file's checksum and shape, in the shape ``workflow.py``'s ``add_file`` already uses.

    ``size_bytes``/``sha256`` and ``rows``/``sequences`` are ``None`` together with ``exists=False``:
    a file that does not exist has nothing to checksum or count, and a placeholder number would be
    read as a checksum of nothing.

    ``path`` is the file's name *relative to the artifact directory* -- always one of the fixed inner
    filenames this module declares, never an absolute path. Two prepares of the same panel bytes into
    two different ``--out-dir`` roots must produce the same manifest apart from ``created_utc``, and
    an absolute path would make the artifact's own reproducibility claim depend on where the run
    happened to be writing (#109).
    """

    path: str = Field(min_length=1)
    exists: bool
    size_bytes: int | None = Field(default=None, ge=0)
    sha256: str | None = Field(default=None)
    rows: int | None = Field(default=None, ge=0, description="Set for a CSV output")
    sequences: int | None = Field(default=None, ge=0, description="Set for a FASTA output")

    model_config = ConfigDict(frozen=True, extra="forbid")

    @model_validator(mode="after")
    def _absent_file_carries_no_measurement(self) -> ManifestOutputEntry:
        if not self.exists and any(
            value is not None for value in (self.size_bytes, self.sha256, self.rows, self.sequences)
        ):
            raise ValueError("exists=False must not carry size_bytes/sha256/rows/sequences")
        return self


class BenchmarkArtifactOutputs(BaseModel):
    """The artifact's four fixed files. The last two are ``None`` until ``benchmark design`` runs."""

    observations_csv: ManifestOutputEntry
    design_inputs_fasta: ManifestOutputEntry
    candidates_all_csv: ManifestOutputEntry | None = Field(default=None, description="Set once design has run")
    accounting_csv: ManifestOutputEntry | None = Field(default=None, description="Set once design has run")

    model_config = ConfigDict(frozen=True, extra="forbid")


class FilterExclusionCounts(BaseModel):
    """How many candidates one gate excluded, under each policy. The shape both count blocks share."""

    default: int = Field(ge=0)
    benchmark: int = Field(ge=0)

    model_config = ConfigDict(frozen=True, extra="forbid")


class BenchmarkArtifactCounts(BaseModel):
    """Observation/mapping/design/exclusion counts, the manifest's audit summary of the whole run.

    ``mapped_native`` is pinned to 0: native transcript mapping is explicitly out of scope for #109
    (only the placeholder fields exist), and #110 is what gets to make that count real. A non-zero
    value here would be a claim this build cannot back up.
    """

    observations_in_source: int = Field(ge=0)
    observations_kept: int = Field(ge=0)
    observations_incompatible: int = Field(ge=0)
    mapped_native: int = Field(ge=0, le=0, description="Always 0 in #109; #110 owns native mapping")
    mapped_panel_local: int = Field(ge=0)
    mapping_unavailable: int = Field(ge=0)
    entered_design: int = Field(ge=0)
    no_candidate: int = Field(ge=0)
    default_pass: int = Field(ge=0)
    benchmark_pass: int = Field(ge=0)
    excluded_by_filter: dict[str, FilterExclusionCounts] = Field(default_factory=dict)

    model_config = ConfigDict(frozen=True, extra="forbid")


class PolynucleotideRunRequirementBlock(BaseModel):
    """The one gate #109 requires to stay active and auditable everywhere: ``max_poly_runs <= 3``.

    Recorded here, in ``BenchmarkAccountingRow.default_filter_reasons``/``benchmark_filter_reasons``
    per row, and in ``BenchmarkArtifactCounts.excluded_by_filter["max_poly_runs"]`` per run --
    three independent traces of the same requirement, which is what "explicitly recorded" means for
    a requirement nothing in the benchmark surface is allowed to widen.
    """

    filter_id: Literal["max_poly_runs"] = "max_poly_runs"
    comparator: FilterComparator
    threshold: int | float = Field(description="Matches FilterDescriptor.threshold's type")
    action: FilterAction
    evaluated: bool
    excluded: FilterExclusionCounts

    model_config = ConfigDict(frozen=True, extra="forbid")


class GCWideningEntry(BaseModel):
    """One GC bound's default value, this run's value, and whether the two differ.

    ``source`` is a :class:`~sirnaforge.models.policy.SettingSource`, not a value comparison: #101
    already established that a value test cannot tell an omitted option from one typed with the
    default, so ``widened`` is read off provenance, not off ``benchmark != default``.
    """

    default: float
    benchmark: float
    source: SettingSource
    widened: bool

    model_config = ConfigDict(frozen=True, extra="forbid")


class GCWideningBlock(BaseModel):
    """Both GC bounds' widening record. The only settings the benchmark surface may move."""

    gc_min: GCWideningEntry
    gc_max: GCWideningEntry

    model_config = ConfigDict(frozen=True, extra="forbid")


class BenchmarkArtifactManifest(BaseModel):
    """``manifest.json``: what was run, on what bytes, under which policy, with what result.

    Reproducibility from checksummed inputs rests on five fields together, none of them prose:
    ``invoked_command`` is ``sys.argv`` verbatim (not reconstructed from options); every entry in
    ``inputs`` and ``outputs`` carries its own sha256; ``panel.descriptor_hash`` covers the registry
    entry so an edited descriptor is as visible as edited data; ``run_policy``/``default_run_policy``
    carry the legacy profile's content hash, covering thresholds and per-filter actions both; and no
    CSV carries a timestamp or a path, so matching input hashes plus the same argv must reproduce
    ``observations.csv``/``accounting.csv`` byte-for-byte. ``created_utc`` is the one field this
    schema deliberately leaves nondeterministic, and it lives only here, never in a CSV.
    """

    schema_version: str = Field(default=BENCHMARK_MANIFEST_SCHEMA_VERSION, description="manifest.json's own shape")
    artifact_schema_version: str = Field(
        default=BENCHMARK_ARTIFACT_SCHEMA_VERSION, description="The two CSV column contracts"
    )
    policy_schema_version: str = Field(default=POLICY_SCHEMA_VERSION, description="Imported, never restated")
    tool: Literal["sirnaforge"] = "sirnaforge"
    tool_version: str = Field(min_length=1, description="sirnaforge.__version__ at run time")
    created_utc: str = Field(min_length=1, description="The only nondeterministic field; not in any CSV")
    invoked_command: tuple[str, ...] = Field(description="sys.argv verbatim")
    panel: ManifestPanelBlock
    paired_length: int = Field(description="Equal to requested_length in #109; #110 separates them")
    requested_length: int
    split_rule: str | None = Field(description="e.g. sha256_accession_parity_v1, or None")
    inputs: tuple[ManifestInputEntry, ...] = Field(default_factory=tuple)
    outputs: BenchmarkArtifactOutputs
    counts: BenchmarkArtifactCounts
    polynucleotide_run_requirement: PolynucleotideRunRequirementBlock
    run_policy: dict[str, Any] = Field(description="ResolvedRunPolicy.as_manifest() for the BENCHMARK policy")
    default_run_policy: dict[str, Any] = Field(description="ResolvedRunPolicy.as_manifest() for the DEFAULT policy")
    gc_widening: GCWideningBlock

    model_config = ConfigDict(frozen=True, extra="forbid")


# --------------------------------------------------------------------------------------------------
# Column order and per-column codecs. Order is declaration order, read straight off the model so it
# cannot drift from the class body above; a codec is declared only for a column whose CSV cell is
# not simply its own string value.
# --------------------------------------------------------------------------------------------------

OBSERVATION_COLUMNS: tuple[str, ...] = tuple(BenchmarkObservation.model_fields)
ACCOUNTING_COLUMNS: tuple[str, ...] = tuple(BenchmarkAccountingRow.model_fields)


def _enc_identity(value: Any) -> str:
    return str(value)


def _dec_identity(cell: str) -> str:
    return cell


def _enc_optional(value: Any) -> str:
    return "" if value is None else str(value)


def _dec_optional(cell: str) -> str | None:
    return None if cell == "" else cell


def _enc_overhang(value: str | None) -> str:
    return _OVERHANG_UNSTATED if value is None else value


def _dec_overhang(cell: str) -> str | None:
    return None if cell == _OVERHANG_UNSTATED else cell


def _enc_int(value: int) -> str:
    return str(value)


def _dec_int(cell: str) -> int:
    return int(cell)


def _enc_optional_int(value: int | None) -> str:
    return "" if value is None else str(value)


def _dec_optional_int(cell: str) -> int | None:
    return None if cell == "" else int(cell)


def _enc_optional_float(value: float | None) -> str:
    return "" if value is None else repr(value)


def _dec_optional_float(cell: str) -> float | None:
    return None if cell == "" else float(cell)


def _enc_bool(value: bool) -> str:
    return "True" if value else "False"


def _dec_bool(cell: str) -> bool:
    if cell == "True":
        return True
    if cell == "False":
        return False
    raise ValueError(f"expected 'True' or 'False', got {cell!r}")


def _enc_evaluation(value: FilterEvaluation) -> str:
    return value.value


def _dec_evaluation(cell: str) -> FilterEvaluation:
    return FilterEvaluation(cell)


_Codec = tuple[Any, Any]
_DEFAULT_CODEC: _Codec = (_enc_identity, _dec_identity)

_OBSERVATION_CODECS: dict[str, _Codec] = {
    "source_row_index": (_enc_int, _dec_int),
    "measured_value": (_enc_optional_float, _dec_optional_float),
    "guide_length": (_enc_int, _dec_int),
    "paired_slice_start_1based": (_enc_int, _dec_int),
    "paired_length": (_enc_int, _dec_int),
    "passenger_sequence": (_enc_optional, _dec_optional),
    "guide_3p_overhang": (_enc_overhang, _dec_overhang),
    "passenger_3p_overhang": (_enc_overhang, _dec_overhang),
    "source_redistribution": (_enc_optional, _dec_optional),
    "target_transcript_id": (_enc_optional, _dec_optional),
    "target_start_1based": (_enc_optional_int, _dec_optional_int),
    "target_end_1based": (_enc_optional_int, _dec_optional_int),
    "target_strand": (_enc_optional, _dec_optional),
    "split": (_enc_optional, _dec_optional),
}

_ACCOUNTING_CODECS: dict[str, _Codec] = {
    "paired_length": (_enc_int, _dec_int),
    "candidate_id": (_enc_optional, _dec_optional),
    "designed_guide_sequence": (_enc_optional, _dec_optional),
    "entered_design": (_enc_bool, _dec_bool),
    "default_filter_status": (_enc_evaluation, _dec_evaluation),
    "benchmark_filter_status": (_enc_evaluation, _dec_evaluation),
}


def _write_rows(
    rows: Sequence[BaseModel],
    path: Path | str,
    *,
    columns: tuple[str, ...],
    codecs: dict[str, _Codec],
    sort_key: Callable[[Any], Any],
) -> Path:
    """Write ``rows`` as a CSV with a declared, fixed column set, sorted for byte-identical reruns."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(rows, key=sort_key)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row in ordered:
            encoded = []
            for name in columns:
                encode, _ = codecs.get(name, _DEFAULT_CODEC)
                encoded.append(encode(getattr(row, name)))
            writer.writerow(encoded)
    return out_path


def _read_rows(
    path: Path | str,
    *,
    model: type[BaseModel],
    columns: tuple[str, ...],
    codecs: dict[str, _Codec],
) -> list[Any]:
    """Read a CSV written by :func:`_write_rows` back into validated model instances.

    Raises :class:`BenchmarkArtifactError` naming the missing column when the header omits one this
    schema declares, rather than constructing a row with a silently defaulted value for it (#109).
    """
    in_path = Path(path)
    if not in_path.is_file():
        raise BenchmarkArtifactError(f"{in_path} does not exist")
    with in_path.open("r", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise BenchmarkArtifactError(f"{in_path} is empty; expected a header row with {list(columns)}") from exc
        missing = [name for name in columns if name not in header]
        if missing:
            raise BenchmarkArtifactError(f"{in_path} header is missing declared column(s) {missing}; found {header}")
        extra = [name for name in header if name not in columns]
        if extra:
            raise BenchmarkArtifactError(
                f"{in_path} header declares column(s) {extra} this schema does not have; expected {list(columns)}"
            )
        index_of = {name: header.index(name) for name in columns}

        out: list[Any] = []
        for line_number, raw_row in enumerate(reader, start=2):
            if len(raw_row) != len(header):
                raise BenchmarkArtifactError(
                    f"{in_path} line {line_number} has {len(raw_row)} field(s), header declares {len(header)}"
                )
            kwargs: dict[str, Any] = {}
            for name in columns:
                _, decode = codecs.get(name, _DEFAULT_CODEC)
                cell = raw_row[index_of[name]]
                try:
                    kwargs[name] = decode(cell)
                except ValueError as exc:
                    raise BenchmarkArtifactError(
                        f"{in_path} line {line_number} column {name!r}={cell!r} is unreadable: {exc}"
                    ) from exc
            try:
                out.append(model(**kwargs))
            except ValidationError as exc:
                raise BenchmarkArtifactError(f"{in_path} line {line_number} failed schema validation: {exc}") from exc
    return out


def write_observations(rows: Sequence[BenchmarkObservation], path: Path | str) -> Path:
    """Write ``BenchmarkObservation`` rows to ``observations.csv``, sorted by ``observation_id``."""
    return _write_rows(
        rows, path, columns=OBSERVATION_COLUMNS, codecs=_OBSERVATION_CODECS, sort_key=lambda row: row.observation_id
    )


def read_observations(path: Path | str) -> list[BenchmarkObservation]:
    """Read ``observations.csv`` back into ``BenchmarkObservation`` rows."""
    return _read_rows(path, model=BenchmarkObservation, columns=OBSERVATION_COLUMNS, codecs=_OBSERVATION_CODECS)


def write_accounting(rows: Sequence[BenchmarkAccountingRow], path: Path | str) -> Path:
    """Write ``BenchmarkAccountingRow`` rows to ``accounting.csv``, sorted by (observation_id, candidate_id)."""
    return _write_rows(
        rows,
        path,
        columns=ACCOUNTING_COLUMNS,
        codecs=_ACCOUNTING_CODECS,
        sort_key=lambda row: (row.observation_id, row.candidate_id or ""),
    )


def read_accounting(path: Path | str) -> list[BenchmarkAccountingRow]:
    """Read ``accounting.csv`` back into ``BenchmarkAccountingRow`` rows."""
    return _read_rows(path, model=BenchmarkAccountingRow, columns=ACCOUNTING_COLUMNS, codecs=_ACCOUNTING_CODECS)


def write_manifest(manifest: BenchmarkArtifactManifest, path: Path | str) -> Path:
    """Write ``manifest.json``. The only file in the artifact allowed to carry a timestamp."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest.model_dump(mode="json"), indent=2, sort_keys=True) + "\n")
    return out_path


def read_manifest(path: Path | str) -> BenchmarkArtifactManifest:
    """Read ``manifest.json`` back into a validated ``BenchmarkArtifactManifest``."""
    in_path = Path(path)
    if not in_path.is_file():
        raise BenchmarkArtifactError(f"{in_path} does not exist")
    try:
        payload = json.loads(in_path.read_text())
    except json.JSONDecodeError as exc:
        raise BenchmarkArtifactError(f"{in_path} is not valid JSON: {exc}") from exc
    try:
        return BenchmarkArtifactManifest.model_validate(payload)
    except ValidationError as exc:
        raise BenchmarkArtifactError(f"{in_path} failed manifest schema validation: {exc}") from exc


__all__ = [
    "ACCOUNTING_COLUMNS",
    "ACCOUNTING_FILENAME",
    "ArchitectureId",
    "BENCHMARK_ARTIFACT_SCHEMA_VERSION",
    "BENCHMARK_MANIFEST_SCHEMA_VERSION",
    "BenchmarkAccountingRow",
    "BenchmarkArtifactCounts",
    "BenchmarkArtifactError",
    "BenchmarkArtifactManifest",
    "BenchmarkArtifactOutputs",
    "BenchmarkObservation",
    "CANDIDATES_ALL_FILENAME",
    "CompatibilityStatus",
    "DESIGN_INPUTS_FASTA_FILENAME",
    "DesignContextSource",
    "DuplexPairingStatus",
    "FilterExclusionCounts",
    "GCWideningBlock",
    "GCWideningEntry",
    "GuideMatch",
    "MANIFEST_FILENAME",
    "ManifestInputEntry",
    "ManifestOutputEntry",
    "ManifestPanelBlock",
    "OBSERVATIONS_FILENAME",
    "OBSERVATION_COLUMNS",
    "PAIRED_LENGTH_BOUNDS",
    "PolynucleotideRunRequirementBlock",
    "SplitLabel",
    "TargetIdentityStatus",
    "artifact_dir_name",
    "read_accounting",
    "read_manifest",
    "read_observations",
    "write_accounting",
    "write_manifest",
    "write_observations",
]
