"""``benchmark prepare``: ingest one panel at one paired length into a checksummed artifact (#109, bm-prepare slice).

``prepare_artifact`` reads a panel's raw table through the registry descriptor
(:mod:`sirnaforge.benchmark.panels`), lets :func:`sirnaforge.benchmark.panels.derive_observation`
decide every architecture-derived field, adds the provenance fields that are this module's own job
(observation id, source file/checksum, citation, design context) and writes the three files
:mod:`sirnaforge.benchmark.artifact` defines: ``observations.csv``, ``design_inputs.fasta`` and
``manifest.json``, into ``<out_dir>/<panel_id>__len<paired_length>/``. It does not decide a slice or
a compatibility verdict itself -- that rule lives in exactly one place, ``panels.py`` -- and it does
not run a design or write ``candidates_all.csv``/``accounting.csv``; those belong to
:func:`sirnaforge.benchmark.design.design_artifact`, a sibling slice's module, which extends this
same ``manifest.json`` in place.

Two facts this module works around, both verified against this repository rather than assumed (full
accounting in ``tests/unit/data/benchmark/README.md`` and ``panels.py``'s own module docstring):

* The issue's cited source, ``docs/prd_benchmark_artifacts_and_variable_length.md``, does not exist
  in the working tree or in git history. The issue body (#109) is the entire specification.
* Five registered panels are vendored and four are not, and this module reads which is which off
  ``PanelDescriptor.vendored_csv`` rather than a second map of its own. A vendored panel refuses an
  explicit ``--panel-csv``; an unvendored one (``huesken_full`` and the three architecture-level
  ``user_supplied_*`` ids) requires one.
* One vendored table is shared: ``tests/data/benchmarks/oligogym/records.csv`` holds all four
  OligoGym datasets, so this module honours ``PanelDescriptor.selects_row`` and
  ``BenchmarkObservation.source_row_index`` counts within the *selected* rows. Ingesting the whole
  file under one panel's declared architecture was measured, and it stamps paired-core-19 on all 356
  asymmetric rows -- the relabelling of a measured sequence #109 exists to prevent.

Because the manifest's design-stage fields (``counts.entered_design``/``no_candidate``/
``default_pass``/``benchmark_pass``, ``polynucleotide_run_requirement``, ``run_policy``,
``default_run_policy``, ``gc_widening``) are declared non-optional on
:class:`~sirnaforge.benchmark.artifact.BenchmarkArtifactManifest`, ``prepare`` writes truthful
zero/unwidened/not-yet-evaluated placeholders for all of them rather than leaving a hole
``design_artifact`` must special-case: a design-only manifest and a not-yet-designed manifest have
the same shape, so ``design_artifact``'s ``model_copy(update=...)`` only ever has to *replace*
fields, never *add* them.
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal, cast

from Bio.Seq import Seq

from sirnaforge import __version__
from sirnaforge.benchmark.artifact import (
    BENCHMARK_ARTIFACT_SCHEMA_VERSION,
    BENCHMARK_MANIFEST_SCHEMA_VERSION,
    DESIGN_INPUTS_FASTA_FILENAME,
    MANIFEST_FILENAME,
    NATIVE_MAPPINGS_FILENAME,
    OBSERVATIONS_FILENAME,
    BenchmarkArtifactCounts,
    BenchmarkArtifactManifest,
    BenchmarkArtifactOutputs,
    BenchmarkNativeMapping,
    BenchmarkObservation,
    FilterExclusionCounts,
    GCWideningBlock,
    GCWideningEntry,
    ManifestInputEntry,
    ManifestOutputEntry,
    ManifestPanelBlock,
    NativeReferenceBlock,
    PolynucleotideRunRequirementBlock,
    artifact_dir_name,
    write_manifest,
    write_native_mappings,
    write_observations,
)
from sirnaforge.benchmark.panels import PanelDescriptor, derive_observation, describe_panel
from sirnaforge.config.run_policy import default_for
from sirnaforge.models.policy import FilterAction, FilterComparator, SettingSource
from sirnaforge.utils.hashing import file_sha256

#: ``prepare.py`` -> ``benchmark`` -> ``sirnaforge`` -> ``src`` -> repo root.
_REPO_ROOT = Path(__file__).resolve().parents[3]


def _checkout_root() -> Path | None:
    """The source checkout the vendored panel tables live in, or ``None`` when there is not one.

    Every vendored panel names a path under ``tests/`` (``panels.py``), and ``pyproject.toml`` ships
    only ``src/sirnaforge`` in the wheel and excludes ``/tests`` from the sdist. So an installed
    sirnaforge -- a wheel, or the release image -- has no vendored bytes at all, and ``_REPO_ROOT``
    resolves to whatever sits four levels above ``site-packages/sirnaforge/benchmark``.

    Left as a checkout-only developer tool deliberately: the tables are third-party redistributions
    whose provenance is documented in ``tests/data/benchmarks/README.md``, and moving 1.8 MB of them
    into the wheel is a redistribution decision, not a packaging tidy-up. What is fixed here is the
    message. ``sirnaforge benchmark prepare --panel huesken_subset`` used to exit 1 with
    ``panel csv not found: /opt/conda/lib/python3.12/tests/unit/data/sirna_efficacy_subset.csv`` --
    a path that reads as a corrupt build rather than as "this panel needs a checkout".
    """
    if (_REPO_ROOT / "pyproject.toml").is_file() and (_REPO_ROOT / "tests").is_dir():
        return _REPO_ROOT
    return None


DEFAULT_ARTIFACT_ROOT = Path("benchmark_artifacts")

#: The one gate #109 requires to stay active and auditable in every artifact, even one no design has
#: touched yet (module docstring). Threshold is read off the shipped default so this placeholder
#: cannot silently drift from the real gate ``design_artifact`` later evaluates against.
_POLYNUCLEOTIDE_FILTER_ID = "max_poly_runs"


class BenchmarkPrepareError(ValueError):
    """A panel cannot be prepared as asked, raised before any file is written.

    A ``ValueError`` so a CLI boundary can catch it alongside argument-parsing errors and
    ``sirnaforge.benchmark.artifact.BenchmarkArtifactError``, matching
    ``config.run_policy.RunPolicyError``'s convention for the same reason.
    """


def _resolve_source_csv(descriptor: PanelDescriptor, panel_csv: Path | None) -> tuple[Path, str]:
    """The table to read, and the path string to record as ``source_file``.

    Refuses exactly two cases, each with the reason in the message: a vendored panel given an
    explicit ``--panel-csv`` anyway (a run must not be able to silently read different bytes than
    the ones its own manifest names as vendored), and an unvendored panel given no override at all.
    """
    if descriptor.data_present:
        if panel_csv is not None:
            raise BenchmarkPrepareError(
                f"panel {descriptor.panel_id!r} ships vendored bytes in this repository; passing "
                "--panel-csv for it would let a run silently read different bytes than the ones its "
                "own manifest names as vendored. Omit --panel-csv to use the vendored table"
            )
        # The descriptor names the path, and its own validator already enforces
        # ``data_present <=> vendored_csv is not None``. This module used to keep a second map keyed by
        # panel_id, on the reasoning that declaring a descriptor and vendoring bytes are different
        # acts. They are -- but two lists of the same fact drift, and this pair did: when the OligoGym
        # panels were re-derived from the bytes f4beab7 vendored, the descriptors gained
        # ``data_present=True`` while the map here still held one entry, so every one of those panels
        # refused with "that is a bug in prepare.py". One owner now.
        relative = descriptor.vendored_csv
        if relative is None:  # pragma: no cover - the descriptor validator forbids this pairing
            raise BenchmarkPrepareError(
                f"panel {descriptor.panel_id!r} is registered data_present=True but names no "
                "vendored_csv; that is a bug in the panel registry, not in your invocation"
            )
        checkout = _checkout_root()
        if checkout is None:
            raise BenchmarkPrepareError(
                f"panel {descriptor.panel_id!r} reads bytes vendored at {relative}, which ship with the "
                "source repository and not with the installed package (pyproject.toml puts only "
                "src/sirnaforge in the wheel and excludes /tests from the sdist). Run this command "
                "from a git checkout of sirnaforge, or use a panel that takes --panel-csv"
            )
        return checkout / relative, relative
    if panel_csv is None:
        # The required columns are named here rather than left to the header check further in: the
        # ``user_supplied_*`` ids read a table this repository has never seen, so "which columns" is
        # the caller's next question and the descriptor is the one thing that can answer it.
        required = ", ".join(descriptor.required_source_columns())
        raise BenchmarkPrepareError(
            f"panel {descriptor.panel_id!r} has no vendored bytes in this repository and no "
            f"--panel-csv was given; supply --panel-csv pointing at a table carrying {required}"
        )
    return panel_csv, str(panel_csv)


def _resolve_paired_length(descriptor: PanelDescriptor, paired_length: int | None) -> int:
    """The paired length this artifact will be tagged with.

    Defaults to the descriptor's declared length. A panel that declares none (only
    :attr:`~sirnaforge.benchmark.panels.PanelArchitecture.ASYMMETRIC` does) has no length #109's
    fixed-length interface could default to, and every asymmetric observation is incompatible
    regardless of which length is picked (``panels.py::_check_compatibility``) -- so rather than
    inventing a number nothing measured, an asymmetric panel with no explicit ``paired_length``
    is refused.
    """
    if paired_length is not None:
        return paired_length
    if descriptor.declared_paired_length is not None:
        return descriptor.declared_paired_length
    raise BenchmarkPrepareError(
        f"panel {descriptor.panel_id!r} is {descriptor.architecture.value} and declares no fixed "
        "paired length; pass paired_length explicitly (19-23) to tag the artifact directory"
    )


def _read_raw_rows(path: Path) -> tuple[tuple[str, ...], tuple[dict[str, str], ...]]:
    """The source table's header and rows. The header is returned so the columns can be checked once."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = tuple(reader.fieldnames or ())
        return fieldnames, tuple(dict(row) for row in reader)


def _load_native_assets(  # noqa: PLR0912 - validates a closed multi-file native import contract
    *,
    observations_csv: Path | str | None,
    mappings_csv: Path | str | None,
    reference_manifest: Path | str | None,
    reference_fasta: Path | str | None,
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, list[BenchmarkNativeMapping]],
    NativeReferenceBlock | None,
    list[tuple[str, Path]],
]:
    supplied = (observations_csv, mappings_csv, reference_manifest)
    if not any(value is not None for value in supplied):
        if reference_fasta is not None:
            raise BenchmarkPrepareError("native_reference_fasta requires the three native import inputs")
        return {}, {}, None, []
    if not all(value is not None for value in supplied):
        raise BenchmarkPrepareError(
            "native import requires native_observations_csv, native_mappings_csv, and native_reference_manifest together"
        )
    assert observations_csv is not None and mappings_csv is not None and reference_manifest is not None
    paths = [
        ("native_observations_csv", Path(observations_csv)),
        ("native_mappings_csv", Path(mappings_csv)),
        ("native_reference_manifest", Path(reference_manifest)),
    ]
    for _, path in paths:
        if not path.is_file():
            raise BenchmarkPrepareError(f"native import input not found: {path}")
    try:
        reference = json.loads(Path(reference_manifest).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise BenchmarkPrepareError(f"invalid native reference manifest: {exc}") from exc
    if not all(reference.get(name) for name in ("species", "assembly", "release", "sha256")):
        raise BenchmarkPrepareError("native reference manifest must declare species, assembly, release, and sha256")
    digest = str(reference["sha256"]).lower().removeprefix("sha256:")
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise BenchmarkPrepareError("native reference manifest sha256 must be 64 hexadecimal characters")
    reference["sha256"] = digest
    if str(reference["species"]).lower() not in {"human", "mouse", "rat", "macaque"}:
        raise BenchmarkPrepareError("native reference species must be human, mouse, rat, or macaque")
    species = cast(Literal["human", "mouse", "rat", "macaque"], str(reference["species"]).lower())
    reference_block = NativeReferenceBlock(
        species=species,
        assembly=str(reference["assembly"]),
        release=str(reference["release"]),
        sha256=digest,
    )
    reference = {
        "species": reference_block.species,
        "assembly": reference_block.assembly,
        "release": reference_block.release,
        "sha256": reference_block.sha256,
    }
    if reference_fasta is not None:
        fasta = Path(reference_fasta)
        if not fasta.is_file() or file_sha256(fasta) != digest:
            raise BenchmarkPrepareError("native reference FASTA does not match the reference manifest sha256")
        paths.append(("native_reference_fasta", fasta))
    native_rows: dict[str, dict[str, str]] = {}
    with Path(observations_csv).open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "benchmark_id",
            "dataset",
            "source_row",
            "full_guide_sequence",
            "paired_guide_sequence",
            "passenger_sequence",
            "paired_length",
            "duplex_pairing_status",
            "guide_overhang_sequence",
            "passenger_overhang_sequence",
            "native_reference_species",
            "native_reference_assembly",
            "native_mapping_count",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise BenchmarkPrepareError(f"native observations csv is missing columns {sorted(missing)}")
        for row in reader:
            key = row["benchmark_id"]
            if not key or key in native_rows:
                raise BenchmarkPrepareError(f"duplicate or empty native observation identity {key!r}")
            if (
                row["native_reference_species"].lower() != reference["species"]
                or row["native_reference_assembly"] != reference["assembly"]
            ):
                raise BenchmarkPrepareError(f"native observation {key!r} reference provenance conflicts with manifest")
            native_rows[key] = row
    native_hits: dict[str, list[BenchmarkNativeMapping]] = {}
    mapping_ids: set[str] = set()
    with Path(mappings_csv).open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "mapping_id",
            "benchmark_id",
            "mapping_status",
            "mapping_transcript_id",
            "mapping_gene_id",
            "mapping_assembly",
            "site_start_1based",
            "site_end_1based",
            "site_sequence",
            "mapping_strand",
            "query_sequence",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise BenchmarkPrepareError(f"native mappings csv is missing columns {sorted(missing)}")
        for row in reader:
            if row["mapping_status"] != "mapped":
                continue
            if not row["mapping_id"] or row["mapping_id"] in mapping_ids:
                raise BenchmarkPrepareError(f"duplicate or empty native mapping identity {row['mapping_id']!r}")
            mapping_ids.add(row["mapping_id"])
            key = row["benchmark_id"]
            if key not in native_rows:
                raise BenchmarkPrepareError(f"native mapping references unknown observation {key!r}")
            strand = row["mapping_strand"]
            if strand not in {"+", "-"}:
                raise BenchmarkPrepareError(f"invalid native mapping strand {strand!r}")
            typed_strand = cast(Literal["+", "-"], strand)
            try:
                hit = BenchmarkNativeMapping(
                    observation_id=key,
                    transcript_id=row["mapping_transcript_id"],
                    gene_id=row["mapping_gene_id"],
                    gene_symbol=row.get("mapping_gene_symbol") or None,
                    start_1based=int(row["site_start_1based"]),
                    end_1based=int(row["site_end_1based"]),
                    strand=typed_strand,
                    site_sequence=row["site_sequence"].upper(),
                    assembly=row["mapping_assembly"],
                    source_mapping_id=row["mapping_id"],
                )
            except (ValueError, TypeError) as exc:
                raise BenchmarkPrepareError(f"invalid native mapping {row.get('mapping_id')!r}: {exc}") from exc
            if hit.assembly != reference["assembly"] or hit.site_sequence != row["query_sequence"].upper():
                raise BenchmarkPrepareError(
                    f"native mapping {hit.source_mapping_id} conflicts with its query/reference"
                )
            native_hits.setdefault(key, []).append(hit)
    for key, native_row in native_rows.items():
        try:
            declared_count = int(native_row["native_mapping_count"])
        except ValueError as exc:
            raise BenchmarkPrepareError(f"native observation {key!r} has an invalid mapping count") from exc
        if declared_count != len(native_hits.get(key, ())):
            raise BenchmarkPrepareError(f"native observation {key!r} mapping count disagrees with all-hit sidecar")
    if reference_fasta is not None:
        native_transcripts = _load_transcripts(Path(reference_fasta))
        for hits in native_hits.values():
            for hit in hits:
                sequence = native_transcripts.get(hit.transcript_id)
                if sequence is None:
                    raise BenchmarkPrepareError(
                        f"native reference FASTA does not contain mapped transcript {hit.transcript_id!r}"
                    )
                span = sequence[hit.start_1based - 1 : hit.end_1based]
                if span != hit.site_sequence and str(Seq(span).reverse_complement()) != hit.site_sequence:
                    raise BenchmarkPrepareError(
                        f"native mapping {hit.source_mapping_id} sequence does not match its reference FASTA coordinates"
                    )
    return native_rows, native_hits, reference_block, paths


def _apply_native_geometry(  # noqa: PLR0912 - validates measured geometry before artifact construction
    observation: BenchmarkObservation,
    native_row: Mapping[str, str],
    requested_length: int,
    *,
    native_key: str,
) -> BenchmarkObservation:
    """Replace panel-wide assumptions with this source row's measured duplex geometry.

    The native OligoGym extract records Martinelli's mixed 19-nt core, fully paired 21-nt, and
    unsupported rows explicitly. Raw strands must still match the panel bytes exactly before its
    geometry is accepted.
    """
    native_guide = native_row.get("full_guide_sequence", "").upper().replace("U", "T")
    native_passenger = native_row.get("passenger_sequence", "").upper().replace("U", "T")
    if native_guide != observation.full_guide_sequence.replace("U", "T"):
        raise BenchmarkPrepareError(f"native source identity {native_key!r} has a different measured guide")
    if native_passenger != (observation.passenger_sequence or "").replace("U", "T"):
        raise BenchmarkPrepareError(f"native source identity {native_key!r} has a different measured passenger")
    try:
        measured_length = int(native_row["paired_length"])
        paired_guide = native_row["paired_guide_sequence"].upper()
    except (KeyError, ValueError) as exc:
        raise BenchmarkPrepareError(f"native source identity {native_key!r} has invalid measured geometry") from exc
    pairing = native_row.get("duplex_pairing_status", "")
    if pairing not in {"fully_complementary", "paired_core_with_overhang"} or measured_length < 19:
        return observation.model_copy(
            update={
                "architecture": "unsupported",
                "paired_guide_sequence": paired_guide if paired_guide else observation.full_guide_sequence,
                "paired_length": measured_length,
                "guide_3p_overhang": native_row.get("guide_3p_overhang", native_row.get("guide_overhang_sequence"))
                or None,
                "passenger_3p_overhang": native_row.get(
                    "passenger_3p_overhang", native_row.get("passenger_overhang_sequence")
                )
                or None,
                "duplex_pairing_status": "non_complementary"
                if pairing == "non_complementary"
                else "partial_unresolved",
                "compatibility_status": "incompatible",
                "compatibility_reason": "review_only_unsupported_or_unresolved_measured_duplex_architecture",
                "native_mapping_review_reason": "review_only_unsupported_or_unresolved_measured_duplex_architecture",
            }
        )
    if len(paired_guide) != measured_length or paired_guide != observation.full_guide_sequence[:measured_length]:
        raise BenchmarkPrepareError(
            f"native source identity {native_key!r} does not preserve the measured paired guide"
        )
    if pairing == "fully_complementary":
        if measured_length == observation.guide_length and len(native_passenger) == measured_length:
            if native_passenger != str(Seq(native_guide).reverse_complement()):
                raise BenchmarkPrepareError(
                    f"native source identity {native_key!r} falsely declares a fully paired duplex"
                )
            architecture = "fully_complementary"
            guide_tail = passenger_tail = ""
        elif measured_length < observation.guide_length and len(native_passenger) == measured_length:
            if str(Seq(paired_guide).reverse_complement()).replace("U", "T") != native_passenger.replace("U", "T"):
                raise BenchmarkPrepareError(
                    f"native source identity {native_key!r} falsely declares a fully paired core"
                )
            architecture = "paired_core_with_overhang"
            guide_tail = observation.full_guide_sequence[measured_length:]
            passenger_tail = ""
            pairing = "paired_core_with_overhang"
        else:
            raise BenchmarkPrepareError(f"native source identity {native_key!r} falsely declares a fully paired duplex")
    else:
        if measured_length >= observation.guide_length:
            raise BenchmarkPrepareError(f"native source identity {native_key!r} has an invalid paired-core length")
        if str(Seq(paired_guide).reverse_complement()).replace("U", "T") != native_passenger[:measured_length].replace(
            "U", "T"
        ):
            raise BenchmarkPrepareError(f"native source identity {native_key!r} paired core does not match passenger")
        architecture = "paired_core_with_overhang"
        guide_tail = observation.full_guide_sequence[measured_length:]
        passenger_tail = (observation.passenger_sequence or "")[measured_length:]
        measured_guide_tail = native_row.get("guide_3p_overhang", native_row.get("guide_overhang_sequence", guide_tail))
        measured_passenger_tail = native_row.get(
            "passenger_3p_overhang", native_row.get("passenger_overhang_sequence", passenger_tail)
        )
        if measured_guide_tail != guide_tail or measured_passenger_tail != passenger_tail:
            raise BenchmarkPrepareError(
                f"native source identity {native_key!r} overhangs disagree with measured strands"
            )
    compatible = measured_length == requested_length
    return observation.model_copy(
        update={
            "architecture": architecture,
            "paired_guide_sequence": paired_guide,
            "paired_length": measured_length,
            "guide_3p_overhang": guide_tail,
            "passenger_3p_overhang": passenger_tail,
            "duplex_pairing_status": pairing,
            "compatibility_status": "compatible" if compatible else "incompatible",
            "compatibility_reason": ""
            if compatible
            else f"measured paired length {measured_length} differs from requested artifact length {requested_length}",
        }
    )


def _check_source_columns(descriptor: PanelDescriptor, fieldnames: Sequence[str], source_file: str) -> None:
    """Refuse a table missing a column the descriptor cannot do without, naming it.

    Checked against the header before any row is read, because the alternative is a ``KeyError``
    from the middle of a file -- and for the ``user_supplied_*`` ids, whose whole input is a table this
    repository has never seen, a mistyped column name is the likeliest way a run goes wrong. The
    descriptor decides which columns are fatal (``PanelDescriptor.required_source_columns``); an
    optional column that is simply absent stays absent, and its field is ``None``.
    """
    missing = [column for column in descriptor.required_source_columns() if column not in fieldnames]
    if missing:
        raise BenchmarkPrepareError(
            f"{source_file} is missing the column(s) {', '.join(missing)} that panel "
            f"{descriptor.panel_id!r} declares. The table's header is "
            f"{', '.join(fieldnames) if fieldnames else '<empty>'}; rename your columns to the ones "
            "the descriptor names (docs/benchmark_artifacts.md lists them per panel)"
        )


def _select_rows(
    descriptor: PanelDescriptor, rows: Sequence[Mapping[str, str]], source_file: str
) -> tuple[Mapping[str, str], ...]:
    """The rows of ``source_file`` that are this panel's, in source order.

    A vendored table may be shared -- ``oligogym/records.csv`` holds four datasets keyed by ``dataset``
    -- so a reader that took every row would derive a sibling panel's measured duplexes under this
    panel's declared architecture. An empty selection is refused rather than written as a zero-row
    artifact, which would look like a panel with nothing in it instead of a table that does not hold
    this panel's rows.
    """
    if descriptor.row_selector is None:
        return tuple(rows)
    selected = tuple(row for row in rows if descriptor.selects_row(row))
    if not selected:
        selector = descriptor.row_selector
        raise BenchmarkPrepareError(
            f"{source_file} has no row whose {selector.column!r} is one of "
            f"{', '.join(selector.values)}, which is how panel {descriptor.panel_id!r} names its own "
            "rows in a shared table; there is nothing for it to ingest"
        )
    return selected


def _load_transcripts(path: Path) -> dict[str, str]:
    """Header id -> upper-cased sequence, for ``--panel-transcripts``."""
    from Bio import SeqIO  # noqa: PLC0415 - only needed on this optional path

    return {record.id: str(record.seq).upper() for record in SeqIO.parse(str(path), "fasta")}


def _design_context(accession: str | None, full_guide_sequence: str, transcripts: Mapping[str, str]) -> tuple[str, str]:
    """The sequence to hand ``sirnaforge design``, and which rule produced it.

    Uses the accession's own transcript when ``--panel-transcripts`` supplied one; falls back to
    the measured target site -- the reverse complement of the guide as measured, never trimmed to
    the requested paired length -- when no transcript is available. Independent of compatibility: a
    target site is knowable for every observation, compatible or not, from the guide alone, which is
    why ``design_context_source`` is a required field on an incompatible row too.
    """
    if accession and accession in transcripts:
        return transcripts[accession], "panel_transcript"
    return str(Seq(full_guide_sequence).reverse_complement()), "measured_target_site"


def _build_observation(
    *,
    descriptor: PanelDescriptor,
    paired_length: int,
    source_row_index: int,
    row: Mapping[str, str],
    source_file: str,
    source_sha256: str,
    transcripts: Mapping[str, str],
) -> tuple[BenchmarkObservation, str]:
    """One provenance-complete observation, plus its design-context sequence.

    Target-locality fields (``target_transcript_id``, ``target_identity_status``,
    ``target_start_1based``/``target_end_1based``/``target_strand``) are always ``None``/
    ``"unavailable"``, which is lossy but never an overclaim. For the three ``user_supplied_*`` ids and
    both Huesken panels it is simply the truth: no coordinate column is mapped, so there is none to
    report. For the OligoGym-derived panels ``derive_observation`` does compute a
    ``synthetic_context_local`` span (1-based 71..) and this writer drops it, because
    ``BenchmarkArtifactCounts`` has no bucket that would keep the tally honest for such a row; carrying
    it through is outstanding work, tracked with the rest of the target-identity story on #110. What is
    structurally impossible either way is a native claim: the vocabulary has no ``confirmed`` member.
    """
    try:
        derived = derive_observation(descriptor, row, requested_paired_length=paired_length)
    except KeyError as exc:
        raise BenchmarkPrepareError(
            f"{source_file} row {source_row_index}: missing declared column {exc}; the panel "
            f"descriptor for {descriptor.panel_id!r} expects it"
        ) from exc

    accession = row.get(descriptor.columns.accession_column) if descriptor.columns.accession_column else None
    design_context_seq, design_context_source = _design_context(accession, derived.full_guide_sequence, transcripts)

    observation_id = f"{descriptor.panel_id}:{source_row_index:06d}"
    observation = BenchmarkObservation(
        observation_id=observation_id,
        panel_id=descriptor.panel_id,
        source_row_index=source_row_index,
        architecture=derived.architecture.value,
        assay_label=derived.assay_label,
        measured_endpoint=derived.measured_endpoint,
        measured_value=derived.measured_value,
        full_guide_sequence=derived.full_guide_sequence,
        guide_length=derived.guide_length,
        paired_guide_sequence=derived.paired_guide_sequence,
        paired_slice_start_1based=derived.paired_slice_start_1based,
        paired_length=derived.paired_length,
        passenger_sequence=derived.passenger_sequence,
        guide_3p_overhang=derived.guide_3p_overhang,
        passenger_3p_overhang=derived.passenger_3p_overhang,
        duplex_pairing_status=derived.duplex_pairing_status.value,
        source_citation=descriptor.citation,
        source_redistribution=descriptor.redistribution,
        source_file=source_file,
        source_sha256=source_sha256,
        target_transcript_id=None,
        target_identity_status="unavailable",
        target_start_1based=None,
        target_end_1based=None,
        target_strand=None,
        design_context_id=observation_id,
        design_context_source=design_context_source,  # type: ignore[arg-type]
        split=derived.split,  # type: ignore[arg-type]
        compatibility_status=derived.compatibility_status.value,
        compatibility_reason=derived.compatibility_reason,
    )
    return observation, design_context_seq


def _write_design_inputs_fasta(
    observations: Sequence[BenchmarkObservation], contexts: Mapping[str, str], path: Path
) -> Path:
    """One FASTA record per *compatible* observation, keyed by ``design_context_id``.

    A plain FASTA: ``sirnaforge design <this file> --length <paired_length>`` and
    ``sirnaforge workflow --input-fasta <this file>`` both consume it with no new code, which is how
    #109's "consumable by the existing fixed-length path" acceptance criterion is met. Sorted by
    ``design_context_id`` and free of any timestamp, so two prepares of the same inputs are
    byte-identical, matching ``write_observations``'s own guarantee for ``observations.csv``.
    """
    ordered = sorted(
        (observation for observation in observations if observation.compatibility_status == "compatible"),
        key=lambda observation: observation.design_context_id,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="\n") as handle:
        for observation in ordered:
            handle.write(f">{observation.design_context_id}\n{contexts[observation.design_context_id]}\n")
    return path


def _count_fasta_records(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.startswith(">"))


def _count_csv_rows(path: Path) -> int:
    with path.open(newline="") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def prepare_artifact(  # noqa: PLR0912 - the optional native import is a separate, bounded input path
    *,
    panel_id: str,
    paired_length: int | None = None,
    panel_csv: Path | str | None = None,
    panel_transcripts: Path | str | None = None,
    native_observations_csv: Path | str | None = None,
    native_mappings_csv: Path | str | None = None,
    native_reference_manifest: Path | str | None = None,
    native_reference_fasta: Path | str | None = None,
    out_dir: Path | str = DEFAULT_ARTIFACT_ROOT,
    overwrite: bool = False,
    invoked_command: Sequence[str] | None = None,
) -> BenchmarkArtifactManifest:
    """Ingest one panel at one paired length into ``<out_dir>/<panel_id>__len<paired_length>/``.

    Args:
        panel_id: A :mod:`sirnaforge.benchmark.panels` registry id.
        paired_length: The fixed design length, 19-23. Defaults to the descriptor's declared
            length; refused if the panel declares none and no length is given (see
            :func:`_resolve_paired_length`).
        panel_csv: The panel's raw table. Required for a registered panel with no vendored bytes;
            refused for one that already ships its own (see :func:`_resolve_source_csv`).
        panel_transcripts: Optional FASTA of full-length transcripts, keyed by the panel's accession
            column. When an observation's accession is found here, its design context is the real
            transcript rather than the bare measured target site.
        native_observations_csv: Optional row-level native geometry and observation table.
        native_mappings_csv: Optional all-hit exact transcript mapping table.
        native_reference_manifest: Optional species/assembly/release/checksum declaration; required
            with both native CSV inputs.
        native_reference_fasta: Optional reference FASTA. Its checksum and every imported mapping
            span are verified when supplied.
        out_dir: Root directory artifacts are written under.
        overwrite: Replace an existing artifact directory rather than refusing it.
        invoked_command: The command that asked for this, verbatim. Defaults to ``sys.argv`` for a
            direct caller; the CLI slice passes its own argv so the manifest never reconstructs the
            command from options.

    Returns:
        The manifest that was also written to ``manifest.json``.

    Raises:
        BenchmarkPrepareError: An aggregate panel id, a panel-csv/vendoring conflict, a paired length
            with no default and none given, a required column the source table does not carry, a
            selector that matches no row, or an existing artifact directory without ``overwrite``.
            Every one of them is raised before the artifact directory is created.
        ValueError: An unregistered ``panel_id`` (:func:`sirnaforge.benchmark.panels.describe_panel`).
        sirnaforge.benchmark.artifact.BenchmarkArtifactError: ``paired_length`` outside 19-23.
    """
    descriptor = describe_panel(panel_id)
    # Refused here rather than by `derive_observation` on the first row, which raises after `mkdir` and
    # so leaves an empty artifact directory behind. Same reason, stated once in `panels.py`.
    if descriptor.aggregate_of is not None:
        raise BenchmarkPrepareError(
            f"panel {descriptor.panel_id!r} is a shared redistribution of "
            f"{', '.join(descriptor.aggregate_of)} spanning more than one duplex architecture, so it "
            "has none of its own to ingest under; prepare one of those panels instead"
        )
    length = _resolve_paired_length(descriptor, paired_length)
    # Validates the 19-23 bound and builds the directory name in the one place that spells the
    # convention, before anything is read from disk.
    dir_name = artifact_dir_name(descriptor.panel_id, length)

    panel_csv_path = Path(panel_csv) if panel_csv is not None else None
    source_path, source_file = _resolve_source_csv(descriptor, panel_csv_path)
    if not source_path.is_file():
        raise BenchmarkPrepareError(f"panel csv not found: {source_path}")
    source_sha256 = file_sha256(source_path)
    native_rows, native_hits, native_reference, native_paths = _load_native_assets(
        observations_csv=native_observations_csv,
        mappings_csv=native_mappings_csv,
        reference_manifest=native_reference_manifest,
        reference_fasta=native_reference_fasta,
    )

    # Read and validated before the artifact directory is created, so a mis-mapped or wrong-panel
    # table is refused without leaving a directory behind for the caller to clean up.
    fieldnames, all_rows = _read_raw_rows(source_path)
    _check_source_columns(descriptor, fieldnames, source_file)
    raw_rows = _select_rows(descriptor, all_rows, source_file)

    transcripts: dict[str, str] = {}
    transcripts_path: Path | None = None
    if panel_transcripts is not None:
        transcripts_path = Path(panel_transcripts)
        if not transcripts_path.is_file():
            raise BenchmarkPrepareError(f"panel transcripts fasta not found: {transcripts_path}")
        transcripts = _load_transcripts(transcripts_path)

    artifact_dir = Path(out_dir) / dir_name
    if artifact_dir.exists() and not overwrite:
        raise BenchmarkPrepareError(
            f"artifact directory {artifact_dir} already exists; pass overwrite=True (--overwrite) "
            "to replace it, or choose a different out_dir"
        )
    artifact_dir.mkdir(parents=True, exist_ok=True)

    observations: list[BenchmarkObservation] = []
    native_mapping_rows: list[BenchmarkNativeMapping] = []
    contexts: dict[str, str] = {}
    for index, row in enumerate(raw_rows):
        observation, context_seq = _build_observation(
            descriptor=descriptor,
            paired_length=length,
            source_row_index=index,
            row=row,
            source_file=source_file,
            source_sha256=source_sha256,
            transcripts=transcripts,
        )
        if native_reference is not None:
            native_key = row.get("benchmark_id", "")
            if not native_key:
                native_key = (
                    f"{row.get('dataset', descriptor.panel_id)}__row_{int(row.get('source_row', index + 1)):06d}"
                )
            native_row = native_rows.get(native_key)
            if native_row is None:
                raise BenchmarkPrepareError(f"native observations are missing source identity {native_key!r}")
            if native_row["dataset"] != row.get("dataset", descriptor.panel_id) or int(native_row["source_row"]) != int(
                row.get("source_row", index + 1)
            ):
                raise BenchmarkPrepareError(f"native source identity {native_key!r} conflicts with panel row order")
            observation = _apply_native_geometry(observation, native_row, length, native_key=native_key)
            hits = native_hits.get(native_key, [])
            expected_target = str(Seq(observation.paired_guide_sequence).reverse_complement()).replace("U", "T")
            geometry_matches = (
                native_row["paired_guide_sequence"].upper().replace("U", "T")
                == observation.paired_guide_sequence.replace("U", "T")
                and int(native_row["paired_length"]) == observation.paired_length
            )
            if geometry_matches:
                projected_hits = []
                for hit in hits:
                    hit_sequence = hit.site_sequence.replace("U", "T")
                    if hit_sequence == expected_target:
                        projected_hits.append(hit)
                    elif (
                        observation.passenger_sequence
                        and hit_sequence == observation.passenger_sequence.replace("U", "T")
                        and hit_sequence.startswith(expected_target)
                    ):
                        # The existing mapper searched the measured passenger strand in full.
                        # Its 5'-anchored prefix is the paired cDNA site; trim the measured 3' tail
                        # from the reported end coordinate while retaining the exact native match.
                        projected_hits.append(
                            hit.model_copy(
                                update={
                                    "end_1based": hit.start_1based + observation.paired_length - 1,
                                    "site_sequence": expected_target,
                                }
                            )
                        )
                    else:
                        raise BenchmarkPrepareError(
                            f"native mapping {hit.source_mapping_id} does not match the measured paired duplex geometry"
                        )
                hits = projected_hits
            elif hits:
                raise BenchmarkPrepareError(f"native source identity {native_key!r} has hits for unsupported geometry")
            if hits:
                representative = hits[0]
                observation = observation.model_copy(
                    update={
                        "target_transcript_id": representative.transcript_id,
                        "target_identity_status": "confirmed" if len(hits) == 1 else "ambiguous",
                        "target_start_1based": representative.start_1based,
                        "target_end_1based": representative.end_1based,
                        "target_strand": representative.strand,
                        "target_gene_id": representative.gene_id,
                        "native_reference_species": native_reference.species,
                        "native_reference_assembly": native_reference.assembly,
                        "native_reference_sha256": native_reference.sha256,
                        "native_mapping_count": len(hits),
                    }
                )
                native_mapping_rows.extend(
                    hit.model_copy(update={"observation_id": observation.observation_id}) for hit in hits
                )
            else:
                observation = observation.model_copy(
                    update={
                        "native_reference_species": native_reference.species,
                        "native_reference_assembly": native_reference.assembly,
                        "native_reference_sha256": native_reference.sha256,
                    }
                )
        observations.append(observation)
        contexts[observation.design_context_id] = context_seq

    observations_path = write_observations(observations, artifact_dir / OBSERVATIONS_FILENAME)
    fasta_path = _write_design_inputs_fasta(observations, contexts, artifact_dir / DESIGN_INPUTS_FASTA_FILENAME)
    native_mappings_path = (
        write_native_mappings(native_mapping_rows, artifact_dir / NATIVE_MAPPINGS_FILENAME)
        if native_reference is not None
        else None
    )

    kept = sum(1 for observation in observations if observation.compatibility_status == "compatible")
    incompatible = len(observations) - kept
    counts = BenchmarkArtifactCounts(
        observations_in_source=len(observations),
        observations_kept=kept,
        observations_incompatible=incompatible,
        mapped_native=sum(1 for observation in observations if observation.target_identity_status == "confirmed"),
        mapped_native_ambiguous=sum(
            1 for observation in observations if observation.target_identity_status == "ambiguous"
        ),
        mapped_panel_local=0,
        mapping_unavailable=sum(
            1 for observation in observations if observation.target_identity_status == "unavailable"
        ),
        entered_design=0,
        no_candidate=0,
        default_pass=0,
        benchmark_pass=0,
        excluded_by_filter={},
    )

    inputs = [
        ManifestInputEntry(
            role="panel_csv", path=source_file, sha256=source_sha256, size_bytes=source_path.stat().st_size
        )
    ]
    if transcripts_path is not None:
        inputs.append(
            ManifestInputEntry(
                role="panel_transcripts_fasta",
                path=str(panel_transcripts),
                sha256=file_sha256(transcripts_path),
                size_bytes=transcripts_path.stat().st_size,
            )
        )
    for role, path in native_paths:
        inputs.append(
            ManifestInputEntry(role=role, path=str(path), sha256=file_sha256(path), size_bytes=path.stat().st_size)
        )

    # Output paths are the fixed inner filenames, not `str(observations_path)`: the manifest must not
    # depend on which --out-dir this run happened to use, or two prepares of the same bytes into two
    # roots would not compare equal (ManifestOutputEntry.path).
    outputs = BenchmarkArtifactOutputs(
        observations_csv=ManifestOutputEntry(
            path=OBSERVATIONS_FILENAME,
            exists=True,
            size_bytes=observations_path.stat().st_size,
            sha256=file_sha256(observations_path),
            rows=_count_csv_rows(observations_path),
        ),
        design_inputs_fasta=ManifestOutputEntry(
            path=DESIGN_INPUTS_FASTA_FILENAME,
            exists=True,
            size_bytes=fasta_path.stat().st_size,
            sha256=file_sha256(fasta_path),
            sequences=_count_fasta_records(fasta_path),
        ),
        native_mappings_csv=(
            ManifestOutputEntry(
                path=NATIVE_MAPPINGS_FILENAME,
                exists=True,
                size_bytes=native_mappings_path.stat().st_size,
                sha256=file_sha256(native_mappings_path),
                rows=_count_csv_rows(native_mappings_path),
            )
            if native_mappings_path is not None
            else None
        ),
    )

    poly_threshold = default_for(_POLYNUCLEOTIDE_FILTER_ID)
    polynucleotide_run_requirement = PolynucleotideRunRequirementBlock(
        comparator=FilterComparator.LE,
        threshold=poly_threshold,
        action=FilterAction.FAIL,
        evaluated=False,
        excluded=FilterExclusionCounts(default=0, benchmark=0),
    )

    def _unwidened(key: str) -> GCWideningEntry:
        value = default_for(key)
        return GCWideningEntry(default=value, benchmark=value, source=SettingSource.BUILTIN_PROFILE, widened=False)

    gc_widening = GCWideningBlock(gc_min=_unwidened("gc_min"), gc_max=_unwidened("gc_max"))

    manifest = BenchmarkArtifactManifest(
        schema_version=BENCHMARK_MANIFEST_SCHEMA_VERSION,
        artifact_schema_version=BENCHMARK_ARTIFACT_SCHEMA_VERSION,
        tool_version=__version__,
        created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        invoked_command=tuple(invoked_command) if invoked_command is not None else tuple(sys.argv),
        panel=ManifestPanelBlock(
            panel_id=descriptor.panel_id,
            display_name=descriptor.display_name,
            architecture=descriptor.architecture.value,
            citation=descriptor.citation,
            redistribution=descriptor.redistribution,
            data_present=descriptor.data_present,
            descriptor_hash=descriptor.content_hash(),
        ),
        native_reference=native_reference,
        paired_length=length,
        requested_length=length,
        split_rule=descriptor.split_rule_id,
        inputs=tuple(inputs),
        outputs=outputs,
        counts=counts,
        polynucleotide_run_requirement=polynucleotide_run_requirement,
        # Filled for real once `design_artifact` runs; `{}` here is schema-valid (both fields are
        # untyped dicts) and is overwritten wholesale, never merged, so an empty placeholder cannot
        # leak a stale value into the extended manifest.
        run_policy={},
        default_run_policy={},
        gc_widening=gc_widening,
    )
    manifest_path = artifact_dir / MANIFEST_FILENAME
    write_manifest(manifest, manifest_path)
    return manifest


__all__ = ["BenchmarkPrepareError", "DEFAULT_ARTIFACT_ROOT", "prepare_artifact"]
