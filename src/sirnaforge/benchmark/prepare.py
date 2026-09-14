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
* Of the panels #109/#110 name, only one is vendored: ``tests/unit/data/sirna_efficacy_subset.csv``,
  a 180-row Huesken subset (issues #95/#102). ``_VENDORED_PANEL_CSV`` below names that one path;
  every other registered panel is ``data_present=False`` and refused without an explicit
  ``--panel-csv``.

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
import hashlib
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

from Bio.Seq import Seq

from sirnaforge import __version__
from sirnaforge.benchmark.artifact import (
    BENCHMARK_ARTIFACT_SCHEMA_VERSION,
    BENCHMARK_MANIFEST_SCHEMA_VERSION,
    DESIGN_INPUTS_FASTA_FILENAME,
    MANIFEST_FILENAME,
    OBSERVATIONS_FILENAME,
    BenchmarkArtifactCounts,
    BenchmarkArtifactManifest,
    BenchmarkArtifactOutputs,
    BenchmarkObservation,
    FilterExclusionCounts,
    GCWideningBlock,
    GCWideningEntry,
    ManifestInputEntry,
    ManifestOutputEntry,
    ManifestPanelBlock,
    PolynucleotideRunRequirementBlock,
    artifact_dir_name,
    write_manifest,
    write_observations,
)
from sirnaforge.benchmark.panels import PanelDescriptor, derive_observation, describe_panel
from sirnaforge.config.run_policy import default_for
from sirnaforge.models.policy import FilterAction, FilterComparator, SettingSource

#: ``prepare.py`` -> ``benchmark`` -> ``sirnaforge`` -> ``src`` -> repo root.
_REPO_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_ARTIFACT_ROOT = Path("benchmark_artifacts")

#: The one place a vendored panel's bytes are located. ``panels.PanelDescriptor`` deliberately
#: carries no path (declaring a descriptor and vendoring bytes are different acts -- a descriptor
#: can exist for a panel this repository never ships), so this module -- the one that actually
#: opens files -- is where a ``data_present=True`` panel resolves to a real path. Exactly one entry
#: today; see the module docstring's second verified fact.
_VENDORED_PANEL_CSV: Mapping[str, str] = {"huesken_subset": "tests/unit/data/sirna_efficacy_subset.csv"}

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
        relative = _VENDORED_PANEL_CSV.get(descriptor.panel_id)
        if relative is None:  # pragma: no cover - registry/vendoring mismatch, not a user error
            raise BenchmarkPrepareError(
                f"panel {descriptor.panel_id!r} is registered data_present=True but this module "
                "declares no vendored path for it; that is a bug in prepare.py, not in your invocation"
            )
        return _REPO_ROOT / relative, relative
    if panel_csv is None:
        raise BenchmarkPrepareError(
            f"panel {descriptor.panel_id!r} has no vendored bytes in this repository and no "
            "--panel-csv was given; supply --panel-csv pointing at your own copy of this panel's table"
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


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_raw_rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return tuple(dict(row) for row in reader)


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
    ``"unavailable"``: no descriptor in :data:`sirnaforge.benchmark.panels.PANEL_REGISTRY` maps a
    per-row local-position column (``PanelColumnMapping`` has none), so there is no local coordinate
    this module could honestly report as ``panel_local`` -- that mapping, like native transcript
    resolution, is #110's, not invented here.
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


def prepare_artifact(
    *,
    panel_id: str,
    paired_length: int | None = None,
    panel_csv: Path | str | None = None,
    panel_transcripts: Path | str | None = None,
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
        out_dir: Root directory artifacts are written under.
        overwrite: Replace an existing artifact directory rather than refusing it.
        invoked_command: The command that asked for this, verbatim. Defaults to ``sys.argv`` for a
            direct caller; the CLI slice passes its own argv so the manifest never reconstructs the
            command from options.

    Returns:
        The manifest that was also written to ``manifest.json``.

    Raises:
        BenchmarkPrepareError: A panel-csv/vendoring conflict, a paired length with no default and
            none given, a missing column, or an existing artifact directory without ``overwrite``.
        ValueError: An unregistered ``panel_id`` (:func:`sirnaforge.benchmark.panels.describe_panel`).
        sirnaforge.benchmark.artifact.BenchmarkArtifactError: ``paired_length`` outside 19-23.
    """
    descriptor = describe_panel(panel_id)
    length = _resolve_paired_length(descriptor, paired_length)
    # Validates the 19-23 bound and builds the directory name in the one place that spells the
    # convention, before anything is read from disk.
    dir_name = artifact_dir_name(descriptor.panel_id, length)

    panel_csv_path = Path(panel_csv) if panel_csv is not None else None
    source_path, source_file = _resolve_source_csv(descriptor, panel_csv_path)
    if not source_path.is_file():
        raise BenchmarkPrepareError(f"panel csv not found: {source_path}")
    source_sha256 = _file_sha256(source_path)

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

    raw_rows = _read_raw_rows(source_path)
    observations: list[BenchmarkObservation] = []
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
        observations.append(observation)
        contexts[observation.design_context_id] = context_seq

    observations_path = write_observations(observations, artifact_dir / OBSERVATIONS_FILENAME)
    fasta_path = _write_design_inputs_fasta(observations, contexts, artifact_dir / DESIGN_INPUTS_FASTA_FILENAME)

    kept = sum(1 for observation in observations if observation.compatibility_status == "compatible")
    incompatible = len(observations) - kept
    counts = BenchmarkArtifactCounts(
        observations_in_source=len(observations),
        observations_kept=kept,
        observations_incompatible=incompatible,
        mapped_native=0,
        # No descriptor maps a per-row local-position column (module docstring), so no observation
        # is ever "panel_local" in #109; every row is honestly "unavailable" until #110.
        mapped_panel_local=0,
        mapping_unavailable=len(observations),
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
                sha256=_file_sha256(transcripts_path),
                size_bytes=transcripts_path.stat().st_size,
            )
        )

    outputs = BenchmarkArtifactOutputs(
        observations_csv=ManifestOutputEntry(
            path=str(observations_path),
            exists=True,
            size_bytes=observations_path.stat().st_size,
            sha256=_file_sha256(observations_path),
            rows=_count_csv_rows(observations_path),
        ),
        design_inputs_fasta=ManifestOutputEntry(
            path=str(fasta_path),
            exists=True,
            size_bytes=fasta_path.stat().st_size,
            sha256=_file_sha256(fasta_path),
            sequences=_count_fasta_records(fasta_path),
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
