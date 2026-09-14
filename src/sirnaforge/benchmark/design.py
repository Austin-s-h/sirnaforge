"""``benchmark design``: run the fixed-length designer once, derive two filter verdict sets (#109, bm-design slice).

One artifact directory in, extended in place: ``design_artifact`` reads ``observations.csv`` and
``design_inputs.fasta`` from a directory ``benchmark prepare`` already wrote, runs
:class:`~sirnaforge.core.design.SiRNADesigner` under a *benchmark* policy that may widen ``gc_min``/
``gc_max``, joins the resulting candidates back onto their observations by RNA-normalised guide
sequence, and writes ``candidates_all.csv`` and ``accounting.csv`` plus an extended
``manifest.json``. It does not read a panel or build the input FASTA -- that is ``prepare``, in a
module this one does not own.

**The run is executed once.** ``accounting.csv`` carries two verdict sets per row
(``default_filter_status``/``default_filter_reasons`` and ``benchmark_filter_status``/
``benchmark_filter_reasons``), but there is exactly one design pass, under the benchmark policy.
The default set is *re-derived* by re-applying the default policy's ``FilterDescriptor`` objects to
the values that one run already recorded in ``SiRNACandidate.filter_observed`` (#100/#103's own
client-side re-application path, :func:`~sirnaforge.config.run_policy.filters_from_manifest`'s
sibling operation performed in-process instead of from a manifest). That re-derivation is *exact*,
not approximate, only because widening is monotone: every design-time observed value
(``gc_content``, ``max_poly_run_length``, ...) is a property of the candidate, not of the threshold,
so a benchmark run enumerates a superset of what a default run would have enumerated. This is why
:func:`design_artifact` **refuses a narrowing** of either GC bound (``gc_min`` raised above, or
``gc_max`` lowered below, the shipped default) before doing any work at all: a narrowed run would
report ``default_filter_status=pass`` for a candidate it never enumerated, which is a claim about a
run that did not happen, not a filter result.

**Why the polynucleotide-run gate cannot leak into the widening.** Only ``gc_min``/``gc_max`` are
ever placed in the ``stated=`` mapping this module builds for
:func:`~sirnaforge.config.run_policy.resolve_run_policy`; ``max_poly_runs`` (``<= 3``, action
``fail``) is resolved identically for the benchmark and the default policy in every call this module
makes, so it fails the same candidates under both -- which is what "kept active and explicitly
recorded" (#109's acceptance criterion) means in practice: the same reason appears in both
``default_filter_reasons`` and ``benchmark_filter_reasons`` for the same row.

**Why ``candidates_all.csv`` is not just ``DesignResult.candidates``.** ``SiRNADesigner`` drops a
candidate outright when an enumeration-time gate (``gc_content_min/max``, ``max_poly_runs``)
rejects it with a ``fail`` action -- it survives only in ``DesignResult.rejected_candidates``. A
homopolymer observation that a *default* run would also have rejected still has to appear in the
accounting with ``entered_design=True`` and a real verdict on both columns (#109's own filter-
accounting contract: "a gate cannot be blamed for a candidate that was never built" only applies
when no candidate was built at all). So this module writes ``candidates.candidates +
candidates.rejected_candidates`` through the unchanged ``DesignResult.save_csv``/
``build_candidate_row`` writer -- "_all" in the filename is the ALL that name promises.
"""

from __future__ import annotations

import hashlib
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from sirnaforge import __version__
from sirnaforge.benchmark.artifact import (
    ACCOUNTING_FILENAME,
    CANDIDATES_ALL_FILENAME,
    DESIGN_INPUTS_FASTA_FILENAME,
    MANIFEST_FILENAME,
    OBSERVATIONS_FILENAME,
    BenchmarkAccountingRow,
    BenchmarkArtifactManifest,
    BenchmarkObservation,
    FilterExclusionCounts,
    GCWideningBlock,
    GCWideningEntry,
    ManifestInputEntry,
    ManifestOutputEntry,
    PolynucleotideRunRequirementBlock,
    read_manifest,
    read_observations,
    write_accounting,
    write_manifest,
)
from sirnaforge.config.run_policy import (
    DEFAULT_PROFILE_NAME,
    EntryPoint,
    ResolvedFilter,
    ResolvedRunPolicy,
    RunPolicyError,
    default_for,
    resolve_run_policy,
)
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.policy import DECLARED_FILTER_IDS, FilterAction, FilterEvaluation, SettingSource
from sirnaforge.models.sirna import DesignResult, SiRNACandidate

#: The one gate this surface must never be able to move (#109 scope: "keep the polynucleotide-run
#: requirement active"). Named once so a comment doesn't have to restate the id.
_POLYNUCLEOTIDE_FILTER_ID = "max_poly_runs"

#: ``FILTER_SPECS`` declares this gate's stage as ``design`` -- always "evaluable" per the resolver
#: -- but the only code that ever records an observed value for it,
#: ``SiRNADesigner.stamp_repeat_verdict``, is called from ``workflow.py``'s post-screening step
#: against a real cDNA reference, never from ``design_from_file``. Native transcript mapping is
#: explicitly out of scope for #109, so a benchmark run has no more evidence for this gate than it
#: has for the nine ``post_screen`` gates it already treats as not evaluated -- treating it as
#: "evaluated, evidence unknown" instead would report ``unknown`` on almost every row of a real
#: panel (confirmed against the vendored Huesken subset: 0 of 180 entered candidates read ``pass``
#: under either policy before this exclusion existed), which drowns every GC/polynucleotide verdict
#: this module exists to report in noise from a channel #109 never claimed to hold.
_UNMEASURED_BY_DESIGN_FROM_FILE: frozenset[str] = frozenset({"max_repeat_transcript_fraction"})


@dataclass(frozen=True)
class DesignArtifactResult:
    """What one ``design_artifact`` call produced, for a thin CLI to summarise.

    Deliberately not the manifest itself: a caller that wants the full record reads
    ``manifest.json`` back with :func:`sirnaforge.benchmark.artifact.read_manifest`, which is the
    one parser this contract has.
    """

    artifact_dir: Path
    candidates_all_csv: Path
    accounting_csv: Path
    manifest_path: Path
    entered_design: int
    no_candidate: int
    default_pass: int
    benchmark_pass: int


def _refuse_narrowing(gc_min: float | None, gc_max: float | None) -> None:
    """Refuse a benchmark run that would narrow either GC bound relative to the shipped default.

    The load-bearing constraint of the whole filter-accounting design (module docstring): the
    default verdict set is re-derived from candidates the *benchmark* run enumerated, which is only
    exact if the benchmark run enumerated a superset of what a default run would have. Raised before
    any file in the artifact directory is even opened, so a caller sees this before any design work.
    """
    if gc_min is not None and gc_min > default_for("gc_min"):
        raise RunPolicyError(
            f"--gc-min {gc_min} narrows the default floor of {default_for('gc_min')}; the benchmark "
            "interface may only widen gc_min, never raise it, because a narrowed run cannot stand in "
            "for the default run its default_filter_status claims to re-derive"
        )
    if gc_max is not None and gc_max < default_for("gc_max"):
        raise RunPolicyError(
            f"--gc-max {gc_max} narrows the default ceiling of {default_for('gc_max')}; the benchmark "
            "interface may only widen gc_max, never lower it, for the same reason gc_min may not be raised"
        )


def _resolve_policies(
    *,
    paired_length: int,
    design_mode: str | None,
    gc_min: float | None,
    gc_max: float | None,
    policy_config: Path | str | None,
    filter_actions: Mapping[str, str] | None,
) -> tuple[ResolvedRunPolicy, ResolvedRunPolicy]:
    """Resolve the benchmark and default policies, identical except for the GC widening keys.

    Two separate :func:`resolve_run_policy` calls rather than one policy mutated after the fact: the
    resolver is pure and its provenance trail is what a manifest reader trusts, so the default policy
    has to come from asking the resolver the same question again with ``gc_min``/``gc_max`` simply
    absent from ``stated`` -- not from copying the benchmark policy and overwriting two numbers.
    """
    common_stated: dict[str, Any] = {"sirna_length": paired_length}
    benchmark_stated = dict(common_stated)
    if gc_min is not None:
        benchmark_stated["gc_min"] = gc_min
    if gc_max is not None:
        benchmark_stated["gc_max"] = gc_max

    benchmark_policy = resolve_run_policy(
        entry_point=EntryPoint.BENCHMARK_COMMAND,
        design_mode=design_mode,
        profile_name=DEFAULT_PROFILE_NAME,
        config_file=policy_config,
        stated=benchmark_stated,
        filter_actions=filter_actions,
    )
    default_policy = resolve_run_policy(
        entry_point=EntryPoint.BENCHMARK_COMMAND,
        design_mode=design_mode,
        profile_name=DEFAULT_PROFILE_NAME,
        config_file=policy_config,
        stated=common_stated,
        filter_actions=filter_actions,
    )
    return benchmark_policy, default_policy


def _rna(sequence: str) -> str:
    """RNA-normalise a guide for comparison: upper-case, T read as U.

    A local twin of ``core.design._as_rna`` rather than an import of it: that name is private to its
    module, and the join here compares *measured* panel sequences (which this repository's fixtures
    write in RNA letters) against *designed* candidate sequences (stored as DNA), so both sides need
    the same normalisation independently of which module last touched either string.
    """
    return sequence.upper().replace("T", "U")


def _count_fasta_records(path: Path) -> int:
    """Count FASTA headers without a Biopython parse, mirroring ``workflow.py``'s own count."""
    if not path.is_file():
        return 0
    with path.open("r") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def _sha256_file(path: Path) -> str:
    """Plain hex SHA-256 of a file's bytes; no ``sha256:`` prefix, matching ``source_sha256``'s convention."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _match_candidate(
    observation: BenchmarkObservation, candidates_by_context: Mapping[str, Sequence[SiRNACandidate]]
) -> tuple[SiRNACandidate | None, str]:
    """Find the one enumerated candidate that reproduces this observation's measured guide, if any.

    ``exact`` when the design produced the observation's full measured guide (only possible when the
    panel's guide length already equals the artifact's paired length); ``paired_core_exact`` when
    only the declared slice matches, which is the ordinary case for a panel whose measured guide is
    longer than the fixed design length. Comparison is RNA-normalised in both directions (#109: T and
    U must read equal, since a panel and the designer do not agree on which letter they store).
    """
    candidates = candidates_by_context.get(observation.design_context_id, ())
    full_rna = _rna(observation.full_guide_sequence)
    paired_rna = _rna(observation.paired_guide_sequence)
    for candidate in candidates:
        if _rna(candidate.guide_sequence) == full_rna:
            return candidate, "exact"
    for candidate in candidates:
        if _rna(candidate.guide_sequence) == paired_rna:
            return candidate, "paired_core_exact"
    return None, "none"


def _evaluate_row(observed: Mapping[str, float | None], filters: Sequence[ResolvedFilter]) -> tuple[str, str]:
    """Re-apply one policy's FAIL-action descriptors to values a run already recorded.

    Only ``FilterAction.FAIL`` gates decide the aggregate status and appear in ``reasons``: a
    ``warn``-action gate can still individually evaluate ``fail`` (its comparator is unaffected by
    its action), but it never rejected the candidate, so it is not a *reason* this row failed. This
    is deliberately not ``SiRNACandidate.FilterStatus``: that enum is single-valued and first-failure-
    wins (#101), which cannot represent two FAIL-action gates failing the same candidate at once --
    exactly what a 4-nt homopolymer with an out-of-range GC content needs to represent.

    Excludes :data:`_UNMEASURED_BY_DESIGN_FROM_FILE` from evaluation entirely, alongside gates the
    resolver itself reports as not evaluated: both are "no claim", just discovered by a different
    route (one from the resolved policy, one from what this run's code path can ever observe).

    Returns:
        ``(status, reasons)`` where ``status`` is a :class:`FilterEvaluation` value and ``reasons``
        is a ``;``-joined, :data:`DECLARED_FILTER_IDS`-ordered list of failing filter ids.
    """
    by_id = {resolved.filter_id: resolved for resolved in filters}
    any_fail = False
    any_unknown = False
    reasons: list[str] = []
    for filter_id in DECLARED_FILTER_IDS:
        if filter_id in _UNMEASURED_BY_DESIGN_FROM_FILE:
            continue
        resolved = by_id.get(filter_id)
        if resolved is None or resolved.descriptor.action is not FilterAction.FAIL or not resolved.is_evaluated:
            continue
        value = observed.get(filter_id)
        if value is None:
            any_unknown = True
            continue
        if not resolved.descriptor.comparator.passes(value, resolved.descriptor.threshold):  # type: ignore[arg-type]
            any_fail = True
            reasons.append(filter_id)
    if any_fail:
        status = FilterEvaluation.FAIL
    elif any_unknown:
        status = FilterEvaluation.UNKNOWN
    else:
        status = FilterEvaluation.PASS
    return status.value, ";".join(reasons)


def _excluded_counts(rows: Sequence[BenchmarkAccountingRow]) -> dict[str, FilterExclusionCounts]:
    """Per declared filter, how many rows that gate's id appears in each reason column."""
    counts: dict[str, dict[str, int]] = {filter_id: {"default": 0, "benchmark": 0} for filter_id in DECLARED_FILTER_IDS}
    for row in rows:
        for filter_id in row.default_filter_reasons.split(";") if row.default_filter_reasons else ():
            counts[filter_id]["default"] += 1
        for filter_id in row.benchmark_filter_reasons.split(";") if row.benchmark_filter_reasons else ():
            counts[filter_id]["benchmark"] += 1
    return {filter_id: FilterExclusionCounts(**tally) for filter_id, tally in counts.items()}


def _gc_widening_entry(policy: ResolvedRunPolicy, key: str, default_value: float) -> GCWideningEntry:
    """One GC bound's widening record: provenance decides ``widened``, not a bare value comparison.

    #101 already established that a value test cannot tell an omitted option from one explicitly
    typed with the default; the same trap applies here, so ``widened`` reads ``source``, not
    ``benchmark != default`` on its own (though the two agree in practice, since narrowing is
    refused and an unstated bound resolves to ``BUILTIN_PROFILE``, never ``EXPLICIT``).
    """
    value = policy.value_of(key)
    source = policy.source_of(key)
    widened = source is SettingSource.EXPLICIT and value != default_value
    return GCWideningEntry(default=default_value, benchmark=value, source=source, widened=widened)


def _empty_design_result(fasta_path: Path, policy: ResolvedRunPolicy) -> DesignResult:
    """A zero-candidate ``DesignResult``, for an artifact whose ``design_inputs.fasta`` has no records.

    ``SiRNADesigner.design_from_file`` raises on an empty FASTA (it has no sequence to report against
    at all); a benchmark artifact with zero compatible observations is a legitimate, if unusual,
    input, and must still produce header-only ``candidates_all.csv``/``accounting.csv`` rather than
    an uncaught exception.
    """
    return DesignResult(
        input_file=str(fasta_path),
        parameters=policy.design_parameters,
        candidates=[],
        top_candidates=[],
        total_sequences=0,
        total_candidates=0,
        filtered_candidates=0,
        processing_time=0.0,
        tool_versions={},
        rejected_candidates=[],
    )


def _prepared_inputs(artifact_path: Path) -> tuple[ManifestInputEntry, ...]:
    """Checksum the two prepared files this design pass actually read.

    Deliberately not a checksum of ``manifest.json``: this function overwrites that file moments
    later, so a hash of it would name bytes the finished artifact no longer holds, and on a second
    ``benchmark design`` over the same directory it would hash an *already designed* manifest under a
    role called "prepared". ``observations.csv`` and ``design_inputs.fasta`` are written by ``prepare``
    and never rewritten here, so hashing them is both stable across re-designs and the honest answer
    to "which bytes did this run consume". Paths are the fixed inner filenames, per
    :class:`~sirnaforge.benchmark.artifact.ManifestOutputEntry`'s own convention.
    """
    entries = []
    for role, filename in (
        ("prepared_observations_csv", OBSERVATIONS_FILENAME),
        ("prepared_design_inputs_fasta", DESIGN_INPUTS_FASTA_FILENAME),
    ):
        path = artifact_path / filename
        entries.append(
            ManifestInputEntry(role=role, path=filename, sha256=_sha256_file(path), size_bytes=path.stat().st_size)
        )
    return tuple(entries)


def _merge_inputs(
    existing: Sequence[ManifestInputEntry], recorded: Sequence[ManifestInputEntry]
) -> tuple[ManifestInputEntry, ...]:
    """Replace, never append, an input entry whose role this run is re-recording.

    ``role`` is a key within ``inputs`` (:class:`~sirnaforge.benchmark.artifact.ManifestInputEntry`).
    Appending instead grew the provenance list by one entry per ``benchmark design`` re-run over the
    same artifact, leaving a reader two same-role entries with different checksums and no rule for
    which one describes the artifact as it now stands.
    """
    by_role = {entry.role: entry for entry in existing}
    order = [entry.role for entry in existing]
    for entry in recorded:
        if entry.role not in by_role:
            order.append(entry.role)
        by_role[entry.role] = entry
    return tuple(by_role[role] for role in order)


def _write_empty_candidates_csv(path: Path) -> None:
    """Write a header-only ``candidates_all.csv`` when a run produced zero candidates.

    ``DesignResult.save_csv`` builds its DataFrame from ``[build_candidate_row(c) for c in
    candidates]``; with zero candidates that list is empty and ``pd.DataFrame([])`` has no columns
    at all, so pandera's ``add_missing_columns`` refuses to fabricate a value for a non-nullable
    column with nothing to infer a dtype from. A benchmark artifact whose FASTA has zero enumerable
    records (or, for a caller that stubs the designer, an artifact with zero *matched* observations
    at all) is a legitimate input, so this writes the same header a non-empty run would, with zero
    data rows, instead of surfacing pandera's schema error as a crash.
    """
    from sirnaforge.models.schemas import SiRNACandidateSchema  # noqa: PLC0415 -- avoid a module-load cycle

    columns = [
        *SiRNACandidateSchema.to_schema().columns,
        *(f"{filter_id}_verdict" for filter_id in DECLARED_FILTER_IDS),
        *(f"{filter_id}_observed" for filter_id in DECLARED_FILTER_IDS),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(",".join(columns) + "\n")


def design_artifact(
    *,
    artifact_dir: Path | str,
    design_mode: str | None = None,
    gc_min: float | None = None,
    gc_max: float | None = None,
    policy_config: Path | str | None = None,
    filter_actions: Mapping[str, str] | None = None,
    invoked_command: Sequence[str] | None = None,
) -> DesignArtifactResult:
    """Run the fixed-length designer over a prepared artifact and extend it with both verdict sets.

    Args:
        artifact_dir: A directory ``benchmark prepare`` already wrote (``manifest.json``,
            ``observations.csv``, ``design_inputs.fasta``).
        design_mode: ``sirna``/``mirna``; forwarded to :func:`resolve_run_policy` unchanged.
        gc_min: Per-run GC floor widening. ``None`` means unstated (the shipped default applies).
            Refused if it narrows the default (see :func:`_refuse_narrowing`).
        gc_max: Per-run GC ceiling widening, same rule in the other direction.
        policy_config: Optional JSON/TOML policy file, forwarded to both policy resolutions.
        filter_actions: ``filter_id -> off|warn|fail`` overrides, forwarded to both resolutions and
            to the designer. Neither GC bound may be switched off through this path (the resolver
            already refuses switching off a threshold with no absent state); nothing here relaxes
            ``max_poly_runs``.
        invoked_command: ``sys.argv`` for this invocation, recorded in the manifest. Defaults to the
            live ``sys.argv`` so a direct API call still records something truthful.

    Returns:
        A summary of what was written, for a thin CLI to print.

    Raises:
        RunPolicyError: A narrowing was requested, or the resolved parameters are invalid.
        sirnaforge.benchmark.artifact.BenchmarkArtifactError: The artifact directory is missing a
            required file or ``manifest.json``/``observations.csv`` fails schema validation.
    """
    _refuse_narrowing(gc_min, gc_max)

    artifact_path = Path(artifact_dir)
    manifest_path = artifact_path / MANIFEST_FILENAME
    manifest = read_manifest(manifest_path)
    observations = read_observations(artifact_path / OBSERVATIONS_FILENAME)
    fasta_path = artifact_path / DESIGN_INPUTS_FASTA_FILENAME

    benchmark_policy, default_policy = _resolve_policies(
        paired_length=manifest.paired_length,
        design_mode=design_mode,
        gc_min=gc_min,
        gc_max=gc_max,
        policy_config=policy_config,
        filter_actions=filter_actions,
    )

    from sirnaforge.workflow import _policy_filter_actions  # noqa: PLC0415 -- avoid a module-load cycle

    designer = SiRNADesigner(
        benchmark_policy.design_parameters, filter_actions=_policy_filter_actions(benchmark_policy)
    )
    if _count_fasta_records(fasta_path) > 0:
        result = designer.design_from_file(str(fasta_path))
    else:
        result = _empty_design_result(fasta_path, benchmark_policy)

    # "_all": every enumerated candidate, including the ones a fail-action gate dropped at
    # enumeration time (module docstring) -- not merely the survivors `DesignResult.candidates` holds.
    all_candidates = list(result.candidates) + list(result.rejected_candidates)
    candidates_all_path = artifact_path / CANDIDATES_ALL_FILENAME
    if all_candidates:
        csv_result = result.model_copy(update={"candidates": all_candidates})
        csv_result.save_csv(str(candidates_all_path))
    else:
        _write_empty_candidates_csv(candidates_all_path)

    candidates_by_context: dict[str, list[SiRNACandidate]] = {}
    for enumerated in all_candidates:
        candidates_by_context.setdefault(enumerated.transcript_id, []).append(enumerated)

    accounting_rows: list[BenchmarkAccountingRow] = []
    for observation in observations:
        candidate, guide_match = _match_candidate(observation, candidates_by_context)
        if candidate is None:
            accounting_rows.append(
                BenchmarkAccountingRow(
                    observation_id=observation.observation_id,
                    panel_id=observation.panel_id,
                    paired_length=observation.paired_length,
                    candidate_id=None,
                    designed_guide_sequence=None,
                    guide_match=guide_match,  # type: ignore[arg-type]
                    entered_design=False,
                    default_filter_status=FilterEvaluation.NOT_EVALUATED,
                    default_filter_reasons="",
                    benchmark_filter_status=FilterEvaluation.NOT_EVALUATED,
                    benchmark_filter_reasons="",
                )
            )
            continue
        default_status, default_reasons = _evaluate_row(candidate.filter_observed, default_policy.filters)
        benchmark_status, benchmark_reasons = _evaluate_row(candidate.filter_observed, benchmark_policy.filters)
        accounting_rows.append(
            BenchmarkAccountingRow(
                observation_id=observation.observation_id,
                panel_id=observation.panel_id,
                paired_length=observation.paired_length,
                candidate_id=candidate.id,
                designed_guide_sequence=candidate.guide_sequence,
                guide_match=guide_match,  # type: ignore[arg-type]
                entered_design=True,
                default_filter_status=FilterEvaluation(default_status),
                default_filter_reasons=default_reasons,
                benchmark_filter_status=FilterEvaluation(benchmark_status),
                benchmark_filter_reasons=benchmark_reasons,
            )
        )

    accounting_path = artifact_path / ACCOUNTING_FILENAME
    write_accounting(accounting_rows, accounting_path)

    entered_design = sum(1 for row in accounting_rows if row.entered_design)
    no_candidate = len(accounting_rows) - entered_design
    default_pass = sum(1 for row in accounting_rows if row.default_filter_status is FilterEvaluation.PASS)
    benchmark_pass = sum(1 for row in accounting_rows if row.benchmark_filter_status is FilterEvaluation.PASS)
    excluded_by_filter = _excluded_counts(accounting_rows)

    updated_counts = manifest.counts.model_copy(
        update={
            "entered_design": entered_design,
            "no_candidate": no_candidate,
            "default_pass": default_pass,
            "benchmark_pass": benchmark_pass,
            "excluded_by_filter": excluded_by_filter,
        }
    )

    poly_descriptor = benchmark_policy.descriptor(_POLYNUCLEOTIDE_FILTER_ID)
    polynucleotide_run_requirement = PolynucleotideRunRequirementBlock(
        comparator=poly_descriptor.comparator,
        threshold=poly_descriptor.threshold,  # type: ignore[arg-type]
        action=poly_descriptor.action,
        evaluated=benchmark_policy.filter(_POLYNUCLEOTIDE_FILTER_ID).is_evaluated,
        excluded=excluded_by_filter[_POLYNUCLEOTIDE_FILTER_ID],
    )

    gc_widening = GCWideningBlock(
        gc_min=_gc_widening_entry(benchmark_policy, "gc_min", default_for("gc_min")),
        gc_max=_gc_widening_entry(benchmark_policy, "gc_max", default_for("gc_max")),
    )

    # Paths are the fixed inner filenames, matching what `prepare` recorded for its own two outputs:
    # a manifest that named an absolute path could not be compared against a second run of the same
    # bytes under a different --out-dir (`ManifestOutputEntry.path`).
    updated_outputs = manifest.outputs.model_copy(
        update={
            "candidates_all_csv": ManifestOutputEntry(
                path=CANDIDATES_ALL_FILENAME,
                exists=True,
                size_bytes=candidates_all_path.stat().st_size,
                sha256=_sha256_file(candidates_all_path),
                rows=len(all_candidates),
            ),
            "accounting_csv": ManifestOutputEntry(
                path=ACCOUNTING_FILENAME,
                exists=True,
                size_bytes=accounting_path.stat().st_size,
                sha256=_sha256_file(accounting_path),
                rows=len(accounting_rows),
            ),
        }
    )

    updated_manifest = manifest.model_copy(
        update={
            "invoked_command": tuple(invoked_command if invoked_command is not None else sys.argv),
            "inputs": _merge_inputs(manifest.inputs, _prepared_inputs(artifact_path)),
            "outputs": updated_outputs,
            "counts": updated_counts,
            "polynucleotide_run_requirement": polynucleotide_run_requirement,
            "run_policy": benchmark_policy.as_manifest(),
            "default_run_policy": default_policy.as_manifest(),
            "gc_widening": gc_widening,
            "tool_version": __version__,
        }
    )
    # BenchmarkArtifactManifest is frozen; model_copy bypasses __init__ validators, so re-validate
    # explicitly rather than trust that every update above produced a jointly-consistent object.
    updated_manifest = BenchmarkArtifactManifest.model_validate(updated_manifest.model_dump(mode="json"))
    write_manifest(updated_manifest, manifest_path)

    return DesignArtifactResult(
        artifact_dir=artifact_path,
        candidates_all_csv=candidates_all_path,
        accounting_csv=accounting_path,
        manifest_path=manifest_path,
        entered_design=entered_design,
        no_candidate=no_candidate,
        default_pass=default_pass,
        benchmark_pass=benchmark_pass,
    )


__all__ = ["DesignArtifactResult", "design_artifact"]
