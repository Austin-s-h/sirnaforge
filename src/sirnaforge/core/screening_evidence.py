"""Evidence producer/parser shared by the Nextflow tasks and workflow.py, per #100.

A screening run's honesty rests on being able to tell "we asked and it ran clean", "we asked and it
failed", "we never asked" and "it ran but we withheld part of the answer" apart, for every
channel x species x guide-set unit. :mod:`sirnaforge.models.evidence` defines those four shapes;
this module is what actually writes them inside a task, reads them back, and reconciles a
:class:`~sirnaforge.models.evidence.ScreeningPlan` (what a run intended) against the envelopes a
run actually produced -- so a species whose reference build crashed produces a FAILED record
instead of silently vanishing from the count. Named ``screening_evidence`` rather than ``evidence``
so it never shadows the pure-type module it builds on.

Requiredness is not decided here: :class:`~sirnaforge.models.policy.EvidenceRequirements` is that
authority, and :func:`completed_pairs` only says what ran to completion, never what was needed.
"""

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from types import MappingProxyType
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from sirnaforge.models.evidence import (
    EVIDENCE_SCHEMA_VERSION,
    EvidenceStatus,
    ObservedCount,
    ObservedCounts,
    ScreeningEvidence,
    ScreeningEvidenceEntry,
    ScreeningPlan,
    ScreeningPlanEntry,
)
from sirnaforge.models.policy import ScreeningChannel

EVIDENCE_FILE_SUFFIX = "_evidence.json"
"""Suffix of one per-unit envelope file, e.g. ``transcriptome_human_evidence.json``."""

RECONCILIATION_FILENAME = "evidence.json"
"""Name of the aggregated reconciliation file. Deliberately does not end in ``_evidence.json`` so
:func:`collect_evidence`'s glob for per-unit envelopes never picks it up as one."""

MISSING_EVIDENCE_DETAIL = "planned but no evidence was published for this unit"
"""Detail attached to a plan entry that reconciliation could not match to any observed envelope."""

#: The 3-tuple every envelope and plan entry is joined on. Not ``reference_id``: an in-container
#: emitter knows only species and an index path, while the plan holds full reference identity, so
#: joining on the 4-tuple ``.key`` would never match. The guide-set digest already carries the
#: identity that matters -- which guides were screened -- and a species holds at most one
#: reference per run by construction.
JoinKey = tuple[str, str, str]


class EvidenceProducer(str, Enum):
    """Which component decided the status carried in one evidence envelope.

    Attributes:
        OFFTARGET_ANALYSIS: The transcriptome alignment task.
        MIRNA_SEED_ANALYSIS: The miRNA seed-matching task.
        AGGREGATE_RESULTS: The aggregation step, for units it can speak to directly.
        WORKFLOW_SYNTHESIS: ``workflow.py``, for an entry it derived rather than observed.
        BASIC_ANALYSIS: The Nextflow-unavailable fallback path.
        STUB: A ``-stub-run`` execution; nothing actually aligned.
    """

    OFFTARGET_ANALYSIS = "offtarget_analysis"
    MIRNA_SEED_ANALYSIS = "mirna_seed_analysis"
    AGGREGATE_RESULTS = "aggregate_results"
    WORKFLOW_SYNTHESIS = "workflow_synthesis"
    BASIC_ANALYSIS = "basic_analysis"
    STUB = "stub"


class EvidenceSource(str, Enum):
    """How one reconciled entry came to exist.

    Attributes:
        ENVELOPE: Read from a per-unit envelope a task actually wrote.
        LEGACY_SUMMARY: Inferred from a pre-#100 summary file that carries no envelope at all.
        SYNTHESIZED: Filled in by reconciliation itself, because nothing was observed.
    """

    ENVELOPE = "envelope"
    LEGACY_SUMMARY = "legacy_summary"
    SYNTHESIZED = "synthesized"


class EvidenceEnvelope(BaseModel):
    """The versioned, on-disk shape of one per-unit evidence file.

    Attributes:
        schema_version: Payload shape version; a reader that sees any other value treats the unit
            as unobserved rather than raising.
        producer: Which component wrote this envelope.
        source: How this envelope came to exist. Defaults to ``ENVELOPE`` because most callers of
            :func:`write_evidence` are the task that produced the result directly.
        entry: The evidence itself, unmodified from :class:`~sirnaforge.models.evidence.ScreeningEvidenceEntry`.
    """

    schema_version: str = Field(default=EVIDENCE_SCHEMA_VERSION, description="Payload shape version")
    producer: EvidenceProducer = Field(description="Component that decided this entry's status")
    source: EvidenceSource = Field(default=EvidenceSource.ENVELOPE, description="How this envelope came to exist")
    entry: ScreeningEvidenceEntry = Field(description="The evidence itself")

    model_config = ConfigDict(frozen=True, extra="forbid")


def guide_set_digest(path: str | Path) -> str:
    """Digest a guide-set FASTA the same way the plan and every envelope key on it.

    Byte-identical to ``workflow.py``'s ``self._file_hash_sha256(path)[:16]``: two computations of
    the same file's digest that disagreed would silently split one guide set into two join keys.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()[:16]


def join_key(entry: ScreeningPlanEntry | ScreeningEvidenceEntry) -> JoinKey:
    """The (channel, species, guide_set_digest) identity a plan entry and its evidence share."""
    return (entry.channel.value, entry.species, entry.guide_set_digest)


def evidence_filename(*, channel: ScreeningChannel, species: str) -> str:
    """Filename one per-unit envelope is written to, e.g. ``transcriptome_human_evidence.json``."""
    return f"{channel.value}_{species}{EVIDENCE_FILE_SUFFIX}"


def write_evidence(
    output_dir: str | Path,
    *,
    producer: EvidenceProducer,
    entry: ScreeningEvidenceEntry,
    source: EvidenceSource = EvidenceSource.ENVELOPE,
) -> Path:
    """Write one per-unit envelope inside a Nextflow task, beside its TSV/summary outputs."""
    envelope = EvidenceEnvelope(producer=producer, source=source, entry=entry)
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / evidence_filename(channel=entry.channel, species=entry.species)
    path.write_text(envelope.model_dump_json(indent=2))
    return path


def read_evidence(path: str | Path) -> EvidenceEnvelope | None:
    """Read one envelope, or ``None`` on anything that makes it unobservable.

    Missing file, unreadable file, malformed JSON and an unrecognised ``schema_version`` are all
    treated identically: the caller must read the unit as unobserved, never as complete. Bumping
    ``EVIDENCE_SCHEMA_VERSION`` is therefore a breaking change by design, not by accident.
    """
    try:
        raw = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None
    if not isinstance(raw, dict) or raw.get("schema_version") != EVIDENCE_SCHEMA_VERSION:
        return None
    try:
        return EvidenceEnvelope.model_validate(raw)
    except ValueError:
        return None


def collect_evidence(root: str | Path) -> tuple[EvidenceEnvelope, ...]:
    """Every readable per-unit envelope under ``root``, sorted by join key for a stable order."""
    envelopes = [
        envelope
        for path in Path(root).rglob(f"*{EVIDENCE_FILE_SUFFIX}")
        if (envelope := read_evidence(path)) is not None
    ]
    envelopes.sort(key=lambda envelope: join_key(envelope.entry))
    return tuple(envelopes)


@dataclass(frozen=True, slots=True)
class Reconciliation:
    """A plan reconciled against what was actually observed.

    Attributes:
        plan: What the run intended to screen.
        evidence: One entry per plan entry, plus any unplanned envelope kept for visibility.
        sources: How each entry, by join key, came to exist.
        unplanned: Join keys present in ``evidence`` but absent from ``plan``.
    """

    plan: ScreeningPlan
    evidence: ScreeningEvidence
    sources: Mapping[JoinKey, EvidenceSource]
    unplanned: tuple[JoinKey, ...]

    def statuses(self) -> Mapping[JoinKey, EvidenceStatus]:
        """Each reconciled entry's status, keyed by join key."""
        return {join_key(entry): entry.status for entry in self.evidence.entries}

    def keys_with_status(self, status: EvidenceStatus) -> tuple[JoinKey, ...]:
        """Join keys whose reconciled entry holds exactly ``status``."""
        return tuple(key for key, entry_status in self.statuses().items() if entry_status is status)


def reconcile(
    plan: ScreeningPlan,
    observed: Sequence[EvidenceEnvelope],
    *,
    missing_detail: str = MISSING_EVIDENCE_DETAIL,
) -> Reconciliation:
    """Reconcile expected screening units against what was actually observed.

    Every plan entry with no matching envelope becomes a FAILED entry with ``missing_detail`` --
    this is how a reference whose index build crashed, or a whole pipeline that aborted before
    aggregation, still produces a failure record instead of vanishing from the count. An observed
    envelope with no plan entry is kept in ``evidence.entries`` and listed in ``unplanned``, never
    dropped. A matched envelope's ``reference_id`` is replaced with the plan's own: an in-container
    emitter knows only species and an index path, so the plan's identity is authoritative.
    """
    observed_by_key: dict[JoinKey, EvidenceEnvelope] = {}
    for obs_envelope in observed:
        observed_by_key[join_key(obs_envelope.entry)] = obs_envelope

    entries: list[ScreeningEvidenceEntry] = []
    sources: dict[JoinKey, EvidenceSource] = {}
    matched: set[JoinKey] = set()

    for plan_entry in plan.entries:
        key = join_key(plan_entry)
        envelope = observed_by_key.get(key)
        if envelope is None:
            entries.append(failed_entry(plan_entry, missing_detail))
            sources[key] = EvidenceSource.SYNTHESIZED
        else:
            entries.append(envelope.entry.model_copy(update={"reference_id": plan_entry.reference_id}))
            sources[key] = envelope.source
            matched.add(key)

    unplanned = tuple(sorted(key for key in observed_by_key if key not in matched))
    for key in unplanned:
        envelope = observed_by_key[key]
        entries.append(envelope.entry)
        sources[key] = envelope.source

    return Reconciliation(
        plan=plan,
        evidence=ScreeningEvidence(entries=tuple(entries)),
        sources=MappingProxyType(sources),
        unplanned=unplanned,
    )


def completed_pairs(
    evidence: ScreeningEvidence,
    *,
    strict: bool = False,
    sources: Mapping[JoinKey, EvidenceSource] | None = None,
) -> frozenset[tuple[str, str]]:
    """(channel, species) pairs whose evidence completed -- the shape the eligibility engine reads.

    Only ``COMPLETE`` counts. ``CENSORED`` never does, even with ``strict=False``: a lower-bound
    count cannot show a ceiling was respected. ``strict=True`` additionally drops any entry whose
    join key maps to ``EvidenceSource.LEGACY_SUMMARY`` in ``sources`` -- the fallback heuristic for
    a pre-#100 run directory that carries no envelope at all. ``sources`` is optional and additive
    to the contract's two-argument call: a caller with no source information gets the lenient
    (``strict=False``-equivalent) answer regardless of what it passes for ``strict``.
    """
    pairs: set[tuple[str, str]] = set()
    for entry in evidence.entries:
        if entry.status is not EvidenceStatus.COMPLETE:
            continue
        if strict and sources is not None and sources.get(join_key(entry)) is EvidenceSource.LEGACY_SUMMARY:
            continue
        pairs.add((entry.channel.value, entry.species))
    return frozenset(pairs)


def build_plan(
    *,
    guide_set_digest: str,
    transcriptome: Sequence[tuple[str, str | None]],
    mirna_species: Sequence[str],
    search_settings: Mapping[str, str | int | float | bool | None] | None = None,
) -> ScreeningPlan:
    """Build the expected plan from what a run requested, before any reference can be dropped.

    Takes requested species, not resolved references: a species whose reference later fails to
    resolve must still appear in the plan, or its absence reconciles as never having been asked
    for rather than as a failure. An unrequested channel gets no plan entry at all, which is what
    lets ``NOT_REQUESTED`` be distinguished from ``FAILED`` downstream.
    """
    settings = dict(search_settings or {})
    entries = tuple(
        ScreeningPlanEntry(
            channel=ScreeningChannel.TRANSCRIPTOME,
            species=species,
            reference_id=reference_id,
            guide_set_digest=guide_set_digest,
            search_settings=settings,
        )
        for species, reference_id in transcriptome
    ) + tuple(
        ScreeningPlanEntry(
            channel=ScreeningChannel.MIRNA_SEED,
            species=species,
            guide_set_digest=guide_set_digest,
            search_settings=settings,
        )
        for species in mirna_species
    )
    return ScreeningPlan(entries=entries)


def censored_counts(
    *,
    retained: int,
    cap: int | None,
    pre_cap: int | None,
    distinct_transcripts: int | None = None,
    distinct_genes: int | None = None,
    unresolved_gene_sites: int | None = None,
) -> ObservedCounts:
    """Counts for a channel whose retention cap makes every reported quantity a lower bound.

    ``sites`` is always ``is_lower_bound=True, truncated=True``: this helper exists only for the
    case a caller has already decided is CENSORED, so an untruncated search must not call it.
    ``pre_cap``, when known, must be at least ``retained`` -- a cap cannot have retained more than
    the search actually produced.
    """
    if pre_cap is not None and pre_cap < retained:
        raise ValueError(f"pre-cap count {pre_cap} is lower than the retained count {retained}")

    def _lower_bound(value: int | None) -> ObservedCount:
        return ObservedCount() if value is None else ObservedCount(value=value, is_lower_bound=True)

    return ObservedCounts(
        sites=ObservedCount(value=retained, is_lower_bound=True, cap=cap, truncated=True),
        distinct_transcripts=_lower_bound(distinct_transcripts),
        distinct_genes=_lower_bound(distinct_genes),
        unresolved_gene_sites=_lower_bound(unresolved_gene_sites),
    )


def failed_entry(plan_entry: ScreeningPlanEntry, detail: str) -> ScreeningEvidenceEntry:
    """A FAILED entry for a planned unit, carrying the plan's own identity forward."""
    return ScreeningEvidenceEntry(
        channel=plan_entry.channel,
        species=plan_entry.species,
        reference_id=plan_entry.reference_id,
        guide_set_digest=plan_entry.guide_set_digest,
        status=EvidenceStatus.FAILED,
        detail=detail,
    )


def not_requested_entry(
    channel: ScreeningChannel, species: str, guide_set_digest: str, detail: str | None = None
) -> ScreeningEvidenceEntry:
    """A NOT_REQUESTED entry for a channel/species the plan deliberately holds no entry for.

    Emitted only by ``workflow_synthesis``: a task that ran was, by definition, requested, so no
    in-container emitter ever produces this status.
    """
    return ScreeningEvidenceEntry(
        channel=channel,
        species=species,
        guide_set_digest=guide_set_digest,
        status=EvidenceStatus.NOT_REQUESTED,
        detail=detail,
    )


def reconciliation_payload(reconciliation: Reconciliation) -> dict[str, Any]:
    """Serialise a reconciliation to the on-disk ``evidence.json`` shape."""
    return {
        "schema_version": EVIDENCE_SCHEMA_VERSION,
        "plan": reconciliation.plan.model_dump(mode="json"),
        "evidence": reconciliation.evidence.model_dump(mode="json"),
        "sources": {"|".join(key): source.value for key, source in reconciliation.sources.items()},
        "unplanned": [list(key) for key in reconciliation.unplanned],
    }


def parse_reconciliation_payload(payload: Mapping[str, Any] | None) -> Reconciliation | None:
    """Parse an ``evidence.json`` payload, or ``None`` on anything that makes it unusable."""
    if payload is None or payload.get("schema_version") != EVIDENCE_SCHEMA_VERSION:
        return None
    try:
        plan = ScreeningPlan.model_validate(payload["plan"])
        evidence = ScreeningEvidence.model_validate(payload["evidence"])
        sources: dict[JoinKey, EvidenceSource] = {}
        for key, value in payload.get("sources", {}).items():
            channel, species, digest = key.split("|", 2)
            sources[(channel, species, digest)] = EvidenceSource(value)
        unplanned = tuple(tuple(key) for key in payload.get("unplanned", []))
    except (KeyError, ValueError):
        return None
    return Reconciliation(plan=plan, evidence=evidence, sources=MappingProxyType(sources), unplanned=unplanned)


def write_reconciliation(output_dir: str | Path, reconciliation: Reconciliation) -> Path:
    """Write the aggregated reconciliation file, unconditionally -- its presence is the contract."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / RECONCILIATION_FILENAME
    path.write_text(json.dumps(reconciliation_payload(reconciliation), indent=2))
    return path
