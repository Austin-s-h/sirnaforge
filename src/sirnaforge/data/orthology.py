"""Orthologue resolution on stable gene IDs, for cross-species hit classification.

Why this module exists: orthology cannot be decided on gene symbols. HGNC and MGI are separate
nomenclature authorities, so the mouse orthologue of human ``TP53`` is ``Trp53`` -- uppercased
symbol equality misses it, and would equally happily match two unrelated genes that share a symbol
across species. Issue #101 requires explicit orthology evidence for any conservation claim, with
symbol equality retained only as a labelled heuristic.

The resolver lives here, in the data layer, because ``core.hit_classification.classify_hit`` is a
pure function and must not perform network I/O. The workflow resolves a mapping once and passes the
resulting gene-ID set into ``ClassificationContext``.

Cost is up to *two* Compara requests per (query gene x target species) -- the gene-ID route, then
the symbol route when the first resolves nothing -- charged once per resolution, never per hit.
An on-disk cache under :data:`ORTHOLOGY_CACHE_SUBDIR` charges that cost once per
(question x target species) rather than once per run: a re-run, a second guide panel against the
same gene, and an offline run behind a warm cache all resolve without a request. What it guarantees
is narrow and stated on every row it serves -- one document per *single* target species, so a
five-species screen that resolved four and failed one re-asks only the fifth; a resolved absence is
cached because it is a real Compara answer, while a failed or partial lookup never is, so a firewall
blip cannot become a 30-day claim that no orthologue exists; and a served answer is published as
:data:`SOURCE_COMPARA_CACHE`, so an offline run cannot claim a REST call it never made. Two guards
keep an unreachable Compara off the critical path: a transport-level failure is never retried, and
:meth:`OrthologueMapping.from_file` resolves from a user-supplied mapping with no network at all
(issue #101: an offline path is required, not optional).

:data:`ORTHOLOGY_BUDGET_SECONDS` is a suggested ceiling for a host that *drops* packets rather than
refusing them, and it is **opt-in** -- the default is no ceiling. Abandoning a slow-but-working
Compara reports a resolvable species as ``unresolved``, and unresolved conservation still publishes
``0.0`` rather than null, so a ceiling trades a stall for a fabricated measurement. Bounding that
case belongs with #101's null-conservation fix.
"""

from __future__ import annotations

import asyncio
import json
import logging
import re
import socket
import ssl
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import aiohttp

from sirnaforge.data.base import ENSEMBL_MAX_ATTEMPTS, ensembl_request_json, ensembl_session
from sirnaforge.data.species_registry import ensembl_species_slug, normalize_species_name
from sirnaforge.utils.cache_utils import (
    is_artifact_stamp_current,
    resolve_cache_subdir,
    stable_cache_key,
    write_artifact_stamp,
)

logger = logging.getLogger(__name__)

ENSEMBL_BASE_URL = "https://rest.ensembl.org"

#: Compara relationship types that are orthologues. Paralogues are deliberately excluded: a
#: within-species duplicate is an off-target liability, not conservation evidence.
ORTHOLOGUE_TYPES = frozenset({"ortholog_one2one", "ortholog_one2many", "ortholog_many2many"})

#: Wall-clock ceiling on one whole resolution, across every species and route. A host that DROPs
#: rather than refuses (the common firewall default) fails only at the per-request timeout, so
#: without a ceiling a multi-species screen stalls for minutes before degrading to the labelled
#: heuristic. Species not reached by then are reported unresolved, which is the honest answer.
ORTHOLOGY_BUDGET_SECONDS = 60.0

#: Provenance labels for :meth:`OrthologueMapping.summary`. ``SOURCE_COMPARA_CACHE`` is a distinct
#: label rather than a flag beside ``SOURCE_COMPARA`` because the two are different claims: one says
#: this run asked Ensembl, the other says this run read an answer a previous run obtained.
SOURCE_COMPARA = "ensembl_compara"
SOURCE_COMPARA_CACHE = "ensembl_compara_cache"
SOURCE_MAPPING_FILE = "ortholog_mapping_file"

#: Cache subdir under the shared cache root, and -- by the convention in ``utils.cache_utils`` that
#: manager-backed classes are named by their subdir -- also the producer-version artifact class, so
#: ``PRODUCER_VERSIONS["orthology"]`` invalidates every cached answer exactly once when this
#: resolver is fixed.
ORTHOLOGY_CACHE_SUBDIR = "orthology"

#: Bumped by hand when the cache *document* shape changes. Part of the cache key rather than a
#: migration, so an old document is simply never addressed again instead of being reinterpreted.
ORTHOLOGY_CACHE_SCHEMA = 1

#: Cache TTL. 30 days matches ``AnnotationManager`` rather than the variant resolver's 90: Compara
#: releases move on a roughly two-month cadence, and an orthologue set that changed between releases
#: is a changed answer, not a stale copy of the same one.
ORTHOLOGY_CACHE_TTL_DAYS = 30

#: A user-supplied mapping, keyed on the uppercased version-stripped query gene ID or symbol.
OrthologTable = dict[str, dict[str, frozenset[str]]]

_VERSION_SUFFIX = re.compile(r"\.\d+$")

#: Transport failures a retry cannot fix: no route, no DNS, or a TLS chain this host will never
#: trust. ``ensembl_request_json`` wraps them in ``DatabaseAccessError``, so the cause chain is
#: what gets inspected.
_UNREACHABLE_ERRORS = (aiohttp.ClientConnectorError, aiohttp.ClientSSLError, ssl.SSLError, socket.gaierror)


def _strip_version(identifier: str) -> str:
    """Drop an Ensembl version suffix so IDs compare equal across releases."""
    return _VERSION_SUFFIX.sub("", identifier.strip())


def _query_key(identifier: str) -> str:
    """Normalise a mapping-file query key so gene IDs and symbols both match case-insensitively."""
    return _strip_version(identifier).upper()


def _budget_spent(deadline: float | None) -> bool:
    """True once a resolution's wall-clock budget is gone."""
    return deadline is not None and time.monotonic() >= deadline


def _is_unreachable(error: BaseException | None) -> bool:
    """True when the failure means "no route from this host", so retrying only adds latency.

    Compara's two documented flaky shapes are HTTP *responses* and keep their retries. A refused
    connection, a DNS failure or a certificate that will not verify is a property of the
    environment (air-gapped host, or a TLS-intercepting proxy whose root CA is missing from the
    trust store), and re-asking after 2s and 4s cannot change it.
    """
    while error is not None:
        if isinstance(error, _UNREACHABLE_ERRORS):
            return True
        error = error.__cause__
    return False


@dataclass(frozen=True)
class OrthologueMapping:
    """Resolved orthologues of one or more query genes, grouped by target species.

    Attributes:
        gene_ids_by_species: Canonical species name -> version-stripped orthologue gene IDs.
        resolved_species: Species for which the lookup completed, even if it found nothing.
        unresolved_species: Species the lookup could not complete for (unregistered species or a
            failed request). Distinguished from "resolved, no orthologue" because only the latter
            is evidence of absence.
        queried_gene_ids: The gene IDs actually sent to Compara. Recorded because a resolved-but-empty
            result is otherwise indistinguishable from a lookup on the wrong identifier -- the first
            mouse run resolved 0 orthologues purely because the "gene IDs" were transcript IDs.
        queried_symbols: The gene symbols used for the fallback route, if any.
        source: Where the orthologues came from -- Compara, or a user-supplied mapping file. An
            offline run must not publish provenance that claims a REST call it never made.
        cached_species: Species answered from the on-disk cache rather than from a request made on
            this run. A set, not a boolean, because a partly-cached resolution is the normal case --
            one added species means one request and N-1 hits -- and one flag over the whole mapping
            would be a lie about most of it.
        provenance_by_species: Resolved species -> :data:`SOURCE_COMPARA`,
            :data:`SOURCE_COMPARA_CACHE` or :data:`SOURCE_MAPPING_FILE`. This is the field that
            answers "did this run make a REST call for this species", which ``source`` alone cannot.
            Unresolved species have no entry: there is no answer to attribute.
        cache_key_by_species: Resolved species -> the cache key of the document written or read, so
            an audit can find the exact file behind a conservation claim.
        cached_at_by_species: Cached species -> the ISO timestamp at which the *request* behind the
            served document was made. Only cache hits have an entry; a species resolved on this run
            is dated by the run itself.
    """

    gene_ids_by_species: dict[str, frozenset[str]]
    resolved_species: frozenset[str]
    unresolved_species: frozenset[str]
    queried_gene_ids: frozenset[str] = frozenset()
    queried_symbols: frozenset[str] = frozenset()
    source: str = SOURCE_COMPARA
    cached_species: frozenset[str] = frozenset()
    provenance_by_species: dict[str, str] = field(default_factory=dict)
    cache_key_by_species: dict[str, str] = field(default_factory=dict)
    cached_at_by_species: dict[str, str] = field(default_factory=dict)

    @property
    def all_gene_ids(self) -> frozenset[str]:
        """Every resolved orthologue gene ID, flattened for ``ClassificationContext``."""
        ids: set[str] = set()
        for gene_ids in self.gene_ids_by_species.values():
            ids |= gene_ids
        return frozenset(ids)

    def summary(self) -> dict[str, Any]:
        """Provenance record for the run summary, so a conservation claim is auditable.

        The cache block is published per species, not per run: with an on-disk cache the run record
        has to state, for each species separately, whether *this* run asked Ensembl. ``source``
        keeps its original meaning (mapping file vs Compara) so existing readers are unaffected.
        """
        return {
            "source": self.source,
            "orthologue_types": sorted(ORTHOLOGUE_TYPES),
            "gene_ids_by_species": {s: sorted(g) for s, g in sorted(self.gene_ids_by_species.items())},
            "resolved_species": sorted(self.resolved_species),
            "unresolved_species": sorted(self.unresolved_species),
            "queried_gene_ids": sorted(self.queried_gene_ids),
            "queried_symbols": sorted(self.queried_symbols),
            "cached_species": sorted(self.cached_species),
            "provenance_by_species": dict(sorted(self.provenance_by_species.items())),
            "cache_keys": dict(sorted(self.cache_key_by_species.items())),
            "cached_at": dict(sorted(self.cached_at_by_species.items())),
        }

    @classmethod
    def empty(cls) -> OrthologueMapping:
        """A mapping that resolved nothing, so every ortholog verdict falls back to the heuristic."""
        return cls(gene_ids_by_species={}, resolved_species=frozenset(), unresolved_species=frozenset())

    @classmethod
    def from_file(
        cls,
        path: Path | str,
        query_gene_ids: frozenset[str] | set[str] | list[str],
        query_species: str,
        target_species: frozenset[str] | set[str] | list[str],
        *,
        query_gene_symbols: frozenset[str] | set[str] | list[str] = frozenset(),
    ) -> OrthologueMapping:
        """Resolve from a user-supplied mapping file instead of Compara. Never touches the network.

        The offline half of issue #101: an air-gapped run, and every deterministic fixture, states
        its orthologues in a file rather than depending on REST being reachable. Same arguments as
        :func:`resolve_orthologues` so the two are interchangeable at the call site.

        Deliberately never reaches the on-disk cache, in either direction. Reading a file is already
        offline and free, so a cache would save nothing; and *writing* one would let a document
        outlive the mapping file it was copied from, so an edited file would be silently ignored for
        the TTL -- exactly the failure the explicit-input-fails-loudly rule in
        :func:`load_ortholog_table` exists to prevent.
        """
        return mapping_from_table(
            load_ortholog_table(path),
            query_gene_ids,
            query_species,
            target_species,
            query_gene_symbols=query_gene_symbols,
        )


def load_ortholog_table(path: Path | str) -> OrthologTable:
    """Parse a mapping file of query gene -> species -> orthologue gene IDs.

    JSON, e.g. ``{"ENSG00000141510": {"mouse": ["ENSMUSG00000059552"], "rat": []}}``. A query key
    may be a gene ID (version optional) or a gene symbol, because an input-FASTA run often has no
    stable gene ID to key on. An explicitly empty list is a *claim* -- "checked, no orthologue" --
    and is kept distinct from a species the file says nothing about.

    Raises:
        FileNotFoundError: The path does not exist.
        ValueError: The file is not JSON, or not of that shape. Unlike a failed network lookup,
            which degrades to the labelled heuristic, a malformed explicit input fails loudly:
            silently ignoring it would screen with weaker evidence than the user asked for.
    """
    mapping_path = Path(path)
    try:
        raw = json.loads(mapping_path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Ortholog mapping file {mapping_path} is not valid JSON: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"Ortholog mapping file {mapping_path} must be an object of query gene -> species -> IDs")

    table: OrthologTable = {}
    for query_key, per_species in raw.items():
        if not isinstance(per_species, dict):
            raise ValueError(f"Ortholog mapping file {mapping_path}: entry {query_key!r} must map species -> gene IDs")
        entry = table.setdefault(_query_key(str(query_key)), {})
        for species, gene_ids in per_species.items():
            if isinstance(gene_ids, str) or not isinstance(gene_ids, (list, tuple, set, frozenset)):
                raise ValueError(
                    f"Ortholog mapping file {mapping_path}: {query_key!r}/{species!r} must be a list of gene IDs"
                )
            canonical = normalize_species_name(str(species))
            resolved = frozenset(_strip_version(str(gene_id)) for gene_id in gene_ids if gene_id)
            entry[canonical] = entry.get(canonical, frozenset()) | resolved
    return table


def mapping_from_table(
    table: OrthologTable,
    query_gene_ids: frozenset[str] | set[str] | list[str],
    query_species: str,
    target_species: frozenset[str] | set[str] | list[str],
    *,
    query_gene_symbols: frozenset[str] | set[str] | list[str] = frozenset(),
) -> OrthologueMapping:
    """Select one query gene's orthologues out of a parsed mapping table.

    A requested species the table says nothing about is ``unresolved``, not resolved-empty: the file
    made no claim about it, so its hits fall back to the labelled symbol heuristic rather than being
    reported as "checked, not an orthologue".
    """
    canonical_query = normalize_species_name(query_species)
    wanted = {normalize_species_name(s) for s in target_species} - {canonical_query}
    genes = {_strip_version(g) for g in query_gene_ids if g}
    symbols = {s.strip() for s in query_gene_symbols if s and s.strip()}
    if (not genes and not symbols) or not wanted:
        return OrthologueMapping.empty()

    stated: dict[str, set[str]] = {}
    for key in {_query_key(identifier) for identifier in genes | symbols}:
        for species, gene_ids in table.get(key, {}).items():
            stated.setdefault(species, set()).update(gene_ids)

    resolved = wanted & set(stated)
    unresolved = wanted - resolved
    if unresolved:
        logger.warning(
            "Ortholog mapping file states nothing for %s; hits in those species fall back to the "
            "labelled symbol heuristic",
            sorted(unresolved),
        )
    return OrthologueMapping(
        gene_ids_by_species={s: frozenset(stated[s]) for s in sorted(resolved) if stated[s]},
        resolved_species=frozenset(resolved),
        unresolved_species=frozenset(unresolved),
        queried_gene_ids=frozenset(genes),
        queried_symbols=frozenset(symbols),
        source=SOURCE_MAPPING_FILE,
        # Per-species provenance is filled in here too, so the field is a complete record for every
        # path rather than one that reads as "unknown" whenever the offline route was taken.
        provenance_by_species=dict.fromkeys(sorted(resolved), SOURCE_MAPPING_FILE),
    )


def orthology_cache_key(
    *,
    query_species: str,
    target_species: str,
    gene_ids: frozenset[str] | set[str],
    symbols: frozenset[str] | set[str],
    base_url: str,
) -> str:
    """Cache key for one (question x *single* target species) Compara answer.

    One species per key, deliberately: a five-species screen that resolves four and fails one must
    record four hits and re-ask only the fifth, which a whole-resolution key could not express.

    Everything that changes what the answer *means* is in the payload. ``target_species`` is in it
    because the mouse document must never be served for rat. ``symbols`` is in it because the symbol
    fallback route is part of the question asked, not an implementation detail: the same gene IDs
    with and without a symbol are two different lookups (#101 -- an input FASTA supplies transcript
    IDs, which only the symbol route rescues). ``base_url`` is in it because a staging or mirrored
    REST host is a different authority. ``orthologue_types`` is in it because widening the accepted
    relationship types must miss every answer computed under the narrower set.

    Deliberately excluded: ``timeout`` and ``budget``, which change how long we were willing to wait
    rather than what the answer means, and the mapping-file path, which is never cached at all.
    """
    return stable_cache_key(
        {
            "schema": ORTHOLOGY_CACHE_SCHEMA,
            "query_species": query_species,
            "target_species": target_species,
            "gene_ids": sorted(gene_ids),
            "symbols": sorted(symbols),
            "base_url": base_url,
            "orthologue_types": sorted(ORTHOLOGUE_TYPES),
        }
    )


def _orthology_cache_root(cache_dir: Path | None) -> Path | None:
    """Resolve the cache directory, or None when it cannot be made writable.

    Uses ``resolve_cache_subdir`` so the orthology cache honours ``SIRNAFORGE_CACHE_DIR`` / XDG /
    ``$HOME/.cache`` / workspace / temp exactly like every other subsystem, instead of inventing a
    second location. A cache root we cannot create is a missing optimisation, never an error.
    """
    try:
        if cache_dir is not None:
            return resolve_cache_subdir(ORTHOLOGY_CACHE_SUBDIR, override=cache_dir)
        return resolve_cache_subdir(ORTHOLOGY_CACHE_SUBDIR)
    except (OSError, RuntimeError) as exc:
        logger.debug("Orthology cache disabled: no writable cache directory (%s)", exc)
        return None


def _read_cached_orthologues(cache_root: Path, cache_key: str) -> tuple[frozenset[str], str] | None:
    """Read one cached answer, or None on any kind of miss.

    Validity is decided entirely by the sidecar stamp -- producer version, TTL and a fingerprint of
    the document's own bytes -- so no second invalidation scheme exists here to disagree with it. The
    output fingerprint is what stops a half-written document being served for its whole TTL.

    Exception-swallowing in the same way as the rest of this module: an unreadable or unparseable
    document is a miss with a debug line, never a raise, because a broken cache must not fail a
    screen. It only costs the request the cache was meant to save.
    """
    document = cache_root / f"{cache_key}.json"
    if not document.exists():
        return None
    if not is_artifact_stamp_current(ORTHOLOGY_CACHE_SUBDIR, document, max_age_days=ORTHOLOGY_CACHE_TTL_DAYS):
        return None
    try:
        payload = json.loads(document.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.debug("Unreadable orthology cache document %s: %s", document, exc)
        return None
    if not isinstance(payload, dict) or payload.get("schema") != ORTHOLOGY_CACHE_SCHEMA:
        logger.debug("Orthology cache document %s has schema %r, not %r", document, payload, ORTHOLOGY_CACHE_SCHEMA)
        return None
    gene_ids = payload.get("orthologue_gene_ids")
    if not isinstance(gene_ids, list):
        logger.debug("Orthology cache document %s carries no orthologue_gene_ids list", document)
        return None
    return frozenset(_strip_version(str(gene_id)) for gene_id in gene_ids if gene_id), str(
        payload.get("resolved_at") or ""
    )


def _write_cached_orthologues(
    cache_root: Path,
    cache_key: str,
    *,
    query_species: str,
    target_species: str,
    gene_ids: frozenset[str] | set[str],
    symbols: frozenset[str] | set[str],
    base_url: str,
    orthologue_gene_ids: set[str],
) -> None:
    """Record one completed Compara answer, then stamp it.

    The question is written into the document alongside the answer even though the key already
    covers it: a bare digest is unreadable, and the point of this cache is that a cached
    conservation claim stays auditable. Stamped *after* the bytes are on disk, so the recorded size
    and digest describe what a later run will read back.
    """
    document = cache_root / f"{cache_key}.json"
    payload = {
        "schema": ORTHOLOGY_CACHE_SCHEMA,
        "cache_key": cache_key,
        "query_species": query_species,
        "target_species": target_species,
        "gene_ids": sorted(gene_ids),
        "symbols": sorted(symbols),
        "base_url": base_url,
        "orthologue_types": sorted(ORTHOLOGUE_TYPES),
        "orthologue_gene_ids": sorted(orthologue_gene_ids),
        "resolved_at": datetime.now().isoformat(),
    }
    try:
        document.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        write_artifact_stamp(
            ORTHOLOGY_CACHE_SUBDIR,
            document,
            extra={"query_species": query_species, "target_species": target_species},
        )
    except OSError as exc:
        logger.debug("Could not cache the %s orthologue answer: %s", target_species, exc)


async def resolve_orthologues(
    query_gene_ids: frozenset[str] | set[str] | list[str],
    query_species: str,
    target_species: frozenset[str] | set[str] | list[str],
    *,
    query_gene_symbols: frozenset[str] | set[str] | list[str] = frozenset(),
    base_url: str = ENSEMBL_BASE_URL,
    timeout: int = 30,
    budget: float | None = None,
    session: aiohttp.ClientSession | None = None,
    cache: bool = True,
    cache_dir: Path | None = None,
) -> OrthologueMapping:
    """Resolve orthologues of the query gene(s) in each target species.

    Two routes, because a run driven by ``input_fasta`` often has no stable gene ID at all. The
    pinned TP53 FASTA in the full-scale kit has headers like ``>ENST00000413465 TP53 ...``, so the
    "gene ID" reaching this function is a *transcript* ID, and
    ``/homology/id/homo_sapiens/ENST00000413465`` answers **HTTP 200 with {"data": []}** -- a
    successful empty result, indistinguishable from "no orthologue exists". So when the ID route
    finds nothing, the symbol route ``/homology/symbol/{species}/{symbol}`` is tried; it resolves
    ``TP53`` to ``ENSMUSG00000059552`` (both verified live on release 116).

    The symbol is only an *input* to the lookup here. What comes back is a stable gene ID from
    Compara, so a verdict reached this way is still ``OrthologEvidence.GENE_ID`` -- unlike the
    symbol-equality heuristic in ``classify_hit``, nothing is being asserted from symbol equality.

    Never raises for a lookup failure: a species that could not be resolved is reported in
    ``unresolved_species`` and its hits fall back to the symbol heuristic. Failing the whole screen
    because Compara was briefly unavailable would be worse than screening with weaker evidence, as
    long as the weaker evidence is labelled -- which ``OrthologEvidence`` does.

    May reach the network. Use :meth:`OrthologueMapping.from_file` for an air-gapped run or any
    deterministic fixture; a unit test must never end up here. A warm cache makes the request
    optional per species, and the returned mapping states per species whether one was made -- so the
    provenance published by an offline run with a warm cache is ``ensembl_compara_cache``, never a
    REST call it did not perform.

    Args:
        query_gene_ids: Stable gene IDs of the query gene (version suffix optional).
        query_species: Species the query gene belongs to.
        target_species: Species to resolve orthologues in; the query species is skipped.
        query_gene_symbols: Gene symbols to fall back to when the ID route resolves nothing.
        base_url: Ensembl REST base URL.
        timeout: Per-request timeout in seconds.
        budget: Optional wall-clock ceiling on the whole resolution; species not reached in time are
            reported unresolved. Defaults to no ceiling, because an abandoned species is
            indistinguishable from an absent orthologue -- see the module docstring.
        session: Session to reuse; one is opened for this call when omitted.
        cache: Read and write the on-disk cache. On by default; ``False`` for a test that must not
            touch the disk, and for a caller that wants to prove Compara answers this question now.
        cache_dir: Cache directory override. Defaults to the shared cache root's ``orthology``
            subdir.

    Returns:
        An :class:`OrthologueMapping`. Empty when there is nothing to resolve.
    """
    canonical_query = normalize_species_name(query_species)
    wanted = {normalize_species_name(s) for s in target_species} - {canonical_query}
    genes = {_strip_version(g) for g in query_gene_ids if g}
    symbols = {s.strip() for s in query_gene_symbols if s and s.strip()}
    if (not genes and not symbols) or not wanted:
        return OrthologueMapping.empty()

    query_slug = ensembl_species_slug(canonical_query)
    if query_slug is None:
        logger.warning("Orthology lookup skipped: %r is not a registered species", query_species)
        return OrthologueMapping(
            gene_ids_by_species={}, resolved_species=frozenset(), unresolved_species=frozenset(wanted)
        )

    by_species: dict[str, set[str]] = {}
    resolved: set[str] = set()
    unresolved: set[str] = set()
    cached: set[str] = set()
    provenance: dict[str, str] = {}
    cache_keys: dict[str, str] = {}
    cached_at: dict[str, str] = {}
    deadline = time.monotonic() + budget if budget else None
    cache_root = _orthology_cache_root(cache_dir) if cache else None

    async with ensembl_session(session, timeout) as active:
        for species in sorted(wanted):
            target_slug = ensembl_species_slug(species)
            if target_slug is None:
                logger.warning("Orthology lookup skipped for unregistered target species %r", species)
                unresolved.add(species)
                continue

            cache_key = (
                orthology_cache_key(
                    query_species=canonical_query,
                    target_species=species,
                    gene_ids=genes,
                    symbols=symbols,
                    base_url=base_url,
                )
                if cache_root is not None
                else None
            )
            # Read before the budget check: serving a cached answer costs no wall clock, so a spent
            # budget must not turn a species we already know about into an unresolved one.
            if cache_root is not None and cache_key is not None:
                hit = _read_cached_orthologues(cache_root, cache_key)
                if hit is not None:
                    cached_ids, resolved_at = hit
                    resolved.add(species)
                    cached.add(species)
                    provenance[species] = SOURCE_COMPARA_CACHE
                    cache_keys[species] = cache_key
                    if resolved_at:
                        cached_at[species] = resolved_at
                    if cached_ids:
                        by_species[species] = set(cached_ids)
                    logger.info(
                        "Served %d %s orthologue gene ID(s) from the on-disk cache; no Compara request was made",
                        len(cached_ids),
                        species,
                    )
                    continue

            if _budget_spent(deadline):
                logger.warning(
                    "Orthology budget of %.0fs is spent; %r is left unresolved rather than stalling the screen",
                    budget,
                    species,
                )
                unresolved.add(species)
                continue

            found, failed = await _lookup_species(
                active,
                genes=genes,
                symbols=symbols,
                query_slug=query_slug,
                target_slug=target_slug,
                base_url=base_url,
                timeout=timeout,
                species=species,
                deadline=deadline,
            )

            if failed and not found:
                unresolved.add(species)
                continue
            resolved.add(species)
            provenance[species] = SOURCE_COMPARA
            if found:
                by_species[species] = found
            # Cached only when *nothing* failed. A resolved absence is cached because it is a real
            # Compara answer, expensive to re-obtain and stable between releases. A failed or
            # partially failed lookup never is: writing one would turn a firewall blip or a single
            # dead identifier into a 30-day claim about what exists, which is the one claim this
            # module makes (#101).
            if cache_root is not None and cache_key is not None and not failed:
                _write_cached_orthologues(
                    cache_root,
                    cache_key,
                    query_species=canonical_query,
                    target_species=species,
                    gene_ids=genes,
                    symbols=symbols,
                    base_url=base_url,
                    orthologue_gene_ids=found,
                )
                cache_keys[species] = cache_key
            logger.info("Resolved %d %s orthologue gene ID(s) for the query gene", len(found), species)

    return OrthologueMapping(
        gene_ids_by_species={s: frozenset(g) for s, g in by_species.items()},
        resolved_species=frozenset(resolved),
        unresolved_species=frozenset(unresolved),
        queried_gene_ids=frozenset(genes),
        queried_symbols=frozenset(symbols),
        cached_species=frozenset(cached),
        provenance_by_species=provenance,
        cache_key_by_species=cache_keys,
        cached_at_by_species=cached_at,
    )


async def _lookup_species(
    session: aiohttp.ClientSession,
    *,
    genes: set[str],
    symbols: set[str],
    query_slug: str,
    target_slug: str,
    base_url: str,
    timeout: int,
    species: str,
    deadline: float | None,
) -> tuple[set[str], bool]:
    """Both Compara routes for one target species: (orthologue gene IDs, any failure).

    Gene IDs first, symbols only when those resolved nothing -- the cheap path stays cheap, and the
    fallback exists because Compara answers a *transcript* ID with a successful empty 200 (#101).
    Extracted from :func:`resolve_orthologues` so the resolution loop reads as one decision per
    species rather than interleaving route selection with cache and budget bookkeeping.
    """
    found, failed = await _lookup_route(
        session, "id", sorted(genes), query_slug, target_slug, base_url, timeout, species, deadline
    )
    if not found and symbols:
        logger.info(
            "No %s orthologue from gene IDs %s; retrying on symbol(s) %s "
            "(an input FASTA supplies transcript IDs, which Compara answers with an empty 200)",
            species,
            sorted(genes) or "<none>",
            sorted(symbols),
        )
        symbol_found, symbol_failed = await _lookup_route(
            session, "symbol", sorted(symbols), query_slug, target_slug, base_url, timeout, species, deadline
        )
        found |= symbol_found
        failed = failed or symbol_failed
    return found, failed


async def _lookup_route(
    session: aiohttp.ClientSession,
    route: str,
    identifiers: list[str],
    query_slug: str,
    target_slug: str,
    base_url: str,
    timeout: int,
    species: str,
    deadline: float | None = None,
) -> tuple[set[str], bool]:
    """Query one Compara route for every identifier, returning (orthologue gene IDs, any failure).

    ``route`` is ``id`` or ``symbol`` -- the two ``/homology/{route}/{species}/{identifier}`` forms.
    ``target_species`` must be the Ensembl slug: ``mouse`` is rejected, ``mus_musculus`` works.
    ``;`` and ``&`` are both accepted as separators (verified live); ``&`` is the conventional form.
    """
    found: set[str] = set()
    failed = False
    for identifier in identifiers:
        if _budget_spent(deadline):
            logger.warning("Orthology budget spent before the %s route reached %s -> %s", route, identifier, species)
            failed = True
            continue
        url = (
            f"{base_url}/homology/{route}/{query_slug}/{identifier}"
            f"?target_species={target_slug}&type=orthologues&format=condensed"
        )
        try:
            payload = await _request_homologies(session, url, timeout, deadline=deadline)
        except Exception as exc:  # noqa: BLE001 - degrade to the heuristic, never fail the screen
            logger.warning("Orthology %s lookup failed for %s -> %s: %s", route, identifier, species, exc)
            failed = True
            continue
        if payload is None:
            failed = True
            continue
        found |= _orthologue_ids(payload)
    return found, failed


async def _request_homologies(
    session: aiohttp.ClientSession,
    url: str,
    timeout: int,
    attempts: int = ENSEMBL_MAX_ATTEMPTS,
    *,
    deadline: float | None = None,
) -> Any | None:
    """One homology request, retrying Compara's two transient failure shapes.

    This endpoint intermittently fails a perfectly valid query, and does so in two ways that the
    shared status-based retry in ``ensembl_request_json`` does not cover. Both were observed live
    against release 116 on a URL that succeeded on the immediately following attempt:

    1. HTTP **200** with ``{"error": "DBD::mysql::st bind_param failed: Illegal parameter number"}``.
       A 200 is never retried, and a naive parser reads the absent ``data`` key as "this gene has no
       orthologue in that species" -- a silent false negative on the one claim this module makes.
    2. HTTP **400** with the same body. 400 is not in ``ENSEMBL_RETRY_STATUSES`` because it normally
       means a genuinely bad request, so ``ensembl_request_json`` raises immediately.

    Measured in-container, unretried resolution succeeded 5 times in 6; both shapes are retried here
    so a flaky lookup does not become a confident wrong answer. Retrying a genuinely malformed
    request costs a bounded ``attempts`` round trips and still ends in the same reported failure.

    Two failures are *not* retried, because the backoff can only make a certain outcome slower: a
    transport-level "no route from this host" (see :func:`_is_unreachable`), and a spent ``deadline``.
    On a TLS-intercepting proxy whose root CA is missing, the three attempts and their 2s + 4s sleeps
    turned an instant certificate rejection into a 6.3s per-route stall.

    Returns:
        The decoded payload, or None when every attempt failed.

    Raises:
        Exception: The last transport error, when every attempt raised.
    """
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        reason: str
        try:
            payload = await ensembl_request_json(
                session, "GET", url, headers={"Accept": "application/json"}, retry_cap=timeout
            )
            if not (isinstance(payload, dict) and payload.get("error")):
                return payload
            reason = f"HTTP 200 with an error body: {payload['error']}"
        except Exception as exc:  # noqa: BLE001 - both shapes are retried; see the docstring
            last_error = exc
            reason = f"{type(exc).__name__}: {exc}"

        logger.warning("Compara lookup attempt %d/%d failed (%s) for %s", attempt, attempts, reason, url)
        if attempt >= attempts:
            break
        if _is_unreachable(last_error):
            logger.warning("Compara is unreachable from this host (%s); not retrying %s", reason, url)
            break
        if _budget_spent(deadline):
            logger.warning("Orthology budget spent after attempt %d; not retrying %s", attempt, url)
            break
        await asyncio.sleep(min(2.0 * attempt, float(timeout)))

    if last_error is not None:
        raise last_error
    return None


def _orthologue_ids(payload: Any) -> set[str]:
    """Extract version-stripped orthologue gene IDs from a condensed homology response.

    The condensed shape is ``{"data": [{"id": ..., "homologies": [{"id": ..., "type": ...}]}]}``.
    Anything unexpected yields nothing rather than raising: a shape change should weaken evidence,
    not crash a screen.
    """
    ids: set[str] = set()
    data = payload.get("data") if isinstance(payload, dict) else None
    for record in data if isinstance(data, list) else []:
        homologies = record.get("homologies") if isinstance(record, dict) else None
        for homology in homologies if isinstance(homologies, list) else []:
            if not isinstance(homology, dict):
                continue
            if str(homology.get("type", "")) not in ORTHOLOGUE_TYPES:
                continue
            gene_id = homology.get("id")
            if gene_id:
                ids.add(_strip_version(str(gene_id)))
    return ids
