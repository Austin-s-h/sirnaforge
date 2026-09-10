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
There is no cache, so a caller that resolves twice pays twice. Two guards keep an unreachable
Compara off the critical path: :data:`ORTHOLOGY_BUDGET_SECONDS` bounds a whole resolution, and
:meth:`OrthologueMapping.from_file` resolves from a user-supplied mapping with no network at all
(issue #101: an offline path is required, not optional).
"""

from __future__ import annotations

import asyncio
import json
import logging
import re
import socket
import ssl
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import aiohttp

from sirnaforge.data.base import ENSEMBL_MAX_ATTEMPTS, ensembl_request_json, ensembl_session
from sirnaforge.data.species_registry import ensembl_species_slug, normalize_species_name

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

#: Provenance labels for :meth:`OrthologueMapping.summary`.
SOURCE_COMPARA = "ensembl_compara"
SOURCE_MAPPING_FILE = "ortholog_mapping_file"

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
    """

    gene_ids_by_species: dict[str, frozenset[str]]
    resolved_species: frozenset[str]
    unresolved_species: frozenset[str]
    queried_gene_ids: frozenset[str] = frozenset()
    queried_symbols: frozenset[str] = frozenset()
    source: str = SOURCE_COMPARA

    @property
    def all_gene_ids(self) -> frozenset[str]:
        """Every resolved orthologue gene ID, flattened for ``ClassificationContext``."""
        ids: set[str] = set()
        for gene_ids in self.gene_ids_by_species.values():
            ids |= gene_ids
        return frozenset(ids)

    def summary(self) -> dict[str, Any]:
        """Provenance record for the run summary, so a conservation claim is auditable."""
        return {
            "source": self.source,
            "orthologue_types": sorted(ORTHOLOGUE_TYPES),
            "gene_ids_by_species": {s: sorted(g) for s, g in sorted(self.gene_ids_by_species.items())},
            "resolved_species": sorted(self.resolved_species),
            "unresolved_species": sorted(self.unresolved_species),
            "queried_gene_ids": sorted(self.queried_gene_ids),
            "queried_symbols": sorted(self.queried_symbols),
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
    )


async def resolve_orthologues(
    query_gene_ids: frozenset[str] | set[str] | list[str],
    query_species: str,
    target_species: frozenset[str] | set[str] | list[str],
    *,
    query_gene_symbols: frozenset[str] | set[str] | list[str] = frozenset(),
    base_url: str = ENSEMBL_BASE_URL,
    timeout: int = 30,
    budget: float | None = ORTHOLOGY_BUDGET_SECONDS,
    session: aiohttp.ClientSession | None = None,
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

    Reaches the network. Use :meth:`OrthologueMapping.from_file` for an air-gapped run or any
    deterministic fixture; a unit test must never end up here.

    Args:
        query_gene_ids: Stable gene IDs of the query gene (version suffix optional).
        query_species: Species the query gene belongs to.
        target_species: Species to resolve orthologues in; the query species is skipped.
        query_gene_symbols: Gene symbols to fall back to when the ID route resolves nothing.
        base_url: Ensembl REST base URL.
        timeout: Per-request timeout in seconds.
        budget: Wall-clock ceiling on the whole resolution; species not reached in time are
            reported unresolved. None removes the ceiling.
        session: Session to reuse; one is opened for this call when omitted.

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
    deadline = time.monotonic() + budget if budget else None

    async with ensembl_session(session, timeout) as active:
        for species in sorted(wanted):
            target_slug = ensembl_species_slug(species)
            if target_slug is None:
                logger.warning("Orthology lookup skipped for unregistered target species %r", species)
                unresolved.add(species)
                continue
            if _budget_spent(deadline):
                logger.warning(
                    "Orthology budget of %.0fs is spent; %r is left unresolved rather than stalling the screen",
                    budget,
                    species,
                )
                unresolved.add(species)
                continue

            # ID route first; the symbol route only if it resolved nothing.
            found, failed = await _lookup_route(
                active, "id", sorted(genes), query_slug, target_slug, base_url, timeout, species, deadline
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
                    active, "symbol", sorted(symbols), query_slug, target_slug, base_url, timeout, species, deadline
                )
                found |= symbol_found
                failed = failed or symbol_failed

            if failed and not found:
                unresolved.add(species)
                continue
            resolved.add(species)
            if found:
                by_species[species] = found
            logger.info("Resolved %d %s orthologue gene ID(s) for the query gene", len(found), species)

    return OrthologueMapping(
        gene_ids_by_species={s: frozenset(g) for s, g in by_species.items()},
        resolved_species=frozenset(resolved),
        unresolved_species=frozenset(unresolved),
        queried_gene_ids=frozenset(genes),
        queried_symbols=frozenset(symbols),
    )


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
