"""The on-disk orthologue cache, and the provenance it is obliged to publish (#101).

``data/orthology.py`` used to say "There is no cache, so a caller that resolves twice pays twice",
and the cost it named is up to two Compara requests per (query gene x target species). This module
pins the cache that replaced that sentence, and every one of its claims is a claim about *honesty*
rather than about speed:

* a served answer is labelled ``ensembl_compara_cache``, so an offline run cannot report a REST call
  it never made;
* one document per single target species, so the mouse answer is never served for rat;
* the producer version invalidates every answer exactly once, which is the axis that matters for a
  resolver that has already shipped a wrong answer class (transcript IDs sent to ``/homology/id/``,
  which Compara answers with a successful empty 200);
* a resolved absence is cached because it is a real Compara answer; a *failed* lookup is not,
  because caching a transport failure would turn a firewall blip into a 30-day claim of absence.

No network: every test drives a monkeypatched ``ensembl_request_json`` and supplies its own
``cache_dir``, so nothing reaches Ensembl and nothing is written to the user's cache root.
"""

import json
import logging
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.data import orthology
from sirnaforge.data.base import DatabaseAccessError
from sirnaforge.data.orthology import (
    ORTHOLOGY_CACHE_SUBDIR,
    SOURCE_COMPARA,
    SOURCE_COMPARA_CACHE,
    OrthologueMapping,
    resolve_orthologues,
)
from sirnaforge.utils import cache_utils
from sirnaforge.utils.cache_utils import artifact_stamp_path

# Public reference genes, used deliberately as a public baseline: TP53 and its rodent orthologues are
# the same worked example the sibling module's tests and docstrings use.
HUMAN_TP53 = "ENSG00000141510"
MOUSE_TRP53 = "ENSMUSG00000059552"
RAT_TP53 = "ENSRNOG00000010756"

#: Slug -> the orthologue that slug's Compara answer carries, so a document served for the wrong
#: species is visible in the gene IDs and not only in the request count.
ORTHOLOGUE_BY_SLUG = {"mus_musculus": MOUSE_TRP53, "rattus_norvegicus": RAT_TP53}


def _condensed(gene_id: str, species: str) -> dict[str, Any]:
    """A condensed homology payload carrying one one2one orthologue."""
    return {
        "data": [
            {
                "id": HUMAN_TP53,
                "homologies": [{"id": gene_id, "species": species, "type": "ortholog_one2one"}],
            }
        ]
    }


async def _no_sleep(_seconds: float) -> None:
    """Collapse retry backoff so the failure tests stay fast."""


@pytest.fixture
def compara(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Patch the Ensembl request helper; the returned list is every URL this run requested.

    Answers per target species, so "how many requests" and "which answer" are two independent
    observations. A test that expects a cache hit asserts on the list *not* growing.
    """
    urls: list[str] = []

    async def fake_request(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        urls.append(url)
        for slug, gene_id in ORTHOLOGUE_BY_SLUG.items():
            if f"target_species={slug}" in url:
                return _condensed(gene_id, slug)
        return {"data": []}

    monkeypatch.setattr(orthology, "ensembl_request_json", fake_request)
    return urls


def _documents(cache_dir: Path) -> list[Path]:
    """Cached answer documents, excluding their sidecar stamps."""
    return sorted(p for p in cache_dir.glob("*.json") if not p.name.endswith(cache_utils.ARTIFACT_STAMP_SUFFIX))


async def _resolve(
    cache_dir: Path, species: frozenset[str] | set[str] = frozenset({"mouse"}), **kwargs: Any
) -> OrthologueMapping:
    """Resolve human TP53 against `species`, always through the supplied cache directory."""
    return await resolve_orthologues({HUMAN_TP53}, "human", set(species), cache_dir=cache_dir, **kwargs)


@pytest.mark.unit
@pytest.mark.asyncio
async def test_a_cached_answer_says_it_was_cached(tmp_path: Path, compara: list[str]) -> None:
    """A cache hit must be published as a cache hit, not as a Compara request.

    This is the whole reason ``provenance_by_species`` exists next to ``source``: the answer is the
    same either way, so without a per-species label an offline run behind a warm cache would report
    ``ensembl_compara`` and claim a REST call it never made. The requirement is provenance, not
    speed, which is why the assertion is on the label and not only on ``len(compara)``.
    """
    first = await _resolve(tmp_path)

    assert len(compara) == 1, "the first resolution pays for the request"
    assert first.cached_species == frozenset(), "nothing could have been cached yet"
    assert first.provenance_by_species == {"mouse": SOURCE_COMPARA}
    assert first.cache_key_by_species["mouse"], "the key of the document just written is published"

    second = await _resolve(tmp_path)

    assert len(compara) == 1, "the second resolution asks Ensembl nothing"
    assert second.gene_ids_by_species == first.gene_ids_by_species == {"mouse": frozenset({MOUSE_TRP53})}
    assert second.resolved_species == frozenset({"mouse"})
    assert second.cached_species == frozenset({"mouse"})
    assert second.provenance_by_species == {"mouse": SOURCE_COMPARA_CACHE}
    assert second.cache_key_by_species == first.cache_key_by_species, "the same document, named"

    summary = second.summary()
    assert summary["cached_species"] == ["mouse"]
    assert summary["provenance_by_species"] == {"mouse": SOURCE_COMPARA_CACHE}
    assert summary["cache_keys"] == {"mouse": first.cache_key_by_species["mouse"]}
    assert summary["cached_at"]["mouse"], "the served answer is dated, so an audit can place it"
    # `source` keeps its original meaning -- mapping file vs Compara -- so existing readers are
    # unaffected by the cache; the per-species block is where "this run" is stated.
    assert summary["source"] == SOURCE_COMPARA


@pytest.mark.unit
@pytest.mark.asyncio
async def test_one_species_hit_does_not_answer_for_another(tmp_path: Path, compara: list[str]) -> None:
    """``target_species`` is part of the cache key, so the mouse document cannot answer for rat.

    A key over the question alone -- gene IDs, symbols, host -- collides across species, and the
    collision is silent in the worst possible way: rat comes back "resolved" carrying the *mouse*
    orthologue gene ID, which is then used as conservation evidence for a rat hit.
    """
    await _resolve(tmp_path, {"mouse"})
    assert len(compara) == 1

    both = await _resolve(tmp_path, {"mouse", "rat"})

    assert both.cached_species == frozenset({"mouse"}), "only the species already asked about is cached"
    assert both.provenance_by_species == {"mouse": SOURCE_COMPARA_CACHE, "rat": SOURCE_COMPARA}
    assert len(compara) == 2, "rat costs exactly one request; mouse costs none"
    assert "target_species=rattus_norvegicus" in compara[1]
    assert both.gene_ids_by_species == {"mouse": frozenset({MOUSE_TRP53}), "rat": frozenset({RAT_TP53})}
    assert both.cache_key_by_species["mouse"] != both.cache_key_by_species["rat"]
    assert len(_documents(tmp_path)) == 2, "two species, two documents"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_bumping_the_producer_version_invalidates_every_answer_once(
    tmp_path: Path, compara: list[str], caplog: pytest.LogCaptureFixture
) -> None:
    """The axis that matters here: a fixed resolver must discard the answers the broken one cached.

    This resolver has already shipped a wrong answer class -- transcript IDs sent to
    ``/homology/id/``, which Compara answers with a successful empty 200 -- so "no mouse orthologue"
    can be a bug rather than biology. TTL and byte fingerprints cannot see that, because the document
    is fresh and intact; only the recorded producer version can. And it must invalidate *once*: a
    permanent miss would be a cache that never caches.
    """
    await _resolve(tmp_path)
    assert len(compara) == 1

    assert ORTHOLOGY_CACHE_SUBDIR in cache_utils.PRODUCER_VERSIONS, "registering the class is what opts it in"

    with pytest.MonkeyPatch.context() as patched:
        patched.setitem(cache_utils.PRODUCER_VERSIONS, "orthology", "99.0")

        with caplog.at_level(logging.WARNING):
            refetched = await _resolve(tmp_path)

        assert len(compara) == 2, "an answer from the superseded producer must be re-asked"
        assert refetched.cached_species == frozenset(), "and reported as a live request, not a hit"
        assert "Discarding stale cached orthology artifact" in caplog.text

        again = await _resolve(tmp_path)
        assert len(compara) == 2, "exactly once: the re-fetched answer is cached under the new version"
        assert again.cached_species == frozenset({"mouse"})


@pytest.mark.unit
@pytest.mark.asyncio
async def test_a_failed_lookup_is_not_cached_as_absence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A transport failure must never become a 30-day claim that no orthologue exists.

    ``unresolved`` and "resolved, no orthologue" are different claims -- the module keeps them in
    separate fields precisely because only the second is evidence of absence. Caching the first would
    collapse them on disk, so a single firewall blip would silently downgrade conservation evidence
    for every later run until the TTL expired.
    """
    attempts: list[str] = []

    async def refused(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        attempts.append(url)
        raise DatabaseAccessError("503 Service Unavailable", "Ensembl")

    monkeypatch.setattr(orthology, "ensembl_request_json", refused)
    monkeypatch.setattr(orthology.asyncio, "sleep", _no_sleep)

    first = await _resolve(tmp_path)
    assert first.unresolved_species == frozenset({"mouse"})
    assert _documents(tmp_path) == [], "a failure writes no document at all"

    second = await _resolve(tmp_path)

    assert len(attempts) == 2 * orthology.ENSEMBL_MAX_ATTEMPTS, "the second run re-asks rather than trusting a failure"
    assert second.unresolved_species == frozenset({"mouse"})
    assert second.cached_species == frozenset()
    assert second.provenance_by_species == {}, "there is no answer to attribute to anything"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_a_resolved_absence_is_cached(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The mirror of the above: "checked, no orthologue" is a real answer and is expensive to re-get.

    Compara charges the same two round trips to say "none" as to say "here it is", and the answer is
    stable between releases, so an empty result must be cached -- while still arriving labelled as a
    cache hit rather than as a request.
    """
    calls: list[str] = []

    async def no_homologies(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        calls.append(url)
        return {"data": [{"id": HUMAN_TP53, "homologies": []}]}

    monkeypatch.setattr(orthology, "ensembl_request_json", no_homologies)

    first = await _resolve(tmp_path)
    assert first.resolved_species == frozenset({"mouse"})
    assert first.all_gene_ids == frozenset()

    second = await _resolve(tmp_path)

    assert len(calls) == 1, "a resolved absence is not re-asked"
    assert second.resolved_species == frozenset({"mouse"}), "and is still a resolved absence, not unresolved"
    assert second.all_gene_ids == frozenset()
    assert second.cached_species == frozenset({"mouse"})
    assert second.provenance_by_species == {"mouse": SOURCE_COMPARA_CACHE}


@pytest.mark.unit
@pytest.mark.asyncio
async def test_a_partly_failed_lookup_is_not_cached_as_the_whole_answer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One dead identifier out of two makes the answer a lower bound, so it is not stored.

    ``resolve_orthologues`` reports a species as resolved when *any* route found something, which is
    right for the run in hand -- the evidence it has is real. It is not right to keep: the cache key
    covers both gene IDs, so storing a partial answer under it would publish a subset as the complete
    orthologue set for that question, indefinitely.
    """
    monkeypatch.setattr(orthology.asyncio, "sleep", _no_sleep)
    calls: list[str] = []

    async def one_id_fails(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        calls.append(url)
        if "ENSG00000012048" in url:
            raise DatabaseAccessError("503 Service Unavailable", "Ensembl")
        return _condensed(MOUSE_TRP53, "mus_musculus")

    monkeypatch.setattr(orthology, "ensembl_request_json", one_id_fails)

    first = await resolve_orthologues({HUMAN_TP53, "ENSG00000012048"}, "human", {"mouse"}, cache_dir=tmp_path)
    assert first.resolved_species == frozenset({"mouse"}), "the evidence actually obtained still counts"
    assert first.all_gene_ids == frozenset({MOUSE_TRP53})
    assert _documents(tmp_path) == [], "but a lower bound is not stored as the answer"

    before = len(calls)
    await resolve_orthologues({HUMAN_TP53, "ENSG00000012048"}, "human", {"mouse"}, cache_dir=tmp_path)
    assert len(calls) > before, "the incomplete question is asked again"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_the_key_covers_the_question_not_just_the_gene(tmp_path: Path, compara: list[str]) -> None:
    """Symbols and the REST host change what was asked, so they cannot share a document.

    The symbol fallback route is part of the question (#101): the same gene IDs with and without a
    symbol are two different lookups, because only the second can rescue a transcript ID. A staging
    or mirrored ``base_url`` is a different authority answering. Both must miss.
    """
    await _resolve(tmp_path)
    assert len(compara) == 1

    await _resolve(tmp_path, query_gene_symbols={"TP53"})
    assert len(compara) == 2, "adding a symbol asks a different question"

    await _resolve(tmp_path, base_url="https://staging.ensembl.example")
    assert len(compara) == 3, "a different REST host is a different authority"

    assert len(_documents(tmp_path)) == 3


@pytest.mark.unit
@pytest.mark.asyncio
async def test_cache_false_neither_reads_nor_writes(tmp_path: Path, compara: list[str]) -> None:
    """``cache=False`` must be a complete opt-out, in both directions.

    Both halves matter: a caller who needs to prove Compara answers this question *now* must not be
    served from disk, and a test that must not touch the disk must not leave a document behind for
    the next one.
    """
    await _resolve(tmp_path, cache=False)
    await _resolve(tmp_path, cache=False)

    assert len(compara) == 2, "no read: both resolutions asked"
    assert list(tmp_path.iterdir()) == [], "no write: nothing was left on disk"


@pytest.mark.unit
def test_the_mapping_file_path_never_reaches_the_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``from_file`` is already offline and free, and a copy of it could outlive an edit.

    Caching a file-sourced mapping would save nothing and would let a document keep answering for the
    TTL after the user edited the file -- silently screening with the evidence they replaced.
    ``SIRNAFORGE_CACHE_DIR`` is pointed here so that "did not reach the cache" is checked against a
    directory we can see, rather than against the user's real cache root.
    """
    monkeypatch.setenv("SIRNAFORGE_CACHE_DIR", str(tmp_path / "root"))
    mapping_file = tmp_path / "orthologs.json"
    mapping_file.write_text(json.dumps({HUMAN_TP53: {"mouse": [MOUSE_TRP53]}}), encoding="utf-8")

    mapping = OrthologueMapping.from_file(mapping_file, {HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert mapping.cached_species == frozenset()
    assert mapping.cache_key_by_species == {}
    assert not (tmp_path / "root").exists(), "the cache directory was never even created"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_a_document_older_than_the_ttl_is_re_asked(tmp_path: Path, compara: list[str]) -> None:
    """Compara releases move on roughly a two-month cadence, so an answer is not kept forever.

    The TTL lives in the sidecar stamp, which is the repo's existing mechanism -- there is no second
    expiry scheme in this module to disagree with it, so ageing the stamp is what ages the answer.
    """
    await _resolve(tmp_path)
    document = _documents(tmp_path)[0]
    stamp_path = artifact_stamp_path(document)
    stamp = json.loads(stamp_path.read_text(encoding="utf-8"))
    stamp["stamped_at"] = "2020-01-01T00:00:00"
    stamp_path.write_text(json.dumps(stamp), encoding="utf-8")

    refreshed = await _resolve(tmp_path)

    assert len(compara) == 2, "an expired answer is re-asked"
    assert refreshed.cached_species == frozenset()
    assert refreshed.all_gene_ids == frozenset({MOUSE_TRP53})


@pytest.mark.unit
@pytest.mark.asyncio
async def test_a_damaged_document_is_a_miss_not_a_failure(tmp_path: Path, compara: list[str]) -> None:
    """A broken cache must cost one request, never a screen.

    Two shapes, and the second is the one a byte fingerprint cannot reach: a half-written document
    (caught by the recorded size, so it is never served for its whole TTL), and a document whose bytes
    are intact and stamped but whose *shape* this release does not recognise.
    """
    await _resolve(tmp_path)
    document = _documents(tmp_path)[0]

    document.write_text("{not json", encoding="utf-8")
    assert (await _resolve(tmp_path)).all_gene_ids == frozenset({MOUSE_TRP53})
    assert len(compara) == 2, "a corrupt document is re-asked rather than raising"

    # Intact, correctly stamped, but written under a schema this release does not read.
    stale_shape = {"schema": orthology.ORTHOLOGY_CACHE_SCHEMA + 1, "orthologue_gene_ids": ["ENSMUSG00000000000"]}
    document.write_text(json.dumps(stale_shape), encoding="utf-8")
    cache_utils.write_artifact_stamp(ORTHOLOGY_CACHE_SUBDIR, document)

    reasked = await _resolve(tmp_path)

    assert len(compara) == 3, "an unreadable shape is a miss, not a reinterpretation"
    assert reasked.all_gene_ids == frozenset({MOUSE_TRP53}), "and never yields the document's contents"
