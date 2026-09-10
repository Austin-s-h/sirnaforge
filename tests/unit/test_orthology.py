"""Unit tests for Compara orthologue resolution.

No network: the Ensembl request function is patched, so these pin response parsing, the
paralogue exclusion, and the degrade-not-raise contract. The live endpoint these fixtures imitate is
``/homology/id/{species}/{gene}?target_species=...&type=orthologues&format=condensed``, whose
condensed payload shape was captured from release 116.

The last group covers the offline path and the two guards that keep an unreachable Compara off the
critical path -- a mapping file (#101 point 4), no retries for a transport failure, and a wall-clock
budget.
"""

import json
import ssl
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.data import orthology
from sirnaforge.data.base import DatabaseAccessError
from sirnaforge.data.ensembl_references import infer_species_from_cdna_headers
from sirnaforge.data.orthology import (
    OrthologueMapping,
    load_ortholog_table,
    resolve_orthologues,
)
from sirnaforge.data.species_registry import ensembl_species_slug

HUMAN_TP53 = "ENSG00000141510"
MOUSE_TRP53 = "ENSMUSG00000059552"


def _condensed(*homologies: dict[str, str]) -> dict[str, Any]:
    """A condensed homology payload for one query gene."""
    return {"data": [{"id": HUMAN_TP53, "homologies": list(homologies)}]}


async def _no_sleep(_seconds: float) -> None:
    """Collapse retry backoff so the retry tests stay fast."""


@pytest.fixture
def captured_urls(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Patch the Ensembl request helper, recording URLs and returning a one2one mouse orthologue."""
    urls: list[str] = []

    async def fake_request(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        urls.append(url)
        return _condensed({"id": MOUSE_TRP53, "species": "mus_musculus", "type": "ortholog_one2one"})

    monkeypatch.setattr(orthology, "ensembl_request_json", fake_request)
    return urls


@pytest.mark.unit
def test_ensembl_species_slug_derives_from_scientific_name():
    """The REST slug is the underscored scientific name, not the canonical common name."""
    assert ensembl_species_slug("mouse") == "mus_musculus"
    assert ensembl_species_slug("hsa") == "homo_sapiens"
    assert ensembl_species_slug("Rattus norvegicus") == "rattus_norvegicus"
    assert ensembl_species_slug("nonesuch") is None


@pytest.mark.unit
@pytest.mark.asyncio
async def test_resolves_mouse_orthologue_gene_id(captured_urls: list[str]):
    """The mouse orthologue arrives as a gene ID -- the evidence symbol equality cannot supply."""
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.gene_ids_by_species == {"mouse": frozenset({MOUSE_TRP53})}
    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert mapping.resolved_species == frozenset({"mouse"})
    assert mapping.unresolved_species == frozenset()
    assert len(captured_urls) == 1
    assert f"/homology/id/homo_sapiens/{HUMAN_TP53}" in captured_urls[0]
    assert "target_species=mus_musculus" in captured_urls[0]
    assert "type=orthologues" in captured_urls[0]
    # The slug, not the common name: target_species=mouse is rejected by the live endpoint.
    assert "target_species=mouse&" not in captured_urls[0]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_persistent_error_body_on_a_200_is_a_failure_not_an_empty_result(
    monkeypatch: pytest.MonkeyPatch,
):
    """Compara answers valid queries with HTTP 200 + {"error": ...} intermittently.

    Status-based retry never sees it, and an absent ``data`` key reads as "no orthologue exists".
    When it persists across attempts the species must land in ``unresolved_species`` -- never be
    reported as a confident absence.
    """
    calls = 0

    async def error_body(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        return {"error": "DBD::mysql::st bind_param failed: Illegal parameter number"}

    monkeypatch.setattr(orthology, "ensembl_request_json", error_body)
    monkeypatch.setattr(orthology.asyncio, "sleep", _no_sleep)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset()
    assert mapping.unresolved_species == frozenset({"mouse"})
    assert mapping.resolved_species == frozenset()
    assert calls == orthology.ENSEMBL_MAX_ATTEMPTS, "the transient error body is retried, not surrendered to"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_transient_error_body_is_retried_then_succeeds(monkeypatch: pytest.MonkeyPatch):
    """The observed live behaviour: the same URL errors, then works. It must recover."""
    calls = 0

    async def flaky(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        if calls == 1:
            return {"error": "DBD::mysql::st bind_param failed: Illegal parameter number"}
        return _condensed({"id": MOUSE_TRP53, "species": "mus_musculus", "type": "ortholog_one2one"})

    monkeypatch.setattr(orthology, "ensembl_request_json", flaky)
    monkeypatch.setattr(orthology.asyncio, "sleep", _no_sleep)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert mapping.unresolved_species == frozenset()
    assert calls == 2


@pytest.mark.unit
@pytest.mark.asyncio
async def test_query_species_is_never_looked_up(captured_urls: list[str]):
    """A same-species screen needs no Compara call at all."""
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"human"})

    assert mapping == OrthologueMapping.empty()
    assert captured_urls == []


@pytest.mark.unit
@pytest.mark.asyncio
async def test_paralogues_are_not_orthologues(monkeypatch: pytest.MonkeyPatch):
    """A within-species duplicate is an off-target liability, never conservation evidence."""

    async def fake_request(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        return _condensed(
            {"id": "ENSMUSG00000000001", "species": "mus_musculus", "type": "within_species_paralog"},
            {"id": MOUSE_TRP53, "species": "mus_musculus", "type": "ortholog_one2many"},
        )

    monkeypatch.setattr(orthology, "ensembl_request_json", fake_request)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53}), "the paralogue must not be admitted"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_versioned_ids_are_stripped(monkeypatch: pytest.MonkeyPatch):
    """Gene IDs compare equal across Ensembl releases, so versions are stripped both ways."""

    async def fake_request(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        assert f"{HUMAN_TP53}?" in url, "the query gene's version suffix is stripped before the request"
        return _condensed({"id": f"{MOUSE_TRP53}.4", "species": "mus_musculus", "type": "ortholog_one2one"})

    monkeypatch.setattr(orthology, "ensembl_request_json", fake_request)
    mapping = await resolve_orthologues({f"{HUMAN_TP53}.18"}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})


@pytest.mark.unit
@pytest.mark.asyncio
async def test_request_failure_degrades_rather_than_raising(monkeypatch: pytest.MonkeyPatch):
    """Compara being unavailable weakens evidence; it must not fail the screen.

    The species lands in ``unresolved_species`` so the caller can say orthology was not decided,
    rather than reporting a confident "no orthologue exists".
    """

    async def boom(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        raise RuntimeError("503 Service Unavailable")

    monkeypatch.setattr(orthology, "ensembl_request_json", boom)
    monkeypatch.setattr(orthology.asyncio, "sleep", _no_sleep)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset()
    assert mapping.unresolved_species == frozenset({"mouse"})
    assert mapping.resolved_species == frozenset()


@pytest.mark.unit
@pytest.mark.asyncio
async def test_resolved_but_empty_is_not_unresolved(monkeypatch: pytest.MonkeyPatch):
    """A resolved absence and an unchecked species are different claims, reported differently."""

    async def no_homologies(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        return {"data": [{"id": HUMAN_TP53, "homologies": []}]}

    monkeypatch.setattr(orthology, "ensembl_request_json", no_homologies)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset()
    assert mapping.resolved_species == frozenset({"mouse"})
    assert mapping.unresolved_species == frozenset()


@pytest.mark.unit
@pytest.mark.asyncio
async def test_unregistered_target_species_is_unresolved(captured_urls: list[str]):
    """A species with no registry entry cannot be addressed on the REST API; say so."""
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"nonesuch"})

    assert mapping.unresolved_species == frozenset({"nonesuch"})
    assert captured_urls == []


@pytest.mark.unit
@pytest.mark.asyncio
async def test_malformed_payload_yields_no_orthologues(monkeypatch: pytest.MonkeyPatch):
    """A shape change should weaken evidence, not crash a screen."""

    async def junk(_session: Any, _method: str, _url: str, **_kwargs: Any) -> Any:
        return {"unexpected": "shape"}

    monkeypatch.setattr(orthology, "ensembl_request_json", junk)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset()
    assert mapping.resolved_species == frozenset({"mouse"})


@pytest.mark.unit
@pytest.mark.asyncio
async def test_summary_records_provenance(captured_urls: list[str]):
    """A conservation claim must be auditable back to its source."""
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})
    summary = mapping.summary()

    assert summary["source"] == "ensembl_compara"
    assert summary["gene_ids_by_species"] == {"mouse": [MOUSE_TRP53]}
    assert summary["resolved_species"] == ["mouse"]
    assert "within_species_paralog" not in summary["orthologue_types"]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_transient_http_error_is_retried_then_succeeds(monkeypatch: pytest.MonkeyPatch):
    """Compara's flakiness also surfaces as HTTP 400, which the shared status retry skips.

    400 is not in ENSEMBL_RETRY_STATUSES because it normally means a bad request, but this endpoint
    returns it for valid queries that succeed on the next attempt (measured 5/6 unretried).
    """
    calls = 0

    async def flaky(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        if calls == 1:
            raise RuntimeError("HTTP 400")
        return _condensed({"id": MOUSE_TRP53, "species": "mus_musculus", "type": "ortholog_one2one"})

    monkeypatch.setattr(orthology, "ensembl_request_json", flaky)
    monkeypatch.setattr(orthology.asyncio, "sleep", _no_sleep)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert mapping.unresolved_species == frozenset()
    assert calls == 2


@pytest.mark.unit
def test_species_inferred_from_cdna_headers(tmp_path):
    """A custom cDNA path must take its species from the reference, not from the parameter (#99)."""
    mouse = tmp_path / "mouse.fa"
    mouse.write_text(
        ">ENSMUST00000108658.11 cdna chromosome:GRCm39:11:69471109:69482701:1 "
        "gene:ENSMUSG00000059552.15 gene_symbol:Trp53\nACGT\n"
    )
    assert infer_species_from_cdna_headers(mouse) == "mouse"

    human = tmp_path / "human.fa"
    human.write_text(
        ">ENST00000269305.9 cdna chromosome:GRCh38:17:7668421:7687550:-1 "
        "gene:ENSG00000141510.18 gene_symbol:TP53\nACGT\n"
    )
    assert infer_species_from_cdna_headers(human) == "human"


@pytest.mark.unit
def test_unidentifiable_and_mixed_headers_infer_nothing(tmp_path):
    """An honest unknown beats inventing a species; a mixed file must not be assigned one label."""
    bare = tmp_path / "bare.fa"
    bare.write_text(">some_contig_1\nACGT\n>some_contig_2\nACGT\n")
    assert infer_species_from_cdna_headers(bare) is None

    mixed = tmp_path / "mixed.fa"
    mixed.write_text(
        ">ENSMUST1 cdna chromosome:GRCm39:11:1:2:1 gene:ENSMUSG1\nACGT\n"
        ">ENST1 cdna chromosome:GRCh38:17:1:2:1 gene:ENSG1\nACGT\n"
    )
    assert infer_species_from_cdna_headers(mixed) is None

    assert infer_species_from_cdna_headers(tmp_path / "absent.fa") is None


@pytest.mark.unit
@pytest.mark.asyncio
async def test_symbol_route_rescues_a_transcript_id_query(monkeypatch: pytest.MonkeyPatch):
    """An input FASTA supplies transcript IDs, which Compara answers with an empty 200.

    /homology/id/homo_sapiens/ENST00000413465 returns HTTP 200 {"data": []} — a successful empty
    result. Without the symbol fallback that is silently read as "TP53 has no mouse orthologue",
    which is what made the first mouse-only run report 0 orthologs against 10,217 off-targets.
    """
    seen: list[str] = []

    async def by_route(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        seen.append(url)
        if "/homology/id/" in url:
            return {"data": []}
        return _condensed({"id": MOUSE_TRP53, "species": "mus_musculus", "type": "ortholog_one2one"})

    monkeypatch.setattr(orthology, "ensembl_request_json", by_route)
    mapping = await resolve_orthologues({"ENST00000413465"}, "human", {"mouse"}, query_gene_symbols={"TP53"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert mapping.resolved_species == frozenset({"mouse"})
    assert any("/homology/id/" in u for u in seen), "the ID route is still tried first"
    assert any("/homology/symbol/homo_sapiens/TP53" in u for u in seen), "the symbol route is the fallback"
    # Provenance must show what was asked, so a zero result is diagnosable.
    assert mapping.summary()["queried_gene_ids"] == ["ENST00000413465"]
    assert mapping.summary()["queried_symbols"] == ["TP53"]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_symbol_route_is_skipped_when_ids_already_resolved(captured_urls: list[str]):
    """The cheap path stays cheap: a working gene ID must not trigger a second lookup."""
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"}, query_gene_symbols={"TP53"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert len(captured_urls) == 1
    assert "/homology/symbol/" not in captured_urls[0]


def _mapping_file(tmp_path: Path, payload: Any) -> Path:
    """Write a user-supplied ortholog mapping file."""
    path = tmp_path / "orthologs.json"
    path.write_text(json.dumps(payload) if not isinstance(payload, str) else payload)
    return path


@pytest.mark.unit
def test_mapping_file_resolves_without_the_network(tmp_path: Path):
    """The offline path #101 requires: orthologues stated in a file, no REST call, honest provenance."""
    path = _mapping_file(tmp_path, {f"{HUMAN_TP53}.18": {"mouse": [f"{MOUSE_TRP53}.4"]}})

    mapping = OrthologueMapping.from_file(path, {HUMAN_TP53}, "human", {"mouse", "human"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53}), "versions are stripped on both sides"
    assert mapping.resolved_species == frozenset({"mouse"}), "the query species is never a target"
    assert mapping.summary()["source"] == "ortholog_mapping_file", "an offline run must not claim a REST call"


@pytest.mark.unit
def test_mapping_file_silence_about_a_species_is_unresolved_not_empty(tmp_path: Path):
    """A species the file omits was never checked; only an explicit empty list claims absence."""
    path = _mapping_file(tmp_path, {HUMAN_TP53: {"mouse": [MOUSE_TRP53], "rhesus": []}})

    mapping = OrthologueMapping.from_file(path, {HUMAN_TP53}, "human", {"mouse", "rat", "rhesus"})

    assert mapping.unresolved_species == frozenset({"rat"}), "silence is not evidence of absence"
    # "rhesus" is canonicalised to "macaque", so a file may spell a species any registered way.
    assert mapping.resolved_species == frozenset({"mouse", "macaque"}), "an empty list is a checked absence"
    assert mapping.gene_ids_by_species == {"mouse": frozenset({MOUSE_TRP53})}


@pytest.mark.unit
def test_mapping_file_may_be_keyed_on_a_symbol(tmp_path: Path):
    """An input-FASTA run has no stable gene ID, so a symbol key must work -- case-insensitively."""
    path = _mapping_file(tmp_path, {"tp53": {"mouse": [MOUSE_TRP53]}})

    mapping = OrthologueMapping.from_file(path, {"ENST00000413465"}, "human", {"mouse"}, query_gene_symbols={"TP53"})

    assert mapping.all_gene_ids == frozenset({MOUSE_TRP53})
    assert mapping.summary()["queried_symbols"] == ["TP53"], "provenance records what was asked"


@pytest.mark.unit
def test_a_malformed_mapping_file_fails_loudly(tmp_path: Path):
    """An explicit user input is not a flaky network call: ignoring it would silently weaken evidence."""
    with pytest.raises(ValueError, match="not valid JSON"):
        load_ortholog_table(_mapping_file(tmp_path, "{not json"))

    with pytest.raises(ValueError, match="must be a list of gene IDs"):
        load_ortholog_table(_mapping_file(tmp_path, {HUMAN_TP53: {"mouse": MOUSE_TRP53}}))

    with pytest.raises(ValueError, match="must map species"):
        load_ortholog_table(_mapping_file(tmp_path, {HUMAN_TP53: [MOUSE_TRP53]}))

    with pytest.raises(FileNotFoundError):
        load_ortholog_table(tmp_path / "absent.json")


@pytest.mark.unit
@pytest.mark.asyncio
async def test_an_unreachable_host_is_not_retried(monkeypatch: pytest.MonkeyPatch):
    """Backoff cannot fix a rejected certificate; it only makes a certain failure slower.

    Behind a TLS-intercepting proxy every attempt died at the handshake in under a second, so the
    three attempts and their 2s + 4s sleeps were pure latency -- 25s per screen, and the reason
    every test in test_hit_class_persistence.py used to take that long.
    """
    attempts = 0
    slept: list[float] = []

    async def refused(_session: Any, _method: str, _url: str, **_kwargs: Any) -> dict[str, Any]:
        nonlocal attempts
        attempts += 1
        try:
            raise ssl.SSLCertVerificationError("unable to get local issuer certificate")
        except ssl.SSLCertVerificationError as exc:
            # The shape ensembl_request_json produces: the cause carries the real failure.
            raise DatabaseAccessError("Connection failed", "Ensembl") from exc

    async def record_sleep(seconds: float) -> None:
        slept.append(seconds)

    monkeypatch.setattr(orthology, "ensembl_request_json", refused)
    monkeypatch.setattr(orthology.asyncio, "sleep", record_sleep)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse"})

    assert attempts == 1, "a transport-level failure is not retried"
    assert slept == [], "and costs no backoff"
    assert mapping.unresolved_species == frozenset({"mouse"}), "still degrades rather than raising"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_the_budget_bounds_a_whole_resolution(monkeypatch: pytest.MonkeyPatch):
    """A host that drops packets fails only at the per-request timeout, so the total needs a ceiling.

    Without one, a three-species screen against a firewalled Compara pays every species' full
    retry ladder before degrading to the labelled heuristic.
    """
    seen: list[str] = []

    async def slow_timeout(_session: Any, _method: str, url: str, **_kwargs: Any) -> dict[str, Any]:
        seen.append(url)
        await orthology.asyncio.sleep(0.05)
        raise DatabaseAccessError("Request timeout", "Ensembl")

    monkeypatch.setattr(orthology, "ensembl_request_json", slow_timeout)
    mapping = await resolve_orthologues({HUMAN_TP53}, "human", {"mouse", "rat"}, budget=0.01)

    assert len(seen) == 1, "the second species is abandoned rather than paying the same ladder"
    assert mapping.unresolved_species == frozenset({"mouse", "rat"})
    assert mapping.resolved_species == frozenset(), "an abandoned species is never reported clean"


@pytest.mark.unit
def test_the_cli_hands_the_mapping_file_down_to_the_workflow(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """An air-gapped user reaches the offline path from the command line, not just from Python."""
    captured: dict[str, Any] = {}

    async def fake_workflow(**kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return {"transcript_summary": {}, "design_summary": {}}

    monkeypatch.setattr("sirnaforge.cli.run_sirna_workflow", fake_workflow)
    fasta = tmp_path / "toy.fa"
    fasta.write_text(">trans1\nATG" + "A" * 300 + "TAA\n")
    path = _mapping_file(tmp_path, {HUMAN_TP53: {"mouse": [MOUSE_TRP53]}})

    result = CliRunner().invoke(
        app,
        [
            "workflow",
            "TP53",
            "--input-fasta",
            str(fasta),
            "--output-dir",
            str(tmp_path / "out"),
            "--species",
            "human,mouse",
            "--ortholog-mapping",
            str(path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["ortholog_mapping_file"] == path
