"""Unit tests for the batched Ensembl REST request paths.

Fetching one sequence per request cost ~350s for a 40-transcript gene and lost several
transcripts to transient 503s, so both the sequence client and the annotation client POST
their identifiers in chunks. These tests pin the batching, the per-id fallback for whatever
a batch omits, and the retry on Ensembl's transient statuses.
"""

from collections.abc import Callable
from typing import Any
from unittest.mock import AsyncMock, patch

import pytest

from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.base import DatabaseAccessError, EnsemblClient, ensembl_request_json
from sirnaforge.data.transcript_annotation import EnsemblTranscriptModelClient


class _FakeResponse:
    """Stand-in for an aiohttp response used as an async context manager."""

    def __init__(self, status: int, payload: Any = None, text: str = "", headers: dict[str, str] | None = None):
        self.status = status
        self.headers = headers or {}
        self._payload = payload
        self._text = text

    async def json(self) -> Any:
        return self._payload

    async def text(self) -> str:
        return self._text

    async def __aenter__(self) -> "_FakeResponse":
        return self

    async def __aexit__(self, *_: object) -> None:
        return None


class _FakeSession:
    """Stand-in for aiohttp.ClientSession that answers from a handler and records calls."""

    def __init__(self, handler: Callable[[str, str], _FakeResponse]):
        self._handler = handler
        self.calls: list[tuple[str, str]] = []

    def request(self, method: str, url: str, **_: Any) -> _FakeResponse:
        self.calls.append((method, url))
        return self._handler(method, url)

    def get(self, url: str, **kwargs: Any) -> _FakeResponse:
        return self.request("GET", url, **kwargs)

    async def __aenter__(self) -> "_FakeSession":
        return self

    async def __aexit__(self, *_: object) -> None:
        return None


def _patch_session(session: _FakeSession) -> Any:
    """Patch aiohttp.ClientSession so every client in the process uses ``session``."""
    return patch("sirnaforge.data.base.aiohttp.ClientSession", return_value=session)


def _sequence_records(identifiers: list[str]) -> list[dict[str, str]]:
    """One record per id, each with a distinct lowercase sequence carrying a newline."""
    return [{"id": i, "query": i, "seq": "acgt\n" * (n + 1)} for n, i in enumerate(identifiers)]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_get_sequences_posts_one_request_per_chunk():
    """120 identifiers become 3 POSTs, not 120 GETs, and sequences are normalised."""
    identifiers = [f"ENST{i:011d}" for i in range(120)]

    def handler(method: str, url: str) -> _FakeResponse:
        assert method == "POST"
        assert url.endswith("/sequence/id?species=homo_sapiens&type=cdna")
        return _FakeResponse(200, payload=_sequence_records(identifiers[:50]))

    session = _FakeSession(handler)
    with _patch_session(session):
        sequences = await EnsemblClient().get_sequences(identifiers)

    assert len(session.calls) == 3
    assert {m for m, _ in session.calls} == {"POST"}
    assert sequences["ENST00000000000"] == "ACGT"
    assert sequences["ENST00000000001"] == "ACGTACGT"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_get_transcripts_retries_per_id_for_what_the_batch_omits():
    """A transcript missing from the POST response is fetched with a single GET."""
    gene_payload = {
        "id": "ENSG1",
        "display_name": "TP53",
        "Transcript": [{"id": "ENST1"}, {"id": "ENST2"}],
    }

    def handler(method: str, url: str) -> _FakeResponse:
        if method == "POST":
            # Ensembl silently omits ids it cannot serve rather than failing the request.
            return _FakeResponse(200, payload=_sequence_records(["ENST1"]))
        return _FakeResponse(200, text="TTTTTTTT")

    session = _FakeSession(handler)
    client = EnsemblClient()
    with _patch_session(session), patch.object(client, "_lookup_gene_data", AsyncMock(return_value=gene_payload)):
        transcripts = await client._get_transcripts("ENSG1", include_sequence=True)

    assert [t.transcript_id for t in transcripts] == ["ENST1", "ENST2"]
    assert transcripts[0].sequence == "ACGT"
    assert transcripts[1].sequence == "TTTTTTTT"
    assert session.calls.count(("POST", session.calls[0][1])) == 1
    assert sum(1 for method, _ in session.calls if method == "GET") == 1


@pytest.mark.unit
@pytest.mark.asyncio
async def test_get_transcripts_keeps_a_transcript_whose_sequence_cannot_be_fetched():
    """An id that fails the batch and the per-id retry is still reported, without a sequence."""
    gene_payload = {"id": "ENSG1", "display_name": "TP53", "Transcript": [{"id": "ENST1"}]}

    def handler(method: str, _: str) -> _FakeResponse:
        return _FakeResponse(200, payload=[]) if method == "POST" else _FakeResponse(404)

    session = _FakeSession(handler)
    client = EnsemblClient()
    with _patch_session(session), patch.object(client, "_lookup_gene_data", AsyncMock(return_value=gene_payload)):
        transcripts = await client._get_transcripts("ENSG1", include_sequence=True)

    assert len(transcripts) == 1
    assert transcripts[0].sequence is None


@pytest.mark.unit
@pytest.mark.asyncio
async def test_request_json_retries_transient_status_then_succeeds():
    """503 is transient: retry rather than dropping the request."""
    statuses = [503, 200]

    def handler(_: str, __: str) -> _FakeResponse:
        status = statuses.pop(0)
        return _FakeResponse(status, payload={"ok": True} if status == 200 else None)

    session = _FakeSession(handler)
    with patch("sirnaforge.data.base.asyncio.sleep", AsyncMock()) as sleep:
        payload = await ensembl_request_json(session, "GET", "https://rest.ensembl.org/lookup/id/X", headers={})

    assert payload == {"ok": True}
    assert len(session.calls) == 2
    sleep.assert_awaited_once()


@pytest.mark.unit
@pytest.mark.asyncio
async def test_request_json_honours_retry_after():
    """A Retry-After header sets the wait, capped by the client timeout."""

    def handler(_: str, __: str) -> _FakeResponse:
        return _FakeResponse(429, headers={"Retry-After": "2"})

    session = _FakeSession(handler)
    with patch("sirnaforge.data.base.asyncio.sleep", AsyncMock()) as sleep, pytest.raises(DatabaseAccessError):
        await ensembl_request_json(session, "GET", "https://rest.ensembl.org/lookup/id/X", headers={})

    assert [call.args[0] for call in sleep.await_args_list] == [2.0, 2.0]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_fetch_by_ids_resolves_annotations_and_symbols_in_two_posts():
    """Transcripts and their parent genes are each looked up in one POST, with no per-id GET."""
    transcript_payload = {
        "ENST1": {
            "id": "ENST1",
            "Parent": "ENSG1",
            "biotype": "protein_coding",
            "seq_region_name": "17",
            "start": 100,
            "end": 200,
            "strand": -1,
        },
        "ENST2": {
            "id": "ENST2",
            "Parent": "ENSG1",
            "biotype": "protein_coding",
            "seq_region_name": "17",
            "start": 300,
            "end": 400,
            "strand": -1,
        },
    }
    gene_payload = {
        "ENSG1": {
            "id": "ENSG1",
            "display_name": "TP53",
            "seq_region_name": "17",
            "start": 100,
            "end": 400,
            "strand": -1,
        }
    }

    def handler(method: str, url: str) -> _FakeResponse:
        assert method == "POST", f"unexpected per-id request: {url}"
        return _FakeResponse(200, payload=gene_payload if "expand=1" not in url else transcript_payload)

    session = _FakeSession(handler)
    with _patch_session(session):
        bundle = await EnsemblTranscriptModelClient().fetch_by_ids(
            ["ENST1", "ENST2"],
            species="human",
            reference=ReferenceChoice.explicit("GRCh38", reason="test"),
        )

    assert len(session.calls) == 2
    assert bundle.resolved_count == 2
    assert bundle.unresolved_count == 0
    assert {a.symbol for a in bundle.transcripts.values()} == {"TP53"}
    assert bundle.transcripts["ENST1"].gene_interval is not None
