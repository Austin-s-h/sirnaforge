"""NCBI ELink ingestion uses nested RNA groups, independently of remote API uptime."""

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from sirnaforge.data.base import DatabaseAccessError, RefSeqClient

pytestmark = [pytest.mark.unit, pytest.mark.asyncio]


@pytest.fixture
def response_factory(monkeypatch):
    """Replace only the HTTP transport; run the real client response parser."""

    def install(payload, status=200, headers=None):
        response = MagicMock(status=status)
        response.json = AsyncMock(return_value=payload)
        response.headers = headers or {}
        session = MagicMock()
        session.__aenter__.return_value = session
        session.get.return_value.__aenter__.return_value = response
        monkeypatch.setattr("sirnaforge.data.base.aiohttp.ClientSession", MagicMock(return_value=session))
        return session

    return install


@pytest.mark.parametrize("database", ["nuccore", "nucleotide"])
async def test_nested_refseq_rna_links_are_retrieved(response_factory, database):
    """A real ELink-shaped response has no dbto on its outer linkset."""
    session = response_factory(
        {
            "linksets": [
                {
                    "dbfrom": "gene",
                    "ids": ["7157"],
                    "linksetdbs": [{"dbto": database, "linkname": "gene_nuccore_refseqrna", "links": ["101", "102"]}],
                }
            ]
        }
    )
    assert await RefSeqClient()._find_linked_transcripts("7157") == ["101", "102"]
    params = session.get.call_args.kwargs["params"]
    assert params["dbfrom"] == "gene"
    assert params["db"] == "nuccore"
    assert params["linkname"] == "gene_nuccore_refseqrna"


async def test_non_rna_groups_are_excluded_and_ids_are_deduplicated(response_factory):
    """Never fetch genomic/protein records as transcript sequences."""
    response_factory(
        {
            "linksets": [
                {
                    "linksetdbs": [
                        {"dbto": "nuccore", "linkname": "gene_nuccore_refseqgene", "links": ["999"]},
                        {"dbto": "protein", "linkname": "gene_protein_refseq", "links": ["888"]},
                        {"dbto": "nuccore", "linkname": "gene_nuccore_refseqrna", "links": ["101", 102]},
                    ]
                },
                {"linksetdbs": [{"dbto": "nuccore", "linkname": "gene_nuccore_refseqrna", "links": ["101"]}]},
            ]
        }
    )
    assert await RefSeqClient()._find_linked_transcripts("7157") == ["101", "102"]


async def test_no_rna_links_remains_an_empty_result(response_factory):
    """A genuinely empty association is not fabricated into a transcript."""
    response_factory({"linksets": [{"dbfrom": "gene", "ids": ["7157"]}]})
    assert await RefSeqClient()._find_linked_transcripts("7157") == []


async def test_http_failures_remain_access_errors(response_factory):
    """Service failures must not silently become empty successful searches."""
    response_factory({}, status=503)
    with pytest.raises(DatabaseAccessError, match="HTTP 503"):
        await RefSeqClient()._find_linked_transcripts("7157")


async def test_rate_limit_retry_honors_retry_after(response_factory):
    """A temporary 429 retries the request without losing its RNA response."""
    payload = {
        "linksets": [{"linksetdbs": [{"dbto": "nuccore", "linkname": "gene_nuccore_refseqrna", "links": ["101"]}]}]
    }
    session = response_factory({}, status=429, headers={"Retry-After": "2"})
    retry_context = session.get.return_value
    success_context = MagicMock()
    success_response = MagicMock(status=200)
    success_response.json = AsyncMock(return_value=payload)
    success_context.__aenter__.return_value = success_response
    session.get.side_effect = [retry_context, success_context]
    client = RefSeqClient()
    with (
        patch.object(client, "_throttle", new_callable=AsyncMock),
        patch("sirnaforge.data.base.asyncio.sleep", new_callable=AsyncMock) as sleep,
    ):
        assert await client._find_linked_transcripts("7157") == ["101"]
        sleep.assert_awaited_once_with(2.0)
    assert session.get.call_count == 2


async def test_rate_limit_retries_are_bounded(response_factory):
    """Persistent throttling stays a visible access error after three attempts."""
    session = response_factory({}, status=429, headers={"Retry-After": "invalid"})
    client = RefSeqClient()
    with (
        patch.object(client, "_throttle", new_callable=AsyncMock),
        patch("sirnaforge.data.base.asyncio.sleep", new_callable=AsyncMock) as sleep,
    ):
        with pytest.raises(DatabaseAccessError, match="HTTP 429"):
            await client._find_linked_transcripts("7157")
        assert [call.args[0] for call in sleep.await_args_list] == [1.0, 2.0]
    assert session.get.call_count == 3


async def test_eutilities_requests_are_paced_without_wall_clock_assertions():
    """Reserve request slots at a safe interval using a deterministic clock."""
    client = RefSeqClient()
    clock = MagicMock()
    clock.time.side_effect = [0.0, 0.0, 0.0, 0.35]
    with (
        patch("sirnaforge.data.base.asyncio.get_running_loop", return_value=clock),
        patch("sirnaforge.data.base.asyncio.sleep", new_callable=AsyncMock) as sleep,
    ):
        await client._throttle()
        await client._throttle()
        sleep.assert_awaited_once_with(0.35)
