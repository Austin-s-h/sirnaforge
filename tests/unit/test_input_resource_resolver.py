"""Tests for the input resource resolver utilities."""

from __future__ import annotations

import threading
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

import pytest

from sirnaforge.utils import resolve_input_source


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
def test_resolve_local_file(tmp_path: Path) -> None:
    content = ">seq1\nAUGCUA\n"
    local = tmp_path / "input.fasta"
    local.write_text(content)

    resolved = resolve_input_source(str(local), tmp_path / "inputs")

    assert resolved.source_type == "local"
    assert resolved.local_path == local.resolve()
    assert resolved.size_bytes == len(content)
    assert resolved.sha256 is not None and len(resolved.sha256) == 64


class _SilentHandler(SimpleHTTPRequestHandler):
    """SimpleHTTPRequestHandler that suppresses logging."""

    def log_message(self, _format: str, *_args):  # noqa: D401
        return


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
def test_resolve_http_download(tmp_path: Path) -> None:
    data_root = tmp_path / "server"
    data_root.mkdir()
    served_file = data_root / "remote.fasta"
    served_content = ">seq\nAUGCUAGCUA\n"
    served_file.write_text(served_content)

    handler = partial(_SilentHandler, directory=str(data_root))
    server = HTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        url = f"http://127.0.0.1:{server.server_port}/{served_file.name}"
        resolved = resolve_input_source(url, tmp_path / "downloads")

        assert resolved.downloaded is True
        assert resolved.source_type == "http"
        assert resolved.local_path.exists()
        assert resolved.local_path.read_text() == served_content
    finally:
        server.shutdown()
        thread.join()


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
def test_resolve_unsupported_scheme(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        resolve_input_source("s3://bucket/key.fasta", tmp_path / "downloads")
