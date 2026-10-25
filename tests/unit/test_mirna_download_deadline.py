"""Actual worker boundaries for slow body transfers and failed refresh preservation."""

import gzip
import threading
import time
from dataclasses import replace
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

from sirnaforge.data.mirna_manager import MiRNADatabaseManager

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("compressed", [False, True])
def test_worker_downloads_and_decodes_local_reference(tmp_path: Path, compressed: bool) -> None:
    """The isolated worker preserves gzip, UTF-8 and HTML-entity decoding behavior."""
    payload = b">hsa-local&#32;fixture<br>ACGUACGU\n"
    path = tmp_path / ("mature.fa.gz" if compressed else "mature.fa")
    path.write_bytes(gzip.compress(payload) if compressed else payload)
    manager = MiRNADatabaseManager(cache_dir=tmp_path / "cache")
    source = replace(manager.get_source_configuration("mirbase", "human"), url=path.as_uri(), compressed=compressed)
    assert manager._download_file(source, timeout=5) == ">hsa-local fixture\nACGUACGU\n"


def test_trickling_body_cannot_extend_deadline_or_replace_valid_cache(tmp_path: Path) -> None:
    """A live local server sends often enough to defeat socket timeouts, but the worker is reaped."""
    stopped = threading.Event()
    started = threading.Event()
    received = threading.Event()

    class TrickleHandler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:  # noqa: N802
            """Send a valid FASTA header followed by a deliberately unfinished body."""
            self.send_response(200)
            self.send_header("Content-Length", "1000000")
            self.end_headers()
            received.set()
            try:
                self.wfile.write(b">hsa-slow\n")
                self.wfile.flush()
                started.set()
                while not stopped.wait(0.02):
                    self.wfile.write(b"A")
                    self.wfile.flush()
            except (BrokenPipeError, ConnectionResetError):
                pass

        def log_message(self, _format: str, *args: object) -> None:
            """Keep HTTP fixture diagnostics out of the application's logs."""

    server = ThreadingHTTPServer(("127.0.0.1", 0), TrickleHandler)
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    manager = MiRNADatabaseManager(
        cache_dir=tmp_path / "cache", download_timeout=1, download_attempts=1, retry_backoff_seconds=0
    )
    original_source = manager.get_source_configuration("mirbase", "human")
    source = replace(original_source, url=f"http://127.0.0.1:{server.server_port}/mature.fa")
    # Prime a valid cache under exactly the key the slow refresh will target.
    cache_file = manager.cache_dir / f"{source.cache_key()}.fa"
    cache_file.write_text(">hsa-valid\nACGUACGU\n")
    manager._record_cache_entry(source.cache_key(), source, cache_file, persist=True)
    before, metadata = cache_file.read_bytes(), manager.metadata_file.read_bytes()
    manager.SOURCES = {**manager.SOURCES, "mirbase": {**manager.SOURCES["mirbase"], "human": source}}
    try:
        start = time.monotonic()
        assert manager.get_database("mirbase", "human", force_refresh=True) is None
        elapsed = time.monotonic() - start
        assert received.is_set() and started.is_set(), "worker must reach the trickling body"
        assert 0.8 <= elapsed < 3.0, elapsed
        assert cache_file.read_bytes() == before
        assert manager.metadata_file.read_bytes() == metadata
        assert manager.get_reference_metadata("mirbase", "human") is None
    finally:
        stopped.set()
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)
