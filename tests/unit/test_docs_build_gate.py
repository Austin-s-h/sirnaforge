"""The documentation gate must fail on documentation defects only (#107).

`make docs` runs Sphinx with `-W`, and `make test-release` runs `make docs` before its first test. So
anything Sphinx warns about blocks the release -- which is the point for a duplicate object or a broken
reference, and wrong for an intersphinx inventory that a proxy would not serve. #107 was filed on 39
duplicate-object warnings; those are fixed, but the same gate still fell over on four unreachable
inventories, so the criterion "`make docs` exits 0 with warnings still treated as errors" could not be
demonstrated on a restricted network at all.

These tests pin the narrow escape hatch: unreachable targets are dropped, everything else still warns.
They are hermetic -- the reachable case is served by a loopback HTTP server, never the real internet.
"""

import importlib.util
import socket
import sys
import threading
from collections.abc import Iterator
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import pytest

CONF_PATH = Path(__file__).resolve().parents[2] / "docs" / "conf.py"


@pytest.fixture(scope="module")
def docs_conf() -> Any:
    """Import ``docs/conf.py`` as a module without running Sphinx."""
    spec = importlib.util.spec_from_file_location("sirnaforge_docs_conf", CONF_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class _Inventory(BaseHTTPRequestHandler):
    def do_HEAD(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler's naming
        self.send_response(200 if self.path.endswith("objects.inv") else 404)
        self.end_headers()

    def log_message(self, *_: object) -> None:
        """Silence the per-request stderr line; the test asserts on our own output."""


@pytest.fixture
def local_inventory() -> Iterator[str]:
    """A loopback base URL whose ``objects.inv`` answers HEAD with 200."""
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    server = ThreadingHTTPServer(("127.0.0.1", port), _Inventory)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{port}/docs/"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def _unreachable_base() -> str:
    """A port nothing is listening on, so the probe fails without leaving the machine."""
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return f"http://127.0.0.1:{probe.getsockname()[1]}/docs/"


def test_the_real_config_declares_every_target_it_intends_to_link(docs_conf: Any) -> None:
    """The targets are declared once, so the probe cannot quietly shrink the intended set."""
    assert set(docs_conf.INTERSPHINX_TARGETS) == {"python", "numpy", "pandas", "biopython"}


def test_an_unreachable_inventory_is_dropped_rather_than_warned_about(
    docs_conf: Any, capsys: pytest.CaptureFixture[str]
) -> None:
    """The #107 gate defect: one unreachable host used to fail the whole build under ``-W``."""
    mapping = docs_conf.reachable_intersphinx_mapping({"nowhere": (_unreachable_base(), None)}, timeout=2)

    assert mapping == {}
    assert "dropping unreachable intersphinx target 'nowhere'" in capsys.readouterr().err


def test_a_reachable_inventory_is_kept(docs_conf: Any, local_inventory: str) -> None:
    """The probe must not be a disguised opt-out: a target that answers is still linked."""
    mapping = docs_conf.reachable_intersphinx_mapping({"local": (local_inventory, None)}, timeout=5)

    assert mapping == {"local": (local_inventory, None)}


def test_one_unreachable_target_does_not_take_the_reachable_ones_with_it(docs_conf: Any, local_inventory: str) -> None:
    """Sphinx reports inventory failures as one warning, so partial reachability must still build."""
    mapping = docs_conf.reachable_intersphinx_mapping(
        {"local": (local_inventory, None), "nowhere": (_unreachable_base(), None)}, timeout=2
    )

    assert mapping == {"local": (local_inventory, None)}


def test_offline_skips_the_probe_entirely(docs_conf: Any, capsys: pytest.CaptureFixture[str]) -> None:
    """A hermetic build asks for no network at all, and says so."""
    mapping = docs_conf.reachable_intersphinx_mapping(
        {"python": ("https://docs.python.org/3", None)}, offline=True, timeout=1
    )

    assert mapping == {}
    assert "SIRNAFORGE_DOCS_OFFLINE" in capsys.readouterr().err


def test_an_explicit_inventory_url_is_probed_instead_of_the_derived_one(docs_conf: Any, local_inventory: str) -> None:
    """Sphinx allows a second element naming the inventory; the probe must honour it, not guess."""
    explicit = f"{local_inventory}objects.inv"

    assert docs_conf.reachable_intersphinx_mapping({"local": ("http://127.0.0.1:1/", explicit)}, timeout=5) == {
        "local": ("http://127.0.0.1:1/", explicit)
    }
