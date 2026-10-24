"""Offline acceptance coverage for pinned mature references and safe refreshes."""

import hashlib
import json
from datetime import datetime, timedelta
from pathlib import Path
from urllib.error import HTTPError

import pytest

from sirnaforge.data.mirna_manager import MiRNADatabaseManager, PinnedMiRNAReference, load_pinned_references

pytestmark = pytest.mark.unit


def pin_file(tmp_path, content=">hsa-test\nACGUACGU\n", **kwargs):
    """Pin file."""
    path = tmp_path / "mature.fa"
    path.write_text(content)
    return PinnedMiRNAReference(
        path=path,
        release="22.1",
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        original_url="https://example.org/releases/22.1/mature.fa",
        **kwargs,
    )


def manager_for(tmp_path, pin, monkeypatch):
    """Manager for."""
    manager = MiRNADatabaseManager(cache_dir=tmp_path / "cache", pinned_references={"mirbase": {"human": pin}})

    def forbidden(*args, **kwargs):
        pytest.fail("Pinned references must never contact upstream")

    monkeypatch.setattr(manager, "_download_file", forbidden)
    return manager


def test_offline_pin_filters_species_and_records_actual_provenance(tmp_path, monkeypatch):
    """Offline pin filters species and records actual provenance."""
    pin = pin_file(tmp_path, ">mmu-other\nCCCCAAAA\n>hsa-test\nACGUACGU\n")
    manager = manager_for(tmp_path, pin, monkeypatch)
    result = manager.get_database("mirbase", "Homo sapiens")
    assert result is not None
    assert result.read_text() == ">hsa-test\nACGUACGU"
    provenance = manager.get_reference_metadata("mirbase", "human")
    assert provenance["release"] == "22.1"
    assert provenance["raw_sha256"] == pin.sha256
    assert provenance["filtered_sha256"] == hashlib.sha256(result.read_bytes()).hexdigest()
    assert provenance["acquisition_url"] == Path(pin.path).as_uri()
    assert provenance["original_url"] == pin.original_url
    assert provenance["species_filtered"] is True
    assert provenance["acquisition_mode"] == "pinned_import"
    provenance["release"] = "mutated"
    assert manager.get_reference_metadata("mirbase", "human")["release"] == "22.1"


def test_pin_cache_survives_ttl_and_missing_original(tmp_path, monkeypatch):
    """Pin cache survives ttl and missing original."""
    pin = pin_file(tmp_path)
    manager = manager_for(tmp_path, pin, monkeypatch)
    result = manager.get_database("mirbase", "human")
    meta = next(iter(manager.metadata.values()))
    meta.downloaded_at = (datetime.now() - timedelta(days=365)).isoformat()
    manager._save_metadata()
    Path(pin.path).unlink()
    reloaded = manager_for(tmp_path, pin, monkeypatch)
    assert reloaded.get_database("mirbase", "human") == result
    assert reloaded.get_reference_metadata("mirbase", "human")["acquisition_mode"] == "pinned_cache"


@pytest.mark.parametrize(
    "payload",
    [
        "",
        "<html>offline</html>",
        ">hsa-test\nINVALID\n",
        ">mmu-test\nACGU\n",
        ">hsa-test\n",
        ">hsa-test\nACGU\n>hsa-test\nUGCA\n",
        ">mmu-test hsa-fake\nACGU\n",
    ],
)
def test_bad_pin_never_downloads_or_installs(tmp_path, monkeypatch, payload):
    """Bad pin never downloads or installs."""
    pin = pin_file(tmp_path, payload)
    manager = manager_for(tmp_path, pin, monkeypatch)
    assert manager.get_database("mirbase", "human") is None
    assert manager.get_reference_metadata("mirbase", "human") is None
    assert not list(manager.cache_dir.glob("*.fa"))


def test_pin_checksum_mismatch_and_failed_reimport_preserve_old_cache(tmp_path, monkeypatch):
    """Pin checksum mismatch and failed reimport preserve old cache."""
    pin = pin_file(tmp_path)
    manager = manager_for(tmp_path, pin, monkeypatch)
    result = manager.get_database("mirbase", "human")
    before = result.read_bytes()
    metadata = manager.metadata_file.read_bytes()
    Path(pin.path).write_text(">hsa-other\nAAAA\n")
    assert manager.get_database("mirbase", "human", force_refresh=True) is None
    assert manager.get_reference_metadata("mirbase", "human") is None
    assert result.read_bytes() == before
    assert manager.metadata_file.read_bytes() == metadata


def test_corrupt_pin_cache_is_not_reused(tmp_path, monkeypatch):
    """Corrupt pin cache is not reused."""
    pin = pin_file(tmp_path)
    manager = manager_for(tmp_path, pin, monkeypatch)
    result = manager.get_database("mirbase", "human")
    result.write_text(">hsa-corrupt\nAAAA\n")
    Path(pin.path).unlink()
    assert manager.get_database("mirbase", "human") is None
    assert manager.get_reference_metadata("mirbase", "human") is None


@pytest.mark.parametrize("failure", [None, "", "<html>down</html>", ">mmu-test\nACGU\n", ">hsa-test\nINVALID\n"])
def test_failed_remote_refresh_preserves_previous_reference(tmp_path, monkeypatch, failure):
    """Failed remote refresh preserves previous reference."""
    manager = MiRNADatabaseManager(cache_dir=tmp_path / "cache")
    monkeypatch.setattr(manager, "_download_file", lambda *_args: ">hsa-test\nACGUACGU\n")
    result = manager.get_database("mirbase", "human")
    before = result.read_bytes()
    metadata = manager.metadata_file.read_bytes()
    monkeypatch.setattr(manager, "_download_file", lambda *_args: failure)
    assert manager.get_database("mirbase", "human", force_refresh=True) is None
    assert result.read_bytes() == before
    assert manager.metadata_file.read_bytes() == metadata
    assert manager.get_reference_metadata("mirbase", "human") is None


def test_interrupted_atomic_refresh_preserves_previous_reference(tmp_path, monkeypatch):
    """Interrupted atomic refresh preserves previous reference."""
    manager = MiRNADatabaseManager(cache_dir=tmp_path / "cache")
    monkeypatch.setattr(manager, "_download_file", lambda *_args: ">hsa-test\nACGUACGU\n")
    result = manager.get_database("mirbase", "human")
    before = result.read_bytes()
    files_before = set(manager.cache_dir.iterdir())

    def interrupted(*args):
        raise OSError("interrupted replacement")

    monkeypatch.setattr(Path, "replace", interrupted)
    assert manager.get_database("mirbase", "human", force_refresh=True) is None
    assert result.read_bytes() == before
    assert set(manager.cache_dir.iterdir()) == files_before


def test_required_union_never_substitutes_missing_pin(tmp_path, monkeypatch):
    """Required union never substitutes missing pin."""
    pin = pin_file(tmp_path, ">hsa-other\nCCCCAAAA\n")
    mirgenedb_path = tmp_path / "mirgenedb.fa"
    mirgenedb_path.write_text(">Hsa-Mir-test\nACGUACGU\n")
    mirgenedb_pin = PinnedMiRNAReference(mirgenedb_path, "3.0", hashlib.sha256(mirgenedb_path.read_bytes()).hexdigest())
    manager = MiRNADatabaseManager(
        cache_dir=tmp_path / "cache",
        pinned_references={"mirbase": {"human": pin}, "mirgenedb": {"human": mirgenedb_pin}},
    )
    monkeypatch.setattr(manager, "_download_file", lambda *_args: pytest.fail("No network for fully pinned union"))
    result = manager.get_database("mirgenedb+mirbase", "human")
    assert result is not None
    provenance = manager.get_reference_metadata("mirgenedb+mirbase", "human")
    assert {entry["source"] for entry in provenance["sources"]} == {"mirbase", "mirgenedb"}
    assert provenance["filtered_sha256"] == hashlib.sha256(result.read_bytes()).hexdigest()
    before = result.read_bytes()
    Path(pin.path).unlink()
    assert manager.get_database("mirgenedb+mirbase", "human", force_refresh=True) is None
    assert manager.get_reference_metadata("mirgenedb+mirbase", "human") is None
    assert result.read_bytes() == before


def test_manifest_relative_paths_and_alias_duplicates(tmp_path):
    """Manifest relative paths and alias duplicates."""
    manifest_path = tmp_path / "references.json"
    entry = {"source": "mirbase", "species": "human", "release": "22.1", "sha256": "a" * 64, "path": "mature.fa"}
    manifest_path.write_text(json.dumps({"schema_version": 1, "references": [entry]}))
    pins = load_pinned_references(manifest_path)
    assert pins["mirbase"]["human"].path == tmp_path / "mature.fa"
    manifest_path.write_text(
        json.dumps({"schema_version": 1, "references": [entry, {**entry, "species": "Homo sapiens"}]})
    )
    with pytest.raises(ValueError, match="Duplicate"):
        load_pinned_references(manifest_path)


@pytest.mark.parametrize(
    "manifest",
    [
        {},
        {"schema_version": True, "references": []},
        {"schema_version": 1, "references": []},
        {"schema_version": 1, "references": [{}]},
        {
            "schema_version": 1,
            "references": [
                {"source": "mirbase", "species": "human", "release": "22", "sha256": "bad", "path": "mature.fa"}
            ],
        },
    ],
)
def test_invalid_manifest_fails_with_actionable_value_error(tmp_path, manifest):
    """Invalid manifest fails with actionable value error."""
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        load_pinned_references(path)


@pytest.mark.parametrize(
    "error",
    [TimeoutError("timed out"), HTTPError("https://www.mirbase.org/download/mature.fa", 503, "unavailable", {}, None)],
)
def test_live_failures_have_bounded_attempts_and_source_diagnostics(tmp_path, monkeypatch, caplog, error):
    """Timeout and HTTP failures respect configured attempt and socket limits."""
    manager = MiRNADatabaseManager(
        cache_dir=tmp_path / "cache", download_timeout=7, download_attempts=3, retry_backoff_seconds=0.5
    )
    attempts = []
    delays = []

    def unavailable(request, *, timeout):
        attempts.append((request.full_url, timeout))
        raise error

    monkeypatch.setattr("urllib.request.urlopen", unavailable)
    monkeypatch.setattr("sirnaforge.data.mirna_manager.time.sleep", delays.append)
    assert manager.get_database("mirbase", "human") is None
    assert attempts == [("https://www.mirbase.org/download/mature.fa", 7)] * 3
    assert delays == [0.5, 1.0]
    assert "after 3 attempts" in caplog.text
    assert "verified pinned mature FASTA" in caplog.text


def test_live_retry_stops_after_success(tmp_path, monkeypatch):
    """A transient failure needs no further attempt after a successful response."""
    manager = MiRNADatabaseManager(cache_dir=tmp_path / "cache", download_attempts=3, retry_backoff_seconds=0)
    responses = iter([None, ">hsa-test\nACGUACGU\n"])
    attempts = []

    def download(source, timeout):
        attempts.append(timeout)
        return next(responses)

    monkeypatch.setattr(manager, "_download_file", download)
    assert manager.get_database("mirbase", "human") is not None
    assert attempts == [30, 30]


@pytest.mark.parametrize(
    "options",
    [
        {"download_timeout": 0},
        {"download_timeout": 301},
        {"download_attempts": 0},
        {"download_attempts": 4},
        {"retry_backoff_seconds": -1},
        {"retry_backoff_seconds": float("inf")},
    ],
)
def test_download_budget_configuration_is_bounded(tmp_path, options):
    """Reject unbounded, invalid, or excessively large retry configuration."""
    with pytest.raises(ValueError):
        MiRNADatabaseManager(cache_dir=tmp_path / "cache", **options)
