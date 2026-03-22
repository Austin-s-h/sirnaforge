"""Regression tests for cache identity and dedupe behavior."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.data.annotation_manager import AnnotationManager
from sirnaforge.data.reference_manager import CacheMetadata
from sirnaforge.data.transcriptome_manager import TranscriptomeManager, TranscriptomeSource
from sirnaforge.utils.cache_utils import resolve_cache_subdir


def test_custom_url_uses_uri_identity_across_names(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Same URI should reuse one cached artifact even with different cache names."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    url = "https://example.org/resources/transcripts.fa.gz"

    download_calls: list[Path] = []

    def _fake_download_to_path(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        download_calls.append(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(">tx1\nAUGC\n", encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _fake_download_to_path)

    first = manager.get_custom_transcriptome(url, build_index=False, cache_name="run_a")
    second = manager.get_custom_transcriptome(url, build_index=False, cache_name="run_b")

    assert first is not None
    assert second is not None
    assert first["fasta"] == second["fasta"]
    assert len(download_calls) == 1


def test_uri_index_backfills_from_legacy_metadata(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Legacy metadata-only files should populate URI index on load."""
    cache_key = "abc123def456"
    cached_file = tmp_path / f"{cache_key}.fa"
    cached_file.write_text(">tx1\nAAAA\n", encoding="utf-8")

    legacy: dict[str, dict[str, Any]] = {
        cache_key: {
            "source": {
                "name": "legacy",
                "url": "https://example.org/legacy.fa.gz",
                "species": "custom",
                "format": "fasta",
                "compressed": True,
                "description": "legacy",
            },
            "downloaded_at": "2099-01-01T00:00:00",
            "file_size": cached_file.stat().st_size,
            "checksum": hashlib.md5(cached_file.read_bytes()).hexdigest(),
            "file_path": str(cached_file),
            "version": "1.0",
            "extra": None,
        }
    }
    (tmp_path / "cache_metadata.json").write_text(json.dumps(legacy), encoding="utf-8")

    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)

    # Ensure we fail if a re-download attempt happens.
    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("download should not be called for legacy cache hit")

    monkeypatch.setattr(manager, "_download_to_path", _fail_download)

    result = manager.get_custom_transcriptome("https://example.org/legacy.fa.gz", build_index=False, cache_name="new")

    assert result is not None
    assert result["fasta"] == cached_file


def test_local_content_dedupe_reuses_identical_files(tmp_path: Path) -> None:
    """Optional content dedupe should reuse the first cached local artifact."""
    input_a = tmp_path / "input_a.fa"
    input_b = tmp_path / "input_b.fa"
    input_a.write_text(">tx1\nCCCC\n", encoding="utf-8")
    input_b.write_text(">tx1\nCCCC\n", encoding="utf-8")

    cache_dir = tmp_path / "cache"
    manager = TranscriptomeManager(cache_dir=cache_dir, auto_build_indices=False, local_content_dedupe=True)

    first = manager.get_custom_transcriptome(input_a, build_index=False, cache_name="a")
    second = manager.get_custom_transcriptome(input_b, build_index=False, cache_name="b")

    assert first is not None
    assert second is not None
    assert first["fasta"] == second["fasta"]

    fasta_files = list(cache_dir.glob("*.fa"))
    assert len(fasta_files) == 1


def test_custom_annotation_uses_uri_identity_across_names(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Same annotation URI should be downloaded once even with different names."""
    manager = AnnotationManager(cache_dir=tmp_path)
    url = "https://example.org/resources/annotation.gtf.gz"

    download_calls: list[Path] = []

    def _fake_download_to_path(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        download_calls.append(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text('chr1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _fake_download_to_path)

    first = manager.get_custom_annotation(url, cache_name="first")
    second = manager.get_custom_annotation(url, cache_name="second")

    assert first is not None
    assert second is not None
    assert first == second
    assert len(download_calls) == 1


def test_remote_transcriptome_recovers_from_existing_file_without_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Existing remote artifact should be reused even if metadata was never persisted."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    source = manager.sources["ensembl_human_cdna"]
    cache_file = tmp_path / f"{source.cache_key()}.fa"
    cache_file.write_text(">tx1\nACGT\n", encoding="utf-8")

    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("download should not be called when cache file exists")

    monkeypatch.setattr(manager, "_download_to_path", _fail_download)

    result = manager.get_transcriptome("ensembl_human_cdna", build_index=False)

    assert result is not None
    assert result["fasta"] == cache_file
    assert source.url in manager.uri_index


def test_remote_custom_annotation_recovers_from_existing_file_without_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Existing annotation artifact should be reused even when metadata is absent."""
    manager = AnnotationManager(cache_dir=tmp_path)
    url = "https://example.org/resources/annotation.gtf.gz"
    source = TranscriptomeSource(name="annotation", url=url, species="custom", format="gtf", compressed=False)
    cache_file = tmp_path / f"{source.cache_key()}.gtf.gz"
    cache_file.write_text('chr1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")

    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("download should not be called when annotation cache file exists")

    monkeypatch.setattr(manager, "_download_to_path", _fail_download)

    result = manager.get_custom_annotation(url, cache_name="annotation")

    assert result == cache_file
    assert url in manager.uri_index


def test_annotation_manager_init_falls_back_from_unwritable_cache_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """AnnotationManager should auto-fallback when explicit cache_dir is not writable."""
    monkeypatch.setattr(AnnotationManager, "_cache_dir_is_writable", staticmethod(lambda _: False))

    manager = AnnotationManager(cache_dir=Path("/cache"))

    assert manager.cache_dir != Path("/cache")
    assert manager.cache_dir.exists()


def test_annotation_manager_rehomes_stale_unwritable_metadata_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Invalid metadata file paths should be re-homed to the active cache directory."""
    manager = AnnotationManager(cache_dir=tmp_path)
    url = "https://example.org/resources/annotation.gtf.gz"
    source = TranscriptomeSource(name="annotation", url=url, species="custom", format="gtf", compressed=False)
    cache_key = source.cache_key()
    manager.metadata[cache_key] = CacheMetadata(
        source=source,
        downloaded_at=datetime.now().isoformat(),
        file_size=1,
        checksum="deadbeef",
        file_path="/cache",
    )

    captured_destinations: list[Path] = []

    def _fake_download_to_path(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        captured_destinations.append(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text('chr1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _fake_download_to_path)
    monkeypatch.setattr(manager, "_cache_dir_is_writable", staticmethod(lambda p: str(p) != "/"))

    result = manager.get_custom_annotation(url)

    assert result is not None
    assert result.parent == manager.cache_dir
    assert captured_destinations[0].parent == manager.cache_dir


def test_resolve_cache_subdir_warns_on_mixed_roots(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Resolving one subdir to multiple roots should emit a warning."""
    subdir = "mixed-root-warning-test"
    root_a = tmp_path / "root_a"
    root_b = tmp_path / "root_b"

    resolve_cache_subdir(subdir, override=root_a / subdir)
    resolve_cache_subdir(subdir, override=root_b / subdir)

    assert any("Multiple cache roots detected" in record.message for record in caplog.records)


def test_force_refresh_bypasses_recovered_remote_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """force_refresh should trigger a new download even when cache file exists."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    source = manager.sources["ensembl_human_cdna"]
    cache_file = tmp_path / f"{source.cache_key()}.fa"
    cache_file.write_text(">tx1\nOLD\n", encoding="utf-8")

    downloads = 0

    def _fake_download_to_path(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        nonlocal downloads
        downloads += 1
        destination.write_text(">tx1\nNEW\n", encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _fake_download_to_path)

    result = manager.get_transcriptome("ensembl_human_cdna", force_refresh=True, build_index=False)

    assert result is not None
    assert downloads == 1
    assert cache_file.read_text(encoding="utf-8") == ">tx1\nNEW\n"


def test_stale_uri_index_is_ignored_and_recomputed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Stale URI index entries should not prevent deterministic key reuse."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    source = manager.sources["ensembl_human_cdna"]
    cache_key = source.cache_key()
    cache_file = tmp_path / f"{cache_key}.fa"
    cache_file.write_text(">tx1\nACGT\n", encoding="utf-8")

    manager.uri_index[source.url] = "missing_key"

    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("download should not be called for existing cached file")

    monkeypatch.setattr(manager, "_download_to_path", _fail_download)

    result = manager.get_transcriptome("ensembl_human_cdna", build_index=False)

    assert result is not None
    assert result["fasta"] == cache_file
    assert manager.uri_index[source.url] == cache_key


def test_cache_persists_across_manager_restarts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A new manager instance should reuse metadata persisted by a prior instance."""
    url = "https://example.org/resources/restart.fa.gz"

    first_manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    first_downloads = 0

    def _first_download(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        nonlocal first_downloads
        first_downloads += 1
        destination.write_text(">tx1\nPERSIST\n", encoding="utf-8")
        return True

    monkeypatch.setattr(first_manager, "_download_to_path", _first_download)
    first = first_manager.get_custom_transcriptome(url, build_index=False, cache_name="persist")

    assert first is not None
    assert first_downloads == 1

    second_manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)

    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("second manager should hit persisted cache")

    monkeypatch.setattr(second_manager, "_download_to_path", _fail_download)
    second = second_manager.get_custom_transcriptome(url, build_index=False, cache_name="persist_again")

    assert second is not None
    assert second["fasta"] == first["fasta"]


def test_malformed_metadata_file_does_not_block_remote_cache_recovery(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Bad metadata JSON should not force redownload when artifact already exists."""
    source = TranscriptomeManager.SOURCES["ensembl_human_cdna"]
    cache_file = tmp_path / f"{source.cache_key()}.fa"
    cache_file.write_text(">tx1\nRECOVER\n", encoding="utf-8")
    (tmp_path / "cache_metadata.json").write_text("{not-json", encoding="utf-8")

    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)

    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("download should not be called when existing artifact is recoverable")

    monkeypatch.setattr(manager, "_download_to_path", _fail_download)

    result = manager.get_transcriptome("ensembl_human_cdna", build_index=False)

    assert result is not None
    assert result["fasta"] == cache_file
    assert source.url in manager.uri_index


def test_remote_uri_identity_is_exact_string_no_normalization(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Different URI strings should not dedupe even if they differ only in query values."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    url_a = "https://example.org/resources/transcripts.fa.gz?version=1"
    url_b = "https://example.org/resources/transcripts.fa.gz?version=2"

    download_calls: list[Path] = []

    def _fake_download_to_path(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        download_calls.append(destination)
        destination.write_text(">tx1\nEXACT\n", encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _fake_download_to_path)

    first = manager.get_custom_transcriptome(url_a, build_index=False, cache_name="v1")
    second = manager.get_custom_transcriptome(url_b, build_index=False, cache_name="v2")

    assert first is not None
    assert second is not None
    assert first["fasta"] != second["fasta"]
    assert len(download_calls) == 2


def test_corrupted_cached_file_triggers_redownload(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Checksum mismatch should invalidate metadata cache and force redownload."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    url = "https://example.org/resources/transcripts.fa.gz"

    first_downloads = 0

    def _first_download(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        nonlocal first_downloads
        first_downloads += 1
        destination.write_text(">tx1\nORIGINAL\n", encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _first_download)
    first = manager.get_custom_transcriptome(url, build_index=False, cache_name="run")
    assert first is not None
    assert first_downloads == 1

    # Corrupt file after metadata/checksum has been persisted in-memory.
    first["fasta"].write_text(">tx1\nCORRUPTED\n", encoding="utf-8")

    second_downloads = 0

    def _second_download(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        nonlocal second_downloads
        second_downloads += 1
        destination.write_text(">tx1\nREPAIRED\n", encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_download_to_path", _second_download)
    second = manager.get_custom_transcriptome(url, build_index=False, cache_name="run")

    assert second is not None
    assert second_downloads == 1
    assert second["fasta"].read_text(encoding="utf-8") == ">tx1\nREPAIRED\n"


def test_clear_cache_clears_uri_index(tmp_path: Path) -> None:
    """Clearing cache should remove URI index entries as well as metadata."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    source = manager.sources["ensembl_human_cdna"]
    cache_file = tmp_path / f"{source.cache_key()}.fa"
    cache_file.write_text(">tx1\nACGT\n", encoding="utf-8")
    manager._record_cache_entry(source.cache_key(), source, cache_file)

    assert manager.uri_index
    assert manager.metadata

    cleared = manager.clear_cache(confirm=True)

    assert cleared["files_deleted"] >= 1
    assert manager.metadata == {}
    assert manager.uri_index == {}
