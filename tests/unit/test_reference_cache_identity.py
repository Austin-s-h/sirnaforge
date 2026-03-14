"""Regression tests for cache identity and dedupe behavior."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.data.annotation_manager import AnnotationManager
from sirnaforge.data.transcriptome_manager import TranscriptomeManager
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


def test_resolve_cache_subdir_warns_on_mixed_roots(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Resolving one subdir to multiple roots should emit a warning."""
    subdir = "mixed-root-warning-test"
    root_a = tmp_path / "root_a"
    root_b = tmp_path / "root_b"

    resolve_cache_subdir(subdir, override=root_a / subdir)
    resolve_cache_subdir(subdir, override=root_b / subdir)

    assert any("Multiple cache roots detected" in record.message for record in caplog.records)
