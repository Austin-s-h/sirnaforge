"""Regression tests for transcriptome cache identity and index invalidation.

Both defects covered here silently corrupt off-target screening:
- a filtered FASTA claiming the base source URI, so later unfiltered requests
  screen against a filtered reference;
- a BWA-MEM2 index surviving a re-download, so hits come from the previous
  reference content.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.data.transcriptome_manager import TranscriptomeManager

# Biotype/canonical mix chosen so `protein_coding` and `canonical_only` each keep a
# different pair of records -- compounded filtering therefore keeps only ENST001.
MIXED_BIOTYPE_FASTA = (
    ">ENST001 gene_biotype:protein_coding canonical:1\nACGTACGTAC\n"
    ">ENST002 gene_biotype:lncRNA canonical:1\nTTTTGGGGCC\n"
    ">ENST003 gene_biotype:protein_coding\nGGCCAATTGG\n"
)

BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt.2bit.64", ".pac")


def _record_count(fasta: Path) -> int:
    return fasta.read_text(encoding="utf-8").count(">")


def _fake_downloader(payload: str, calls: list[Path]) -> Any:
    def _download(source: Any, destination: Path, timeout: int = 600) -> bool:  # noqa: ARG001
        calls.append(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(payload, encoding="utf-8")
        return True

    return _download


@pytest.mark.unit
def test_filtered_transcriptome_does_not_hijack_unfiltered_cache_entry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An unfiltered request must keep returning the full transcriptome after a filtered run."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    calls: list[Path] = []
    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))

    unfiltered = manager.get_transcriptome("ensembl_human_cdna", build_index=False)
    assert unfiltered is not None
    assert _record_count(unfiltered["fasta"]) == 3

    filtered = manager.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"], build_index=False)
    assert filtered is not None
    assert _record_count(filtered["fasta"]) == 2

    again = manager.get_transcriptome("ensembl_human_cdna", build_index=False)
    assert again is not None
    assert again["fasta"] == unfiltered["fasta"]
    assert _record_count(again["fasta"]) == 3


@pytest.mark.unit
def test_filtered_cache_poisoning_does_not_persist_across_processes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A fresh manager must not inherit a filtered entry as the base source's cache entry."""
    first = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    calls: list[Path] = []
    monkeypatch.setattr(first, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))

    assert first.get_transcriptome("ensembl_human_cdna", build_index=False) is not None
    assert first.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"], build_index=False) is not None
    first._save_metadata()

    second = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    monkeypatch.setattr(second, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))

    result = second.get_transcriptome("ensembl_human_cdna", build_index=False)

    assert result is not None
    assert _record_count(result["fasta"]) == 3


@pytest.mark.unit
def test_repeated_filtering_does_not_compound(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Filtering twice must filter the full reference again, not the previous filtered output."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    calls: list[Path] = []
    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))

    protein_coding = manager.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"], build_index=False)
    assert protein_coding is not None
    assert _record_count(protein_coding["fasta"]) == 2

    canonical = manager.get_filtered_transcriptome("ensembl_human_cdna", ["canonical_only"], build_index=False)

    # canonical_only must be applied to the full reference, so the lncRNA survives.
    assert canonical is not None
    assert _record_count(canonical["fasta"]) == 2
    assert "ENST002" in canonical["fasta"].read_text(encoding="utf-8")


@pytest.mark.unit
def test_legacy_poisoned_metadata_is_repaired_on_load(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Metadata written before the fix must not keep serving a filtered FASTA."""
    source = TranscriptomeManager.SOURCES["ensembl_human_cdna"]
    base_key = source.cache_key()
    base_fasta = tmp_path / f"{base_key}.fa"
    base_fasta.write_text(MIXED_BIOTYPE_FASTA, encoding="utf-8")
    filtered_key = f"{base_key}_protein_coding"
    filtered_fasta = tmp_path / f"{filtered_key}.fa"
    filtered_fasta.write_text(
        ">ENST001 gene_biotype:protein_coding canonical:1\nACGTACGTAC\n"
        ">ENST003 gene_biotype:protein_coding\nGGCCAATTGG\n",
        encoding="utf-8",
    )

    def _entry(path: Path, extra: dict[str, Any] | None) -> dict[str, Any]:
        return {
            "source": {
                "name": source.name,
                # The pre-fix bug: the filtered entry carries the bare source URL.
                "url": source.url,
                "species": source.species,
                "format": "fasta",
                "compressed": True,
                "description": source.description,
            },
            "downloaded_at": datetime.now().isoformat(),
            "file_size": path.stat().st_size,
            "checksum": hashlib.md5(path.read_bytes()).hexdigest(),
            "file_path": str(path),
            "version": "1.0",
            "extra": extra,
        }

    (tmp_path / "cache_metadata.json").write_text(
        json.dumps(
            {
                "version": "2.0",
                "entries": {
                    base_key: _entry(base_fasta, None),
                    filtered_key: _entry(filtered_fasta, {"filters": ["protein_coding"], "kept_count": 2}),
                },
                "uri_index": {source.url: filtered_key},
            }
        ),
        encoding="utf-8",
    )

    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)

    def _fail_download(*args: Any, **kwargs: Any) -> bool:  # noqa: ARG001
        raise AssertionError("download should not be needed; the base artifact is cached")

    monkeypatch.setattr(manager, "_download_to_path", _fail_download)

    result = manager.get_transcriptome("ensembl_human_cdna", build_index=False)

    assert result is not None
    assert result["fasta"] == base_fasta
    assert _record_count(result["fasta"]) == 3


def _fake_index_builder(manager: TranscriptomeManager, monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace real bwa-mem2 index building with marker files stamped from the FASTA."""

    def _build(fasta_path: Path, index_prefix: Path) -> bool:
        stamp = fasta_path.read_text(encoding="utf-8")
        index_prefix.parent.mkdir(parents=True, exist_ok=True)
        for suffix in BWA_INDEX_SUFFIXES:
            (index_prefix.parent / f"{index_prefix.name}{suffix}").write_text(stamp, encoding="utf-8")
        return True

    monkeypatch.setattr(manager, "_build_index", _build)


@pytest.mark.unit
def test_refreshed_transcriptome_rebuilds_stale_index(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A re-download must invalidate the index built from the previous content."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    old_fasta = ">ENST001 release_113\nACGTACGTAC\n"
    new_fasta = ">ENST999 release_114\nTTTTGGGGCC\n"

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(old_fasta, calls))
    first = manager.get_transcriptome("ensembl_human_cdna")
    assert first is not None
    assert Path(f"{first['index']}.amb").read_text(encoding="utf-8") == old_fasta

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(new_fasta, calls))
    second = manager.get_transcriptome("ensembl_human_cdna", force_refresh=True)

    assert second is not None
    assert second["fasta"].read_text(encoding="utf-8") == new_fasta
    assert Path(f"{second['index']}.amb").read_text(encoding="utf-8") == new_fasta


@pytest.mark.unit
def test_expired_cache_refresh_rebuilds_stale_index(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """TTL-driven re-download must invalidate the index too, without force_refresh."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    old_fasta = ">ENST001 release_113\nACGTACGTAC\n"
    new_fasta = ">ENST999 release_114\nTTTTGGGGCC\n"

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(old_fasta, calls))
    assert manager.get_transcriptome("ensembl_human_cdna") is not None

    cache_key = manager.sources["ensembl_human_cdna"].cache_key()
    manager.metadata[cache_key].downloaded_at = "2000-01-01T00:00:00"

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(new_fasta, calls))
    refreshed = manager.get_transcriptome("ensembl_human_cdna")

    assert refreshed is not None
    assert Path(f"{refreshed['index']}.amb").read_text(encoding="utf-8") == new_fasta


@pytest.mark.unit
def test_unchanged_content_keeps_existing_index(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Re-downloading identical bytes must reuse the index instead of rebuilding it."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []
    payload = ">ENST001 release_113\nACGTACGTAC\n"

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(payload, calls))
    first = manager.get_transcriptome("ensembl_human_cdna")
    assert first is not None

    builds: list[Path] = []
    real_build = manager._build_index

    def _counting_build(fasta_path: Path, index_prefix: Path) -> bool:
        builds.append(index_prefix)
        return real_build(fasta_path, index_prefix)

    monkeypatch.setattr(manager, "_build_index", _counting_build)
    second = manager.get_transcriptome("ensembl_human_cdna", force_refresh=True)

    assert second is not None
    assert second["index"] == first["index"]
    assert builds == []
