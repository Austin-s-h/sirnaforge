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
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.data.transcriptome_manager import TranscriptomeManager, TranscriptomeSource

# Biotype/canonical mix chosen so `protein_coding` and `canonical_only` each keep a
# different pair of records -- compounded filtering therefore keeps only ENST001.
MIXED_BIOTYPE_FASTA = (
    ">ENST001 gene_biotype:protein_coding canonical:1\nACGTACGTAC\n"
    ">ENST002 gene_biotype:lncRNA canonical:1\nTTTTGGGGCC\n"
    ">ENST003 gene_biotype:protein_coding\nGGCCAATTGG\n"
)

# Same release-agnostic URL, next Ensembl release: ENST003 retired, ENST004 added.
NEXT_RELEASE_FASTA = (
    ">ENST001 gene_biotype:protein_coding canonical:1\nACGTACGTAC\n"
    ">ENST002 gene_biotype:lncRNA canonical:1\nTTTTGGGGCC\n"
    ">ENST004 gene_biotype:protein_coding\nAATTCCGGAA\n"
)

# Every file bwa-mem2 index writes for a prefix, including the .0123 suffix that
# validate_index_files() does not require but the tool (and our own Nextflow stub)
# does produce.
BWA_INDEX_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")

OLD_CONTENT = ">ENST001 release_113\nACGTACGTAC\n"
NEW_CONTENT = ">ENST999 release_114\nTTTTGGGGCC\n"


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

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(OLD_CONTENT, calls))
    first = manager.get_transcriptome("ensembl_human_cdna")
    assert first is not None
    assert Path(f"{first['index']}.amb").read_text(encoding="utf-8") == OLD_CONTENT

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))
    second = manager.get_transcriptome("ensembl_human_cdna", force_refresh=True)

    assert second is not None
    assert second["fasta"].read_text(encoding="utf-8") == NEW_CONTENT
    assert Path(f"{second['index']}.amb").read_text(encoding="utf-8") == NEW_CONTENT


@pytest.mark.unit
def test_expired_cache_refresh_rebuilds_stale_index(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """TTL-driven re-download must invalidate the index too, without force_refresh."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(OLD_CONTENT, calls))
    assert manager.get_transcriptome("ensembl_human_cdna") is not None

    cache_key = manager.sources["ensembl_human_cdna"].cache_key()
    manager.metadata[cache_key].downloaded_at = "2000-01-01T00:00:00"

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))
    refreshed = manager.get_transcriptome("ensembl_human_cdna")

    assert refreshed is not None
    assert Path(f"{refreshed['index']}.amb").read_text(encoding="utf-8") == NEW_CONTENT


@pytest.mark.unit
def test_unchanged_content_keeps_existing_index(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Re-downloading identical bytes must reuse the index instead of rebuilding it."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(OLD_CONTENT, calls))
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

    # The bookkeeping, not just the files, has to survive the re-record: a metadata
    # entry that forgot its index_path falls back to the cache-key-derived prefix,
    # which is how a stale index gets adopted in the first place.
    cache_key = manager.sources["ensembl_human_cdna"].cache_key()
    assert (manager.metadata[cache_key].extra or {}).get("index_path") == str(first["index"])


@pytest.mark.unit
def test_filtered_transcriptome_refilters_when_base_release_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A new base release must invalidate the filtered FASTA *and* its index.

    The filtered entry has no TTL and no checksum of its own to compare, so without a
    recorded link to the base bytes it is served forever - screening against a
    transcript set that no longer exists.
    """
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))
    first = manager.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"])
    assert first is not None
    assert "ENST003" in first["fasta"].read_text(encoding="utf-8")

    # Expire the base entry so the TTL path (not force_refresh) re-downloads it.
    base_key = manager.sources["ensembl_human_cdna"].cache_key()
    manager.metadata[base_key].downloaded_at = "2000-01-01T00:00:00"
    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(NEXT_RELEASE_FASTA, calls))

    second = manager.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"])

    assert second is not None
    filtered_text = second["fasta"].read_text(encoding="utf-8")
    assert "ENST004" in filtered_text
    assert "ENST003" not in filtered_text
    assert Path(f"{second['index']}.amb").read_text(encoding="utf-8") == filtered_text


@pytest.mark.unit
def test_filtered_transcriptome_reuses_cache_when_base_is_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An unchanged base must not trigger re-filtering or an index rebuild."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []
    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))

    first = manager.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"])
    assert first is not None

    builds: list[Path] = []
    real_build = manager._build_index

    def _counting_build(fasta_path: Path, index_prefix: Path) -> bool:
        builds.append(index_prefix)
        return real_build(fasta_path, index_prefix)

    monkeypatch.setattr(manager, "_build_index", _counting_build)
    second = manager.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"])

    assert second is not None
    assert second["fasta"] == first["fasta"]
    assert second["index"] == first["index"]
    assert builds == []


@pytest.mark.unit
def test_filter_qualified_uri_survives_a_reload(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The persisted URI index must keep filtered and base identities separate.

    The filtered entry owns its own filter-qualified URI and the base entry owns the
    bare download URI; if a reload drops the filtered mapping while
    `_load_metadata` keeps backfilling it, the separation only holds in-process.
    """
    source = TranscriptomeManager.SOURCES["ensembl_human_cdna"]
    first = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    calls: list[Path] = []
    monkeypatch.setattr(first, "_download_to_path", _fake_downloader(MIXED_BIOTYPE_FASTA, calls))
    assert first.get_filtered_transcriptome("ensembl_human_cdna", ["protein_coding"], build_index=False) is not None
    first._save_metadata()

    second = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)

    base_key = source.cache_key()
    filtered_key = f"{base_key}_protein_coding"
    assert second.uri_index.get(source.url) == base_key
    assert second.uri_index.get(f"{source.url}#filters=protein_coding") == filtered_key


@pytest.mark.unit
@pytest.mark.parametrize("dedupe", [True, False])
def test_file_already_in_cache_dir_rebuilds_index_when_content_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dedupe: bool
) -> None:
    """A rewritten FASTA inside the cache dir must not inherit the old index.

    `_handle_cached_file` names the index from the file *stem* while the cache key is
    content-derived, so the new content lands on the previous content's index prefix.
    """
    manager = TranscriptomeManager(cache_dir=tmp_path, local_content_dedupe=dedupe)
    _fake_index_builder(manager, monkeypatch)
    in_cache = manager.cache_dir / "transcripts.fa"
    in_cache.write_text(OLD_CONTENT, encoding="utf-8")

    first = manager.get_custom_transcriptome(in_cache)
    assert first is not None
    assert Path(f"{first['index']}.amb").read_text(encoding="utf-8") == OLD_CONTENT

    in_cache.write_text(NEW_CONTENT, encoding="utf-8")
    second = manager.get_custom_transcriptome(in_cache)

    assert second is not None
    assert Path(f"{second['index']}.amb").read_text(encoding="utf-8") == NEW_CONTENT


@pytest.mark.unit
def test_edited_local_fasta_refreshes_cached_copy_and_index(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Without content dedupe the cache key is path-derived, so edits must still land.

    `_cache_local_file` only copied when the destination was missing, so a user editing
    their own FASTA silently kept designing against the first version they ever ran.
    """
    manager = TranscriptomeManager(cache_dir=tmp_path / "cache", local_content_dedupe=False)
    _fake_index_builder(manager, monkeypatch)
    local_fasta = tmp_path / "mine.fa"
    local_fasta.write_text(OLD_CONTENT, encoding="utf-8")

    first = manager.get_custom_transcriptome(local_fasta)
    assert first is not None
    assert first["fasta"].read_text(encoding="utf-8") == OLD_CONTENT

    local_fasta.write_text(NEW_CONTENT, encoding="utf-8")
    second = manager.get_custom_transcriptome(local_fasta)

    assert second is not None
    assert second["fasta"] == first["fasta"]  # path-derived key: same cache file
    assert second["fasta"].read_text(encoding="utf-8") == NEW_CONTENT
    assert Path(f"{second['index']}.amb").read_text(encoding="utf-8") == NEW_CONTENT


@pytest.mark.unit
def test_clean_cache_removes_every_bwa_mem2_index_file(tmp_path: Path) -> None:
    """Index cleanup must cover .0123, which bwa-mem2 writes but validation ignores."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    fasta = manager.cache_dir / "retired.fa"
    fasta.write_text(OLD_CONTENT, encoding="utf-8")
    index_prefix = manager.cache_dir / "retired_index"
    index_prefix.touch()
    for suffix in BWA_INDEX_SUFFIXES:
        (manager.cache_dir / f"{index_prefix.name}{suffix}").write_text("x", encoding="utf-8")

    manager._record_cache_entry(
        "retiredkey",
        TranscriptomeSource(name="retired", url=str(fasta), species="custom"),
        fasta,
        extra={"index_path": str(index_prefix)},
        downloaded_at="2000-01-01T00:00:00",
    )
    manager.clean_cache(older_than_days=1)

    assert sorted(p.name for p in manager.cache_dir.glob("retired_index*")) == []


@pytest.mark.unit
def test_index_removal_survives_unwritable_cache_parent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A read-only cache parent must not turn index cleanup into a crashed run."""
    manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=False)
    index_prefix = manager.cache_dir / "locked_index"
    index_prefix.touch()
    for suffix in BWA_INDEX_SUFFIXES:
        (manager.cache_dir / f"{index_prefix.name}{suffix}").write_text("x", encoding="utf-8")

    real_unlink = Path.unlink

    def _refuse(self: Path, missing_ok: bool = False) -> None:
        if self.name.startswith("locked_index"):
            raise PermissionError(f"read-only cache: {self}")
        real_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", _refuse)

    with caplog.at_level(logging.WARNING, logger="sirnaforge.data.transcriptome_manager"):
        manager._remove_index_files(index_prefix)

    assert "locked_index" in caplog.text


def _refuse_to_unlink(monkeypatch: pytest.MonkeyPatch, index_prefix: Path) -> None:
    """Make every file sharing `index_prefix` undeletable, as a read-only cache dir would."""
    real_unlink = Path.unlink

    def _refuse(self: Path, missing_ok: bool = False) -> None:
        if self.name.startswith(index_prefix.name):
            raise PermissionError(f"read-only cache: {self}")
        real_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", _refuse)


@pytest.mark.unit
def test_undeletable_stale_index_is_not_served(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A stale index we failed to delete must not come back as a cached index.

    Reached through the public `get_transcriptome`, because the escalation and the
    "using cached index" decision live in different methods: logging the failure and
    then handing the same files out anyway is the corruption this branch exists to stop.
    """
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(OLD_CONTENT, calls))
    first = manager.get_transcriptome("ensembl_human_cdna")
    assert first is not None
    index_prefix = first["index"]
    assert Path(f"{index_prefix}.amb").read_text(encoding="utf-8") == OLD_CONTENT

    # The cache directory has become read-only (shared/bind-mounted cache), so neither
    # deleting the stale index nor writing a replacement can succeed.
    _refuse_to_unlink(monkeypatch, index_prefix)
    monkeypatch.setattr(manager, "_build_index", lambda fasta_path, prefix: False)  # noqa: ARG005
    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))

    with caplog.at_level(logging.INFO, logger="sirnaforge.data.transcriptome_manager"):
        second = manager.get_transcriptome("ensembl_human_cdna", force_refresh=True)

    assert second is not None
    assert second["fasta"].read_text(encoding="utf-8") == NEW_CONTENT
    # Scenario sanity check: the stale files really did survive, still holding OLD_CONTENT.
    assert Path(f"{index_prefix}.amb").read_text(encoding="utf-8") == OLD_CONTENT
    # The point of the test: an index built from replaced content is never handed back...
    assert "index" not in second
    # ...and the logs stop contradicting themselves about it.
    assert "Using cached BWA-MEM2 index" not in caplog.text


@pytest.mark.unit
def test_undeletable_stale_index_stays_refused_in_a_later_process(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The refusal must be persisted, not just remembered by the manager that saw it."""
    first_manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(first_manager, monkeypatch)
    calls: list[Path] = []

    monkeypatch.setattr(first_manager, "_download_to_path", _fake_downloader(OLD_CONTENT, calls))
    first = first_manager.get_transcriptome("ensembl_human_cdna")
    assert first is not None
    index_prefix = first["index"]

    _refuse_to_unlink(monkeypatch, index_prefix)
    monkeypatch.setattr(first_manager, "_build_index", lambda fasta_path, prefix: False)  # noqa: ARG005
    monkeypatch.setattr(first_manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))
    assert first_manager.get_transcriptome("ensembl_human_cdna", force_refresh=True) is not None

    # A later run finds the cache entry valid and in-TTL, so it never re-reconciles;
    # only the persisted refusal keeps it from adopting the leftovers.
    second_manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(second_manager, monkeypatch)
    monkeypatch.setattr(second_manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))

    reused = second_manager.get_transcriptome("ensembl_human_cdna")

    assert reused is not None
    assert reused["fasta"].read_text(encoding="utf-8") == NEW_CONTENT
    assert "index" not in reused


@pytest.mark.unit
def test_index_is_rebuilt_once_the_undeletable_leftovers_are_gone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Refusing the prefix must be recoverable, not a permanent loss of indexing."""
    manager = TranscriptomeManager(cache_dir=tmp_path)
    _fake_index_builder(manager, monkeypatch)
    calls: list[Path] = []

    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(OLD_CONTENT, calls))
    first = manager.get_transcriptome("ensembl_human_cdna")
    assert first is not None
    index_prefix = first["index"]

    with monkeypatch.context() as locked:
        _refuse_to_unlink(locked, index_prefix)
        locked.setattr(manager, "_build_index", lambda fasta_path, prefix: False)  # noqa: ARG005
        locked.setattr(manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))
        assert "index" not in (manager.get_transcriptome("ensembl_human_cdna", force_refresh=True) or {})

    # The operator follows the escalation and deletes the leftovers by hand.
    for stale in tmp_path.glob(f"{index_prefix.name}*"):
        stale.unlink()

    _fake_index_builder(manager, monkeypatch)
    monkeypatch.setattr(manager, "_download_to_path", _fake_downloader(NEW_CONTENT, calls))
    recovered = manager.get_transcriptome("ensembl_human_cdna", force_refresh=True)

    assert recovered is not None
    assert Path(f"{recovered['index']}.amb").read_text(encoding="utf-8") == NEW_CONTENT
