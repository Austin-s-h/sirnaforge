"""Producer-version (schema) invalidation for on-disk cache artifacts.

The recorded md5 of a cache file is computed from that same file when it is
written, so it can only ever prove the bytes have not rotted since - never that a
buggy writer produced the right bytes in the first place. These tests cover the
one mechanism that can: the producer version recorded on write and compared on
load, both for files tracked in `cache_metadata.json` and for derived artifacts
(combined FASTAs, BWA-MEM2 indices) tracked by sidecar stamps.
"""

import json
import logging
import os
from datetime import datetime, timedelta
from pathlib import Path

import pytest
from Bio import SeqIO

from sirnaforge.data.mirna_manager import MiRNADatabaseManager, MiRNASource
from sirnaforge.data.reference_manager import CacheMetadata
from sirnaforge.data.transcriptome_manager import TranscriptomeManager, TranscriptomeSource
from sirnaforge.utils import cache_utils
from sirnaforge.utils.cache_utils import (
    ARTIFACT_MIRNA_COMBINED,
    ARTIFACT_TRANSCRIPTOME_INDEX,
    DIGEST_MODE_FULL,
    DIGEST_MODE_SAMPLED,
    LEGACY_PRODUCER_VERSION,
    artifact_stamp_path,
    current_producer_version,
    fingerprint_outputs,
    is_artifact_stamp_current,
    is_producer_version_current,
    read_artifact_stamp,
    write_artifact_stamp,
)

HSA_FASTA = ">hsa-miR-1-3p MIMAT0000416\nUGGAAUGUAAAGAAGUAUGUAU\n"
MMU_FASTA = ">mmu-miR-1a-3p MIMAT0000123\nUGGAAUGUAAAGAAGUAUGUAG\n"


def _mirna_entry(manager: MiRNADatabaseManager, cache_file: Path, version: str) -> str:
    """Register `cache_file` as a cache hit written by producer `version`."""
    source = MiRNASource(name="test", url="https://example.org/mature.fa", species="human")
    cache_key = source.cache_key()
    manager.metadata[cache_key] = CacheMetadata(
        source=source,
        downloaded_at=datetime.now().isoformat(),
        file_size=cache_file.stat().st_size,
        checksum=manager._compute_file_checksum(cache_file),
        file_path=str(cache_file),
        version=version,
    )
    return cache_key


class TestProducerVersionScope:
    """Invalidation is scoped to the artifact class whose producer moved."""

    def test_mirna_is_versioned_past_the_pre_versioning_writer(self) -> None:
        """MiRNA artifacts must be invalidated: their writer mis-parsed miRBase."""
        assert current_producer_version("mirna") != LEGACY_PRODUCER_VERSION
        assert current_producer_version(ARTIFACT_MIRNA_COMBINED) != LEGACY_PRODUCER_VERSION

    def test_bulk_downloads_are_not_invalidated_by_other_subsystems(self) -> None:
        """A user with a multi-GB transcriptome must not re-download it for a miRNA fix."""
        for artifact_class in ("transcriptomes", "genomes", "annotations"):
            assert current_producer_version(artifact_class) == LEGACY_PRODUCER_VERSION
            assert is_producer_version_current(artifact_class, None)

    def test_unregistered_class_never_self_invalidates(self) -> None:
        """Registering a class is what opts it in; an unknown name must be inert."""
        assert is_producer_version_current("some_future_cache", LEGACY_PRODUCER_VERSION)

    def test_missing_recorded_version_reads_as_pre_versioning(self) -> None:
        """An absent version (None, or NaN from a column that did not exist) is stale."""
        assert not is_producer_version_current(ARTIFACT_MIRNA_COMBINED, None)
        assert not is_producer_version_current(ARTIFACT_MIRNA_COMBINED, float("nan"))
        assert not is_producer_version_current(ARTIFACT_MIRNA_COMBINED, "")
        assert is_producer_version_current(ARTIFACT_MIRNA_COMBINED, current_producer_version(ARTIFACT_MIRNA_COMBINED))

    def test_only_classes_that_are_actually_checked_are_registered(self) -> None:
        """A registered class nobody validates reads as a guarantee we do not give.

        Every entry in the table must be reachable from a call site: either a
        `ReferenceManager.artifact_class` (its cache subdir) or an artifact-stamp
        constant. `variants` and `transcriptome_filtered` were registered without one
        and are gone; re-registering either without wiring it up must fail here.
        """
        checked = {"mirna", "transcriptomes", "genomes", "annotations"} | {
            ARTIFACT_MIRNA_COMBINED,
            ARTIFACT_TRANSCRIPTOME_INDEX,
        }
        assert set(cache_utils.PRODUCER_VERSIONS) == checked


class TestReferenceMetadataVersion:
    """`_is_cache_valid` must consult the recorded producer version."""

    def test_artifact_from_old_producer_is_a_miss(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """An artifact left by the pre-versioning writer must not be served."""
        manager = MiRNADatabaseManager(cache_dir=tmp_path)
        cache_file = tmp_path / "old.fa"
        cache_file.write_text(HSA_FASTA, encoding="utf-8")
        cache_key = _mirna_entry(manager, cache_file, LEGACY_PRODUCER_VERSION)

        with caplog.at_level(logging.WARNING):
            assert manager._is_cache_valid(cache_key) is False

        # The user has to be able to see why a re-download is happening.
        assert "Discarding stale cached mirna artifact" in caplog.text

    def test_ttl_valid_and_checksum_valid_but_wrong_version_is_not_served(self, tmp_path: Path) -> None:
        """The two pre-existing checks passing must not be enough on their own."""
        manager = MiRNADatabaseManager(cache_dir=tmp_path, cache_ttl_days=30)
        cache_file = tmp_path / "fresh.fa"
        cache_file.write_text(HSA_FASTA, encoding="utf-8")
        cache_key = _mirna_entry(manager, cache_file, "0.9")
        meta = manager.metadata[cache_key]

        # Both of the checks that existed before are satisfied.
        assert manager._compute_file_checksum(cache_file) == meta.checksum
        assert datetime.now() - datetime.fromisoformat(meta.downloaded_at) < manager.cache_ttl

        assert manager._is_cache_valid(cache_key) is False
        assert manager.stale_producer_entries() == [cache_key]

    def test_artifact_from_current_producer_is_reused(self, tmp_path: Path) -> None:
        """Versioning must not turn every load into a re-download."""
        manager = MiRNADatabaseManager(cache_dir=tmp_path)
        cache_file = tmp_path / "current.fa"
        cache_file.write_text(HSA_FASTA, encoding="utf-8")
        cache_key = _mirna_entry(manager, cache_file, current_producer_version("mirna"))

        assert manager._is_cache_valid(cache_key) is True
        assert manager.stale_producer_entries() == []

    def test_recorded_entries_carry_the_current_version(self, tmp_path: Path) -> None:
        """Writing an artifact records the producer version, without callers asking."""
        manager = MiRNADatabaseManager(cache_dir=tmp_path)
        cache_file = tmp_path / "recorded.fa"
        cache_file.write_text(HSA_FASTA, encoding="utf-8")
        source = MiRNASource(name="test", url="https://example.org/x.fa", species="human")

        meta = manager._record_cache_entry(source.cache_key(), source, cache_file)

        assert meta.version == current_producer_version("mirna")
        assert manager._is_cache_valid(source.cache_key()) is True

    def test_stale_artifact_is_regenerated_end_to_end(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A stale artifact is re-downloaded once and then served from cache again."""
        manager = MiRNADatabaseManager(cache_dir=tmp_path)
        downloads: list[str] = []

        def fake_download(source: MiRNASource, timeout: int = 600) -> str:  # noqa: ARG001
            downloads.append(source.url)
            return HSA_FASTA

        monkeypatch.setattr(manager, "_download_file", fake_download)

        first = manager.get_database("mirgenedb", "human")
        assert first is not None
        assert len(downloads) == 1

        # Rewrite the entry as if the file on disk had been produced by the buggy
        # writer: same bytes, same checksum, still inside TTL.
        cache_key = manager.SOURCES["mirgenedb"]["hsa"].cache_key()
        manager.metadata[cache_key].version = LEGACY_PRODUCER_VERSION

        with caplog.at_level(logging.WARNING):
            second = manager.get_database("mirgenedb", "human")

        assert second == first
        assert len(downloads) == 2, "stale-producer artifact should have been re-downloaded"
        assert manager.metadata[cache_key].version == current_producer_version("mirna")
        assert "Discarding stale cached mirna artifact" in caplog.text

        third = manager.get_database("mirgenedb", "human")
        assert third == first
        assert len(downloads) == 2, "a current-version artifact must be served from cache"


class TestCombinedDatabaseStamp:
    """`combined_*.fa` is derived and previously reused on mtime alone."""

    @staticmethod
    def _manager_with_sources(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> tuple[MiRNADatabaseManager, dict[str, Path]]:
        manager = MiRNADatabaseManager(cache_dir=tmp_path)
        files = {
            "mirgenedb": tmp_path / "mirgenedb_hsa.fa",
            "mirbase": tmp_path / "mirbase_hsa.fa",
        }
        files["mirgenedb"].write_text(HSA_FASTA, encoding="utf-8")
        files["mirbase"].write_text(MMU_FASTA, encoding="utf-8")

        def fake_get_database(source_name: str, species: str, force_refresh: bool = False) -> Path:  # noqa: ARG001
            return files[source_name]

        monkeypatch.setattr(manager, "get_database", fake_get_database)
        return manager, files

    @staticmethod
    def _rebuilds(caplog: pytest.LogCaptureFixture) -> int:
        """How many times the combining step ran.

        The obvious probe - append a marker to the artifact and see whether it
        survives - cannot be used any more, because appending to a stamped artifact is
        exactly the post-stamp rewrite the output fingerprint now rejects. It would
        make every rebuild assertion below pass for the wrong reason, so reuse is
        measured by whether the producer ran at all.
        """
        return caplog.text.count("🔄 Combining")

    def test_combined_database_is_reused_when_stamp_matches(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """An unchanged combined FASTA is still served from cache."""
        caplog.set_level(logging.INFO)
        manager, _ = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None
        assert artifact_stamp_path(combined).exists()
        assert self._rebuilds(caplog) == 1

        assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined
        assert self._rebuilds(caplog) == 1, "an untouched artifact must be reused, not rebuilt"
        assert "Using existing combined database" in caplog.text

    def test_combined_database_is_rebuilt_when_a_source_changes_behind_its_mtime(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Rewritten source content with an unchanged mtime fooled the old check."""
        caplog.set_level(logging.INFO)
        manager, files = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None

        source = files["mirgenedb"]
        stat = source.stat()
        source.write_text(HSA_FASTA.replace("UGGAAUGUAAAGAAGUAUGUAU", "AAAAAAAAAAAAAAAAAAAAAA"), encoding="utf-8")
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        assert source.stat().st_mtime_ns <= combined.stat().st_mtime_ns

        assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined
        assert self._rebuilds(caplog) == 2
        assert "AAAAAAAAAAAAAAAAAAAAAA" in combined.read_text(encoding="utf-8")

    def test_unstamped_combined_database_is_rebuilt(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A combined FASTA left by an older release cannot vouch for itself."""
        caplog.set_level(logging.INFO)
        manager, _ = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None
        artifact_stamp_path(combined).unlink()

        assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined
        assert self._rebuilds(caplog) == 2

    def test_combined_database_is_rebuilt_when_producer_version_moves(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Bumping the combined-FASTA producer regenerates the artifact."""
        caplog.set_level(logging.INFO)
        manager, _ = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None

        monkeypatch.setitem(cache_utils.PRODUCER_VERSIONS, ARTIFACT_MIRNA_COMBINED, "99.0")

        assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined
        assert self._rebuilds(caplog) == 2

    def test_truncated_combined_database_is_regenerated(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A file truncated *after* stamping was served behind a "current" stamp.

        This is the failure input-only fingerprints cannot see: the producer version,
        the source digests and the TTL are all still correct, because none of them
        describe the artifact's own bytes.
        """
        manager, _ = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None
        assert len(list(SeqIO.parse(combined, "fasta"))) == 2

        text = combined.read_text(encoding="utf-8")
        stat = combined.stat()
        combined.write_text(text[: text.index("\n>") + 1], encoding="utf-8")
        os.utime(combined, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        assert len(list(SeqIO.parse(combined, "fasta"))) == 1

        with caplog.at_level(logging.WARNING):
            assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined

        assert len(list(SeqIO.parse(combined, "fasta"))) == 2, "a truncated artifact must be rebuilt, not served"
        assert "truncated or rewritten" in caplog.text

    def test_combined_database_corrupted_without_changing_its_length_is_regenerated(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Same byte count, same mtime: only a digest of the output can catch this."""
        manager, _ = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None
        original = combined.read_text(encoding="utf-8")

        stat = combined.stat()
        combined.write_text(original.replace("UGGAAUGUAAAGAAGUAUGUAU", "NNNNNNNNNNNNNNNNNNNNNN"), encoding="utf-8")
        os.utime(combined, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        assert combined.stat().st_size == stat.st_size
        assert combined.stat().st_mtime_ns == stat.st_mtime_ns

        with caplog.at_level(logging.WARNING):
            assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined

        assert combined.read_text(encoding="utf-8") == original
        assert "corrupted after it was stamped" in caplog.text

    def test_expired_combined_database_is_rebuilt(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The old mtime comparison had no TTL at all."""
        caplog.set_level(logging.INFO)
        manager, _ = self._manager_with_sources(tmp_path, monkeypatch)
        combined = manager.get_combined_database(["mirgenedb", "mirbase"], "human")
        assert combined is not None
        stamp = artifact_stamp_path(combined)
        aged = stamp.read_text(encoding="utf-8").replace(
            datetime.now().isoformat()[:10], (datetime.now() - timedelta(days=365)).isoformat()[:10]
        )
        stamp.write_text(aged, encoding="utf-8")

        assert manager.get_combined_database(["mirgenedb", "mirbase"], "human") == combined
        assert self._rebuilds(caplog) == 2


class TestTranscriptomeIndexStamp:
    """A BWA-MEM2 index must not outlive the content it was built from."""

    # What bwa-mem2 writes beside the prefix; the real `_is_index_complete` requires
    # all four to exist and be non-empty.
    MEMBER_EXTENSIONS = (".amb", ".ann", ".bwt.2bit.64", ".pac")

    @classmethod
    def _manager(cls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[TranscriptomeManager, list[Path]]:
        """Manager whose index build writes real member files.

        The stand-in has to produce files, not just a flag: an index is a multi-file
        artifact, and the point of the output fingerprint is that each member's bytes
        are checked. `_is_index_complete` is left unpatched so completeness and
        integrity are exercised as the two separate checks they are.
        """
        manager = TranscriptomeManager(cache_dir=tmp_path, auto_build_indices=True, local_content_dedupe=False)
        builds: list[Path] = []

        def fake_build(fasta: Path, prefix: Path) -> bool:
            builds.append(prefix)
            payload = f"bwa-mem2 index of {fasta.read_text(encoding='utf-8')}"
            for ext in cls.MEMBER_EXTENSIONS:
                (prefix.parent / f"{prefix.name}{ext}").write_text(payload, encoding="utf-8")
            return True

        monkeypatch.setattr(manager, "_build_index", fake_build)
        return manager, builds

    @staticmethod
    def _record(manager: TranscriptomeManager, cached: Path) -> CacheMetadata:
        source = TranscriptomeSource(name="ref", url=str(cached), species="custom")
        return manager._record_cache_entry("refkey", source, cached)

    def test_index_is_reused_while_content_is_unchanged(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Stamping must not cause an index rebuild on every run."""
        manager, builds = self._manager(tmp_path, monkeypatch)
        cached = manager.cache_dir / "ref.fa"
        cached.write_text(">t1\nACGT\n", encoding="utf-8")
        self._record(manager, cached)
        prefix = manager.cache_dir / "refkey_index"

        manager._prepare_result_with_index(cached, prefix, "refkey", True)
        manager._prepare_result_with_index(cached, prefix, "refkey", True)

        assert len(builds) == 1

    def test_index_is_discarded_when_the_content_it_indexed_is_replaced(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """An index whose source FASTA was replaced must be rebuilt."""
        manager, builds = self._manager(tmp_path, monkeypatch)
        cached = manager.cache_dir / "ref.fa"
        cached.write_text(">t1\nACGT\n", encoding="utf-8")
        self._record(manager, cached)
        prefix = manager.cache_dir / "refkey_index"
        manager._prepare_result_with_index(cached, prefix, "refkey", True)

        # New upstream release lands on the same cache key and the same index prefix.
        cached.write_text(">t2\nTTTTTTTT\n", encoding="utf-8")
        self._record(manager, cached)

        with caplog.at_level(logging.WARNING):
            manager._prepare_result_with_index(cached, prefix, "refkey", True)

        assert len(builds) == 2, "index built from replaced content must not be reused"
        assert "transcriptome_index" in caplog.text

    def test_legacy_index_built_after_the_cached_content_is_adopted(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Guard against over-invalidation: adopting costs nothing, rebuilding costs an hour."""
        manager, builds = self._manager(tmp_path, monkeypatch)
        cached = manager.cache_dir / "ref.fa"
        cached.write_text(">t1\nACGT\n", encoding="utf-8")
        meta = self._record(manager, cached)
        prefix = manager.cache_dir / "refkey_index"
        manager._prepare_result_with_index(cached, prefix, "refkey", True)
        assert len(builds) == 1

        # Simulate an index from a release that stamped nothing.
        artifact_stamp_path(prefix).unlink()
        meta.extra = {"index_path": str(prefix), "index_built_at": datetime.now().isoformat()}

        manager._prepare_result_with_index(cached, prefix, "refkey", True)

        assert len(builds) == 1
        assert artifact_stamp_path(prefix).exists(), "adopted index should be stamped for next time"

    def test_legacy_index_that_predates_the_cached_content_is_rebuilt(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """An undatable or older index cannot be shown to match, so it is rebuilt."""
        manager, builds = self._manager(tmp_path, monkeypatch)
        cached = manager.cache_dir / "ref.fa"
        cached.write_text(">t1\nACGT\n", encoding="utf-8")
        meta = self._record(manager, cached)
        prefix = manager.cache_dir / "refkey_index"
        manager._prepare_result_with_index(cached, prefix, "refkey", True)

        artifact_stamp_path(prefix).unlink()
        meta.extra = {
            "index_path": str(prefix),
            "index_built_at": (datetime.now() - timedelta(days=7)).isoformat(),
        }

        manager._prepare_result_with_index(cached, prefix, "refkey", True)

        assert len(builds) == 2

    def test_index_is_rebuilt_when_a_member_file_is_truncated_after_stamping(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A half-written index member must be caught, not aligned against.

        `_is_index_complete` only asks whether the members exist and are non-empty, so
        an interrupted copy or a full disk leaves an index it happily accepts while the
        input fingerprint still matches.
        """
        manager, builds = self._manager(tmp_path, monkeypatch)
        cached = manager.cache_dir / "ref.fa"
        cached.write_text(">t1\nACGT\n", encoding="utf-8")
        self._record(manager, cached)
        prefix = manager.cache_dir / "refkey_index"
        manager._prepare_result_with_index(cached, prefix, "refkey", True)
        assert len(builds) == 1

        member = prefix.parent / f"{prefix.name}.bwt.2bit.64"
        intact = member.read_text(encoding="utf-8")
        stat = member.stat()
        member.write_text(intact[: len(intact) // 2], encoding="utf-8")
        os.utime(member, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        assert manager._is_index_complete(prefix), "the completeness helper cannot see this"

        with caplog.at_level(logging.WARNING):
            manager._prepare_result_with_index(cached, prefix, "refkey", True)

        assert len(builds) == 2, "a truncated index member must force a rebuild"
        assert "truncated or rewritten" in caplog.text
        assert member.read_text(encoding="utf-8") == intact

    def test_stamp_alone_rejects_an_incomplete_index(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Every member is recorded, so the stamp no longer needs to be taken on trust.

        `_is_index_complete` would also catch a missing member; the point here is that
        the stamp is asked directly and answers "no" on its own, which is what lets
        `is_artifact_stamp_current` be documented as a sufficient reuse gate.
        """
        manager, _ = self._manager(tmp_path, monkeypatch)
        cached = manager.cache_dir / "ref.fa"
        cached.write_text(">t1\nACGT\n", encoding="utf-8")
        meta = self._record(manager, cached)
        prefix = manager.cache_dir / "refkey_index"
        manager._prepare_result_with_index(cached, prefix, "refkey", True)
        assert manager._index_is_trustworthy(meta, prefix)

        stamp = read_artifact_stamp(prefix)
        assert stamp is not None
        assert set(stamp["outputs"]) == {f"{prefix.name}{ext}" for ext in self.MEMBER_EXTENSIONS}

        (prefix.parent / f"{prefix.name}.pac").unlink()
        assert not manager._index_is_trustworthy(meta, prefix)


class TestArtifactOutputFingerprint:
    """A stamp must fingerprint the artifact's own bytes, not only its inputs."""

    def test_stamp_records_size_and_digest_of_the_artifact(self, tmp_path: Path) -> None:
        """Writing a stamp records the output's own size and digest, unasked."""
        artifact = tmp_path / "derived.fa"
        artifact.write_text(">x\nACGT\n", encoding="utf-8")
        write_artifact_stamp(ARTIFACT_MIRNA_COMBINED, artifact)

        stamp = read_artifact_stamp(artifact)
        assert stamp is not None
        recorded = stamp["outputs"]["derived.fa"]
        assert recorded["size"] == artifact.stat().st_size
        assert recorded["digest_mode"] == DIGEST_MODE_FULL
        assert recorded["digest"]

    def test_truncation_is_rejected_by_size_without_reading_the_bytes(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Size is the cheap pre-check; a mismatch must not pay for a digest at all."""
        artifact = tmp_path / "derived.fa"
        artifact.write_text(">x\nACGTACGT\n", encoding="utf-8")
        write_artifact_stamp(ARTIFACT_MIRNA_COMBINED, artifact)
        artifact.write_text(">x\n", encoding="utf-8")

        def explode(*_args: object, **_kwargs: object) -> str:
            raise AssertionError("size mismatch must short-circuit before any digest")

        monkeypatch.setattr(cache_utils, "_digest_output", explode)
        assert not is_artifact_stamp_current(ARTIFACT_MIRNA_COMBINED, artifact)

    def test_stamp_with_no_output_fingerprint_is_not_reusable(self, tmp_path: Path) -> None:
        """An artifact whose bytes were never fingerprinted cannot be vouched for."""
        artifact = tmp_path / "derived.fa"
        artifact.write_text(">x\nACGT\n", encoding="utf-8")
        write_artifact_stamp(ARTIFACT_MIRNA_COMBINED, artifact)

        stamp_path = artifact_stamp_path(artifact)
        stamp = json.loads(stamp_path.read_text(encoding="utf-8"))
        del stamp["outputs"]
        stamp_path.write_text(json.dumps(stamp), encoding="utf-8")

        assert not is_artifact_stamp_current(ARTIFACT_MIRNA_COMBINED, artifact)

    def test_large_outputs_are_sampled_instead_of_read_whole(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A multi-GB index must not be hashed end-to-end on every cache hit.

        With the budget scaled down to bytes, this pins the documented trade-off: size
        plus a head/tail window, so truncation and head/tail damage are caught and a
        same-length rewrite of the middle is knowingly not.
        """
        monkeypatch.setattr(cache_utils, "STAMP_FULL_DIGEST_MAX_BYTES", 16)
        monkeypatch.setattr(cache_utils, "STAMP_SAMPLE_BYTES", 4)
        artifact = tmp_path / "big.idx"
        artifact.write_bytes(b"HEAD" + b"m" * 32 + b"TAIL")
        write_artifact_stamp(ARTIFACT_TRANSCRIPTOME_INDEX, artifact)

        stamp = read_artifact_stamp(artifact)
        assert stamp is not None
        assert stamp["outputs"]["big.idx"]["digest_mode"] == DIGEST_MODE_SAMPLED
        assert is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, artifact)

        # Caught: truncation (size) and damage inside the sampled windows.
        artifact.write_bytes(b"HEAD" + b"m" * 32)
        assert not is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, artifact)
        artifact.write_bytes(b"head" + b"m" * 32 + b"TAIL")
        assert not is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, artifact)

        # Known gap, stated in the module docs: an equal-length rewrite between the
        # windows is invisible without a full read we refuse to pay for.
        artifact.write_bytes(b"HEAD" + b"X" * 32 + b"TAIL")
        assert is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, artifact)

    def test_fingerprint_skips_outputs_it_cannot_read(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """Bookkeeping must never abort a build; an unreadable output is just absent."""
        present = tmp_path / "present.fa"
        present.write_text(">x\nAC\n", encoding="utf-8")

        with caplog.at_level(logging.WARNING):
            fingerprints = fingerprint_outputs([present, tmp_path / "gone.fa"])

        assert set(fingerprints) == {"present.fa"}
        assert "Could not fingerprint cache output" in caplog.text


class TestArtifactStampPrimitives:
    """Direct coverage of the stamp helpers the sibling fixes call."""

    def test_stamp_round_trip(self, tmp_path: Path) -> None:
        """A stamp vouches for its own class and inputs, and nothing else."""
        artifact = tmp_path / "derived.fa"
        artifact.write_text(">x\nAC\n", encoding="utf-8")
        write_artifact_stamp(ARTIFACT_TRANSCRIPTOME_INDEX, artifact, inputs={"fasta": "abc"})

        assert is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, artifact, inputs={"fasta": "abc"})
        assert not is_artifact_stamp_current(ARTIFACT_TRANSCRIPTOME_INDEX, artifact, inputs={"fasta": "def"})
        # A stamp for one artifact class must not vouch for another.
        assert not is_artifact_stamp_current(ARTIFACT_MIRNA_COMBINED, artifact, inputs={"fasta": "abc"})

    def test_absent_artifact_is_quietly_a_miss(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        """A first run has nothing to discard, so it must not warn."""
        with caplog.at_level(logging.WARNING):
            assert not is_artifact_stamp_current(ARTIFACT_MIRNA_COMBINED, tmp_path / "never_built.fa")
        assert caplog.text == ""

    def test_stamp_survives_a_relocated_cache_directory(self, tmp_path: Path) -> None:
        """Container runs mount the cache elsewhere; that must not invalidate anything."""
        first = tmp_path / "a"
        first.mkdir()
        artifact = first / "derived.fa"
        artifact.write_text(">x\nAC\n", encoding="utf-8")
        write_artifact_stamp(ARTIFACT_MIRNA_COMBINED, artifact, inputs={"0:src.fa": "abc"})

        second = tmp_path / "b"
        second.mkdir()
        moved = second / "derived.fa"
        moved.write_text(artifact.read_text(encoding="utf-8"), encoding="utf-8")
        artifact_stamp_path(moved).write_text(artifact_stamp_path(artifact).read_text(encoding="utf-8"), "utf-8")

        assert is_artifact_stamp_current(ARTIFACT_MIRNA_COMBINED, moved, inputs={"0:src.fa": "abc"})
