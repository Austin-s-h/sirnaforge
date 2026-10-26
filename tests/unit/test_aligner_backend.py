"""Backend selection and index separation for transcriptome screening."""

from pathlib import Path

import pytest

from sirnaforge.core.aligner_backend import normalize_aligner
from sirnaforge.core.off_target import BwaAnalyzer, validate_index_files
from sirnaforge.data.reference_manager import CacheMetadata
from sirnaforge.data.transcriptome_manager import TranscriptomeManager


@pytest.mark.unit
def test_aligner_selector_is_explicit() -> None:
    """Aligner selector is explicit."""
    assert normalize_aligner("bwa-mem2") == "bwa-mem2"
    assert normalize_aligner("bwa-mem3") == "bwa-mem3"
    with pytest.raises(ValueError, match="Unsupported aligner"):
        normalize_aligner("bwa")


@pytest.mark.unit
def test_mem3_has_backend_specific_cache_prefix_and_ignores_legacy_mem2_index(tmp_path: Path) -> None:
    """Mem3 has backend specific cache prefix and ignores legacy mem2 index."""
    manager2 = TranscriptomeManager(cache_dir=tmp_path, aligner="bwa-mem2")
    manager3 = TranscriptomeManager(cache_dir=tmp_path, aligner="bwa-mem3")
    metadata = CacheMetadata(
        file_path=str(tmp_path / "reference.fa"),
        file_size=1,
        checksum="abc",
        downloaded_at="2026-01-01T00:00:00",
        source=manager2.sources["ensembl_human_cdna"],
        extra={"index_path": str(tmp_path / "reference_index")},
    )

    assert manager2._index_prefix_for(tmp_path / "reference") == tmp_path / "reference_index"
    assert manager3._index_prefix_for(tmp_path / "reference") == tmp_path / "reference_bwa_mem3_index"
    assert manager2._get_index_path(metadata) == tmp_path / "reference_index"
    assert manager3._get_index_path(metadata) is None


@pytest.mark.unit
def test_backend_index_validation_requires_complete_nonempty_files(tmp_path: Path) -> None:
    """Backend index validation requires complete nonempty files."""
    prefix = tmp_path / "sample_bwa_mem3_index"
    for suffix in (".amb", ".ann", ".bwt.2bit.64", ".pac"):
        (tmp_path / f"{prefix.name}{suffix}").write_bytes(b"index")

    assert validate_index_files(prefix, tool="bwa-mem3")
    assert validate_index_files(prefix, tool="bwa-mem2")
    (tmp_path / f"{prefix.name}.pac").write_bytes(b"")
    assert not validate_index_files(prefix, tool="bwa-mem3")


@pytest.mark.unit
def test_analyzer_builds_selected_backend_command() -> None:
    """Analyzer builds selected backend command."""
    analyzer = BwaAnalyzer("/tmp/mem3_index", mode="transcriptome", aligner="bwa-mem3")
    command = analyzer._build_bwa_command("/usr/bin/bwa-mem3", "/tmp/query.fa")

    assert command[0:2] == ["/usr/bin/bwa-mem3", "mem"]
    assert command[-2:] == ["/tmp/mem3_index", "/tmp/query.fa"]
