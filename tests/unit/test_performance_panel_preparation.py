"""Tests for deterministic OligoGym performance panel inputs."""

from pathlib import Path

from scripts.performance.prepare_oligogym_panel import OBSERVATIONS, prepare


def test_oligogym_panel_fasta_is_pinned_and_reproducible(tmp_path: Path) -> None:
    """The tracked source data produce the same pinned FASTA on repeated runs."""
    first = prepare(OBSERVATIONS, tmp_path / "first.fa")
    second = prepare(OBSERVATIONS, tmp_path / "second.fa")

    assert first["observation_rows"] == 4_113
    assert first["distinct_guides"] == 3_395
    assert first["observations_sha256"] == "05a16b6a2ae4e025fcbd408ebbc5d547f85c72daf3764b7baef4f96ca32293ed"
    assert first["guide_fasta_sha256"] == "e80d5cae150fd38526f4eb41cc1225c7d91fc43fed54484a97d24d608cbd0a2b"
    assert second["guide_fasta_sha256"] == first["guide_fasta_sha256"]
    assert (tmp_path / "first.fa").read_bytes() == (tmp_path / "second.fa").read_bytes()
