"""Tests for packaging the four-species benchmark artifacts."""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest

from scripts import package_four_species_mappings


@pytest.mark.unit
def test_split_assay_manifest_references_only_existing_parts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Split assay archives replace the deleted CSV path with their part list."""
    monkeypatch.setattr(package_four_species_mappings, "_MAX_PART_BYTES", 8)
    assay = tmp_path / "assay_observations.csv"
    original = b"guide,score\nACGU,0.7\n"
    assay.write_bytes(original)
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps({"assay_observations_csv": {"path": assay.name}, "species": {}}))

    package_four_species_mappings._package_manifest(tmp_path, manifest_path, include_assay=True)

    artifact = json.loads(manifest_path.read_text())["assay_observations_csv"]
    assert "path" not in artifact
    assert artifact["parts"]
    compressed = b"".join((tmp_path / part).read_bytes() for part in artifact["parts"])
    assert hashlib.sha256(compressed).hexdigest() == artifact["sha256"]
    assert hashlib.sha256(gzip.decompress(compressed)).hexdigest() == artifact["uncompressed_sha256"]


@pytest.mark.unit
def test_packaging_removes_only_manifest_csvs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Unrelated CSV evidence in the artifact directory is preserved."""
    species_dir = tmp_path / "species" / "human"
    species_dir.mkdir(parents=True)
    (tmp_path / "assay.csv").write_text("guide,score\nACGU,0.7\n")
    (species_dir / "mappings.csv").write_text("guide,transcript\nACGU,tx1\n")
    (species_dir / "observations.csv").write_text("guide,score\nACGU,0.7\n")
    unrelated = tmp_path / "review_notes.csv"
    unrelated.write_text("note\nretain me\n")
    species = {
        "human": {
            "mappings_csv": "human/mappings.csv",
            "observations_csv": "human/observations.csv",
        }
    }
    (tmp_path / "manifest.json").write_text(
        json.dumps({"assay_observations_csv": {"path": "assay.csv"}, "species": species})
    )
    (tmp_path / "species" / "manifest.json").write_text(json.dumps({"species": species}))
    monkeypatch.setattr(sys, "argv", ["package_four_species_mappings", "--artifact-dir", str(tmp_path)])

    package_four_species_mappings.main()

    assert unrelated.read_text() == "note\nretain me\n"
    assert not (tmp_path / "assay.csv").exists()
    assert not (species_dir / "mappings.csv").exists()
    assert not (species_dir / "observations.csv").exists()
    assert (tmp_path / "assay.csv.gz").exists()
    assert (species_dir / "mappings.csv.gz").exists()
