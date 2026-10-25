"""Tests for packaging the four-species benchmark artifacts."""

from __future__ import annotations

import gzip
import hashlib
import json
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
