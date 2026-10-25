#!/usr/bin/env python3
"""Deterministically gzip release-116 mapping CSVs and update their manifests."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

_MAX_PART_BYTES = 480 * 1024


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _gzip_csv(path: Path) -> tuple[Path, str, str]:
    if path.suffix != ".csv":
        raise ValueError(f"expected an uncompressed CSV, got {path}")
    compressed = path.with_suffix(".csv.gz")
    with (
        path.open("rb") as source,
        compressed.open("wb") as raw_output,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as output,
    ):
        while chunk := source.read(1024 * 1024):
            output.write(chunk)
    return compressed, _sha256(compressed), _sha256(path)


def _split_if_large(compressed: Path) -> list[Path] | None:
    """Split a gzip byte stream into concatenable pieces below the repository file limit."""
    if compressed.stat().st_size <= _MAX_PART_BYTES:
        return None
    parts: list[Path] = []
    with compressed.open("rb") as stream:
        index = 1
        while chunk := stream.read(_MAX_PART_BYTES):
            part = compressed.with_name(f"{compressed.name}.part-{index:03d}")
            part.write_bytes(chunk)
            parts.append(part)
            index += 1
    compressed.unlink()
    return parts


def _package_manifest(
    artifact_root: Path,
    manifest_path: Path,
    *,
    include_assay: bool,
) -> None:
    manifest: dict[str, Any] = json.loads(manifest_path.read_text())
    if include_assay:
        artifact = manifest["assay_observations_csv"]
        source = artifact_root / artifact["path"]
        compressed, compressed_sha, csv_sha = _gzip_csv(source)
        parts = _split_if_large(compressed)
        artifact["sha256"] = compressed_sha
        artifact["uncompressed_sha256"] = csv_sha
        if parts:
            artifact.pop("path", None)
            artifact["parts"] = [part.relative_to(artifact_root).as_posix() for part in parts]
        else:
            artifact["path"] = compressed.relative_to(artifact_root).as_posix()
    for species in manifest["species"].values():
        for name in ("mappings", "observations"):
            path_key, sha_key = f"{name}_csv", f"{name}_sha256"
            source = artifact_root / "species" / species[path_key]
            compressed, compressed_sha, csv_sha = _gzip_csv(source)
            species[sha_key] = compressed_sha
            species[f"{name}_uncompressed_sha256"] = csv_sha
            parts = _split_if_large(compressed)
            if parts:
                species.pop(path_key)
                species[f"{path_key}_parts"] = [
                    part.relative_to(artifact_root / "species").as_posix() for part in parts
                ]
            else:
                species[path_key] = compressed.relative_to(artifact_root / "species").as_posix()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    """Package the generator's plain CSV outputs in place and preserve both digests."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.artifact_dir.resolve()
    _package_manifest(root, root / "manifest.json", include_assay=True)
    _package_manifest(root, root / "species" / "manifest.json", include_assay=False)
    for path in root.rglob("*.csv"):
        path.unlink()


if __name__ == "__main__":
    main()
