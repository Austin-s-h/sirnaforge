#!/usr/bin/env python3
"""Fetch pinned public benchmark review data and verify every file against its recorded SHA-256."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "tests/data/benchmarks/public_sources.json"


def fetch_file(entry: dict[str, str], destination: Path, *, verify_only: bool = False) -> None:
    """Install a verified download atomically, or verify an existing local copy."""
    target = destination / entry["path"]
    if not target.resolve().is_relative_to(destination.resolve()):
        raise ValueError(f"Path escapes destination: {entry['path']}")
    if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == entry["sha256"]:
        return
    if verify_only:
        raise ValueError(f"Missing or changed file: {target}")
    if not entry["url"].startswith("https://"):
        raise ValueError("Public benchmark downloads must use HTTPS")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with (
            urllib.request.urlopen(entry["url"], timeout=60) as response,
            tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle,
        ):
            temporary = Path(handle.name)
            digest = hashlib.sha256()
            while chunk := response.read(1024 * 1024):
                digest.update(chunk)
                handle.write(chunk)
        if digest.hexdigest() != entry["sha256"]:
            raise ValueError(f"Checksum mismatch: {entry['path']}; existing file retained")
        temporary.replace(target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    """Download the selected collections or verify an existing download directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "tests/data/external")
    parser.add_argument("--dataset", choices=["oligogym", "aso_atlas_2", "txbench_oligo_discovery"], action="append")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    selected = [entry for entry in manifest["files"] if not args.dataset or entry["path"].split("/")[0] in args.dataset]
    for entry in selected:
        fetch_file(entry, args.output_dir, verify_only=args.verify_only)
    print(f"Verified {len(selected)} public benchmark files in {args.output_dir}")


if __name__ == "__main__":
    main()
