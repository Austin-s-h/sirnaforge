"""Download pinned Ensembl cDNA releases for the four-species native mapping review."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

SPECIES = {"human": "homo_sapiens", "mouse": "mus_musculus", "rat": "rattus_norvegicus", "macaque": "macaca_mulatta"}


def fetch(release: int, species: str, output: Path) -> Path:
    """Verify upstream Unix checksums and record SHA-256 for both downloaded and expanded bytes."""
    output.mkdir(parents=True, exist_ok=True)
    base = f"https://ftp.ensembl.org/pub/release-{release}/fasta/{SPECIES[species]}/cdna/"

    def download(url: str, target: Path) -> None:
        subprocess.run(
            [
                "curl",
                "--fail",
                "--location",
                "--max-time",
                "300",
                "--retry",
                "2",
                "--silent",
                "--show-error",
                url,
                "-o",
                str(target),
            ],
            check=True,
        )

    listing = output / f"{species}-listing.html"
    download(base, listing)
    names = re.findall(r'href="([A-Za-z0-9_.-]+\.cdna\.all\.fa\.gz)"', listing.read_text())
    if len(set(names)) != 1:
        raise ValueError(f"Expected one all-cDNA FASTA in pinned release directory {base}")
    name = names[0]
    checksums = output / f"{species}-CHECKSUMS"
    download(base + "CHECKSUMS", checksums)
    declared = [line.split() for line in checksums.read_text().splitlines() if line.split()[-1:] == [name]]
    if len(declared) != 1:
        raise ValueError(f"No unique upstream CHECKSUMS entry for {name}")
    compressed = output / name
    if not compressed.exists():
        download(base + name, compressed)
    observed = subprocess.run(["sum", str(compressed)], check=True, text=True, capture_output=True).stdout.split()
    if observed[:2] != declared[0][:2]:
        raise ValueError(f"Upstream Unix checksum mismatch for {name}; remove incomplete download and retry")
    fasta = output / f"{species}.fa"
    digest = hashlib.sha256()
    with gzip.open(compressed, "rb") as source, fasta.open("wb") as target:
        while chunk := source.read(1024 * 1024):
            target.write(chunk)
            digest.update(chunk)
    assembly = name.split(".", 1)[1].removesuffix(".cdna.all.fa.gz")
    record = {
        "species": species,
        "assembly": assembly,
        "release": str(release),
        "fasta": str(fasta.resolve()),
        "sha256": digest.hexdigest(),
        "source_url": base + name,
        "download_sha256": hashlib.sha256(compressed.read_bytes()).hexdigest(),
        "upstream_unix_sum": declared[0][:2],
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }
    manifest = output / f"{species}_manifest.json"
    manifest.write_text(json.dumps(record, indent=2))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=int, required=True)
    parser.add_argument("--species", choices=SPECIES, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.release < 1:
        parser.error("--release must be positive")
    print(fetch(args.release, args.species, args.output))
