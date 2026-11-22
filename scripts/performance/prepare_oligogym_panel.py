"""Build a deterministic BWA query FASTA from the pinned OligoGym observations."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OBSERVATIONS = ROOT / "tests/data/benchmarks/oligogym_four_species_release116/assay_observations.csv.gz"
EXPECTED_ROWS = 4_113
EXPECTED_UNIQUE_GUIDES = 3_395
EXPECTED_OBSERVATIONS_SHA256 = "05a16b6a2ae4e025fcbd408ebbc5d547f85c72daf3764b7baef4f96ca32293ed"
EXPECTED_GUIDE_FASTA_SHA256 = "e80d5cae150fd38526f4eb41cc1225c7d91fc43fed54484a97d24d608cbd0a2b"


def prepare(observations: Path, output: Path) -> dict:
    """Deduplicate full guide sequences in source order and record input/output hashes."""
    input_hash = hashlib.sha256(observations.read_bytes()).hexdigest()
    if input_hash != EXPECTED_OBSERVATIONS_SHA256:
        raise ValueError(f"Pinned OligoGym observation checksum mismatch: {input_hash}")
    with gzip.open(observations, "rt", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    guides = list(dict.fromkeys(row["full_guide_sequence"].strip().upper().replace("U", "T") for row in rows))
    if len(rows) != EXPECTED_ROWS or len(guides) != EXPECTED_UNIQUE_GUIDES:
        raise ValueError(f"Expected {EXPECTED_ROWS} observations and {EXPECTED_UNIQUE_GUIDES} unique guides")
    if any(len(guide) not in (20, 21) or set(guide) - set("ACGT") for guide in guides):
        raise ValueError("Panel contains a guide outside the measured 20–21 nt unambiguous sequence set")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii", newline="\n") as handle:
        for index, guide in enumerate(guides, start=1):
            handle.write(f">og_{index:05d}\n{guide}\n")

    output_hash = hashlib.sha256(output.read_bytes()).hexdigest()
    if output_hash != EXPECTED_GUIDE_FASTA_SHA256:
        raise ValueError(f"Generated guide FASTA checksum mismatch: {output_hash}")
    return {
        "schema_version": 1,
        "observations": str(observations),
        "observations_sha256": input_hash,
        "observation_rows": len(rows),
        "distinct_guides": len(guides),
        "guide_lengths": [20, 21],
        "guide_fasta": str(output),
        "guide_fasta_sha256": output_hash,
        "deduplication": "full_guide_sequence, first occurrence order; U converted to T",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=Path, default=OBSERVATIONS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.observations, args.output), indent=2))
