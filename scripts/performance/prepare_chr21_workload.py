"""Create a reproducible short-guide workload from the pinned chromosome-21 fixture."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def prepare(output: Path, transcript_count: int = 64) -> dict:
    """Verify the reference bytes, then select exact complementary 19–23 nt guides."""
    fixture = ROOT / "tests/container/data/human_grch38_chr21.cdna.fa.gz"
    provenance = json.loads(fixture.with_suffix(".json").read_text())
    compressed = fixture.read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == provenance["gzip_sha256"]
    reference_bytes = gzip.decompress(compressed)
    assert hashlib.sha256(reference_bytes).hexdigest() == provenance["fasta_sha256"]
    output.mkdir(parents=True, exist_ok=True)
    reference = output / "reference.fa"
    reference.write_bytes(reference_bytes)
    records = []
    name = ""
    chunks: list[str] = []
    for line in reference_bytes.decode().splitlines() + [">end"]:
        if line.startswith(">"):
            sequence = "".join(chunks).upper().replace("U", "T")
            if name and len(sequence) >= 46 and set(sequence) <= set("ACGT"):
                records.append((name, sequence))
                if len(records) == transcript_count:
                    break
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if len(records) != transcript_count:
        raise ValueError("Not enough unambiguous transcripts in the pinned fixture")
    queries = output / "guides.fa"
    with queries.open("w") as handle:
        for name, sequence in records:
            start = len(sequence) // 2 - 12
            for length in range(19, 24):
                site = sequence[start : start + length]
                guide = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                handle.write(f">{name}_site{start + 1}_len{length}\n{guide}\n")
    result = {
        "schema_version": 1,
        "reference": {"fixture": str(fixture.relative_to(ROOT)), "sha256": provenance["fasta_sha256"]},
        "queries": {
            "transcripts": len(records),
            "count": len(records) * 5,
            "lengths": list(range(19, 24)),
            "sha256": hashlib.sha256(queries.read_bytes()).hexdigest(),
        },
        "limitations": [
            "Chromosome 21 only; not a full transcriptome or a biological efficacy benchmark.",
            "Queries are exact antisense sites; mismatch/repeat stress cases use the separate synthetic fixture.",
        ],
    }
    (output / "workload.json").write_text(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--transcripts", type=int, default=64)
    args = parser.parse_args()
    if args.transcripts < 1:
        parser.error("--transcripts must be positive")
    print(json.dumps(prepare(args.output, args.transcripts), indent=2))
