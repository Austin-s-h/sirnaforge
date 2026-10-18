#!/usr/bin/env python3
"""Prepare OligoGym siRNA observations as runnable siRNAforge design inputs.

OligoGym stores oligo duplexes, not full target transcripts. siRNAforge's ``design`` command
expects transcript FASTA and enumerates candidate sites from each transcript. This adapter therefore
creates an explicitly synthetic, neutral-context transcript per observation:

    A * 70 + reverse_complement(guide) + A * 70

The target site starts at position 71 (1-based), safely inside the default RNAplfold context. The
synthetic nature is recorded in every output row and must not be mistaken for native transcript
accessibility evidence. The guide is the second dot-separated FASTA component, corresponding to
the antisense strand in the OligoGym HELM records; this agrees with the existing Huesken fixture.

The generated FASTA files are useful for parameter smoke tests, ranking/filter comparisons, and
chemistry-independent sequence scoring. They are not a replacement for native transcript
sequences when evaluating target accessibility.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DIR = ROOT / "tests" / "data" / "external" / "oligogym"
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "data" / "benchmarks" / "oligogym"
FLANK_NT = 70
SITE_START_1BASED = FLANK_NT + 1

DATASETS = {
    "ichihara_2007_1": {
        "modality": "siRNA",
        "label_semantics": "percentage inhibition relative to control; higher is better",
        "direction": "higher_is_better",
        "source_url": "https://doi.org/10.1093/nar/gkm699",
    },
    "ichihara_2007_2": {
        "modality": "siRNA",
        "label_semantics": "percentage inhibition relative to control; higher is better",
        "direction": "higher_is_better",
        "source_url": "https://doi.org/10.1093/nar/gkm699",
    },
    "martinelli_2023_1": {
        "modality": "siRNA",
        "label_semantics": "percentage knockdown of target mRNA; higher is better",
        "direction": "higher_is_better",
        "source_url": "https://doi.org/10.1016/j.ygeno.2024.110815",
    },
    "shmushkovich_2018_1": {
        "modality": "siRNA",
        "label_semantics": "percentage remaining relative to control; lower is better",
        "direction": "lower_is_better",
        "source_url": "https://doi.org/10.1093/nar/gky745",
    },
}


def _dna(sequence: str) -> str:
    """Normalize OligoGym RNA spelling to the DNA spelling accepted by siRNAforge."""
    return "".join(sequence.upper().replace("U", "T").split())


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _clean_optional(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.lower() in {"nan", "none", "undefined", "null"} else value


def _number(value: str | None) -> float | None:
    value = _clean_optional(value)
    return float(value) if value else None


def _read_dataset(source_dir: Path, dataset: str) -> list[dict[str, str | float | None]]:
    path = source_dir / f"{dataset}.csv.gz"
    if not path.exists():
        raise FileNotFoundError(path)

    rows: list[dict[str, str | float | None]] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for source_row, raw in enumerate(reader, start=1):
            pair = [_dna(part) for part in str(raw.get("fasta", "")).split(".")]
            if len(pair) != 2 or not pair[0] or not pair[1]:
                raise ValueError(f"{dataset} row {source_row}: expected two FASTA strands")
            if any(set(strand) - set("ACGT") for strand in pair):
                raise ValueError(f"{dataset} row {source_row}: non-canonical base in FASTA field")

            # OligoGym's HELM marks RNA1 as ss and RNA2 as as. For siRNAforge's input convention,
            # use the antisense/guide component and construct its target site by reverse complement.
            passenger, guide = pair
            processed = _number(raw.get("y"))
            raw_label = _number(raw.get("y_raw"))
            if processed is None:
                raise ValueError(f"{dataset} row {source_row}: missing numeric y label")
            higher_is_better = processed if DATASETS[dataset]["direction"] == "higher_is_better" else 100.0 - processed
            source_file = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()
            rows.append(
                {
                    "source_row": source_row,
                    "target": _clean_optional(raw.get("targets")),
                    "passenger_sequence": passenger,
                    "guide_sequence": guide,
                    "guide_length": len(guide),
                    "label_processed": processed,
                    "label_raw": raw_label,
                    "efficacy_higher_is_better": higher_is_better,
                    "source_file": source_file,
                }
            )
    return rows


def _relative_or_absolute(path: Path) -> str:
    """Use repository-relative paths when possible, otherwise retain the absolute path."""
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


def prepare(source_dir: Path, output_dir: Path, datasets: tuple[str, ...]) -> dict[str, Any]:
    """Normalize source observations and write benchmark inputs with a manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = output_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    counts: Counter[str] = Counter()

    for dataset in datasets:
        config = DATASETS[dataset]
        for row in _read_dataset(source_dir, dataset):
            source_row = int(str(row["source_row"]))
            length = int(str(row["guide_length"]))
            transcript_id = f"bench_{dataset}_row{source_row:06d}"
            benchmark_id = f"{dataset}__row_{source_row:06d}"
            guide = str(row["guide_sequence"])
            target_site = _reverse_complement(guide)
            transcript = "A" * FLANK_NT + target_site + "A" * FLANK_NT
            record = {
                "benchmark_id": benchmark_id,
                "dataset": dataset,
                "modality": config["modality"],
                "source_row": source_row,
                "target": row["target"],
                "transcript_id": transcript_id,
                "input_fasta": f"inputs/{dataset}_len{length}.fasta",
                "site_start": SITE_START_1BASED,
                "site_end": SITE_START_1BASED + length - 1,
                "guide_sequence": guide,
                "passenger_sequence": row["passenger_sequence"],
                "guide_length": length,
                "label_processed": row["label_processed"],
                "label_raw": row["label_raw"],
                "efficacy_higher_is_better": row["efficacy_higher_is_better"],
                "label_direction": config["direction"],
                "label_semantics": config["label_semantics"],
                "context_type": "synthetic_neutral_flanks",
                "guide_source": "OligoGym fasta second component / HELM antisense strand",
                "source_file": row["source_file"],
                "source_url": config["source_url"],
            }
            records.append(record)
            grouped[(dataset, length)].append({"transcript_id": transcript_id, "sequence": transcript})
            counts[dataset] += 1

    records.sort(key=lambda row: (str(row["dataset"]), int(str(row["source_row"]))))
    fieldnames = list(records[0]) if records else []
    with (output_dir / "records.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    input_counts: dict[str, int] = {}
    commands: dict[str, str] = {}
    for (dataset, length), entries in sorted(grouped.items()):
        filename = f"{dataset}_len{length}.fasta"
        path = input_dir / filename
        with path.open("w", encoding="utf-8") as handle:
            for entry in entries:
                handle.write(f">{entry['transcript_id']}\n{entry['sequence']}\n")
        key = f"{dataset}_len{length}"
        input_counts[key] = len(entries)
        commands[key] = (
            f"uv run sirnaforge design {_relative_or_absolute(path)} --length {length} --output benchmark_results/{key}.csv"
        )

    manifest = {
        "description": "Runnable synthetic-context design inputs derived from OligoGym siRNA rows",
        "source_directory": _relative_or_absolute(source_dir),
        "output_directory": _relative_or_absolute(output_dir),
        "flank_nt": FLANK_NT,
        "site_start_1based": SITE_START_1BASED,
        "context_warning": (
            "Synthetic neutral flanks make design enumeration and parameter benchmarking runnable; "
            "they do not validate native transcript accessibility."
        ),
        "guide_convention": "second dot-separated fasta strand / HELM RNA2 antisense strand",
        "datasets": dict(counts),
        "total_records": len(records),
        "input_files": input_counts,
        "commands": commands,
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return manifest


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--dataset", action="append", choices=sorted(DATASETS), dest="datasets")
    args = parser.parse_args()
    selected = tuple(args.datasets or DATASETS)
    manifest = prepare(args.source_dir, args.output_dir, selected)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
