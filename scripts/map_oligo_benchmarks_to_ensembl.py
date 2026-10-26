#!/usr/bin/env python3
"""Map OligoGym duplex observations to native Ensembl cDNA sequences.

OligoGym's ``fasta`` column contains the two oligo strands, not an mRNA.  The first strand is
the passenger/sense strand and the second is the antisense guide.  Only the complementary region
belongs in the mRNA target site; unequal strand lengths represent asymmetric architecture or an
overhang.  This script therefore searches the passenger strand in native cDNA and preserves the
full guide and its unpaired tail separately.

The output has two levels:

* ``observations.csv``: exactly one row per source observation, including mapping status/count.
* ``mappings.csv``: one normalized row per unique query-sequence/native transcript hit with
  coordinates and Ensembl header metadata; join observations on ``paired_target_sequence``.

Compatible panels also receive a native target FASTA containing each matching transcript once.  The
current siRNAforge designer supports 19-23 nt target sites; the 15 nt/20 nt Shmushkovich hsiRNA
panel is mapped and documented but is not emitted as a runnable design input.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sqlite3
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any, TextIO

import ahocorasick

from scripts.native_mapping_generator import (
    NativeMappingInputError,
    file_sha256,
    load_reference_manifests,
    write_species_sidecars,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DIR = ROOT / "tests" / "data" / "external" / "oligogym"
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "data" / "benchmarks" / "oligogym_native"

DATASETS = {
    "ichihara_2007_1": {
        "source_url": "https://doi.org/10.1093/nar/gkm699",
        "label_semantics": "percentage inhibition relative to control; higher is better",
        "label_direction": "higher_is_better",
    },
    "ichihara_2007_2": {
        "source_url": "https://doi.org/10.1093/nar/gkm699",
        "label_semantics": "percentage inhibition relative to control; higher is better",
        "label_direction": "higher_is_better",
    },
    "martinelli_2023_1": {
        "source_url": "https://doi.org/10.1016/j.ygeno.2024.110815",
        "label_semantics": "percentage knockdown of target mRNA; higher is better",
        "label_direction": "higher_is_better",
    },
    "shmushkovich_2018_1": {
        "source_url": "https://doi.org/10.1093/nar/gky745",
        "label_semantics": "percentage remaining relative to control; lower is better",
        "label_direction": "lower_is_better",
    },
}


def _dna(sequence: str) -> str:
    return "".join(sequence.upper().replace("U", "T").split())


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _clean(value: str | None) -> str:
    value = (value or "").strip()
    return "" if value.lower() in {"nan", "none", "undefined", "null"} else value


def _number(value: str | None) -> float | None:
    value = _clean(value)
    return float(value) if value else None


def _paired_prefix(passenger: str, full_guide: str) -> int:
    """Return the longest prefix that is exactly complementary in the two reported strands."""
    paired_length = 0
    for length in range(1, min(len(passenger), len(full_guide)) + 1):
        if _reverse_complement(full_guide[:length]) == passenger[:length]:
            paired_length = length
    return paired_length


def _read_observations(source_dir: Path, dataset: str) -> list[dict[str, Any]]:
    path = source_dir / f"{dataset}.csv.gz"
    if not path.exists():
        raise FileNotFoundError(path)

    rows: list[dict[str, Any]] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for source_row, raw in enumerate(csv.DictReader(handle), start=1):
            strands = [_dna(part) for part in str(raw.get("fasta", "")).split(".")]
            if len(strands) != 2 or not all(strands):
                raise ValueError(f"{dataset} row {source_row}: expected two canonical FASTA strands")
            passenger, full_guide = strands
            paired_length = _paired_prefix(passenger, full_guide)
            paired_guide = full_guide[:paired_length]
            passenger_overhang = passenger[paired_length:]
            guide_overhang = full_guide[paired_length:]
            if paired_length == 0:
                pairing_status = "non_complementary"
            elif paired_length == min(len(passenger), len(full_guide)):
                pairing_status = "fully_complementary"
            else:
                pairing_status = "paired_core_with_overhang"
            processed = _number(raw.get("y"))
            raw_label = _number(raw.get("y_raw"))
            if processed is None:
                raise ValueError(f"{dataset} row {source_row}: missing numeric y label")
            direction = DATASETS[dataset]["label_direction"]
            higher_is_better = processed if direction == "higher_is_better" else 100.0 - processed
            source_file = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()
            rows.append(
                {
                    "benchmark_id": f"{dataset}__row_{source_row:06d}",
                    "dataset": dataset,
                    "source_row": source_row,
                    "target": _clean(raw.get("targets")),
                    "full_guide_sequence": full_guide,
                    "paired_guide_sequence": paired_guide,
                    "passenger_sequence": passenger,
                    "paired_target_sequence": passenger[:paired_length],
                    "guide_length": len(full_guide),
                    "paired_length": paired_length,
                    "guide_overhang_sequence": guide_overhang,
                    "passenger_overhang_sequence": passenger_overhang,
                    "duplex_pairing_status": pairing_status,
                    "label_processed": processed,
                    "label_raw": raw_label,
                    "efficacy_higher_is_better": higher_is_better,
                    "label_direction": direction,
                    "label_semantics": DATASETS[dataset]["label_semantics"],
                    "source_file": source_file,
                    "source_url": DATASETS[dataset]["source_url"],
                }
            )
    return rows


def _header_metadata(description: str) -> dict[str, str]:
    assembly_match = re.search(r"(?:chromosome|scaffold):([^: ]+):", description)
    gene_match = re.search(r"\bgene:([^ ]+)", description)
    return {
        "gene_id": gene_match.group(1) if gene_match else "",
        "assembly": assembly_match.group(1) if assembly_match else "",
    }


def _build_automaton(patterns: list[str]) -> ahocorasick.Automaton:
    automaton = ahocorasick.Automaton()
    for pattern in sorted(set(patterns)):
        automaton.add_word(pattern, pattern)
    automaton.make_automaton()
    return automaton


def _iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Stream a plain FASTA without materializing the native transcriptome."""
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header.split()[0], header, "".join(chunks)


def _relative_or_absolute(path: Path) -> str:
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


# Keep mapping and deduplicated FASTA emission in a single streaming pass.
def map_datasets(  # noqa: PLR0912
    dataset_rows: dict[str, list[dict[str, Any]]],
    cdna_fasta: Path,
    mapping_writer: csv.DictWriter,
    native_input_dir: Path,
) -> tuple[dict[str, dict[str, int]], dict[tuple[str, int], int], int]:
    """Map all selected datasets in one pass, streaming normalized results to disk."""
    by_sequence: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for rows in dataset_rows.values():
        for row in rows:
            paired_target = str(row["paired_target_sequence"])
            if paired_target:
                by_sequence[paired_target].append(row)
    automaton = _build_automaton(list(by_sequence))
    pattern_counts: dict[str, int] = defaultdict(int)
    mapping_ranks: dict[str, int] = defaultdict(int)
    mapping_serial = 0
    transcript_counts: dict[tuple[str, int], int] = defaultdict(int)
    fasta_handles: dict[tuple[str, int], TextIO] = {}
    state_path = native_input_dir.parent / ".native_fasta_seen.sqlite3"
    state = sqlite3.connect(state_path)
    state.execute(
        "CREATE TABLE IF NOT EXISTS emitted_transcripts "
        "(dataset TEXT NOT NULL, paired_length INTEGER NOT NULL, transcript_id TEXT NOT NULL, "
        "PRIMARY KEY (dataset, paired_length, transcript_id))"
    )

    for record_number, (record_id, description, raw_sequence) in enumerate(_iter_fasta(cdna_fasta), start=1):
        sequence = raw_sequence.upper().replace("U", "T")
        metadata = _header_metadata(description)
        if record_number % 100_000 == 0:
            print(f"scanned {record_number:,} native transcripts", flush=True)
        for end_index, pattern in automaton.iter(sequence):
            start_zero = end_index - len(pattern) + 1
            pattern_counts[pattern] += 1
            mapping_ranks[pattern] += 1
            mapping_serial += 1
            mapping_writer.writerow(
                {
                    "mapping_id": f"mapping_{mapping_serial:012d}",
                    "query_sequence": pattern,
                    "mapping_status": "mapped",
                    "mapping_transcript_id": record_id,
                    "mapping_transcript_header": description,
                    "mapping_gene_id": metadata["gene_id"],
                    "mapping_assembly": metadata["assembly"],
                    "site_start_1based": start_zero + 1,
                    "site_end_1based": start_zero + len(pattern),
                    "site_sequence": pattern,
                    "mapping_strand": "+",
                    "mapping_rank": mapping_ranks[pattern],
                }
            )
            if mapping_serial % 10_000 == 0:
                state.commit()
            for row in by_sequence[pattern]:
                dataset = str(row["dataset"])
                paired_length = int(row["paired_length"])
                group = (dataset, paired_length)
                if 19 <= paired_length <= 23:
                    inserted = state.execute(
                        "INSERT OR IGNORE INTO emitted_transcripts "
                        "(dataset, paired_length, transcript_id) VALUES (?, ?, ?)",
                        (dataset, paired_length, record_id),
                    ).rowcount
                else:
                    inserted = 0
                if inserted:
                    fasta_path = native_input_dir / f"{dataset}_native_len{paired_length}.fasta"
                    handle = fasta_handles.get(group)
                    if handle is None:
                        handle = fasta_path.open("w", encoding="utf-8")
                        fasta_handles[group] = handle
                    handle.write(f">{description}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")
                    transcript_counts[group] += 1

    state.commit()
    state.close()
    for handle in fasta_handles.values():
        handle.close()
    mapping_counts: dict[str, dict[str, int]] = {
        dataset: {str(row["benchmark_id"]): pattern_counts[str(row["paired_target_sequence"])] for row in rows}
        for dataset, rows in dataset_rows.items()
    }
    return mapping_counts, transcript_counts, mapping_serial


def prepare(source_dir: Path, output_dir: Path, cdna_fasta: Path, datasets: tuple[str, ...]) -> dict[str, Any]:
    """Write native transcript mappings, observation tables, and workflow commands."""
    output_dir.mkdir(parents=True, exist_ok=True)
    native_input_dir = output_dir / "inputs"
    native_input_dir.mkdir(parents=True, exist_ok=True)
    all_observations: list[dict[str, Any]] = []
    dataset_summary: dict[str, Any] = {}

    rows_by_dataset = {dataset: _read_observations(source_dir, dataset) for dataset in datasets}
    mapping_fields = [
        "mapping_id",
        "query_sequence",
        "mapping_status",
        "mapping_transcript_id",
        "mapping_transcript_header",
        "mapping_gene_id",
        "mapping_assembly",
        "site_start_1based",
        "site_end_1based",
        "site_sequence",
        "mapping_strand",
        "mapping_rank",
    ]
    with (output_dir / "mappings.csv").open("w", encoding="utf-8", newline="") as mapping_handle:
        mapping_writer = csv.DictWriter(mapping_handle, fieldnames=mapping_fields)
        mapping_writer.writeheader()
        mapping_counts_by_dataset, transcript_counts, mapping_rows_written = map_datasets(
            rows_by_dataset, cdna_fasta, mapping_writer, native_input_dir
        )

    for dataset in datasets:
        rows = rows_by_dataset[dataset]
        mapping_count_by_id = mapping_counts_by_dataset[dataset]
        for row in rows:
            benchmark_id = str(row["benchmark_id"])
            count = mapping_count_by_id.get(benchmark_id, 0)
            compatible = 19 <= int(row["paired_length"]) <= 23
            row.update(
                {
                    "native_reference_species": "human",
                    "native_reference_assembly": "GRCh38",
                    "native_reference_fasta": _relative_or_absolute(cdna_fasta),
                    "native_mapping_count": count,
                    "native_mapping_status": "mapped" if count else "unmapped",
                    "sirnaforge_compatible": compatible,
                    "sirnaforge_length": int(row["paired_length"]) if compatible else "",
                    "workflow_input_fasta": "",
                }
            )
            all_observations.append(row)

        length_summary: dict[str, Any] = {}
        for paired_length in sorted({int(row["paired_length"]) for row in rows}):
            length_rows = [row for row in rows if int(row["paired_length"]) == paired_length]
            compatible_length = 19 <= paired_length <= 23
            fasta_path = native_input_dir / f"{dataset}_native_len{paired_length}.fasta"
            transcript_count = transcript_counts.get((dataset, paired_length), 0)
            if compatible_length and transcript_count and fasta_path.exists():
                workflow_input = _relative_or_absolute(fasta_path)
            else:
                workflow_input = ""
            for row in length_rows:
                if mapping_count_by_id.get(str(row["benchmark_id"]), 0) and workflow_input:
                    row["workflow_input_fasta"] = workflow_input
            mapped_length_observations = sum(
                1 for row in length_rows if mapping_count_by_id.get(str(row["benchmark_id"]), 0)
            )
            length_summary[str(paired_length)] = {
                "observations": len(length_rows),
                "mapped_observations": mapped_length_observations,
                "unmapped_observations": len(length_rows) - mapped_length_observations,
                "native_transcripts": transcript_count,
                "observation_mapping_rows": sum(
                    mapping_count_by_id.get(str(row["benchmark_id"]), 0) for row in length_rows
                ),
                "sirnaforge_compatible": compatible_length,
                "native_input_fasta": workflow_input,
            }

        mapped_observations = sum(1 for row in rows if mapping_count_by_id.get(str(row["benchmark_id"]), 0))
        dataset_summary[dataset] = {
            "observations": len(rows),
            "paired_lengths": sorted({int(row["paired_length"]) for row in rows}),
            "mapped_observations": mapped_observations,
            "unmapped_observations": len(rows) - mapped_observations,
            "native_transcripts": sum(
                transcript_counts.get((dataset, paired_length), 0)
                for paired_length in {int(row["paired_length"]) for row in rows}
            ),
            "observation_mapping_rows": sum(mapping_count_by_id.values()),
            "sirnaforge_compatible": any(item["sirnaforge_compatible"] for item in length_summary.values()),
            "native_input_fastas": [
                item["native_input_fasta"] for item in length_summary.values() if item["native_input_fasta"]
            ],
            "paired_length_summary": length_summary,
        }

    all_observations.sort(key=lambda row: (str(row["dataset"]), int(row["source_row"])))
    observation_fields = list(all_observations[0]) if all_observations else []
    with (output_dir / "observations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=observation_fields)
        writer.writeheader()
        writer.writerows(all_observations)

    commands: dict[str, str] = {}
    for dataset, summary in dataset_summary.items():
        for length_key, length_summary in summary.get("paired_length_summary", {}).items():
            if not length_summary["native_input_fasta"]:
                continue
            length = int(length_key)
            commands[f"{dataset}_len{length}"] = (
                f"uv run sirnaforge workflow benchmark_oligogym_{dataset} "
                f"--input-fasta {length_summary['native_input_fasta']} --length {length} "
                f"--run-mode exploratory --design-mode sirna --database ensembl "
                f"--species human,mouse,rat,macaque --query-species human "
                f"--transcriptome-fasta ensembl_human_cdna,ensembl_mouse_cdna,"
                f"ensembl_rat_cdna,ensembl_macaque_cdna --gc-min 0 --gc-max 100 "
                f"--modifications none --output-dir workflow_results/oligogym_native/"
                f"{dataset}_len{length}"
            )

    manifest = {
        "description": "Native Ensembl cDNA mapping of OligoGym duplex observations",
        "source_directory": _relative_or_absolute(source_dir),
        "native_reference": {
            "species": "human",
            "assembly": "GRCh38",
            "fasta": _relative_or_absolute(cdna_fasta),
        },
        "guide_convention": (
            "first fasta component is passenger/sense; second is full antisense guide; "
            "the paired region is the prefix of the guide complementary to the passenger"
        ),
        "dataset_summary": dataset_summary,
        "total_observations": len(all_observations),
        "total_mapping_rows": mapping_rows_written,
        "commands": commands,
        "shmushkovich_note": (
            "The 15 nt passenger / 20 nt guide hsiRNA panel is mapped for review but is not emitted "
            "as a runnable design input because siRNAforge requires 19-23 nt target sites."
        ),
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return manifest


def prepare_four_species(
    source_dir: Path,
    output_dir: Path,
    reference_manifests: tuple[Path, ...],
    datasets: tuple[str, ...],
) -> dict[str, Any]:
    """Generate exact native/projection mappings from four caller-pinned references."""
    references = load_reference_manifests(reference_manifests)
    rows_by_dataset = {dataset: _read_observations(source_dir, dataset) for dataset in datasets}
    source_rows = [row for dataset in datasets for row in rows_by_dataset[dataset]]
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "assay_observations.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(source_rows[0]) if source_rows else []
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(source_rows)
    species_summary = write_species_sidecars(source_rows, references, output_dir / "species")
    manifest = json.loads((output_dir / "species" / "manifest.json").read_text())
    if not isinstance(manifest, dict):
        raise ValueError("multi-species sidecar writer returned a malformed manifest")
    manifest.update(
        {
            "description": "Exact paired-target mapping of human OligoGym assays against four caller-pinned cDNA references",
            "source_directory": _relative_or_absolute(source_dir),
            "datasets": {
                dataset: {
                    "observation_count": len(rows),
                    "paired_lengths": sorted({int(row["paired_length"]) for row in rows}),
                    "source_file": _relative_or_absolute(source_dir / f"{dataset}.csv.gz"),
                    "source_sha256": file_sha256(source_dir / f"{dataset}.csv.gz"),
                }
                for dataset, rows in rows_by_dataset.items()
            },
            "assay_observations_csv": {
                "path": "assay_observations.csv",
                "sha256": file_sha256(output_dir / "assay_observations.csv"),
                "rows": len(source_rows),
            },
            "species": species_summary,
        }
    )
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cdna-fasta", type=Path, help="Legacy single human Ensembl cDNA FASTA to search")
    parser.add_argument(
        "--species-reference-manifest",
        type=Path,
        action="append",
        default=[],
        help="Pinned JSON reference declaration; pass once for each human, mouse, rat, and macaque FASTA",
    )
    parser.add_argument("--dataset", action="append", choices=sorted(DATASETS), dest="datasets")
    args = parser.parse_args()
    selected = tuple(args.datasets or DATASETS)
    if args.species_reference_manifest:
        if args.cdna_fasta is not None:
            parser.error("--species-reference-manifest cannot be combined with --cdna-fasta")
        try:
            result = prepare_four_species(
                args.source_dir, args.output_dir, tuple(args.species_reference_manifest), selected
            )
        except NativeMappingInputError as exc:
            parser.error(str(exc))
    else:
        if args.cdna_fasta is None:
            parser.error("provide --cdna-fasta or all four --species-reference-manifest declarations")
        result = prepare(args.source_dir, args.output_dir, args.cdna_fasta, selected)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
