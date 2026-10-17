#!/usr/bin/env python3
"""Prepare target-aware, version-compatible native cDNA design inputs.

This is deliberately narrower than an exhaustive off-target catalog.  It creates the native
transcript FASTA and join tables needed to run the current siRNAforge design/workflow commands:

* Ichihara/Huesken: map the 19-nt paired target within the recorded source gene.
* Martinelli: use the nominal 21-nt passenger/target strand for the design artifact, while
  preserving the inferred pairing status (the source has mixed 19-nt cores and invalid pairs).
* Shmushkovich: preserve source metadata as review-only because its 15/20-nt architecture is not
  supported by the current designer.

Off-target screening against the complete four-species transcriptome is intentionally not done by
this preparer; it is a separate siRNAforge operation.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, TextIO, TypedDict

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import map_oligo_benchmarks_to_ensembl as source  # noqa: E402

ROOT = source.ROOT
DEFAULT_SOURCE_DIR = source.DEFAULT_SOURCE_DIR
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "data" / "benchmarks" / "oligogym_native_design"


class DesignConfig(TypedDict):
    """Describe a panel's native design settings."""

    design_basis: str
    design_length: int
    target_filter: bool
    workflow_status: str


DESIGN_CONFIG: dict[str, DesignConfig] = {
    "ichihara_2007_1": {
        "design_basis": "paired_target_sequence",
        "design_length": 19,
        "target_filter": True,
        "workflow_status": "runnable",
    },
    "ichihara_2007_2": {
        "design_basis": "paired_target_sequence",
        "design_length": 19,
        "target_filter": True,
        "workflow_status": "runnable",
    },
    "martinelli_2023_1": {
        "design_basis": "full_passenger_sequence_nominal_21nt",
        "design_length": 21,
        "target_filter": False,
        "workflow_status": "runnable_unqualified_target_identity",
    },
    "shmushkovich_2018_1": {
        "design_basis": "paired_target_sequence",
        "design_length": 15,
        "target_filter": False,
        "workflow_status": "review_only_unsupported_length",
    },
}


def _header_metadata(description: str) -> dict[str, str]:
    assembly_match = re.search(r"(?:chromosome|scaffold):([^: ]+):", description)
    gene_match = re.search(r"\bgene:([^ ]+)", description)
    symbol_match = re.search(r"\bgene_symbol:([^ ]+)", description)
    return {
        "gene_id": gene_match.group(1) if gene_match else "",
        "gene_symbol": symbol_match.group(1) if symbol_match else "",
        "assembly": assembly_match.group(1) if assembly_match else "",
    }


def _query_sequence(row: dict[str, Any], dataset: str) -> str:
    if dataset == "martinelli_2023_1":
        return str(row["passenger_sequence"])
    return str(row["paired_target_sequence"])


def _write_fasta_record(handle: TextIO, description: str, sequence: str) -> None:
    handle.write(f">{description}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


# Keep panel preparation and manifest assembly together for this standalone data conversion.
def prepare(source_dir: Path, output_dir: Path, cdna_fasta: Path, datasets: tuple[str, ...]) -> dict[str, Any]:  # noqa: PLR0912
    """Map panel targets to native transcripts and write design inputs with provenance."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = output_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    rows_by_dataset = {dataset: source._read_observations(source_dir, dataset) for dataset in datasets}

    by_sequence: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for dataset, rows in rows_by_dataset.items():
        config = DESIGN_CONFIG[dataset]
        if config["workflow_status"].startswith("review_only"):
            continue
        for row in rows:
            query = _query_sequence(row, dataset)
            if query:
                row["design_target_sequence"] = query
                row["design_target_length"] = config["design_length"]
                row["design_target_basis"] = config["design_basis"]
                by_sequence[query].append(row)
            else:
                row["design_target_sequence"] = ""
                row["design_target_length"] = config["design_length"]
                row["design_target_basis"] = config["design_basis"]

    mapping_fields = [
        "mapping_id",
        "benchmark_id",
        "dataset",
        "query_sequence",
        "design_target_length",
        "target_symbol",
        "mapping_status",
        "mapping_transcript_id",
        "mapping_transcript_header",
        "mapping_gene_id",
        "mapping_gene_symbol",
        "mapping_assembly",
        "site_start_1based",
        "site_end_1based",
        "site_sequence",
        "mapping_strand",
        "mapping_rank",
    ]
    mapping_counts: dict[str, int] = defaultdict(int)
    mapping_ranks: dict[str, int] = defaultdict(int)
    transcript_counts: dict[tuple[str, int], int] = defaultdict(int)
    fasta_handles: dict[tuple[str, int], TextIO] = {}
    state_path = output_dir / ".native_design_seen.sqlite3"
    state_path.unlink(missing_ok=True)
    state = sqlite3.connect(state_path)
    state.execute(
        "CREATE TABLE emitted (dataset TEXT, design_length INTEGER, transcript_id TEXT, "
        "PRIMARY KEY(dataset, design_length, transcript_id))"
    )
    mapping_id = 0
    with (output_dir / "mappings.csv").open("w", encoding="utf-8", newline="") as mapping_handle:
        mapping_writer = csv.DictWriter(mapping_handle, fieldnames=mapping_fields)
        mapping_writer.writeheader()
        automaton = source._build_automaton(list(by_sequence))
        for record_number, (record_id, description, raw_sequence) in enumerate(source._iter_fasta(cdna_fasta), start=1):
            sequence = raw_sequence.upper().replace("U", "T")
            metadata = _header_metadata(description)
            if record_number % 100_000 == 0:
                print(f"scanned {record_number:,} native transcripts", flush=True)
            for end_index, pattern in automaton.iter(sequence):
                start_zero = end_index - len(pattern) + 1
                for row in by_sequence[pattern]:
                    dataset = str(row["dataset"])
                    config = DESIGN_CONFIG[dataset]
                    target_symbol = str(row["target"])
                    if config["target_filter"] and target_symbol and metadata["gene_symbol"] != target_symbol:
                        continue
                    benchmark_id = str(row["benchmark_id"])
                    mapping_id += 1
                    mapping_counts[benchmark_id] += 1
                    mapping_ranks[benchmark_id] += 1
                    mapping_writer.writerow(
                        {
                            "mapping_id": f"mapping_{mapping_id:012d}",
                            "benchmark_id": benchmark_id,
                            "dataset": dataset,
                            "query_sequence": pattern,
                            "design_target_length": config["design_length"],
                            "target_symbol": target_symbol,
                            "mapping_status": "mapped",
                            "mapping_transcript_id": record_id,
                            "mapping_transcript_header": description,
                            "mapping_gene_id": metadata["gene_id"],
                            "mapping_gene_symbol": metadata["gene_symbol"],
                            "mapping_assembly": metadata["assembly"],
                            "site_start_1based": start_zero + 1,
                            "site_end_1based": start_zero + len(pattern),
                            "site_sequence": pattern,
                            "mapping_strand": "+",
                            "mapping_rank": mapping_ranks[benchmark_id],
                        }
                    )
                    group = (dataset, config["design_length"])
                    if config["workflow_status"].startswith("runnable"):
                        inserted = state.execute(
                            "INSERT OR IGNORE INTO emitted VALUES (?, ?, ?)",
                            (dataset, config["design_length"], record_id),
                        ).rowcount
                        if inserted:
                            fasta_path = input_dir / f"{dataset}_native_len{config['design_length']}.fasta"
                            handle = fasta_handles.get(group)
                            if handle is None:
                                handle = fasta_path.open("w", encoding="utf-8")
                                fasta_handles[group] = handle
                            _write_fasta_record(handle, description, sequence)
                            transcript_counts[group] += 1
                    if mapping_id % 10_000 == 0:
                        state.commit()
    state.commit()
    state.close()
    for handle in fasta_handles.values():
        handle.close()
    state_path.unlink(missing_ok=True)

    all_observations: list[dict[str, Any]] = []
    dataset_summary: dict[str, Any] = {}
    for dataset, rows in rows_by_dataset.items():
        config = DESIGN_CONFIG[dataset]
        for row in rows:
            benchmark_id = str(row["benchmark_id"])
            count = mapping_counts.get(benchmark_id, 0)
            row.update(
                {
                    "native_reference_species": "human",
                    "native_reference_assembly": "GRCh38",
                    "native_reference_fasta": source._relative_or_absolute(cdna_fasta),
                    "native_mapping_count": count,
                    "native_mapping_status": (
                        "review_only_not_mapped"
                        if config["workflow_status"].startswith("review_only")
                        else "mapped"
                        if count
                        else "unmapped"
                    ),
                    "target_identity_status": ("confirmed" if row["target"] else "unavailable"),
                    "workflow_input_fasta": (
                        f"inputs/{dataset}_native_len{config['design_length']}.fasta"
                        if count and config["workflow_status"].startswith("runnable")
                        else ""
                    ),
                    "workflow_status": config["workflow_status"],
                }
            )
            all_observations.append(row)
        length = int(config["design_length"])
        group = (dataset, length)
        dataset_summary[dataset] = {
            "observations": len(rows),
            "design_basis": config["design_basis"],
            "design_length": length,
            "paired_lengths": sorted({int(row["paired_length"]) for row in rows}),
            "mapped_observations": sum(mapping_counts.get(str(row["benchmark_id"]), 0) > 0 for row in rows),
            "unmapped_observations": sum(mapping_counts.get(str(row["benchmark_id"]), 0) == 0 for row in rows),
            "mapping_rows": sum(mapping_counts.get(str(row["benchmark_id"]), 0) for row in rows),
            "native_transcripts": transcript_counts.get(group, 0),
            "native_input_fasta": (
                source._relative_or_absolute(input_dir / f"{dataset}_native_len{length}.fasta")
                if transcript_counts.get(group, 0)
                else ""
            ),
            "workflow_status": config["workflow_status"],
        }

    all_observations.sort(key=lambda row: (str(row["dataset"]), int(row["source_row"])))
    observation_fields = list(all_observations[0]) if all_observations else []
    with (output_dir / "observations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=observation_fields)
        writer.writeheader()
        writer.writerows(all_observations)

    manifest = {
        "description": "Target-aware native Ensembl design artifacts for OligoGym siRNA panels",
        "source_directory": source._relative_or_absolute(source_dir),
        "native_reference": {
            "species": "human",
            "assembly": "GRCh38",
            "fasta": source._relative_or_absolute(cdna_fasta),
        },
        "species_for_separate_offtarget_screen": ["human", "mouse", "rat", "macaque"],
        "dataset_summary": dataset_summary,
        "total_observations": len(all_observations),
        "total_mapping_rows": mapping_id,
        "mapping_join": "mappings.csv is keyed directly by benchmark_id; native FASTA headers retain Ensembl transcript/gene metadata",
        "offtarget_note": "This preparer creates design inputs only; run the separate siRNAforge off-target screen against all four species.",
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
    parser.add_argument("--cdna-fasta", type=Path, required=True)
    parser.add_argument("--dataset", action="append", choices=sorted(DESIGN_CONFIG), dest="datasets")
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(args.source_dir, args.output_dir, args.cdna_fasta, tuple(args.datasets or DESIGN_CONFIG)), indent=2
        )
    )


if __name__ == "__main__":
    main()
