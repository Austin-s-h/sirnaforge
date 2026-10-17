#!/usr/bin/env python3
"""Link siRNAforge design candidates back to native benchmark observations.

The native design workflow intentionally operates on a panel-level FASTA containing every
transcript with at least one measured-site mapping.  This utility adds the missing review join:
candidate transcript/coordinate/passenger sequence -> native mapping -> benchmark row.

It reads the standard ``candidates_all.csv`` files and writes compact trace tables next to the
native design manifest.  No candidate is changed and no off-target result is inferred.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NATIVE_DIR = ROOT / "tests" / "data" / "benchmarks" / "oligogym_native_design"
DEFAULT_WORKFLOW_ROOT = ROOT / "workflow_results" / "oligogym_native_design"

RUNS = {
    "ichihara_2007_1": "ichihara_2007_1_len19",
    "ichihara_2007_2": "ichihara_2007_2_len19",
    "martinelli_2023_1": "martinelli_2023_1_len21",
}


MappingKey = tuple[str, str, str, str]
MappingIndex = dict[MappingKey, list[dict[str, str]]]
TranscriptIndex = dict[tuple[str, str], dict[str, str]]


def _read_mapping_indexes(mapping_path: Path) -> tuple[MappingIndex, MappingIndex, TranscriptIndex]:
    by_one_based: MappingIndex = defaultdict(list)
    by_zero_based: MappingIndex = defaultdict(list)
    transcripts: TranscriptIndex = {}
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["dataset"], row["mapping_transcript_id"], row["site_start_1based"], row["site_sequence"])
            by_one_based[key].append(row)
            zero_key = (
                row["dataset"],
                row["mapping_transcript_id"],
                str(int(row["site_start_1based"]) - 1),
                row["site_sequence"],
            )
            by_zero_based[zero_key].append(row)
            transcript_key = (row["dataset"], row["mapping_transcript_id"])
            transcripts.setdefault(
                transcript_key,
                {
                    "dataset": row["dataset"],
                    "mapping_transcript_id": row["mapping_transcript_id"],
                    "mapping_transcript_header": row["mapping_transcript_header"],
                    "mapping_gene_id": row["mapping_gene_id"],
                    "mapping_gene_symbol": row["mapping_gene_symbol"],
                    "mapping_assembly": row["mapping_assembly"],
                },
            )
    return by_one_based, by_zero_based, transcripts


def _write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
            count += 1
    return count


def link(native_dir: Path, workflow_root: Path) -> dict[str, object]:
    """Link workflow candidates to native benchmark mappings and write join tables."""
    mapping_path = native_dir / "mappings.csv"
    by_one_based, by_zero_based, transcript_metadata = _read_mapping_indexes(mapping_path)
    links_path = native_dir / "candidate_benchmark_links.csv"
    transcript_path = native_dir / "candidate_transcripts.csv"
    link_fields = [
        "dataset",
        "workflow_run",
        "candidate_source",
        "candidate_id",
        "transcript_id",
        "position",
        "coordinate_convention",
        "guide_sequence",
        "passenger_sequence",
        "gc_content",
        "design_score",
        "composite_score",
        "passes_filters",
        "selection_state",
        "quality_issues",
        "max_poly_run_length",
        "gc_content_min_verdict",
        "gc_content_max_verdict",
        "max_poly_runs_verdict",
        "benchmark_id",
        "source_row",
        "target",
        "label_processed",
        "duplex_pairing_status",
        "design_target_basis",
        "mapping_transcript_id",
        "mapping_gene_id",
        "mapping_gene_symbol",
        "mapping_assembly",
        "site_start_1based",
        "site_end_1based",
        "site_sequence",
        "mapping_strand",
        "mapping_rank",
    ]
    transcript_fields = [
        "dataset",
        "workflow_run",
        "transcript_id",
        "mapping_gene_id",
        "mapping_gene_symbol",
        "mapping_assembly",
        "mapping_transcript_header",
        "candidate_count",
        "candidate_min_position",
        "candidate_max_position",
        "measured_site_candidate_count",
    ]
    observation_by_id = {}
    with (native_dir / "observations.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            observation_by_id[row["benchmark_id"]] = row

    link_rows = []
    transcript_stats = {}
    for dataset, run_name in RUNS.items():
        candidate_path = workflow_root / run_name / "sirnaforge" / "candidates_all.csv"
        candidate_count_by_transcript: dict[str, int] = defaultdict(int)
        positions_by_transcript = defaultdict(list)
        measured_count_by_transcript: dict[str, int] = defaultdict(int)
        with candidate_path.open(encoding="utf-8", newline="") as handle:
            for candidate in csv.DictReader(handle):
                transcript_id = candidate["transcript_id"]
                position = candidate["position"]
                sequence = candidate["passenger_sequence"]
                candidate_count_by_transcript[transcript_id] += 1
                positions_by_transcript[transcript_id].append(int(position))
                key = (dataset, transcript_id, position, sequence)
                matches = by_one_based.get(key)
                convention = "1based" if matches else ""
                if not matches:
                    matches = by_zero_based.get(key, [])
                    convention = "0based" if matches else ""
                for mapping in matches:
                    observation = observation_by_id[mapping["benchmark_id"]]
                    measured_count_by_transcript[transcript_id] += 1
                    link_rows.append(
                        {
                            "dataset": dataset,
                            "workflow_run": run_name,
                            "candidate_source": "candidates_all.csv",
                            "candidate_id": candidate["id"],
                            "transcript_id": transcript_id,
                            "position": position,
                            "coordinate_convention": convention,
                            "guide_sequence": candidate["guide_sequence"],
                            "passenger_sequence": sequence,
                            "gc_content": candidate["gc_content"],
                            "design_score": candidate["design_score"],
                            "composite_score": candidate["composite_score"],
                            "passes_filters": candidate["passes_filters"],
                            "selection_state": candidate["selection_state"],
                            "quality_issues": candidate["quality_issues"],
                            "max_poly_run_length": candidate["max_poly_run_length"],
                            "gc_content_min_verdict": candidate["gc_content_min_verdict"],
                            "gc_content_max_verdict": candidate["gc_content_max_verdict"],
                            "max_poly_runs_verdict": candidate["max_poly_runs_verdict"],
                            "benchmark_id": mapping["benchmark_id"],
                            "source_row": observation["source_row"],
                            "target": observation["target"],
                            "label_processed": observation["label_processed"],
                            "duplex_pairing_status": observation["duplex_pairing_status"],
                            "design_target_basis": observation["design_target_basis"],
                            "mapping_transcript_id": mapping["mapping_transcript_id"],
                            "mapping_gene_id": mapping["mapping_gene_id"],
                            "mapping_gene_symbol": mapping["mapping_gene_symbol"],
                            "mapping_assembly": mapping["mapping_assembly"],
                            "site_start_1based": mapping["site_start_1based"],
                            "site_end_1based": mapping["site_end_1based"],
                            "site_sequence": mapping["site_sequence"],
                            "mapping_strand": mapping["mapping_strand"],
                            "mapping_rank": mapping["mapping_rank"],
                        }
                    )
        for transcript_id, count in candidate_count_by_transcript.items():
            meta = transcript_metadata[(dataset, transcript_id)]
            positions = positions_by_transcript[transcript_id]
            transcript_stats[(dataset, transcript_id)] = {
                **meta,
                "workflow_run": run_name,
                "transcript_id": transcript_id,
                "candidate_count": count,
                "candidate_min_position": min(positions),
                "candidate_max_position": max(positions),
                "measured_site_candidate_count": measured_count_by_transcript[transcript_id],
            }

    link_rows.sort(key=lambda row: (row["dataset"], int(row["source_row"]), row["candidate_id"], row["mapping_rank"]))
    transcript_rows = sorted(transcript_stats.values(), key=lambda row: (row["dataset"], row["transcript_id"]))
    link_count = _write_csv(links_path, link_fields, link_rows)
    transcript_count = _write_csv(transcript_path, transcript_fields, transcript_rows)
    summary = {
        "description": "Trace links from siRNAforge design candidates to native OligoGym observations",
        "candidate_source": "candidates_all.csv",
        "workflow_root": str(workflow_root),
        "mapping_source": str(mapping_path),
        "coordinate_join": "candidate position is matched as 1-based first, with a 0-based fallback",
        "offtarget_note": "No off-target result is inferred; all workflow runs were design_only.",
        "candidate_benchmark_link_rows": link_count,
        "candidate_native_transcript_rows": transcript_count,
        "datasets": {
            dataset: {
                "workflow_run": run_name,
                "candidate_benchmark_link_rows": sum(row["dataset"] == dataset for row in link_rows),
                "candidate_native_transcript_rows": sum(row["dataset"] == dataset for row in transcript_rows),
            }
            for dataset, run_name in RUNS.items()
        },
    }
    with (native_dir / "design_trace_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    return summary


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-dir", type=Path, default=DEFAULT_NATIVE_DIR)
    parser.add_argument("--workflow-root", type=Path, default=DEFAULT_WORKFLOW_ROOT)
    args = parser.parse_args()
    print(json.dumps(link(args.native_dir, args.workflow_root), indent=2))


if __name__ == "__main__":
    main()
