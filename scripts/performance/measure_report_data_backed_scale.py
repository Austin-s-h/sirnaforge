#!/usr/bin/env python3
"""Measure data-backed report generation at the historical large-report row counts.

The guide sequences come from the committed native OligoGym observations. Candidate expansion,
off-target hits, and miRNA rows are synthetic resource stress data, not biological observations.
The run directory, report, Parquet tables, and detail shards live only in a temporary directory;
the optional JSON output records compact measurement metadata.

Run from the repository root with::

    uv run python scripts/performance/measure_report_data_backed_scale.py \
      --output scripts/performance/measurements/report-data-backed-large.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS
from sirnaforge.reporting import write_report_from_run

ROOT = Path(__file__).resolve().parents[2]
OBSERVATIONS = ROOT / "tests/data/benchmarks/oligogym_native_design/observations.csv"
N_CANDIDATES = 29_605
N_HITS = 254_487
N_MIRNA = 118_559
N_NATIVE_OBSERVATIONS = 4_113

CANDIDATE_COLUMNS = [
    "id",
    "guide_sequence",
    "passenger_sequence",
    "transcript_id",
    "position",
    "gc_content",
    "asymmetry_score",
    "empirical_score",
    "paired_fraction",
    "design_score",
    "composite_score",
    "weight_vector",
    "passes_filters",
    "off_target_screened",
]
HIT_COLUMNS = [
    "qseq",
    "species",
    "rname",
    "strand",
    "cigar",
    "nm",
    "seed_mismatches",
    *CLASSIFICATION_COLUMNS,
]
MIRNA_COLUMNS = [
    "qseq",
    "mirna_id",
    "rname",
    "mirna",
    "name",
    "nm",
    "seed_mismatches",
    "coord",
    "species",
    "source",
    "database",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _guides() -> tuple[list[str], int]:
    with OBSERVATIONS.open(encoding="utf-8", newline="") as handle:
        observed = [str(row["full_guide_sequence"]).replace("U", "T") for row in csv.DictReader(handle)]
    guides = list(dict.fromkeys(observed))
    if len(observed) != N_NATIVE_OBSERVATIONS or len(guides) != 3_395:
        raise RuntimeError("native OligoGym observation fixture changed; review the scale inputs")
    return guides, len(observed)


def _write_candidate_table(path: Path, guides: list[str]) -> None:
    columns = CANDIDATE_COLUMNS + [f"unused_candidate_{i:02d}" for i in range(52)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for index in range(N_CANDIDATES):
            guide = guides[index % len(guides)]
            row: dict[str, Any] = {
                "id": f"candidate-{index:06d}",
                "guide_sequence": guide,
                "passenger_sequence": "",
                "transcript_id": f"ENST_SYNTH_{index % 20:04d}",
                "position": index % 2_200 + 1,
                "gc_content": "45.0",
                "asymmetry_score": "0.5",
                "empirical_score": "0.5",
                "paired_fraction": "1.0",
                "design_score": "50.0",
                "composite_score": "50.0",
                "weight_vector": "synthetic_resource_fixture",
                "passes_filters": "PASS",
                "off_target_screened": "True",
            }
            writer.writerow([row.get(column, "x") for column in columns])


def _write_hit_table(path: Path, guides: list[str]) -> None:
    columns = HIT_COLUMNS + [f"unused_hit_{i:02d}" for i in range(37)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for index in range(N_HITS):
            guide = guides[index % len(guides)]
            row: dict[str, Any] = {
                "qseq": guide,
                "species": "human",
                "rname": f"ENST_SYNTH_HIT_{index:07d}",
                "strand": "+" if index % 2 == 0 else "-",
                "cigar": f"{len(guide)}M",
                "nm": index % 3,
                "seed_mismatches": index % 2,
                "hit_class": "off_target",
                "matched_symbol": "SYNTHETIC",
                "symbol_lookup_missing": "False",
                "hit_symbol": "SYNTHETIC",
                "hit_symbol_missing": "False",
                "species_index_missing": "False",
                "ortholog_evidence": "not_applicable",
            }
            writer.writerow([row.get(column, "x") for column in columns])


def _write_mirna_table(path: Path, guides: list[str]) -> None:
    columns = MIRNA_COLUMNS + [f"unused_mirna_{i:02d}" for i in range(10)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for index in range(N_MIRNA):
            guide = guides[index % len(guides)]
            row: dict[str, Any] = {
                "qseq": guide,
                "mirna_id": f"SYNTHETIC_MIRNA_{index % 97:03d}",
                "rname": f"SYNTHETIC_MIRNA_{index % 97:03d}",
                "mirna": f"SYNTHETIC_MIRNA_{index % 97:03d}",
                "name": f"SYNTHETIC_MIRNA_{index % 97:03d}",
                "nm": index % 2,
                "seed_mismatches": index % 2,
                "coord": index + 1,
                "species": "human",
                "source": "synthetic resource fixture",
                "database": "synthetic",
            }
            writer.writerow([row.get(column, "x") for column in columns])


def _directory_bytes(path: Path) -> int:
    return sum(file.stat().st_size for file in path.rglob("*") if file.is_file())


def measure() -> dict[str, Any]:
    """Generate a temporary large fixture and return the resulting report measurements."""
    guides, native_observations = _guides()
    with tempfile.TemporaryDirectory(prefix="sirnaforge-report-scale-") as temporary:
        root = Path(temporary)
        run = root / "run"
        aggregate = run / "off_target/results/aggregated"
        (run / "sirnaforge").mkdir(parents=True)
        aggregate.mkdir(parents=True)
        _write_candidate_table(run / "sirnaforge/candidates_all.csv", guides)
        _write_hit_table(aggregate / "combined_offtargets.tsv", guides)
        _write_mirna_table(aggregate / "combined_mirna_hits.tsv", guides)
        (run / "sirnaforge/manifest.json").write_text(
            json.dumps({"gene_query": "synthetic OligoGym scale fixture", "tool_version": "resource-only"})
        )

        report = root / "report.html"
        started = time.perf_counter()
        result = write_report_from_run(run, report, mode="data_backed")
        wall_seconds = time.perf_counter() - started
        manifest = result.manifest
        detail_dir = report.with_name(f"{report.stem}.details")
        artifact_bytes = _directory_bytes(result.artifact_dir) if result.artifact_dir else 0
        sidecar_bytes = _directory_bytes(detail_dir)
        inputs = {
            "candidate_rows": manifest["tables"]["candidates"]["rows"],
            "off_target_rows": manifest["tables"]["hits"]["rows"],
            "mirna_rows": manifest["tables"]["mirna"]["rows"],
            "distinct_guides": manifest["tables"]["candidates"]["distinct_guides"],
            "source_input_bytes": sum(
                int(manifest["tables"][name]["source_bytes"]) for name in ("candidates", "hits", "mirna")
            ),
        }
        if inputs != {
            "candidate_rows": N_CANDIDATES,
            "off_target_rows": N_HITS,
            "mirna_rows": N_MIRNA,
            "distinct_guides": len(guides),
            "source_input_bytes": inputs["source_input_bytes"],
        }:
            raise RuntimeError(f"staged input counts differ from fixture declaration: {inputs}")
        if report.stat().st_size >= 10_000_000 or sidecar_bytes <= report.stat().st_size:
            raise RuntimeError("large report missed the expected bounded-entrypoint/sidecar budget")
        return {
            "scope": "synthetic_resource_stress_only",
            "interpretation": "Native OligoGym guide strings; candidate, hit, and miRNA rows are synthetic.",
            "machine": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "git_revision": subprocess.run(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
                ).stdout.strip(),
            },
            "guide_source": {
                "path": OBSERVATIONS.relative_to(ROOT).as_posix(),
                "sha256": _sha256(OBSERVATIONS),
                "observation_rows": native_observations,
                "distinct_guides": len(guides),
            },
            "inputs": inputs,
            "outputs": {
                "html_bytes": report.stat().st_size,
                "sidecar_bytes": sidecar_bytes,
                "normalized_artifact_bytes": artifact_bytes,
                "detail_shards": len(manifest["detail_shards"]),
                "guide_index_bytes": manifest["guide_index"]["bytes"],
            },
            "diagnostics": {
                "wall_seconds": wall_seconds,
                "peak_rss_bytes": result.diagnostics["peak_rss_bytes"],
                "peak_memory_scope": result.diagnostics["peak_memory_scope"],
                "phases_seconds": result.diagnostics["phases_seconds"],
            },
        }


def main() -> int:
    """Run the resource measurement and emit its compact JSON summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write compact measurement JSON here (default: stdout)")
    args = parser.parse_args()
    result = measure()
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        sys.stdout.write(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
