#!/usr/bin/env python3
"""Exploratory, frozen GC-feature check on Shabalina et al.'s additional 32-guide set.

This is not a profile validation or promotion gate. The records lack verified target identity,
full duplex geometry, chemistry, and liability measurements. It computes only the default GC
feature's univariate association with the reported residual-expression endpoint and checks that the
guide sequences do not occur in the paper's 653-guide training extract.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

from sirnaforge.core.scoring import SCORING_WEIGHT_SET_VERSION
from sirnaforge.models.scoring_profile import SHIPPED_PROFILE
from sirnaforge.models.sirna import ScoringWeights

ROOT = Path(__file__).resolve().parents[1]
PANEL_CSV = ROOT / "tests/data/benchmarks/ml_reviews/built/records.csv"
TRAINING_JSONL = ROOT / "tests/data/benchmarks/ml_reviews/source_extracts/shabalina_2006.jsonl"
DOI = "10.1186/1471-2105-7-65"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = ((start + 1) + end) / 2
        for position in range(start, end):
            ranks[order[position]] = average_rank
        start = end
    return ranks


def _pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) != len(y) or not x:
        return None
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    centered_x = [value - mean_x for value in x]
    centered_y = [value - mean_y for value in y]
    norm = math.sqrt(sum(value**2 for value in centered_x) * sum(value**2 for value in centered_y))
    return sum(a * b for a, b in zip(centered_x, centered_y, strict=True)) / norm if norm else None


def _read_panel(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["collection_id"] == "shabalina_independent"]
    if not rows:
        raise ValueError(f"No shabalina_independent rows found in {path}")
    return rows


def analyze(panel_path: Path = PANEL_CSV, training_path: Path = TRAINING_JSONL) -> dict[str, Any]:
    """Compute frozen univariate metrics and provenance for the independent subset."""
    rows = _read_panel(panel_path)
    training_sequences = {
        json.loads(line)["sequence"].upper().replace("U", "T")
        for line in training_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }

    design_vector = ScoringWeights().design
    gc_weight = float(design_vector.as_mapping()["gc_content"])
    x: list[float] = []
    y: list[float] = []
    groups: dict[str, list[int]] = defaultdict(list)
    sequences: list[str] = []
    for row in rows:
        sequence = row["paired_guide_sequence"].upper().replace("U", "T")
        if len(sequence) != 19 or set(sequence) - set("ACGT"):
            raise ValueError(f"Expected a 19-nt A/C/G/T guide core: {row['record_id']}")
        if row["endpoint"] != "residual_expression" or row["direction"] != "lower_is_better":
            raise ValueError(f"Unexpected endpoint/direction for {row['record_id']}")
        accession = row["target_accession"]
        if not accession:
            raise ValueError(f"Missing source accession for {row['record_id']}")

        gc_percent = 100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence)
        gc_score = math.exp(-(((gc_percent - 40.0) / 10.0) ** 2))
        sequences.append(sequence)
        x.append(gc_score)
        y.append(float(row["measurement_value"]))
        groups[accession].append(len(x) - 1)

    centered_x: list[float] = []
    centered_y: list[float] = []
    for indices in groups.values():
        mean_x = sum(x[index] for index in indices) / len(indices)
        mean_y = sum(y[index] for index in indices) / len(indices)
        centered_x.extend(x[index] - mean_x for index in indices)
        centered_y.extend(y[index] - mean_y for index in indices)

    slope_denom = sum(value**2 for value in centered_x)
    fixed_effect_slope = (
        sum(a * b for a, b in zip(centered_x, centered_y, strict=True)) / slope_denom if slope_denom else None
    )
    return {
        "analysis_id": "shabalina-independent-gc-exploratory-v1",
        "purpose": "exploratory single-feature efficacy check; not profile validation or promotion evidence",
        "source": {"doi": DOI, "dataset": "shabalina_independent", "endpoint": "residual_expression_percent"},
        "profile": {
            "profile_id": SHIPPED_PROFILE.profile_id,
            "weight_set_version": SCORING_WEIGHT_SET_VERSION,
            "feature": "gc_content",
            "feature_formula": "exp(-((GC_percent - 40) / 10)^2)",
            "design_vector_weight": gc_weight,
        },
        "inputs": {
            "panel_csv_sha256": _sha256(panel_path),
            "training_jsonl_sha256": _sha256(training_path),
        },
        "counts": {
            "observations": len(rows),
            "distinct_guides": len(set(sequences)),
            "exact_guide_overlaps_with_653_training_records": sum(
                sequence in training_sequences for sequence in sequences
            ),
            "source_accession_groups": len(groups),
            "group_sizes": {key: len(groups[key]) for key in sorted(groups)},
        },
        "outcome_direction": "lower residual expression indicates greater efficacy",
        "metrics": {
            "pooled_spearman_gc_score_vs_residual_expression": _pearson(_average_ranks(x), _average_ranks(y)),
            "accession_fixed_effect_gc_score_slope_residual_percent_per_score_unit": fixed_effect_slope,
            "accession_centered_pearson_gc_score_vs_residual_expression": _pearson(centered_x, centered_y),
        },
        "limitations": [
            "32 observations and 10 source accession groups; point estimates only, no inferential claim",
            "all records are paired-core-only; full duplex geometry and chemistry are not stated",
            "target identity is unestablished; accessions are source grouping labels, not verified gene clusters",
            "no target accessibility or unintended-repression endpoint is available",
            "the 2006 paper's own reported model metrics are not siRNAforge validation results",
        ],
    }


def main() -> None:
    """Run the analysis and print or write its compact JSON report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional path for the JSON summary")
    args = parser.parse_args()
    result = analyze()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
