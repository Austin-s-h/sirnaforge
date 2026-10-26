#!/usr/bin/env python
"""Run a fixed, exploratory feature check on Wang et al.'s validation panel.

This reports associations only. The panel has no unintended-repression endpoint and this script
does not fit or promote scoring weights. Repeated doses are reported separately; the descriptive
target/algorithm-adjusted correlations are OLS residual correlations, not held-out estimates.

Run with the project environment:

    .venv/bin/python scripts/validate_wang2009_independent.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from sirnaforge.core.design import SiRNADesigner, au_content_5p_score
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.scoring_profile import SHIPPED_PROFILE
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "tests/data/benchmarks/ml_reviews/source_extracts/wang_2009.jsonl"
TRAINING = ROOT / "tests/data/benchmarks/ml_reviews/source_extracts/shabalina_2006.jsonl"
FEATURE_CODE = (
    Path(__file__).resolve(),
    ROOT / "src/sirnaforge/core/design.py",
    ROOT / "src/sirnaforge/core/thermodynamics.py",
    ROOT / "src/sirnaforge/models/scoring_profile.py",
)


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_jsonl(path: Path) -> list[dict[str, object]]:
    """Read one JSON object per line."""
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def ranks(values: np.ndarray) -> np.ndarray:
    """Compute average ranks for ties, using one-based ranks."""
    order = np.argsort(values, kind="mergesort")
    result = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        result[order[start:stop]] = (start + stop - 1) / 2 + 1
        start = stop
    return result


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    """Pearson correlation, returning NaN for a constant input."""
    if np.ptp(left) == 0 or np.ptp(right) == 0:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    """Spearman correlation with average ranks for tied observations."""
    return correlation(ranks(left), ranks(right))


def residualize(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Remove the least-squares intercept, target-gene and algorithm-group effects."""
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def make_features(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Calculate the shipped sequence features on the reported 19-nt paired cores."""
    thermodynamics = ThermodynamicCalculator()
    designer = SiRNADesigner(DesignParameters())
    result: list[dict[str, object]] = []

    for row in rows:
        guide = str(row["guide_sequence"])[:19].upper().replace("T", "U")
        passenger = str(row["passenger_sequence"])[:19].upper().replace("T", "U")
        if len(guide) != 19 or len(passenger) != 19:
            raise ValueError(f"Expected 19-nt paired cores: {row['source_record_id']}")
        gc_percent = sum(base in "GC" for base in guide) / len(guide) * 100.0
        candidate = SiRNACandidate(
            id=str(row["source_record_id"]),
            transcript_id=str(row["target_gene"]),
            position=1,
            guide_sequence=guide,
            passenger_sequence=passenger,
            gc_content=gc_percent,
            length=len(guide),
            asymmetry_score=0.0,
        )
        result.append(
            {
                "guide": guide,
                "target": str(row["target_gene"]),
                "algorithm": str(row["design_algorithm"]),
                "dose": str(row["concentration"]),
                "residual_expression": float(row["measurement_raw"]),
                "au_1_5": float(au_content_5p_score(guide)),
                "asymmetry": float(thermodynamics.calculate_asymmetry_score(candidate)[2]),
                "gc_feature": float(designer._calculate_gc_score(gc_percent)),
            }
        )
    return result


def report_dose(rows: list[dict[str, object]], dose: str) -> None:
    """Print pooled rank associations and descriptive adjusted correlations for one dose."""
    panel = [row for row in rows if row["dose"] == dose]
    y = np.asarray([row["residual_expression"] for row in panel], dtype=float)
    targets = sorted({str(row["target"]) for row in panel})
    algorithms = sorted({str(row["algorithm"]) for row in panel})
    design = np.column_stack(
        [np.ones(len(panel))]
        + [[float(row["target"] == target) for row in panel] for target in targets[1:]]
        + [[float(row["algorithm"] == algorithm) for row in panel] for algorithm in algorithms[1:]]
    )
    y_residual = residualize(y, design)

    print(f"\n{dose}: observations={len(panel)}, unique_guides={len({row['guide'] for row in panel})}")
    print(f"  residual expression range: {min(y):g}–{max(y):g}% (lower means more knockdown)")
    for feature in ("au_1_5", "asymmetry", "gc_feature"):
        x = np.asarray([row[feature] for row in panel], dtype=float)
        x_residual = residualize(x, design)
        print(
            f"  {feature}: pooled Spearman rho={spearman(x, y):+.3f}; "
            f"target+algorithm-adjusted Pearson r={correlation(x_residual, y_residual):+.3f}"
        )


def main() -> None:
    """Validate source shape and print reproducible exploratory statistics."""
    source_rows = read_jsonl(PANEL)
    train_rows = read_jsonl(TRAINING)
    if any(row.get("collection_id") != "wang_2009" for row in source_rows):
        raise ValueError("Unexpected collection in Wang panel")
    if any(row.get("endpoint") != "residual_expression" for row in source_rows):
        raise ValueError("Unexpected endpoint in Wang panel")
    if any(row.get("chemistry") != "RNA_core_DNA_overhangs" for row in source_rows):
        raise ValueError("Unexpected chemistry in Wang panel")

    features = make_features(source_rows)
    panel_guides = {str(row["guide"]) for row in features}
    training_guides: set[str] = set()
    for row in train_rows:
        sequence = str(row.get("sequence", "")).upper().replace("T", "U")
        if len(sequence) != 19 or set(sequence) - set("ACGU"):
            raise ValueError(f"Invalid Shabalina training guide sequence: {row.get('source_record_id')}")
        training_guides.add(sequence)
    overlap = panel_guides & training_guides

    profile_payload = {
        "profile": SHIPPED_PROFILE.model_dump(mode="json"),
        "term_records": [term.model_dump(mode="json") for term in SHIPPED_PROFILE.term_records()],
    }
    profile_json = json.dumps(profile_payload, sort_keys=True, separators=(",", ":"))
    profile_hash = hashlib.sha256(profile_json.encode()).hexdigest()
    print("Source: Wang et al., Nucleic Acids Research (2009), doi:10.1093/nar/gkp864")
    print(f"Panel SHA-256: {sha256(PANEL)}")
    print(f"Shabalina 653-guide extract SHA-256: {sha256(TRAINING)}")
    print(f"Profile: {SHIPPED_PROFILE.profile_id}; canonical profile SHA-256: {profile_hash}")
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()
    print(f"Git revision: {revision}")
    for path in FEATURE_CODE:
        print(f"Implementation {path.relative_to(ROOT)} SHA-256: {sha256(path)}")
    print(f"Observations: {len(features)}; unique paired-core guides: {len(panel_guides)}")
    print(f"Exact guide overlap with Shabalina 653-guide extract: {len(overlap)}")
    print("All observations use the same RNA-core/DNA-overhang chemistry and residual-expression endpoint.")
    print("Features: frozen guide A/U(1–5), shipped asymmetry, and shipped GC transform.")
    for dose in ("3nM", "30nM"):
        report_dose(features, dose)
    print(
        "\nLimits: repeated-dose observations are correlated; these are point associations without "
        "clustered uncertainty or held-out prediction. The panel has no unintended-repression "
        "endpoint and cannot validate the complete default score or justify weight promotion."
    )


if __name__ == "__main__":
    main()
