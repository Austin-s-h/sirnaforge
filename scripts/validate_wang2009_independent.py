#!/usr/bin/env python
"""Run a fixed, exploratory feature check on Wang et al.'s validation panel.

This reports exploratory feature associations and target-held-out ranking. The panel has no
unintended-repression endpoint. Fold-local coefficients are not fitted scoring-profile weights and
are never promoted. Repeated doses are reported separately; the pooled target/algorithm-adjusted
correlations are OLS residual correlations.

Run with the project environment:

    .venv/bin/python scripts/validate_wang2009_independent.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import cast

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
    return cast(np.ndarray, values - design @ np.linalg.lstsq(design, values, rcond=None)[0])


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
        au_score = au_content_5p_score(guide)
        if au_score is None:
            raise ValueError(f"Guide too short for A/U(1–5) score: {row['source_record_id']}")
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
                "residual_expression": float(str(row["measurement_raw"])),
                "au_1_5": au_score,
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


def unique_guide_rows(rows: list[dict[str, object]], dose: str) -> list[dict[str, object]]:
    """Collapse repeated design-group entries to one target/guide at the selected dose."""
    by_guide: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        if row["dose"] == dose:
            key = (str(row["target"]), str(row["guide"]))
            by_guide.setdefault(key, []).append(row)

    collapsed: list[dict[str, object]] = []
    for (target, guide), matches in sorted(by_guide.items()):
        first = matches[0]
        collapsed.append(
            {
                "target": target,
                "guide": guide,
                "residual_expression": float(np.mean([float(str(row["residual_expression"])) for row in matches])),
                "au_1_5": float(str(first["au_1_5"])),
                "asymmetry": float(str(first["asymmetry"])),
                "gc_feature": float(str(first["gc_feature"])),
                # If a guide was selected by multiple algorithms, retain its fractional group mix.
                "algorithm_mix": {
                    algorithm: sum(str(row["algorithm"]) == algorithm for row in matches) / len(matches)
                    for algorithm in ("New algorithm", "Previous algorithm", "Random")
                },
            }
        )
    return collapsed


def bootstrap_mean(values: np.ndarray, seed: int = 97) -> tuple[float, float, float]:
    """Mean and percentile interval from a paired, target-cluster bootstrap."""
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(values), size=(10_000, len(values)))
    means = values[draws].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(values.mean()), float(low), float(high)


def report_leave_one_target_out(rows: list[dict[str, object]], dose: str) -> None:
    """Compare A/U and asymmetry in target-held-out ranking with target-cluster uncertainty."""
    panel = unique_guide_rows(rows, dose)
    targets = sorted({str(row["target"]) for row in panel})
    sequence_targets: dict[str, set[str]] = {}
    for row in panel:
        sequence_targets.setdefault(str(row["guide"]), set()).add(str(row["target"]))
    leaked = {sequence: groups for sequence, groups in sequence_targets.items() if len(groups) > 1}
    if leaked:
        raise ValueError(f"Guide sequences occur across target folds: {sorted(leaked)[:5]}")

    models = {
        "AU only": ("au_1_5",),
        "asymmetry only": ("asymmetry",),
        "AU + asymmetry": ("au_1_5", "asymmetry"),
    }
    observed_by_target: dict[str, np.ndarray] = {}
    feature_by_target: dict[str, dict[str, np.ndarray]] = {}
    mix_by_target: dict[str, np.ndarray] = {}
    for target in targets:
        group = [row for row in panel if row["target"] == target]
        observed_by_target[target] = np.asarray([row["residual_expression"] for row in group], dtype=float)
        feature_by_target[target] = {
            feature: np.asarray([row[feature] for row in group], dtype=float) for feature in ("au_1_5", "asymmetry")
        }
        mix_by_target[target] = np.asarray(
            [
                [
                    cast(dict[str, float], row["algorithm_mix"])[algorithm]
                    for algorithm in ("New algorithm", "Previous algorithm")
                ]
                for row in group
            ],
            dtype=float,
        )

    per_target_rho: dict[str, dict[str, float]] = {name: {} for name in models}
    for held_out in targets:
        train_targets = [target for target in targets if target != held_out]
        for model_name, selected_features in models.items():
            train_x: list[np.ndarray] = []
            train_y: list[np.ndarray] = []
            for target in train_targets:
                x = np.column_stack(
                    [feature_by_target[target][feature] for feature in selected_features] + [mix_by_target[target]]
                )
                train_x.append(x - x.mean(axis=0))
                y = observed_by_target[target]
                train_y.append(y - y.mean())
            x_fit = np.vstack(train_x)
            y_fit = np.concatenate(train_y)
            coefficients = np.linalg.lstsq(x_fit, y_fit, rcond=None)[0]

            test_x = np.column_stack(
                [feature_by_target[held_out][feature] for feature in selected_features] + [mix_by_target[held_out]]
            )
            test_prediction = (test_x - test_x.mean(axis=0)) @ coefficients
            per_target_rho[model_name][held_out] = spearman(test_prediction, observed_by_target[held_out])

    print(f"\n{dose}: leave-one-target-out validation on {len(targets)} targets")
    print(
        f"  one row per target/guide ({len(panel)} guides); repeated doses kept separate; "
        "algorithm selection mix used as a nuisance covariate"
    )
    print("  per-target held-out Spearman: target | n guides | AU | asymmetry | both")
    for target in targets:
        n_guides = len(observed_by_target[target])
        target_values = [per_target_rho[name][target] for name in models]
        print(f"    {target} | {n_guides} | " + " | ".join(f"{value:+.3f}" for value in target_values))
    for name in models:
        metric_values = np.asarray([per_target_rho[name][target] for target in targets])
        mean, low, high = bootstrap_mean(metric_values)
        print(f"  {name}: mean target Spearman rho={mean:+.3f} (95% target-bootstrap CI {low:+.3f}, {high:+.3f})")
    for comparator in ("AU only", "asymmetry only"):
        diff = np.asarray(
            [per_target_rho["AU + asymmetry"][target] - per_target_rho[comparator][target] for target in targets]
        )
        mean, low, high = bootstrap_mean(diff)
        print(
            f"  AU + asymmetry minus {comparator}: Δrho={mean:+.3f} (95% paired target-bootstrap CI {low:+.3f}, {high:+.3f})"
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
        report_leave_one_target_out(features, dose)
    print(
        "\nLimits: leave-one-target-out prediction and intervals remain exploratory, with only "
        "14 target clusters and this one assay. The panel has no unintended-repression endpoint "
        "and cannot validate the complete default score or justify weight promotion."
    )


if __name__ == "__main__":
    main()
