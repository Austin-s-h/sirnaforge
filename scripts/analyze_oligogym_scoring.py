#!/usr/bin/env python3
"""Explore OligoGym activity associations and alternative siRNA score weights.

This is an exploratory analysis, not a training step for the production profile. It joins native
benchmark mappings to the completed full-workflow candidate tables, collapses isoform-level rows to
one row per measured observation, and reports panel-stratified associations, simple calibration,
GC-bound sensitivity, and a leave-one-panel-out weight stress test.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "tests/data/benchmarks/oligogym_native_design"
RUNS = {
    "ichihara_2007_1": "ichihara_2007_1_qualified_all_representative_orthologs",
    "ichihara_2007_2": "ichihara_2007_2_qualified_all_representative_orthologs",
    "martinelli_2023_1": "martinelli_2023_1_qualified_all_representative_unqualified_target",
}
TERMS = ["off_target", "target_accessibility", "asymmetry", "gc_content"]
BASELINE_WEIGHTS = np.array([0.25, 0.30, 0.25, 0.20])
COMPONENT_COLUMNS = [
    "score_off_target",
    "score_target_accessibility",
    "score_asymmetry",
    "score_gc_content",
]


def spearman_value(x: object, y: object) -> float:
    """Calculate Spearman rho with pandas ranks, without an optional SciPy dependency."""
    x_series = pd.Series(x).reset_index(drop=True)
    y_series = pd.Series(y).reset_index(drop=True)
    if len(x_series) < 3 or x_series.nunique() < 2 or y_series.nunique() < 2:
        return float("nan")
    ranked_x = x_series.rank(method="average")
    ranked_y = y_series.rank(method="average")
    return float(ranked_x.corr(ranked_y))


def bootstrap_corr(frame: pd.DataFrame, column: str, iterations: int = 2000, seed: int = 20260913) -> dict:
    """Bootstrap a Spearman interval at the observation level."""
    frame = frame[[column, "label_processed"]].dropna().reset_index(drop=True)
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(iterations):
        sample = frame.iloc[rng.integers(0, len(frame), len(frame))]
        value = spearman_value(sample[column], sample["label_processed"])
        if np.isfinite(value):
            values.append(value)
    if not values:
        return {"iterations": iterations, "valid_iterations": 0, "rho_median": None, "rho_p025": None, "rho_p975": None}
    quantiles = np.quantile(values, [0.025, 0.5, 0.975])
    return {
        "iterations": iterations,
        "valid_iterations": len(values),
        "rho_median": float(quantiles[1]),
        "rho_p025": float(quantiles[0]),
        "rho_p975": float(quantiles[2]),
    }


def corr(x: pd.Series, y: pd.Series) -> dict[str, float | int | None]:
    """Calculate paired Pearson and Spearman correlations after dropping missing values."""
    frame = pd.concat([x, y], axis=1).dropna()
    if len(frame) < 3 or frame.iloc[:, 0].nunique() < 2 or frame.iloc[:, 1].nunique() < 2:
        return {"n": int(len(frame)), "spearman_rho": None, "pearson_r": None}
    ranked_x = frame.iloc[:, 0].rank(method="average")
    ranked_y = frame.iloc[:, 1].rank(method="average")
    return {
        "n": int(len(frame)),
        "spearman_rho": float(ranked_x.corr(ranked_y)),
        "pearson_r": float(frame.iloc[:, 0].corr(frame.iloc[:, 1])),
    }


def load_linked(
    native_dir: Path,
    workflow_root: Path,
    runs: dict[str, str],
    source_weights: np.ndarray,
) -> pd.DataFrame:
    """Join measured observations to candidates and normalize score components."""
    links = pd.read_csv(native_dir / "candidate_benchmark_links.csv")
    observations = pd.read_csv(native_dir / "observations.csv", usecols=["benchmark_id", "dataset", "label_processed"])
    rows = []
    candidate_columns = [
        "id",
        "gc_content",
        "max_poly_run_length",
        "composite_score",
        *COMPONENT_COLUMNS,
    ]
    for dataset, run in runs.items():
        candidate_path = workflow_root / run / "sirnaforge/candidates_all.csv"
        candidates = pd.read_csv(candidate_path, usecols=candidate_columns)
        panel_links = links.loc[links["dataset"].eq(dataset), ["dataset", "benchmark_id", "candidate_id"]]
        joined = panel_links.merge(observations, on=["benchmark_id", "dataset"], how="left")
        joined = joined.merge(candidates, left_on="candidate_id", right_on="id", how="inner")
        joined["label_processed"] = pd.to_numeric(joined["label_processed"], errors="coerce")
        rows.append(joined)
    linked: pd.DataFrame = pd.concat(rows, ignore_index=True)
    for term, column, weight in zip(TERMS, COMPONENT_COLUMNS, source_weights, strict=True):
        linked[f"component_{term}"] = pd.to_numeric(linked[column], errors="coerce") / weight
    linked["current_score"] = linked[[f"component_{term}" for term in TERMS]].mul(source_weights).sum(axis=1)
    return linked


def aggregate(linked: pd.DataFrame, method: str = "median") -> pd.DataFrame:
    """Collapse candidate scores to one row per measured observation."""
    score_columns = [f"component_{term}" for term in TERMS]
    aggregations: dict[str, str] = {"label_processed": "first", "gc_content": method, "max_poly_run_length": "max"}
    aggregations.update(dict.fromkeys(score_columns, method))
    aggregations["current_score"] = method
    result = linked.groupby(["dataset", "benchmark_id"], as_index=False).agg(aggregations)
    result["candidate_rows"] = linked.groupby(["dataset", "benchmark_id"]).size().to_numpy()
    result["alternative_score"] = result[score_columns].mul(BASELINE_WEIGHTS).sum(axis=1)
    return result


def panel_correlations(aggregated: pd.DataFrame) -> dict[str, dict[str, dict]]:
    """Report score and component correlations per panel and across pooled data."""
    output: dict[str, dict[str, dict]] = {}
    for dataset, frame in aggregated.groupby("dataset"):
        metrics = {"current_score": corr(frame["current_score"], frame["label_processed"])}
        for term in TERMS:
            metrics[term] = corr(frame[f"component_{term}"], frame["label_processed"])
        output[str(dataset)] = metrics
    pooled = {"current_score": corr(aggregated["current_score"], aggregated["label_processed"])}
    for term in TERMS:
        pooled[term] = corr(aggregated[f"component_{term}"], aggregated["label_processed"])
    output["pooled"] = pooled
    return output


def component_variation(aggregated: pd.DataFrame) -> dict[str, dict[str, dict[str, float | int]]]:
    """Summarize the diversity and spread of component scores per panel."""
    output: dict[str, dict[str, dict[str, float | int]]] = {}
    for dataset, frame in aggregated.groupby("dataset"):
        output[str(dataset)] = {}
        for term in TERMS:
            series = frame[f"component_{term}"]
            output[str(dataset)][term] = {
                "n": int(series.notna().sum()),
                "unique_values": int(series.nunique()),
                "sd": float(series.std(ddof=1)) if len(series) > 1 else 0.0,
            }
    return output


def calibrated(aggregated: pd.DataFrame) -> dict[str, dict]:
    """Fit a linear score-to-activity calibration separately for each panel."""
    output: dict[str, dict] = {}
    for dataset, frame in aggregated.groupby("dataset"):
        x = frame["current_score"].to_numpy(float)
        y = frame["label_processed"].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 3:
            continue
        slope, intercept = np.polyfit(x, y, 1)
        predicted = intercept + slope * x
        ss_res = float(np.sum((y - predicted) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        output[str(dataset)] = {
            "n": int(len(x)),
            "intercept": float(intercept),
            "slope": float(slope),
            "rmse": float(np.sqrt(np.mean((y - predicted) ** 2))),
            "r_squared": float(1 - ss_res / ss_tot) if ss_tot else None,
            "activity_mean": float(y.mean()),
            "activity_sd": float(y.std(ddof=1)),
            "score_mean": float(x.mean()),
            "score_sd": float(x.std(ddof=1)),
        }
    return output


def score_for_weights(frame: pd.DataFrame, weights: np.ndarray) -> np.ndarray:
    """Combine normalized score components using the supplied weights."""
    columns = [f"component_{term}" for term in TERMS]
    return np.asarray(frame[columns].to_numpy(float) @ weights)


def mean_panel_rho(frame: pd.DataFrame, weights: np.ndarray, datasets: list[str]) -> float:
    """Average finite Spearman correlations across the selected panels."""
    scores = score_for_weights(frame, weights)
    values = []
    for dataset in datasets:
        mask = frame["dataset"].to_numpy() == dataset
        if mask.sum() < 3:
            continue
        value = spearman_value(scores[mask], frame.loc[mask, "label_processed"])
        if np.isfinite(value):
            values.append(float(value))
    return float(np.mean(values)) if values else float("nan")


def weight_grid(frame: pd.DataFrame, datasets: list[str], step: float = 0.02) -> dict:
    """Search nonnegative weights summing to one for the best panel correlation."""
    best: dict[str, Any] = {"objective_mean_panel_spearman": -np.inf, "weights": None, "panel_rho": {}}
    values = np.arange(0.0, 1.0 + step / 2, step)
    for a in values:
        for b in values:
            for c in values:
                d = 1.0 - a - b - c
                if d < -1e-9:
                    continue
                weights = np.array([a, b, c, max(0.0, d)])
                objective = mean_panel_rho(frame, weights, datasets)
                if objective > best["objective_mean_panel_spearman"]:
                    scores = score_for_weights(frame, weights)
                    panel_rho = {}
                    for dataset in datasets:
                        mask = frame["dataset"].to_numpy() == dataset
                        panel_rho[dataset] = spearman_value(scores[mask], frame.loc[mask, "label_processed"])
                    best = {
                        "objective_mean_panel_spearman": objective,
                        "weights": {term: float(weight) for term, weight in zip(TERMS, weights, strict=False)},
                        "panel_rho": panel_rho,
                        "grid_step": step,
                    }
    return best


def leave_one_panel_out(frame: pd.DataFrame, fit_datasets: list[str]) -> dict:
    """Fit weights on selected panels and evaluate the remaining panel."""
    holdout = next(dataset for dataset in RUNS if dataset not in fit_datasets)
    fit = weight_grid(frame.loc[frame["dataset"].isin(fit_datasets)], fit_datasets, step=0.02)
    weights = np.array([fit["weights"][term] for term in TERMS])
    mask = frame["dataset"].eq(holdout)
    score = score_for_weights(frame.loc[mask], weights)
    rho = spearman_value(score, frame.loc[mask, "label_processed"])
    return {"fit_panels": fit_datasets, "holdout_panel": holdout, "fit": fit, "holdout_spearman_rho": float(rho)}


def fixed_off_target_grid(frame: pd.DataFrame, datasets: list[str], fixed: float = 0.25, step: float = 0.02) -> dict:
    """Search efficacy terms while preserving the declared specificity allocation."""
    best: dict[str, Any] = {"objective_mean_panel_spearman": -np.inf, "weights": None, "panel_rho": {}}
    remaining = 1.0 - fixed
    values = np.arange(0.0, remaining + step / 2, step)
    for accessibility in values:
        for asymmetry in values:
            gc = remaining - accessibility - asymmetry
            if gc < -1e-9:
                continue
            weights = np.array([fixed, accessibility, asymmetry, max(0.0, gc)])
            objective = mean_panel_rho(frame, weights, datasets)
            if objective > best["objective_mean_panel_spearman"]:
                scores = score_for_weights(frame, weights)
                panel_rho = {
                    dataset: spearman_value(
                        scores[frame["dataset"].to_numpy() == dataset],
                        frame.loc[frame["dataset"].eq(dataset), "label_processed"],
                    )
                    for dataset in datasets
                }
                best = {
                    "objective_mean_panel_spearman": objective,
                    "weights": {term: float(weight) for term, weight in zip(TERMS, weights, strict=False)},
                    "panel_rho": panel_rho,
                    "grid_step": step,
                    "fixed_off_target": fixed,
                }
    return best


def profile_result(frame: pd.DataFrame, weights: list[float]) -> dict:
    """Evaluate a fixed scoring profile on each benchmark panel."""
    scores = score_for_weights(frame, np.array(weights))
    panel_rho = {
        dataset: spearman_value(
            scores[frame["dataset"].to_numpy() == dataset], frame.loc[frame["dataset"].eq(dataset), "label_processed"]
        )
        for dataset in RUNS
    }
    return {"weights": dict(zip(TERMS, weights, strict=False)), "panel_rho": panel_rho}


def gc_sensitivity(linked: pd.DataFrame, aggregated: pd.DataFrame) -> dict:
    """Measure retention and score correlations under alternative GC bounds."""
    output: dict[str, dict] = {}
    for lower, upper in [(0, 100), (25, 70), (30, 65), (35, 60)]:
        key = f"gc_{lower:g}_{upper:g}"
        panel_result = {}
        filtered = linked.loc[
            linked["gc_content"].between(lower, upper, inclusive="both") & (linked["max_poly_run_length"] <= 3)
        ]
        filtered_obs = aggregate(filtered, "median")
        for dataset in RUNS:
            original = aggregated.loc[aggregated["dataset"].eq(dataset), "benchmark_id"]
            current = filtered_obs.loc[filtered_obs["dataset"].eq(dataset)]
            panel_result[dataset] = {
                "candidate_rows": int(filtered.loc[filtered["dataset"].eq(dataset)].shape[0]),
                "observations_with_candidate": int(current.shape[0]),
                "observations_total_mapped": int(original.nunique()),
                "observation_retention_fraction": float(current.shape[0] / original.nunique())
                if original.nunique()
                else None,
                "median_candidates_per_observation": float(current["candidate_rows"].median())
                if len(current)
                else None,
                "current_score_spearman": corr(current["current_score"], current["label_processed"])
                if len(current)
                else None,
            }
        output[key] = {"lower": lower, "upper": upper, "panels": panel_result}
    return output


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "benchmark_results/oligogym_scoring_analysis.json")
    parser.add_argument("--native-dir", type=Path, default=NATIVE)
    parser.add_argument("--workflow-root", type=Path, default=ROOT / "workflow_results/oligogym_native_design")
    parser.add_argument("--run", action="append", default=[], metavar="DATASET=DIRECTORY")
    parser.add_argument(
        "--component-weights",
        type=float,
        nargs=4,
        default=BASELINE_WEIGHTS.tolist(),
        help="Weights used in source CSVs, ordered off-target, accessibility, asymmetry, GC.",
    )
    args = parser.parse_args()
    runs = dict(RUNS)
    for assignment in args.run:
        dataset, separator, directory = assignment.partition("=")
        if not separator or dataset not in runs or not directory:
            parser.error("--run must name a known DATASET=DIRECTORY")
        runs[dataset] = directory
    source_weights = np.asarray(args.component_weights)
    if (
        not np.all(np.isfinite(source_weights))
        or np.any(source_weights <= 0)
        or not np.isclose(source_weights.sum(), 1)
    ):
        parser.error("--component-weights must be positive finite weights summing to one")
    native_dir = args.native_dir
    workflow_root = args.workflow_root
    linked = load_linked(native_dir, workflow_root, runs, source_weights)
    median = aggregate(linked, "median")
    maximum = aggregate(linked, "max")
    fit_panels = ["ichihara_2007_1", "ichihara_2007_2"]
    pooled = median.loc[median["dataset"].isin(fit_panels)]
    result = {
        "analysis": "Exploratory OligoGym scoring analysis",
        "source": {
            "candidate_links": str(native_dir / "candidate_benchmark_links.csv"),
            "observations": str(native_dir / "observations.csv"),
            "completed_runs": {dataset: str(workflow_root / run) for dataset, run in runs.items()},
        },
        "observation_counts": {
            dataset: int(frame["benchmark_id"].nunique()) for dataset, frame in median.groupby("dataset")
        },
        "aggregation": {
            "primary": "median across exact native transcript-site candidate rows per benchmark observation",
            "sensitivity": "maximum across exact native transcript-site candidate rows",
            "candidate_rows_linked": int(len(linked)),
        },
        "source_component_weights": dict(zip(TERMS, source_weights.tolist(), strict=True)),
        "baseline_weights": {term: float(weight) for term, weight in zip(TERMS, BASELINE_WEIGHTS, strict=False)},
        "panel_correlations_median": panel_correlations(median),
        "panel_correlations_maximum": panel_correlations(maximum),
        "component_variation_median": component_variation(median),
        "bootstrap_current_score_median": {
            dataset: bootstrap_corr(frame, "current_score") for dataset, frame in median.groupby("dataset")
        },
        "calibration_median_current_score": calibrated(median),
        "weight_optimization": {
            "fit_panels": fit_panels,
            "objective": "equal-weight mean of panel Spearman correlations; nonnegative weights summing to 1",
            "grid": weight_grid(pooled, fit_panels, step=0.02),
            "fixed_off_target_grid": fixed_off_target_grid(pooled, fit_panels, fixed=0.25, step=0.02),
            "leave_one_panel_out": [
                leave_one_panel_out(median, ["ichihara_2007_1"]),
                leave_one_panel_out(median, ["ichihara_2007_2"]),
            ],
        },
        "candidate_profiles": {
            "current": profile_result(median, [0.25, 0.30, 0.25, 0.20]),
            "fixed_off_target_rounded": profile_result(median, [0.25, 0.10, 0.30, 0.35]),
            "unconstrained_rounded": profile_result(median, [0.10, 0.20, 0.30, 0.40]),
        },
        "gc_sensitivity": gc_sensitivity(linked, median),
        "expanded_reference_comparison": {
            "status": "not_estimable",
            "reason": "No successful full-expanded-isoform off-target run was retained; the final runs used deterministic representative references.",
        },
        "caveats": [
            "The analysis is observation-level after transcript-site aggregation; it is not a candidate-row test.",
            "The two Ichihara panels are the primary fitting data because Martinelli target identity is unavailable and its off-target orthology is intentionally unqualified.",
            "The weight search is exploratory and uses a coarse 0.02 grid; it is not a validated production fit.",
            "Calibration is descriptive in-sample linear calibration, not an independently validated activity predictor.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
