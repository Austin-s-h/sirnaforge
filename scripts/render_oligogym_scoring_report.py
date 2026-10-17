#!/usr/bin/env python3
"""Render reproducible SVG visuals for the OligoGym scoring analysis.

The benchmark image intentionally has pandas but not a plotting stack.  This script therefore
emits small, dependency-light SVG figures directly from the analysis JSON and observation tables.
It is intended for documentation, not for changing production scoring behavior.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd

PANEL_ORDER = ["ichihara_2007_1", "ichihara_2007_2", "martinelli_2023_1"]
PANEL_LABELS = {
    "ichihara_2007_1": "Ichihara 1 (2007)",
    "ichihara_2007_2": "Ichihara 2 (2007)",
    "martinelli_2023_1": "Martinelli (2023)",
}
COMPONENTS = ["off_target", "target_accessibility", "asymmetry", "gc_content"]
COMPONENT_LABELS = {
    "off_target": "Off-target",
    "target_accessibility": "Accessibility",
    "asymmetry": "Asymmetry",
    "gc_content": "GC content",
}
PROFILE_ORDER = ["current", "fixed_off_target_rounded", "unconstrained_rounded"]
PROFILE_LABELS = {
    "current": "Historical baseline",
    "fixed_off_target_rounded": "Pre-production",
    "unconstrained_rounded": "Unconstrained",
}
COLORS = {
    "ink": "#243447",
    "muted": "#5c6b73",
    "grid": "#d8e0e5",
    "blue": "#2f6f9f",
    "teal": "#2a9d8f",
    "orange": "#e58f3d",
    "red": "#c95c54",
    "purple": "#7461a8",
    "green": "#5f9e6e",
    "white": "#ffffff",
}


def esc(value: object) -> str:
    """Escape a value for use in SVG text or attributes."""
    return html.escape(str(value), quote=True)


def svg_start(width: int, height: int, title: str, description: str) -> list[str]:
    """Create an SVG header with an accessible title, description, and styles."""
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f'<title id="title">{esc(title)}</title>',
        f'<desc id="desc">{esc(description)}</desc>',
        "<style>"
        "text{font-family:Arial,Helvetica,sans-serif;fill:#243447}"
        ".title{font-size:24px;font-weight:700}.subtitle{font-size:13px;fill:#5c6b73}"
        ".axis{font-size:12px;fill:#5c6b73}.label{font-size:13px}"
        ".small{font-size:11px;fill:#5c6b73}.value{font-size:12px;font-weight:700}"
        "</style>",
    ]


def text(parts: list[str], x: float, y: float, value: object, cls: str = "label", anchor: str = "start") -> None:
    """Append an escaped text element to the SVG."""
    parts.append(f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}">{esc(value)}</text>')


def rect(
    parts: list[str], x: float, y: float, width: float, height: float, fill: str, stroke: str = "none", rx: int = 0
) -> None:
    """Append a rectangle with the requested fill, stroke, and corner radius."""
    parts.append(
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" '
        f'fill="{fill}" stroke="{stroke}" rx="{rx}" />'
    )


def line(
    parts: list[str],
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    stroke: str = COLORS["grid"],
    width: int = 1,
    dash: str = "",
) -> None:
    """Append a line segment to the SVG."""
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    parts.append(
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{width}"{extra} />'
    )


def finish(parts: list[str], path: Path) -> None:
    """Close the SVG document and write it to disk."""
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def scale(value: float, low: float, high: float, start: float, end: float) -> float:
    """Map a value linearly from a data interval to a display interval."""
    return start + (value - low) / (high - low) * (end - start)


def render_overview(analysis: dict, observations: pd.DataFrame, links: pd.DataFrame, out: Path) -> None:
    """Render panel sample counts, mapping coverage, and measured activity."""
    width, height = 1200, 425
    parts = svg_start(width, height, "OligoGym panel overview", "Panel-level sample, mapping, and activity summary.")
    text(parts, 40, 42, "OligoGym panel overview", "title")
    text(parts, 40, 67, "Full native-reference workflow; observation-level joins used for scoring analysis", "subtitle")
    cards = [(40, 95, 350), (425, 95, 350), (810, 95, 350)]
    colors = [COLORS["blue"], COLORS["teal"], COLORS["orange"]]
    for panel, (x, y, w), color in zip(PANEL_ORDER, cards, colors, strict=False):
        obs = observations[observations["dataset"] == panel]
        lnk = links[links["dataset"] == panel]
        mapped = analysis["observation_counts"][panel]
        cal = analysis["calibration_median_current_score"][panel]
        rect(parts, x, y, w, 280, COLORS["white"], COLORS["grid"], 10)
        rect(parts, x, y, w, 48, color, color, 10)
        rect(parts, x, y + 38, w, 10, color, color)
        text(parts, x + 18, y + 31, PANEL_LABELS[panel], "label")
        parts[-1] = parts[-1].replace('fill="#243447"', 'fill="#ffffff"')
        rows = [
            ("Benchmark observations", f"{len(obs):,}"),
            ("Mapped observations used", f"{mapped:,}"),
            ("Linked candidate rows", f"{len(lnk):,}"),
            ("Mean measured activity", f"{obs['label_processed'].mean():.1f}"),
            ("Activity SD", f"{obs['label_processed'].std():.1f}"),
            ("Current score R²", f"{cal['r_squared']:.3f}"),
        ]
        yy = y + 82
        for label, value in rows:
            text(parts, x + 18, yy, label, "small")
            text(parts, x + w - 18, yy, value, "value", "end")
            yy += 31
    text(
        parts,
        40,
        408,
        "Ichihara panels are the primary fitting set; Martinelli is an external diagnostic pending target-identity resolution.",
        "small",
    )
    finish(parts, out / "panel_overview.svg")


def render_score_performance(analysis: dict, out: Path) -> None:
    """Render panel correlations for baseline and alternative score weights."""
    width, height = 1200, 500
    parts = svg_start(
        width,
        height,
        "Score performance by panel",
        "Spearman rank correlation for current and experimental score profiles by benchmark panel.",
    )
    text(parts, 55, 42, "Score performance by panel", "title")
    text(
        parts,
        55,
        67,
        "Observation-level Spearman rho; historical baseline versus pre-production candidates",
        "subtitle",
    )
    left, right, top, bottom = 95, 1145, 110, 405
    ymin, ymax = -0.5, 0.4
    for tick in [-0.4, -0.2, 0.0, 0.2, 0.4]:
        yy = scale(tick, ymin, ymax, bottom, top)
        line(parts, left, yy, right, yy, COLORS["grid"], 1, "4 4" if tick else "")
        text(parts, left - 12, yy + 4, f"{tick:.1f}", "axis", "end")
    line(parts, left, bottom, right, bottom, COLORS["muted"])
    line(parts, left, top, left, bottom, COLORS["muted"])
    colors = [COLORS["blue"], COLORS["teal"], COLORS["orange"]]
    group_width = (right - left) / len(PANEL_ORDER)
    bar_width = 48
    for idx, panel in enumerate(PANEL_ORDER):
        center = left + group_width * (idx + 0.5)
        text(parts, center, bottom + 34, PANEL_LABELS[panel], "label", "middle")
        text(parts, center, bottom + 52, "n=" + str(analysis["observation_counts"][panel]), "small", "middle")
        for j, profile in enumerate(PROFILE_ORDER):
            value = analysis["candidate_profiles"][profile]["panel_rho"][panel]
            x = center + (j - 1) * 62 - bar_width / 2
            zero = scale(0, ymin, ymax, bottom, top)
            yy = scale(value, ymin, ymax, bottom, top)
            rect(parts, x, min(zero, yy), bar_width, abs(zero - yy), colors[j], colors[j], 3)
            text(
                parts,
                x + bar_width / 2,
                min(zero, yy) - 7 if value >= 0 else max(zero, yy) + 17,
                f"{value:.3f}",
                "value",
                "middle",
            )
    legend_x = 95
    for j, profile in enumerate(PROFILE_ORDER):
        x = legend_x + j * 205
        rect(parts, x, 450, 16, 16, colors[j], colors[j], 2)
        text(parts, x + 24, 463, PROFILE_LABELS[profile], "small")
    finish(parts, out / "score_performance.svg")


def heat_color(value: float | None) -> str:
    """Map a correlation to a diverging color, using gray for missing values."""
    if value is None:
        return "#edf1f3"
    value = max(-0.35, min(0.35, value))
    if value >= 0:
        alpha = value / 0.35
        r = int(235 - 150 * alpha)
        g = int(242 - 55 * alpha)
        b = int(245 - 20 * alpha)
    else:
        alpha = abs(value) / 0.35
        r = int(245 - 10 * alpha)
        g = int(238 - 110 * alpha)
        b = int(238 - 95 * alpha)
    return f"#{r:02x}{g:02x}{b:02x}"


def render_component_signals(analysis: dict, out: Path) -> None:
    """Render a heatmap of component correlations with measured activity."""
    width, height = 1200, 445
    parts = svg_start(
        width,
        height,
        "Component signals by panel",
        "Component Spearman correlations with measured activity by benchmark panel.",
    )
    text(parts, 55, 42, "Component signals by panel", "title")
    text(
        parts,
        55,
        67,
        "Positive values indicate higher component values tend to rank more active observations",
        "subtitle",
    )
    x0, y0, cell_w, cell_h = 310, 125, 190, 72
    text(parts, x0 + 2 * cell_w, 102, "Spearman rho", "axis", "middle")
    for j, component in enumerate(COMPONENTS):
        text(parts, x0 + j * cell_w + cell_w / 2, y0 - 17, COMPONENT_LABELS[component], "label", "middle")
    for i, panel in enumerate(PANEL_ORDER):
        yy = y0 + i * cell_h
        text(parts, x0 - 20, yy + cell_h / 2 + 5, PANEL_LABELS[panel], "label", "end")
        for j, component in enumerate(COMPONENTS):
            raw = analysis["panel_correlations_median"][panel][component]["spearman_rho"]
            value = None if raw is None else float(raw)
            xx = x0 + j * cell_w
            rect(parts, xx + 4, yy + 4, cell_w - 8, cell_h - 8, heat_color(value), COLORS["white"], 5)
            label = "constant" if value is None else f"{value:+.3f}"
            text(parts, xx + cell_w / 2, yy + cell_h / 2 + 5, label, "value", "middle")
    # Scale legend.
    legend_y = 370
    for k in range(15):
        value = -0.35 + k * 0.05
        rect(parts, 310 + k * 26, legend_y, 27, 17, heat_color(value), COLORS["white"])
    text(parts, 300, legend_y + 13, "−0.35", "small", "end")
    text(parts, 710, legend_y + 13, "+0.35", "small")
    text(
        parts,
        55,
        425,
        "Off-target is constant in Ichihara 2 and Martinelli, so those panels cannot identify its efficacy weight.",
        "small",
    )
    finish(parts, out / "component_signals.svg")


def render_gc_sensitivity(analysis: dict, out: Path) -> None:
    """Render retention and score performance under alternative GC bounds."""
    width, height = 1200, 540
    parts = svg_start(
        width,
        height,
        "GC-bound sensitivity",
        "Observation retention and current-score correlation under alternative GC bounds.",
    )
    text(parts, 55, 42, "GC-bound sensitivity by panel", "title")
    text(
        parts,
        55,
        67,
        "Existing historical candidate rows were re-filtered; these are sensitivity analyses, not new design runs",
        "subtitle",
    )
    conditions = ["0–100", "25–70", "30–65", "35–60"]
    keys = ["gc_0_100", "gc_25_70", "gc_30_65", "gc_35_60"]
    colors = [COLORS["blue"], COLORS["teal"], COLORS["orange"]]
    charts = [(95, 395, "Observation retention (%)", 0, 105), (675, 395, "Current score rho", -0.5, 0.4)]
    top = 125
    for chart_idx, (left, right, title, ymin, ymax) in enumerate(charts):
        text(parts, (left + right) / 2, 102, title, "label", "middle")
        for tick in [0, 25, 50, 75, 100] if chart_idx == 0 else [-0.4, -0.2, 0, 0.2, 0.4]:
            yy = scale(tick, ymin, ymax, 395, top)
            line(parts, left, yy, right, yy, COLORS["grid"], 1, "4 4" if tick != 0 else "")
            text(parts, left - 10, yy + 4, f"{tick:.0f}" if chart_idx == 0 else f"{tick:.1f}", "axis", "end")
        line(parts, left, top, left, 395, COLORS["muted"])
        line(parts, left, 395, right, 395, COLORS["muted"])
        for i, cond in enumerate(conditions):
            xx = left + i * (right - left) / 3
            text(parts, xx, 421, cond, "axis", "middle")
        for i, panel in enumerate(PANEL_ORDER):
            points = []
            for j, key in enumerate(keys):
                entry = analysis["gc_sensitivity"][key]["panels"][panel]
                value = (
                    100 * entry["observation_retention_fraction"]
                    if chart_idx == 0
                    else entry["current_score_spearman"]["spearman_rho"]
                )
                xx = left + j * (right - left) / 3
                yy = scale(value, ymin, ymax, 395, top)
                points.append((xx, yy))
            poly = " ".join(f"{xx:.1f},{yy:.1f}" for xx, yy in points)
            parts.append(f'<polyline points="{poly}" fill="none" stroke="{colors[i]}" stroke-width="3" />')
            for xx, yy in points:
                parts.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="5" fill="{colors[i]}" />')
    for i, panel in enumerate(PANEL_ORDER):
        x = 95 + i * 205
        rect(parts, x, 472, 16, 16, colors[i], colors[i], 2)
        text(parts, x + 24, 485, PANEL_LABELS[panel], "small")
    finish(parts, out / "gc_sensitivity.svg")


def render_calibration(analysis: dict, out: Path) -> None:
    """Render panel calibration slopes and score-to-activity fit metrics."""
    width, height = 1200, 420
    parts = svg_start(
        width,
        height,
        "Current score calibration",
        "Descriptive calibration of current score against measured activity by panel.",
    )
    text(parts, 55, 42, "Current score calibration", "title")
    text(
        parts,
        55,
        67,
        "Historical baseline descriptive fit; R² is not an independently validated efficacy predictor",
        "subtitle",
    )
    left, right, top, bottom = 110, 1145, 125, 325
    for tick in [0, 0.02, 0.04, 0.06, 0.08]:
        yy = scale(tick, 0, 0.08, bottom, top)
        line(parts, left, yy, right, yy, COLORS["grid"], 1, "4 4")
        text(parts, left - 12, yy + 4, f"{tick:.2f}", "axis", "end")
    line(parts, left, top, left, bottom, COLORS["muted"])
    line(parts, left, bottom, right, bottom, COLORS["muted"])
    group_width = (right - left) / 3
    colors = [COLORS["blue"], COLORS["teal"], COLORS["orange"]]
    for i, panel in enumerate(PANEL_ORDER):
        cal = analysis["calibration_median_current_score"][panel]
        center = left + group_width * (i + 0.5)
        yy = scale(cal["r_squared"], 0, 0.08, bottom, top)
        rect(parts, center - 65, yy, 130, bottom - yy, colors[i], colors[i], 4)
        text(parts, center, yy - 9, f"R²={cal['r_squared']:.3f}", "value", "middle")
        text(parts, center, bottom + 31, PANEL_LABELS[panel], "label", "middle")
        text(parts, center, bottom + 50, f"slope={cal['slope']:+.2f}  RMSE={cal['rmse']:.1f}", "small", "middle")
    text(
        parts,
        55,
        385,
        "The score is suitable for relative ranking experiments, not for direct conversion to knockdown percentage.",
        "small",
    )
    finish(parts, out / "calibration.svg")


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", default="benchmark_results/oligogym_scoring_analysis.json")
    parser.add_argument("--observations", default="tests/data/benchmarks/oligogym_native_design/observations.csv")
    parser.add_argument("--links", default="tests/data/benchmarks/oligogym_native_design/candidate_benchmark_links.csv")
    parser.add_argument("--output-dir", default="docs/_static/oligogym_scoring")
    args = parser.parse_args()
    with Path(args.analysis).open(encoding="utf-8") as handle:
        analysis = json.load(handle)
    observations = pd.read_csv(args.observations)
    links = pd.read_csv(args.links)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    render_overview(analysis, observations, links, out)
    render_score_performance(analysis, out)
    render_component_signals(analysis, out)
    render_gc_sensitivity(analysis, out)
    render_calibration(analysis, out)
    print(f"Rendered five SVG figures to {out}")


if __name__ == "__main__":
    main()
