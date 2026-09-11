"""Transcript-position tracks as inline SVG, independent of any one report.

A design's position on its transcript is the axis every other panel is missing: whether a guide sits
in the 5'UTR, the CDS or the 3'UTR, where the enumerated windows cluster, and which stretches the
design gates emptied. This module draws that, from plain numbers, and returns an SVG string.

Reusable on purpose. It takes coordinates and series rather than a run directory or a DataFrame, so
the HTML report, a notebook and a standalone figure all draw the same map:

    from sirnaforge.reporting.tracks import PointSeries, transcript_map_svg, transcript_regions

    regions = transcript_regions(run_dir)["ENST00000265081"]
    svg, geometry = transcript_map_svg(
        regions,
        [PointSeries("passes every gate", "pass", passing_points)],
        gaps=not_enumerated_stretches(positions, regions.length, window=23),
    )

Constraints inherited from #103: no external URLs, no fonts, no script, no raster. Everything here is
shapes and text in one string that works over ``file://`` and inside a sandboxed iframe.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

#: Verdict colours, matching the report's own status pills so a point and a pill cannot disagree.
SERIES_FILL = {
    "pass": "#15803d",
    "warn": "#b45309",
    "fail": "#cbd5e1",
    "unknown": "#94a3b8",
    "reference": "#7c3aed",
    "selected": "#1d4ed8",
}

_REGION_FILL = {"utr5": "#cbd5e1", "cds": "#94a3b8", "utr3": "#cbd5e1"}
_AXIS = "#6b7280"
_GRID = "#e5e7eb"
_GAP_FILL = "#f1f5f9"


#: The SVG namespace, required for a standalone ``.svg`` file and pointless inside an HTML document,
#: where the parser already knows the element. The report omits it: its own contract is that the file
#: contains no URL at all, and a namespace declaration would satisfy the letter and not the spirit.
_XMLNS = 'xmlns="http://www.w3.org/2000/svg"'


@dataclass(frozen=True)
class TranscriptRegions:
    """One transcript's length and the CDS it carries, in 1-based inclusive coordinates.

    ``cds_start``/``cds_end`` are None for a transcript with no called ORF; the map then draws one
    undifferentiated bar rather than inventing a boundary.
    """

    transcript_id: str
    length: int
    cds_start: int | None = None
    cds_end: int | None = None

    @property
    def spans(self) -> list[tuple[str, int, int]]:
        """``(region, start, end)`` in transcript order, 1-based inclusive."""
        if self.cds_start is None or self.cds_end is None:
            return [("cds", 1, self.length)]
        out = []
        if self.cds_start > 1:
            out.append(("utr5", 1, self.cds_start - 1))
        out.append(("cds", self.cds_start, min(self.cds_end, self.length)))
        if self.cds_end < self.length:
            out.append(("utr3", self.cds_end + 1, self.length))
        return out

    def region_of(self, position: int) -> str:
        """Which region a 1-based position falls in."""
        for region, start, end in self.spans:
            if start <= position <= end:
                return region
        return "cds"


@dataclass(frozen=True)
class PointSeries:
    """One class of points: a label for the legend, a colour key, and ``(position, value)`` pairs."""

    label: str
    fill_key: str
    points: Sequence[tuple[int, float]]


@dataclass(frozen=True)
class TickSeries:
    """A one-dimensional lane of positions under the plot -- a reference panel, a shortlist."""

    label: str
    fill_key: str
    positions: Sequence[int]
    note: str | None = None


@dataclass(frozen=True)
class PlotGeometry:
    """Where the plot area landed, so a client can place a marker without re-deriving the scale."""

    x0: float
    x1: float
    y0: float
    y1: float
    length: int
    value_min: float
    value_max: float

    def as_dict(self) -> dict[str, float | int]:
        """A JSON-ready mapping, for embedding beside the SVG."""
        return {
            "x0": self.x0,
            "x1": self.x1,
            "y0": self.y0,
            "y1": self.y1,
            "length": self.length,
            "value_min": self.value_min,
            "value_max": self.value_max,
        }


def transcript_regions(run_dir: Path | str) -> dict[str, TranscriptRegions]:
    """Read CDS boundaries from a run's ORF validation report.

    The report already publishes ``longest_orf_start``/``longest_orf_end`` per transcript, so the map
    annotates regions from the run's own call rather than re-deriving an ORF. Returns an empty mapping
    when the run published no report, which is a map without a region bar, not an error.
    """
    report = next(iter(sorted(Path(run_dir).glob("**/orf_validation.txt"))), None)
    if report is None:
        return {}
    out: dict[str, TranscriptRegions] = {}
    with report.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tid = (row.get("transcript_id") or "").strip()
            if not tid:
                continue
            out[tid] = TranscriptRegions(
                transcript_id=tid,
                length=_int(row.get("sequence_length")) or 0,
                cds_start=_int(row.get("longest_orf_start")),
                cds_end=_int(row.get("longest_orf_end")),
            )
    return out


def not_enumerated_stretches(
    positions: Iterable[int], length: int, *, window: int, min_nt: int = 40
) -> list[tuple[int, int]]:
    """Stretches of the transcript no enumerated window covers, at least ``min_nt`` long.

    Drawn because an empty region of the map is otherwise ambiguous: a gate rejecting every window and
    no window existing look identical, and only one of them is a design decision.
    """
    covered = bytearray(length + 2)
    for start in positions:
        for i in range(max(1, start), min(length, start + window - 1) + 1):
            covered[i] = 1
    gaps: list[tuple[int, int]] = []
    run_start = None
    for i in range(1, length + 1):
        if not covered[i]:
            run_start = i if run_start is None else run_start
        elif run_start is not None:
            if i - run_start >= min_nt:
                gaps.append((run_start, i - 1))
            run_start = None
    if run_start is not None and length + 1 - run_start >= min_nt:
        gaps.append((run_start, length))
    return gaps


def transcript_map_svg(
    regions: TranscriptRegions,
    series: Sequence[PointSeries],
    *,
    title: str | None = None,
    value_label: str = "composite score",
    value_range: tuple[float, float] | None = None,
    gaps: Sequence[tuple[int, int]] = (),
    ticks: Sequence[TickSeries] = (),
    width: int = 1120,
    plot_height: int = 260,
    standalone: bool = True,
) -> tuple[str, PlotGeometry]:
    """One transcript, its regions, and every design plotted along it.

    Args:
        regions: The transcript's length and CDS boundaries.
        series: Point classes to plot, drawn in order so the first sits underneath.
        title: Heading above the region bar.
        value_label: Y-axis label.
        value_range: Y bounds; taken from the data when omitted.
        gaps: Stretches to shade as not enumerated, from :func:`not_enumerated_stretches`.
        ticks: One-dimensional lanes drawn beneath the plot.
        width: Total SVG width in px.
        plot_height: Height of the scatter area alone.
        standalone: Emit the SVG namespace. True for a ``.svg`` file, False when embedding in HTML.

    Returns:
        ``(svg, geometry)``. The geometry lets a caller place its own marker in the same coordinates.
    """
    pad_left, pad_right, pad_top = 64.0, 18.0, 14.0
    bar_h, bar_gap, axis_h = 22.0, 10.0, 46.0
    lane_h = 34.0

    x0, x1 = pad_left, float(width) - pad_right
    y_title = pad_top + (14.0 if title else 0.0)
    y_bar = y_title + (10.0 if title else 0.0)
    y_plot0 = y_bar + bar_h + bar_gap + 16.0
    y_plot1 = y_plot0 + plot_height
    total_h = y_plot1 + axis_h + lane_h * len(ticks) + 26.0

    values = [v for s in series for _, v in s.points]
    if value_range:
        lo, hi = value_range
    elif values:
        lo, hi = _round_bounds(min(values), max(values))
    else:
        lo, hi = 0.0, 1.0
    if hi <= lo:
        hi = lo + 1.0
    span = max(regions.length, 1)

    def sx(position: float) -> float:
        return x0 + (x1 - x0) * (min(max(position, 1.0), span) - 1.0) / max(span - 1, 1)

    def sy(value: float) -> float:
        return y_plot1 - (y_plot1 - y_plot0) * (min(max(value, lo), hi) - lo) / (hi - lo)

    parts: list[str] = [
        f'<svg viewBox="0 0 {width:g} {total_h:g}" width="100%" '
        f'style="max-width:{width:g}px;font:12px ui-sans-serif,system-ui,sans-serif" '
        f'{_XMLNS if standalone else ""} role="img">'
    ]
    if title:
        parts.append(
            f'<text x="{x0:g}" y="{y_title:g}" font-size="12.5" font-weight="650" fill="#1a1d21">{_esc(title)}</text>'
        )

    parts.extend(_region_bar(regions, sx, y_bar, bar_h))

    # Not-enumerated shading sits behind the points, so a point always wins the pixel.
    for start, end in gaps:
        gx, gw = sx(start), max(sx(end) - sx(start), 1.0)
        parts.append(f'<rect x="{gx:g}" y="{y_plot0:g}" width="{gw:g}" height="{plot_height:g}" fill="{_GAP_FILL}"/>')

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        gy = y_plot1 - (y_plot1 - y_plot0) * frac
        value = lo + (hi - lo) * frac
        parts.append(f'<line x1="{x0:g}" y1="{gy:g}" x2="{x1:g}" y2="{gy:g}" stroke="{_GRID}" stroke-width="1"/>')
        text = f"{value:,.0f}" if abs(hi - lo) >= 4 else f"{value:,.2f}"
        parts.append(
            f'<text x="{x0 - 8:g}" y="{gy + 4:g}" text-anchor="end" font-size="10" fill="{_AXIS}">{text}</text>'
        )
    parts.append(
        f'<text x="{14:g}" y="{(y_plot0 + y_plot1) / 2:g}" font-size="10.5" fill="{_AXIS}" '
        f'transform="rotate(-90 14 {(y_plot0 + y_plot1) / 2:g})" text-anchor="middle">{_esc(value_label)}</text>'
    )

    # One path per series rather than one element per point: 40,000 <circle> nodes cost megabytes and
    # seconds of layout, while the same marks as path subpaths cost neither.
    for s in series:
        if not s.points:
            continue
        marks = "".join(f"M{sx(p):.1f} {sy(v):.1f}h2v2h-2z" for p, v in s.points)
        parts.append(f'<path d="{marks}" fill="{SERIES_FILL.get(s.fill_key, _AXIS)}" fill-opacity="0.85"/>')

    parts.extend(_axis(sx, span, x0, x1, y_plot1))

    parts.extend(_tick_lanes(ticks, sx, x0, x1, y_plot1 + axis_h + 6.0, lane_h))

    parts.append("</svg>")
    geometry = PlotGeometry(x0=x0, x1=x1, y0=y_plot0, y1=y_plot1, length=span, value_min=lo, value_max=hi)
    return "".join(parts), geometry


def legend_html(series: Sequence[PointSeries], *, gaps: bool = False) -> str:
    """A legend as HTML rather than SVG text, so it wraps with the card instead of clipping."""
    items = [
        f'<span style="margin-right:14px;white-space:nowrap"><span style="display:inline-block;width:9px;height:9px;'
        f'border-radius:2px;background:{SERIES_FILL.get(s.fill_key, _AXIS)};margin-right:5px"></span>'
        f"{_esc(s.label)} ({len(s.points):,})</span>"
        for s in series
    ]
    if gaps:
        items.append(
            f'<span style="white-space:nowrap"><span style="display:inline-block;width:9px;height:9px;border-radius:2px;'
            f'background:{_GAP_FILL};border:1px solid {_GRID};margin-right:5px"></span>no enumerated window</span>'
        )
    return f'<div style="font-size:11.5px;color:#6b7280;margin:6px 0 2px">{"".join(items)}</div>'


def _region_bar(regions: TranscriptRegions, sx: Callable[[float], float], y: float, height: float) -> list[str]:
    """The 5'UTR / CDS / 3'UTR bar, labelled where a label fits."""
    out: list[str] = []
    for region, start, end in regions.spans:
        rx, rw = sx(start), max(sx(end) - sx(start), 1.0)
        out.append(f'<rect x="{rx:g}" y="{y:g}" width="{rw:g}" height="{height:g}" fill="{_REGION_FILL[region]}"/>')
        label = {"utr5": "5'UTR", "cds": f"CDS {start:,}-{end:,}", "utr3": "3'UTR"}[region]
        if rw > 8 * len(label) * 0.55:
            out.append(
                f'<text x="{rx + rw / 2:g}" y="{y + height / 2 + 4:g}" text-anchor="middle" '
                f'font-size="11" fill="#1f2937">{_esc(label)}</text>'
            )
    return out


def _tick_lanes(
    ticks: Sequence[TickSeries],
    sx: Callable[[float], float],
    x0: float,
    x1: float,
    y: float,
    lane_h: float,
) -> list[str]:
    """One labelled lane per tick series, under the plot."""
    out: list[str] = []
    for lane in ticks:
        out.append(
            f'<text x="{x0 - 8:g}" y="{y + 4:g}" text-anchor="end" font-size="10.5" fill="{_AXIS}">{_esc(lane.label)}</text>'
        )
        out.append(f'<line x1="{x0:g}" y1="{y:g}" x2="{x1:g}" y2="{y:g}" stroke="{_GRID}" stroke-width="1"/>')
        if lane.positions:
            marks = "".join(f"M{sx(p):.1f} {y - 7:.1f}h1.6v14h-1.6z" for p in lane.positions)
            out.append(f'<path d="{marks}" fill="{SERIES_FILL.get(lane.fill_key, _AXIS)}"/>')
        if lane.note:
            out.append(
                f'<text x="{x1:g}" y="{y + 18:g}" text-anchor="end" font-size="10" fill="{_AXIS}">{_esc(lane.note)}</text>'
            )
        y += lane_h
    return out


def _axis(sx: Callable[[float], float], span: int, x0: float, x1: float, y: float) -> list[str]:
    """A position axis with round-number ticks."""
    step = _tick_step(span)
    out = [f'<line x1="{x0:g}" y1="{y:g}" x2="{x1:g}" y2="{y:g}" stroke="{_AXIS}" stroke-width="1"/>']
    position = step
    while position < span:
        tx = sx(position)
        out.append(f'<line x1="{tx:g}" y1="{y:g}" x2="{tx:g}" y2="{y + 4:g}" stroke="{_AXIS}" stroke-width="1"/>')
        out.append(
            f'<text x="{tx:g}" y="{y + 16:g}" text-anchor="middle" font-size="10" fill="{_AXIS}">{position:,}</text>'
        )
        position += step
    out.append(
        f'<text x="{(x0 + x1) / 2:g}" y="{y + 31:g}" text-anchor="middle" font-size="10.5" fill="{_AXIS}">position on transcript (nt)</text>'
    )
    return out


def _round_bounds(lo: float, hi: float) -> tuple[float, float]:
    """Bounds on a round step containing the data, so gridline labels read as numbers not extremes."""
    span = max(hi - lo, 1e-9)
    step = 10.0 ** math.floor(math.log10(span))
    for factor in (1.0, 2.0, 2.5, 5.0, 10.0):
        if span / (step * factor) <= 4:
            step *= factor
            break
    return math.floor(lo / step) * step, math.ceil(hi / step) * step


def _tick_step(span: int) -> int:
    """A round tick interval giving roughly eight labels."""
    for step in (50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000):
        if span / step <= 10:
            return step
    return 50000


def _int(value: object) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _esc(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )
