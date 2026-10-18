"""Contract tests for the reusable transcript-position tracks.

The map exists to answer where designs sit and where none could be enumerated, so the tests pin the
region call, the gap detection and the coordinate mapping a caller places its own markers with. The
SVG assertions are deliberately structural rather than pixel-level: what matters is that the shapes
exist, that a series is one path rather than thousands of elements, and that nothing reaches outward.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from sirnaforge.reporting.tracks import (
    PointSeries,
    TickSeries,
    TranscriptRegions,
    legend_html,
    transcript_map_svg,
    transcript_regions,
    uncovered_stretches,
)

_ORF_HEADER = "transcript_id\tsequence_length\tutr5_length\tlongest_orf_start\tlongest_orf_end\tutr3_length"


def _write_orf_report(tmp_path: Path, rows: list[str]) -> Path:
    run = tmp_path / "run" / "orf_reports"
    run.mkdir(parents=True)
    (run / "orf_validation.txt").write_text("\n".join([_ORF_HEADER, *rows]) + "\n")
    return tmp_path / "run"


@pytest.mark.unit
def test_regions_come_from_the_runs_own_orf_call(tmp_path: Path) -> None:
    """The run already published the ORF, so the map annotates it rather than re-deriving one."""
    run = _write_orf_report(tmp_path, ["ENST1\t4443\t76\t76\t3490\t953"])
    regions = transcript_regions(run)

    assert list(regions) == ["ENST1"]
    assert regions["ENST1"].spans == [("utr5", 1, 75), ("cds", 76, 3490), ("utr3", 3491, 4443)]
    assert regions["ENST1"].region_of(10) == "utr5"
    assert regions["ENST1"].region_of(2000) == "cds"
    assert regions["ENST1"].region_of(4000) == "utr3"


@pytest.mark.unit
def test_a_run_without_an_orf_report_yields_no_regions(tmp_path: Path) -> None:
    """A map without a region bar is a lesser map, not an error."""
    assert transcript_regions(tmp_path) == {}


@pytest.mark.unit
def test_a_transcript_with_no_called_orf_is_not_labelled_coding() -> None:
    """Inventing a CDS is the one thing the region bar must never do.

    Returning ``("cds", 1, length)`` drew every lncRNA, NMD isoform and failed ORF call as
    ``CDS 1-N``: fully coding, on no evidence.
    """
    regions = TranscriptRegions("ENST1", 1000)
    assert regions.spans == [("unknown", 1, 1000)]
    assert regions.region_of(500) == "unknown"

    svg, _ = transcript_map_svg(regions, [PointSeries("p", "pass", [(1, 1.0)])])
    assert "no ORF called" in svg
    assert "CDS" not in svg


@pytest.mark.unit
def test_unenumerated_stretches_are_found_and_short_ones_are_not_called_out() -> None:
    """An empty region of the map is ambiguous unless the run says no window existed there.

    The tail of a transcript is always uncovered -- a 23-mer cannot start in the last 22 nt -- so a
    minimum length keeps arithmetic from being reported as a design decision.
    """
    positions = [1, 24, 47, 400, 423]
    gaps = uncovered_stretches(positions, 500, window=23, min_nt=40)

    assert (70, 399) in gaps, "the 330 nt with no candidate is the point of the panel"
    assert all(end - start + 1 >= 40 for start, end in gaps)
    assert uncovered_stretches([1], 30, window=23, min_nt=40) == [], "an 8 nt tail is arithmetic"


@pytest.mark.unit
def test_the_map_draws_one_path_per_series_not_one_element_per_point() -> None:
    """40,000 <circle> nodes cost megabytes and seconds of layout; the same marks as a path cost neither."""
    regions = TranscriptRegions("ENST1", 1000, 100, 900)
    series = [
        PointSeries("rejected", "fail", [(i, 40.0) for i in range(1, 400)]),
        PointSeries("passes", "pass", [(500, 80.0), (600, 70.0)]),
    ]
    svg, _ = transcript_map_svg(regions, series)

    assert svg.count("<circle") == 0
    assert svg.count("<path") == 2, "one path per non-empty series"
    assert "CDS 100-900" in svg and "5&#39;UTR" in svg


@pytest.mark.unit
def test_the_geometry_lets_a_caller_place_a_marker_in_the_same_coordinates() -> None:
    """The report's selection marker is a second overlay, so it must share the base map's scale."""
    regions = TranscriptRegions("ENST1", 1000)
    svg, geometry = transcript_map_svg(regions, [PointSeries("p", "pass", [(1, 1.0), (1000, 2.0)])])

    assert geometry.length == 1000
    assert geometry.x0 < geometry.x1 and geometry.y0 < geometry.y1
    assert f'viewBox="0 0 {int(geometry.x1 + 18)}' in svg, "the viewBox and the geometry agree on width"
    assert set(geometry.as_dict()) == {"x0", "x1", "y0", "y1", "length", "value_min", "value_max"}


@pytest.mark.unit
def test_value_bounds_are_round_numbers() -> None:
    """Gridlines labelled with the data's extremes read as noise; a round step reads as a scale."""
    regions = TranscriptRegions("ENST1", 100)
    _, geometry = transcript_map_svg(regions, [PointSeries("p", "pass", [(1, 7.08), (50, 93.66)])])

    assert geometry.value_min == 0.0
    assert geometry.value_max == 100.0


@pytest.mark.unit
def test_an_embedded_map_declares_no_namespace_and_reaches_nothing() -> None:
    """The report's contract is that the file contains no URL at all, namespace declarations included."""
    regions = TranscriptRegions("ENST1", 1000, 100, 900)
    series = [PointSeries("passes", "pass", [(10, 1.0)])]
    embedded, _ = transcript_map_svg(regions, series, standalone=False, ticks=[TickSeries("panel", "reference", [5])])
    standalone, _ = transcript_map_svg(regions, series, standalone=True)

    assert not re.findall(r"https?://", embedded)
    assert "xmlns" not in embedded
    assert 'xmlns="http://www.w3.org/2000/svg"' in standalone, "a .svg file on disk needs it"
    assert "<script" not in embedded and "fetch(" not in embedded


@pytest.mark.unit
def test_the_legend_reports_the_count_it_drew() -> None:
    """A legend whose counts come from anywhere but the drawn points can disagree with the picture."""
    series = [
        PointSeries("passes every gate", "pass", [(1, 1.0), (2, 2.0), (3, 3.0)]),
        PointSeries("not established", "unknown", []),
    ]
    html = legend_html(series, gaps=True)

    assert "passes every gate (3)" in html
    assert "no candidate in this table" in html
    assert "not established" not in html, "a class the map drew nothing for is not a legend key"
