"""Contract tests for the guide-structure rendering.

The load-bearing property is that this module lays out a structure and never predicts one: the picture
has to be of the same fold the gates were decided on, or ``paired_fraction`` in the caption and the
shape above it can disagree. The rest pins the degenerate fold as a real answer rather than a failure,
because it is 34.8% of the passing pool on one MSH3 run.
"""

from __future__ import annotations

import builtins
import re

import pytest

from sirnaforge.reporting.structure import (
    StructureError,
    is_degenerate,
    layouts_for,
    pair_table,
    structure_svg,
    summarise,
)

GUIDE = "UGACGUUUCGAAUCUUCAGAAGA"
FOLDED = ".........((((....)))).."
FLAT = "." * 23


@pytest.mark.unit
def test_pair_table_names_both_partners() -> None:
    """A rung needs both ends, so the table is symmetric by construction."""
    partners = pair_table(FOLDED)

    assert len(partners) == len(FOLDED)
    assert partners[9] == 20 and partners[20] == 9
    assert partners[0] == -1
    for i, j in enumerate(partners):
        if j >= 0:
            assert partners[j] == i


@pytest.mark.unit
@pytest.mark.parametrize("bad", ["((.", ".)", "..x.."])
def test_a_dot_bracket_that_cannot_describe_a_structure_is_refused(bad: str) -> None:
    """Silently dropping an unmatched bracket would draw a structure the run never published."""
    with pytest.raises(StructureError):
        pair_table(bad)


@pytest.mark.unit
def test_the_open_chain_is_a_real_answer_not_a_failure() -> None:
    """All dots with mfe 0 is the physical floor of the MFE, and a third of one run's passing pool."""
    assert is_degenerate(FLAT) is True
    assert is_degenerate(FOLDED) is False

    summary = summarise(FLAT)
    assert summary.pairs == 0 and summary.paired_fraction == 0.0
    assert "no pairs predicted" in summary.caption

    svg = structure_svg(GUIDE, FLAT, mfe=0.0)
    assert svg.count("<line") == 0, "no pairs means no rungs"
    assert svg.count("<text") == len(GUIDE) + 2, "every base, plus the seed label and the caption"
    assert "no pairs predicted" in svg


@pytest.mark.unit
def test_the_caption_counts_what_the_picture_draws() -> None:
    """A caption sourced from anywhere but the drawn dot-bracket can contradict the drawing."""
    summary = summarise(FOLDED)
    assert summary.pairs == 4
    assert summary.paired_fraction == pytest.approx(8 / 23, abs=1e-3)

    svg = structure_svg(GUIDE, FOLDED, mfe=-1.9)
    assert svg.count("<line") == 4, "one rung per pair"
    assert "4 pairs, 8 of 23 nt paired" in svg
    assert "mfe -1.90 kcal/mol" in svg


@pytest.mark.unit
def test_a_structure_that_does_not_match_its_sequence_is_refused() -> None:
    """Drawing bases against the wrong dot-bracket would mislabel which positions are paired."""
    with pytest.raises(StructureError, match="23 nt and dot-bracket is 5"):
        structure_svg(GUIDE, ".....")


@pytest.mark.unit
def test_layouts_are_computed_once_per_distinct_structure() -> None:
    """Deduplicating is what makes the coordinates small enough to embed: 40,079 rows, 1,333 shapes."""
    layouts = layouts_for([FOLDED, FOLDED, FLAT, None, "", FOLDED])

    assert set(layouts) == {FOLDED, FLAT}
    assert len(layouts[FOLDED]) == len(FOLDED), "one coordinate pair per base, sentinel dropped"
    assert all(len(xy) == 2 for xy in layouts[FOLDED])


@pytest.mark.unit
def test_the_arc_diagram_carries_the_structure_when_no_layout_is_available(monkeypatch) -> None:
    """ViennaRNA is optional, so its absence must cost the layout and not the panel."""
    real_import = builtins.__import__

    def no_rna(name: str, *args: object, **kwargs: object) -> object:
        if name == "RNA":
            raise ImportError("ViennaRNA is not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_rna)

    svg = structure_svg(GUIDE, FOLDED, mfe=-1.9)
    assert svg.count("<path") == 4, "each pair becomes one arc"
    assert svg.count("<text") == len(GUIDE) + 1
    assert layouts_for([FOLDED]) == {}, "and no coordinates are published"


@pytest.mark.unit
def test_an_embedded_structure_declares_no_namespace_and_reaches_nothing() -> None:
    """Same contract as the rest of the report: one file, no URL, no script."""
    embedded = structure_svg(GUIDE, FOLDED, standalone=False)
    standalone = structure_svg(GUIDE, FOLDED, standalone=True)

    assert not re.findall(r"https?://", embedded)
    assert 'xmlns="http://www.w3.org/2000/svg"' in standalone
    assert "<script" not in embedded
