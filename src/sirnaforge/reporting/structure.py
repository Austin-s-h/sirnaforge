"""Lay out and draw a guide's published secondary structure.

**Layout, not folding.** The dot-bracket comes from the run's own ``structure`` column and is never
recomputed: refolding here would let the picture disagree with the ``paired_fraction`` the gates were
decided on. ViennaRNA supplies the naview coordinates when it is installed, which is a layout of a
structure, not a prediction of one.

Two renderings, because a 23-mer guide has two honest pictures:

* ``layout_xy`` -- naview coordinates, drawn as the backbone with rungs between paired bases. This is
  the recognisable hairpin.
* ``pair_table`` -- partners only, drawn as an arc diagram over a straight sequence. The fallback when
  ViennaRNA is absent, and the clearer picture when a caller wants the seed register to stay readable.

A degenerate fold -- all dots, ``mfe`` 0 -- is a real and common answer, not a failure: 34.8% of the
passing pool on one MSH3 run. It renders as an open chain and says so, rather than as an empty panel.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

#: Guide positions carrying the seed, 1-based inclusive. The seed decides off-target reach, so it is
#: shaded on every rendering.
SEED_SPAN = (2, 8)

_BACKBONE = "#9ca3af"
_RUNG = "#1d4ed8"
_SEED = "#fef3c7"
_BASE_TEXT = "#1a1d21"


#: The SVG namespace, required for a standalone ``.svg`` file and pointless inside an HTML document,
#: where the parser already knows the element. The report omits it: its own contract is that the file
#: contains no URL at all, and a namespace declaration would satisfy the letter and not the spirit.
_XMLNS = 'xmlns="http://www.w3.org/2000/svg"'


class StructureError(ValueError):
    """A dot-bracket string that cannot describe a structure."""


def pair_table(structure: str) -> list[int]:
    """Partner index per position, 0-based, ``-1`` for unpaired.

    Args:
        structure: Dot-bracket string. ``.``, ``(`` and ``)`` only.

    Returns:
        A list as long as ``structure``.

    Raises:
        StructureError: Unbalanced brackets or an unexpected character.
    """
    partners = [-1] * len(structure)
    stack: list[int] = []
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise StructureError(f"unbalanced dot-bracket at position {i + 1}: {structure!r}")
            j = stack.pop()
            partners[i], partners[j] = j, i
        elif ch != ".":
            raise StructureError(f"unexpected character {ch!r} in dot-bracket {structure!r}")
    if stack:
        raise StructureError(f"{len(stack)} unclosed pair(s) in dot-bracket {structure!r}")
    return partners


def is_degenerate(structure: str) -> bool:
    """True when no pair was predicted -- the open chain, the physical floor of the MFE."""
    return bool(structure) and set(structure) <= {"."}


def layout_xy(structure: str) -> list[tuple[float, float]] | None:
    """Naview coordinates for each base, or None when there is no layout to be had.

    ``None`` covers three cases: ViennaRNA is not installed, the dot-bracket does not describe a
    structure, or the layout call failed. All three fall back to the arc diagram.

    **The validity check is not defensive tidiness -- it is load-bearing.**
    ``RNA.naview_xy_coordinates`` *segfaults* on an unbalanced dot-bracket (verified on ViennaRNA
    2.7.2 for both ``..((..`` and ``..)(..``), which no ``except`` can catch and which would take the
    whole report process down over one malformed cell in a ``structure`` column. ``pair_table`` rejects
    the same strings in Python first, so ViennaRNA only ever sees input it can lay out.

    ViennaRNA returns one sentinel entry past the end; it is dropped here so the result is exactly as
    long as the structure.
    """
    try:
        pair_table(structure)
    except StructureError:
        return None
    try:
        import RNA  # noqa: PLC0415  # optional at runtime; the arc fallback covers its absence
    except ImportError:
        return None
    try:
        coordinates = RNA.naview_xy_coordinates(structure)
    except Exception:  # a layout failure costs the picture, never the report
        return None
    return [(round(float(c.X), 2), round(float(c.Y), 2)) for c in coordinates][: len(structure)]


def layouts_for(structures: Iterable[object]) -> dict[str, list[list[float]]]:
    """Coordinates keyed by dot-bracket, computed once per distinct structure.

    Deduplicating is what makes this affordable to embed: one MSH3 run has 40,079 candidates and only
    1,333 distinct structures among them.
    """
    out: dict[str, list[list[float]]] = {}
    for structure in structures:
        text = str(structure or "")
        if not text or text in out:
            continue
        coordinates = layout_xy(text)
        if coordinates is not None:
            out[text] = [[x, y] for x, y in coordinates]
    return out


@dataclass(frozen=True)
class StructureSummary:
    """What the picture is claiming, in numbers, so a caption cannot drift from the drawing."""

    length: int
    pairs: int
    paired_fraction: float
    degenerate: bool

    @property
    def caption(self) -> str:
        """A one-line description of the structure as drawn."""
        if self.degenerate:
            return f"no pairs predicted over {self.length} nt - open chain, the physical floor of the MFE"
        return f"{self.pairs} pairs, {2 * self.pairs} of {self.length} nt paired ({self.paired_fraction:.3f})"


def summarise(structure: str) -> StructureSummary:
    """Count the pairs a dot-bracket describes."""
    pairs = structure.count("(")
    length = len(structure)
    return StructureSummary(
        length=length,
        pairs=pairs,
        paired_fraction=round(2 * pairs / length, 3) if length else 0.0,
        degenerate=is_degenerate(structure),
    )


def structure_svg(
    sequence: str,
    structure: str,
    *,
    mfe: float | None = None,
    seed: tuple[int, int] | None = SEED_SPAN,
    width: int = 460,
    height: int = 260,
    standalone: bool = True,
) -> str:
    """The guide's structure as inline SVG: naview layout when available, else an arc diagram.

    Args:
        sequence: Guide sequence, drawn base by base. Spelled as given.
        structure: The run's published dot-bracket.
        mfe: Published minimum free energy, for the caption.
        seed: 1-based inclusive span to shade, or None.
        width: SVG width in px.
        height: SVG height in px.
        standalone: Emit the SVG namespace. True for a ``.svg`` file, False when embedding in HTML.

    Returns:
        One ``<svg>`` element. No external references.

    Raises:
        StructureError: The dot-bracket does not describe a structure, or does not match the sequence.
    """
    if len(sequence) != len(structure):
        raise StructureError(f"sequence is {len(sequence)} nt and dot-bracket is {len(structure)}")
    partners = pair_table(structure)
    coordinates = layout_xy(structure)
    body = (
        _draw_layout(sequence, partners, coordinates, seed, width, height)
        if coordinates
        else _draw_arcs(sequence, partners, seed, width, height)
    )
    summary = summarise(structure)
    caption = summary.caption + (f" - mfe {mfe:.2f} kcal/mol" if mfe is not None else "")
    return (
        f'<svg viewBox="0 0 {width:g} {height:g}" width="100%" style="max-width:{width:g}px;'
        'font:11px ui-monospace,SFMono-Regular,Menlo,monospace" '
        f'{_XMLNS if standalone else ""} role="img">'
        f"{body}"
        f'<text x="{width / 2:g}" y="{height - 6:g}" text-anchor="middle" font-size="10.5" '
        f'fill="#6b7280" font-family="ui-sans-serif,system-ui,sans-serif">{_esc(caption)}</text>'
        "</svg>"
    )


def _draw_layout(
    sequence: str,
    partners: list[int],
    coordinates: list[tuple[float, float]],
    seed: tuple[int, int] | None,
    width: int,
    height: int,
) -> str:
    """Backbone plus rungs, in naview coordinates scaled to the viewport."""
    xs = [x for x, _ in coordinates]
    ys = [y for _, y in coordinates]
    pad, foot = 20.0, 26.0
    sx, sy = _fit(min(xs), max(xs), min(ys), max(ys), pad, width - pad, pad, height - foot)

    parts = [
        f'<path d="{"".join(("M" if i == 0 else "L") + f"{sx(x):.1f} {sy(y):.1f}" for i, (x, y) in enumerate(coordinates))}" fill="none" stroke="{_BACKBONE}" stroke-width="1.6"/>'
    ]
    for i, j in enumerate(partners):
        if j > i:
            xi, yi = coordinates[i]
            xj, yj = coordinates[j]
            parts.append(
                f'<line x1="{sx(xi):.1f}" y1="{sy(yi):.1f}" x2="{sx(xj):.1f}" y2="{sy(yj):.1f}" '
                f'stroke="{_RUNG}" stroke-width="1.2" stroke-opacity="0.55"/>'
            )
    for i, (x, y) in enumerate(coordinates):
        px, py = sx(x), sy(y)
        if seed and seed[0] <= i + 1 <= seed[1]:
            parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="7" fill="{_SEED}"/>')
        parts.append(
            f'<text x="{px:.1f}" y="{py + 3.5:.1f}" text-anchor="middle" font-size="9.5" '
            f'fill="{_BASE_TEXT}">{_esc(sequence[i])}</text>'
        )
    if seed:
        parts.append(
            f'<text x="{width - 20:g}" y="{20:g}" text-anchor="end" font-size="10" fill="#92400e" '
            f'font-family="ui-sans-serif,system-ui,sans-serif">seed {seed[0]}-{seed[1]} shaded</text>'
        )
    return "".join(parts)


def _draw_arcs(
    sequence: str,
    partners: list[int],
    seed: tuple[int, int] | None,
    width: int,
    height: int,
) -> str:
    """Sequence on a line with semicircles joining partners: the layout-free fallback."""
    n = len(sequence)
    pad, foot = 18.0, 30.0
    baseline = height - foot - 12.0
    step = (width - 2 * pad) / max(n - 1, 1)
    px = [pad + i * step for i in range(n)]
    parts: list[str] = []
    if seed:
        left, right = px[seed[0] - 1] - step / 2, px[min(seed[1], n) - 1] + step / 2
        parts.append(
            f'<rect x="{left:.1f}" y="{baseline - 12:.1f}" width="{right - left:.1f}" height="20" fill="{_SEED}"/>'
        )
    parts.append(
        f'<line x1="{pad:g}" y1="{baseline:g}" x2="{width - pad:g}" y2="{baseline:g}" stroke="{_BACKBONE}" stroke-width="1.2"/>'
    )
    for i, j in enumerate(partners):
        if j > i:
            radius = (px[j] - px[i]) / 2
            top = max(baseline - radius, 14.0)
            parts.append(
                f'<path d="M{px[i]:.1f} {baseline:.1f}Q{(px[i] + px[j]) / 2:.1f} {2 * top - baseline:.1f} '
                f'{px[j]:.1f} {baseline:.1f}" fill="none" stroke="{_RUNG}" stroke-width="1.2" stroke-opacity="0.55"/>'
            )
    for i, base in enumerate(sequence):
        parts.append(
            f'<text x="{px[i]:.1f}" y="{baseline + 14:.1f}" text-anchor="middle" font-size="9.5" fill="{_BASE_TEXT}">{_esc(base)}</text>'
        )
    return "".join(parts)


def _fit(
    x_lo: float, x_hi: float, y_lo: float, y_hi: float, left: float, right: float, top: float, bottom: float
) -> tuple[Callable[[float], float], Callable[[float], float]]:
    """Scale factories preserving aspect ratio, so a hairpin is not drawn stretched."""
    scale = min((right - left) / max(x_hi - x_lo, 1e-6), (bottom - top) / max(y_hi - y_lo, 1e-6))
    dx = left + ((right - left) - (x_hi - x_lo) * scale) / 2
    dy = top + ((bottom - top) - (y_hi - y_lo) * scale) / 2
    return (lambda x: dx + (x - x_lo) * scale), (lambda y: dy + (y_hi - y) * scale)


def _esc(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
