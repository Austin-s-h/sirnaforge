"""What the two SVG builders share: the namespace declaration and text escaping.

structure.py and tracks.py both emit standalone-or-embedded SVG, and neither imports the other.
"""

from __future__ import annotations

__all__ = ["XMLNS", "escape_text"]

#: The SVG namespace, required for a standalone ``.svg`` file and pointless inside an HTML document,
#: where the parser already knows the element. The report omits it: its own contract is that the file
#: contains no URL at all, and a namespace declaration would satisfy the letter and not the spirit.
XMLNS = 'xmlns="http://www.w3.org/2000/svg"'


def escape_text(text: str) -> str:
    """Escape a string for an SVG text node or attribute value."""
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )
