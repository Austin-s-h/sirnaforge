"""This repository is public, so an internal programme target may not name itself in it.

Twenty-four comments, docstrings and fixture values justified a threshold or a design decision by
citing measurements from an internal drug-programme run and naming that programme's target. The
measurements are the evidence and stay; the identifier does not, because a public reader learns the
programme from it. The established phrasing is "one internal run" / "one internal reference set".

The token is assembled at runtime rather than written out, so a plain ``grep`` over the tree stays
clean and this guard is not its own only hit.
"""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Split so the file does not itself contain the identifier it forbids.
PRIVATE_TARGET = "MSH" + "3"

#: Source trees that ship publicly. Comments and test fixtures are the surface that leaked.
SCANNED_TREES = ("src", "tests")


@pytest.mark.unit
def test_no_private_programme_target_in_source_or_tests() -> None:
    """The programme target appears nowhere under ``src/`` or ``tests/``, in any case.

    Case-insensitive because the species convention writes the mouse symbol in title case, and both
    spellings were present: the prose carried the upper-case human form and a conservation fixture's
    symbol column carried the title-case rodent one. A public reference gene used deliberately as a
    public baseline (``TP53``) is a different thing and is untouched.
    """
    needle = PRIVATE_TARGET.lower()
    offenders = [
        f"{path.relative_to(REPO_ROOT)}:{n}"
        for tree in SCANNED_TREES
        for path in sorted((REPO_ROOT / tree).rglob("*.py"))
        for n, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1)
        if needle in line.lower()
    ]

    assert not offenders, (
        f"{len(offenders)} site(s) name the internal programme target. Keep the measurement, drop the "
        f"identifier -- say 'one internal run' instead: {offenders}"
    )
