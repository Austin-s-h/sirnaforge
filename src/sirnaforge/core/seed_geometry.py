"""Guide/transcript pairing geometry: the one place the antiparallel map is written down (#101).

Issue #101 adds a transcript-seed liability channel and an intent-aware coverage numerator, and both
ask the same question: *where in a transcript does a guide pair?* Getting that backwards is the
defect #101 names explicitly, because it does not change any count. Searching a cDNA for the guide's
seed window **forward** (which is what ``off_target._prepare_seed_queries`` does, correctly, against
a miRNA database) finds guide-*identical* windows -- passenger-orientation sites. In a random
transcriptome those occur at the same expected rate as real target sites, so a count-based null is
indistinguishable and only the identity is wrong. Hence this module, and hence tests that assert the
matched **string and coordinates** rather than the number of hits.

Nothing here is a scanner. This is a leaf of pure arithmetic and string algebra: no I/O, no
alignment, no configuration. It deliberately imports
:func:`sirnaforge.core.repeat_detection.normalize_guide_sequence` instead of adding a third
normaliser beside that one and ``off_target._normalize_nucleotide_sequence``; the three must stay
byte-equivalent, and ``tests/unit/test_target_intent.py`` pins that they do without either of those
modules being edited.

DEFINITIONS
-----------
- Guide positions are 1-based from the guide's 5' end: ``g[1], g[2], ... g[21]``.
- Transcript positions are 1-based on the cDNA **as written in the reference FASTA**, i.e. the
  transcript's own sense strand. Spelled ``coordinate_system = "transcript_cdna_1based"`` wherever a
  site is published, so it is stated rather than implied.
- Both sequences are normalised to DNA (uppercase, ``U`` -> ``T``) before any comparison.

THE MAP
-------
Pairing is ANTIPARALLEL, so the transcript position paired with guide position ``i`` *decreases* as
``i`` increases::

    p(i) = anchor + 2 - i,   where anchor := p(2)

``anchor`` is the transcript position paired with guide position 2, and is the site's 3'-most
complementary base. It is the site identity (see the transcript-seed dedup rule) precisely because
all four seed classes are nested windows sharing that one base and nothing else.

THE EQUATION
------------
A site of the guide window ``[2..b]`` exists at ``anchor`` iff::

    transcript[p(i)] == complement(guide[i])   for every i in 2..b

equivalently, as one string comparison::

    transcript[anchor-(b-2) .. anchor] == reverse_complement(guide[2..b])

Class assignment at one ``anchor``::

    6mer     : b = 7
    7mer-m8  : b = 8
    7mer-A1  : b = 7 AND transcript[anchor + 1] == 'A'      (the A1 base sits at p(1) = anchor+1)
    8mer     : b = 8 AND transcript[anchor + 1] == 'A'

The A1 base is NOT required to be complementary to ``g[1]``: it is an unpaired target adenosine,
which is why it is tested as a literal ``'A'`` and never through :func:`complement`.

WORKED EXAMPLE (let-7a-like guide, so the answer is independently checkable)
---------------------------------------------------------------------------
::

    guide       = 5'-TGAGGTAGTAGGTTGTATAGT-3'   (normalised DNA)
    guide[2..8] = "GAGGTAG"                     (m8 window)
    reverse_complement("GAGGTAG") = "CTACCTC"
    guide[2..7] = "GAGGTA" ; reverse_complement = "TACCTC"

Take a transcript whose bases 101..108 are ``"CTACCTCA"``:

- ``"CTACCTC"`` occupies 101..107, so ``p(8) = 101`` and ``anchor = p(2) = 107``.
- ``p(1) = 108`` holds ``'A'`` -> 8mer, ``site_start = 101``, ``site_end = anchor = 107``.
- If position 108 held ``'G'`` instead -> 7mer-m8, same anchor 107.
- If only 102..107 matched ``"TACCTC"`` -> 6mer, anchor still 107; with 108 == ``'A'`` -> 7mer-A1.

All four cases share ``anchor = 107`` and therefore one dedup identity, with the maximal class
published.
"""

from __future__ import annotations

from sirnaforge.core.repeat_detection import normalize_guide_sequence

__all__ = [
    "A1_POSITION",
    "M8_POSITION",
    "SEED_START_POSITION",
    "complement",
    "normalize_guide_sequence",
    "paired_transcript_position",
    "reverse_complement",
    "seed_window",
]

# Guide position whose paired transcript position, p(1) = anchor + 1, carries the unpaired target
# adenosine of a 7mer-A1/8mer site. Named because "position 1" appears in three unrelated senses in
# this codebase (guide 5' base, 1-based coordinate origin, first array element) and only this one
# means the A1 base.
A1_POSITION = 1

# Guide position 8, the "m8" base. Its pairing is a real base pair, unlike the A1 adenosine, which
# is why 7mer-m8 outranks 7mer-A1 wherever the classes are ordered.
M8_POSITION = 8

# Guide position 2: the 5'-most base of the canonical seed and the anchor's partner.
SEED_START_POSITION = 2

# Same table as ``repeat_detection._reverse_complement`` uses, deliberately: a non-ACGT base (``N``,
# an IUPAC ambiguity code) passes through unchanged rather than raising, so an ambiguous reference
# base simply fails to match instead of aborting a scan. Byte-equivalence with that function is
# pinned by test, because a private divergence here would silently halve a site count.
_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")


def complement(sequence: str) -> str:
    """Base-wise complement, in the same 5'->3' reading order as the input.

    Not the reverse complement: this is the per-position partner used when the pairing map is
    applied one position at a time (``transcript[p(i)] == complement(guide[i])``). #101 keeps both
    spellings because the per-position form is what the tests assert the map with, and the string
    form is what a scanner searches with; deriving one from the other in the test would let a
    reversed map agree with itself.

    Args:
        sequence: Nucleotide sequence. Normalised to DNA first, so RNA input is accepted.

    Returns:
        The complement, un-reversed. Non-ACGT characters are returned unchanged.
    """
    return normalize_guide_sequence(sequence).translate(_COMPLEMENT_TABLE)


def reverse_complement(sequence: str) -> str:
    """Reverse complement, i.e. the sense-strand string a guide window pairs with.

    This is the search string for a transcript-seed site and for :func:`observe_coverage`'s
    full-site coverage test: a guide is antisense to its target, so what occurs verbatim in the
    cDNA is the guide's reverse complement and never the guide itself.

    Args:
        sequence: Nucleotide sequence. Normalised to DNA first, so RNA input is accepted.

    Returns:
        The reverse complement. Non-ACGT characters are complemented to themselves.
    """
    return normalize_guide_sequence(sequence).translate(_COMPLEMENT_TABLE)[::-1]


def seed_window(guide: str, first: int, last: int) -> str:
    """The guide substring ``guide[first..last]``, 1-based and inclusive on both ends.

    1-based inclusive because every coordinate #101 publishes is, and because the off-by-one
    between ``guide[2..8]`` and ``guide[1:8]`` is exactly the class-boundary error the 6mer/7mer/8mer
    fixtures exist to catch. A window running past the guide's end raises rather than silently
    truncating: ``off_target._prepare_seed_queries`` truncates with a warning, which is right for a
    best-effort miRNA query but would here publish a 6mer as though it were the requested 8mer.

    Args:
        guide: Guide sequence, RNA or DNA, any case.
        first: 1-based inclusive start position in the guide.
        last: 1-based inclusive end position in the guide.

    Returns:
        The normalised DNA window.

    Raises:
        ValueError: If the window is empty, starts before position 1, or ends past the guide.
    """
    normalized = normalize_guide_sequence(guide)
    if first < 1 or last < first:
        raise ValueError(f"seed window [{first}..{last}] is not a 1-based inclusive range")
    if last > len(normalized):
        raise ValueError(f"seed window [{first}..{last}] runs past the {len(normalized)}-base guide {normalized!r}")
    return normalized[first - 1 : last]


def paired_transcript_position(guide_position: int, anchor: int) -> int:
    """Transcript position paired with ``guide_position``, given the site's anchor.

    The whole antiparallel map, in one line: ``p(i) = anchor + 2 - i``. It decreases as the guide
    position increases, which is the property a reversed implementation (``anchor - 2 + i``) breaks
    while leaving ``p(2) == anchor`` intact -- so the tests assert the map for ``i`` in 2..8, not
    just at the anchor.

    Args:
        guide_position: 1-based guide position ``i``. ``A1_POSITION`` is legitimate and returns
            ``anchor + 1``, the unpaired target adenosine's position; it is a *coordinate*, not a
            claim that ``g[1]`` pairs there.
        anchor: Transcript position paired with guide position 2.

    Returns:
        The 1-based transcript (cDNA sense-strand) position paired with that guide position. May
        fall outside the transcript; a caller scanning near either end must bounds-check, because
        this is arithmetic and knows nothing about transcript length.

    Raises:
        ValueError: If ``guide_position`` is not a 1-based position.
    """
    if guide_position < 1:
        raise ValueError(f"guide position {guide_position} is not 1-based")
    return anchor + SEED_START_POSITION - guide_position
