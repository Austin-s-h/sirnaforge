"""Ranking utilities for ZFN off-target sites.

This module centralizes deterministic site ordering for reporting/regression use.
The tie-break behavior follows the PROGNOS-inspired policy used in tests:

1. Higher score first.
2. Region priority for equal scores: Exon > Promoter > Intron > Intergenic > Unknown.
3. Chromosomal location ordering.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

from sirnaforge.models.zfn import ZFNDesignParameters, ZFNOffTargetSite

_REGION_PRIORITY: dict[str, int] = {
    "unknown": 0,
    "intergenic": 1,
    "intron": 2,
    "promoter": 3,
    "exon": 4,
}


def _chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    """Return a sortable chromosome key with numeric chromosomes first."""
    normalized = chrom.lower()
    if normalized.startswith("chr"):
        normalized = normalized[3:]

    if normalized.isdigit():
        return (0, int(normalized))

    if normalized in {"x", "y", "m", "mt"}:
        special_order = {"x": 23, "y": 24, "m": 25, "mt": 25}
        return (0, special_order[normalized])

    match = re.match(r"^(\d+)(.*)$", normalized)
    if match:
        suffix = match.group(2)
        return (1, f"{int(match.group(1)):08d}{suffix}")

    return (2, normalized)


def rank_sites(
    sites: Sequence[ZFNOffTargetSite],
    params: ZFNDesignParameters | None = None,
) -> list[ZFNOffTargetSite]:
    """Rank off-target sites with score-first ordering and deterministic tie-breaks.

    Args:
        sites: Candidate sites to rank.
        params: Optional design parameters. Currently only used to keep a typed,
            future-proof API for algorithm-specific ranking extensions.

    Returns:
        A new sorted list of sites.
    """
    algorithm_label = params.algorithm.value if params is not None else ""

    return sorted(
        sites,
        key=lambda site: (
            -site.score,
            -_REGION_PRIORITY.get(site.region, 0),
            _chrom_sort_key(site.chrom),
            site.start_1based,
            site.end_1based,
            algorithm_label,
            site.site_id,
        ),
    )


__all__ = ["rank_sites"]
