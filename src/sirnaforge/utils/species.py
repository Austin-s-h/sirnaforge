"""Shared helpers for reasoning about species labels.

These utilities centralize the normalization and categorization rules we use
when aggregating off-target hits across genome/transcriptome and miRNA
pipelines. The helpers intentionally accept very loose inputs (common names,
scientific labels, assembly identifiers) so downstream callers can rely on a
single implementation when deciding whether a hit should be treated as
human-specific.
"""

from __future__ import annotations

import collections.abc as cabc
import re

# Common aliases that should be treated as human, regardless of formatting.
HUMAN_SPECIES_ALIASES: set[str] = {
    "human",
    "homo_sapiens",
    "hsapiens",
    "h_sapiens",
    "hsa",
    "grch37",
    "grch38",
    "hg19",
    "hg38",
}


def normalize_species_label(label: str | None) -> str:
    """Normalize a species label to a lowercase slug.

    The transformation removes punctuation, collapses whitespace, and replaces
    runs of non-alphanumeric characters with a single underscore so strings
    like ``"Homo sapiens (GRCh38)"`` become ``"homo_sapiens_grch38"``.
    """
    if label is None:
        return ""
    normalized = str(label).strip().lower()
    normalized = re.sub(r"[^a-z0-9]+", "_", normalized)
    return normalized.strip("_")


def is_human_species(label: str | None) -> bool:
    """Return ``True`` when the provided label should be bucketed as human."""
    normalized = normalize_species_label(label)
    if not normalized:
        return False
    if normalized in HUMAN_SPECIES_ALIASES:
        return True
    if normalized.startswith("human"):
        return True
    if normalized.startswith("homo_sapiens"):
        return True
    return normalized.startswith("grch") or normalized.startswith("hg")


def bucket_species(label: str | None) -> str:
    """Map an arbitrary label to ``"human"`` or ``"other"``."""
    return "human" if is_human_species(label) else "other"


def human_vs_other_totals(counts: cabc.Mapping[str, int]) -> tuple[int, int]:
    """Collapse per-species hit counts into ``(human, other)`` totals."""
    human_total = 0
    other_total = 0
    for species_label, value in counts.items():
        if is_human_species(species_label):
            human_total += value
        else:
            other_total += value
    return human_total, other_total


__all__ = [
    "HUMAN_SPECIES_ALIASES",
    "bucket_species",
    "human_vs_other_totals",
    "is_human_species",
    "normalize_species_label",
]
