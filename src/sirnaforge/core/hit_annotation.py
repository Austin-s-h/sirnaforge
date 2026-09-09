"""Persist per-hit classification onto the aggregated off-target table.

``hit_classification`` decides a hit's class once. This module writes that verdict onto the hit
row and derives the per-candidate counters *from the written row*, so the row-level TSV and the
candidate-level columns are one computation rather than two that can drift apart. It adds no
classification logic of its own.

The three columns are ``hit_class``, ``matched_symbol`` and ``symbol_lookup_missing``. A symbol
the index could not resolve is the literal string ``unknown``: an empty cell renders as "no gene"
and hides how much of a run had no annotation behind it.
"""

from __future__ import annotations

import csv
import os
from collections.abc import Iterable, Mapping, MutableMapping, Sequence
from pathlib import Path
from typing import Any

from sirnaforge.core.hit_classification import HitClass, HitClassCounts, HitClassification

HIT_CLASS_COLUMN = "hit_class"
MATCHED_SYMBOL_COLUMN = "matched_symbol"
SYMBOL_LOOKUP_MISSING_COLUMN = "symbol_lookup_missing"
CLASSIFICATION_COLUMNS: tuple[str, ...] = (
    HIT_CLASS_COLUMN,
    MATCHED_SYMBOL_COLUMN,
    SYMBOL_LOOKUP_MISSING_COLUMN,
)

UNKNOWN_SYMBOL = "unknown"


def annotate_hit_row(row: MutableMapping[str, Any], classification: HitClassification) -> HitClass:
    """Write a classification onto its hit row and return the class as persisted.

    The return value is read back out of the row rather than taken from ``classification``, so a
    caller cannot count one class and publish another.
    """
    row[HIT_CLASS_COLUMN] = classification.hit_class.value
    row[MATCHED_SYMBOL_COLUMN] = classification.matched_symbol or UNKNOWN_SYMBOL
    row[SYMBOL_LOOKUP_MISSING_COLUMN] = classification.symbol_lookup_missing
    return hit_class_of(row)


def hit_class_of(row: Mapping[str, Any]) -> HitClass:
    """Read the persisted class off a hit row.

    Raises:
        KeyError: The row was never annotated.
        ValueError: The row carries a value outside the four-way taxonomy.
    """
    return HitClass(str(row[HIT_CLASS_COLUMN]))


def is_annotated(row: Mapping[str, Any]) -> bool:
    """True when this row already carries a persisted class."""
    return HIT_CLASS_COLUMN in row


def accumulate_hit_class(
    row: Mapping[str, Any],
    counts: HitClassCounts,
    species_bucket: MutableMapping[str, int],
    hit_species: str,
) -> HitClass:
    """Add one persisted row to the candidate-level and per-species counters.

    Every quantity here is derived from the row's own columns, which is what keeps the per-hit
    table and the per-candidate totals from disagreeing. The counter field names are the
    ``HitClass`` values verbatim, so there is no class-to-counter mapping to get wrong.
    """
    hit_class = hit_class_of(row)
    field = hit_class.value
    setattr(counts, field, getattr(counts, field) + 1)
    species_bucket[field] = species_bucket.get(field, 0) + 1

    if _as_bool(row.get(SYMBOL_LOOKUP_MISSING_COLUMN)):
        counts.symbol_lookup_missing += 1
        species_bucket["symbol_lookup_missing"] = species_bucket.get("symbol_lookup_missing", 0) + 1

    if hit_class is HitClass.ORTHOLOG and str(row.get(MATCHED_SYMBOL_COLUMN)) != UNKNOWN_SYMBOL:
        counts.ortholog_species = frozenset(counts.ortholog_species | {hit_species})

    return hit_class


def count_persisted_classes(rows: Iterable[Mapping[str, Any]]) -> dict[str, int]:
    """Tally the persisted class column across rows, for logging and reconciliation."""
    tally = dict.fromkeys((member.value for member in HitClass), 0)
    for row in rows:
        if is_annotated(row):
            tally[hit_class_of(row).value] += 1
    return tally


def write_classified_hits(
    tsv_path: Path,
    rows: Sequence[Mapping[str, Any]],
    fieldnames: Sequence[str],
) -> int:
    """Rewrite an aggregated off-target TSV with the classification columns appended.

    Values are written back as the strings they were read as, so only the three new columns
    change. The replacement is atomic and replaces the directory entry, which matters because
    Nextflow's publishDir may have left a symlink into the work directory here.

    Args:
        tsv_path: The aggregated TSV to rewrite in place.
        rows: The parsed rows, annotated, in file order.
        fieldnames: The file's original header, in order.

    Returns:
        Number of rows written.
    """
    header = list(fieldnames) + [column for column in CLASSIFICATION_COLUMNS if column not in fieldnames]
    tmp_path = tsv_path.with_name(tsv_path.name + ".tmp")
    with tmp_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in header})
    os.replace(tmp_path, tsv_path)  # noqa: PTH105 - replaces the symlink itself, not its target
    return len(rows)


def _as_bool(value: Any) -> bool:
    """Interpret a persisted flag, which may be a real bool or the string a TSV round-tripped."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}
