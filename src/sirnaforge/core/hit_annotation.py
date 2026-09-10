"""Persist per-hit classification onto the aggregated off-target table.

Owned by #100; extended in place by #101 (new ``HitClass`` members). #103 consumes only and must not
edit it. ``hit_classification`` decides a hit's class once; this module writes that verdict onto the row
and derives the per-candidate counters *from the written row*, so no hit row can reach a candidate
counter without reaching a published table. That is a row-count guarantee, reported as a run warning
when it is violated — it is not per-candidate attribution, which can still be wrong in ways that net to
zero across the table.

Seven columns are written. ``hit_class`` is the persisted class; ``matched_symbol`` carries the
symbol that *established* the class (ortholog match, or a symbol-recognised on-target) and is the
literal string ``unknown`` otherwise — it is not a per-hit gene name. ``hit_symbol`` is the
per-row gene name: the symbol the transcript index resolves for ``rname``, independent of class,
with ``hit_symbol_missing`` flagging the rows it could not resolve. ``unknown`` rather than an
empty cell because a blank renders as "no gene", and on a real reference 14.6% of transcripts carry
no symbol at all, so ``unknown`` is a common and honest state.

The one decision this module does make is UNDETERMINED. A hit species with no transcript index
cannot be checked for orthology or for the query gene, so its alignments used to fall through to an
unqualified ``off_target``; ``species_index_missing`` is persisted alongside so the reason is on the
row.

``ortholog_evidence`` records *how* an ORTHOLOG verdict was reached — a resolved Compara gene-ID
mapping, or gene-symbol equality — and is the literal ``not_applicable``, never an empty cell, on
every other class. Without it a symbol-heuristic ortholog and a validated one are the same cell,
which #101 forbids: HGNC and MGI are different nomenclature authorities, so symbol equality both
misses real orthologues (mouse TP53 is ``Trp53``) and can assert false ones.
"""

from __future__ import annotations

import csv
import os
from collections.abc import Iterable, Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from sirnaforge.core.hit_classification import HitClass, HitClassCounts, HitClassification, OrthologEvidence
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcript_index import TranscriptGeneIndex

HIT_CLASS_COLUMN = "hit_class"
MATCHED_SYMBOL_COLUMN = "matched_symbol"
SYMBOL_LOOKUP_MISSING_COLUMN = "symbol_lookup_missing"
HIT_SYMBOL_COLUMN = "hit_symbol"
HIT_SYMBOL_MISSING_COLUMN = "hit_symbol_missing"
SPECIES_INDEX_MISSING_COLUMN = "species_index_missing"
ORTHOLOG_EVIDENCE_COLUMN = "ortholog_evidence"
CLASSIFICATION_COLUMNS: tuple[str, ...] = (
    HIT_CLASS_COLUMN,
    MATCHED_SYMBOL_COLUMN,
    SYMBOL_LOOKUP_MISSING_COLUMN,
    HIT_SYMBOL_COLUMN,
    HIT_SYMBOL_MISSING_COLUMN,
    SPECIES_INDEX_MISSING_COLUMN,
    ORTHOLOG_EVIDENCE_COLUMN,
)

UNKNOWN_SYMBOL = "unknown"

#: What the producer writes in every classification column before anything has classified the row.
#: The columns are written by ``aggregate_offtarget_results`` so the published column set does not
#: depend on which entry point produced the table, and this cell says the verdict has not been
#: reached yet -- distinct from ``unknown`` (looked up and not found) and never blank, because a
#: blank cell reads as missing data and ``is_annotated`` cannot tell it from a dropped column.
UNCLASSIFIED_CELL = "not_classified"

#: ``ortholog_evidence`` for a row that made no orthology claim. Spelled out rather than left blank
#: because the published table never carries an empty cell (a blank renders as missing data), and
#: distinct from ``unknown`` because nothing was looked up and failed -- there was nothing to look up.
ORTHOLOG_EVIDENCE_NOT_APPLICABLE = "not_applicable"

#: Every value the ``ortholog_evidence`` cell may hold. Must stay in step with
#: :data:`~sirnaforge.models.schemas.ORTHOLOG_EVIDENCE_VALUES`, which the schema enforces; kept here
#: rather than imported to avoid a core -> models dependency.
ORTHOLOG_EVIDENCE_CELL_VALUES: frozenset[str] = frozenset(
    {member.value for member in OrthologEvidence} | {ORTHOLOG_EVIDENCE_NOT_APPLICABLE}
)

# What the gates count. UNDETERMINED is absent evidence, not innocence: excluding it here would
# make a run with no transcript index pass candidates that today fail, i.e. a missing reference
# would loosen the screen. Until a resolved policy says what an unknown costs
# (EvidenceRequirements.unknown_evidence_action), an undecidable alignment is gated as a liability
# and reported separately, so the counted total and the qualified total are both visible.
LIABILITY_CLASSES: frozenset[HitClass] = frozenset({HitClass.OFF_TARGET, HitClass.UNDETERMINED})


def liabilities_counted(counts: HitClassCounts) -> int:
    """Hits gated as liabilities: genuine off-targets plus those whose class could not be decided."""
    return counts.off_target + counts.undetermined


@dataclass(frozen=True)
class HitAnnotator:
    """Per-row reference lookups: the gene name for a hit, and whether its species has an index.

    Separate from ``ClassificationContext`` because these are properties of the reference
    inventory, not of the query gene: the same lookups answer "what gene did this land on" for
    every row regardless of class.

    Attributes:
        index: Multi-species transcript→gene index.
        query_species: Canonical species a blank ``species`` label belongs to.
    """

    index: TranscriptGeneIndex
    query_species: str

    def species_of(self, row: Mapping[str, Any]) -> str:
        """Canonical species of a hit row; a blank label is the query species (see classify_hit)."""
        label = row.get("species")
        if label is None or str(label).strip() == "":
            return self.query_species
        return normalize_species_name(str(label))

    def has_index_for(self, row: Mapping[str, Any]) -> bool:
        """Whether any transcript index exists for this row's species."""
        return self.index.for_species(self.species_of(row)) is not None

    def symbol_for(self, row: Mapping[str, Any]) -> str | None:
        """Gene symbol the index resolves for this row's ``rname``, or None."""
        species_index = self.index.for_species(self.species_of(row))
        if species_index is None:
            return None
        rname = str(row.get("rname") or "")
        return species_index.symbol_for(rname) if rname else None


def unclassified_cells() -> dict[str, str]:
    """The classification columns as the producer writes them: present, and explicitly undecided.

    ``is_annotated`` is False for these rows, so the workflow's classifier fills them in place and
    the column set never changes between the two entry points.
    """
    return dict.fromkeys(CLASSIFICATION_COLUMNS, UNCLASSIFIED_CELL)


def annotate_hit_row(
    row: MutableMapping[str, Any],
    classification: HitClassification,
    annotator: HitAnnotator,
) -> HitClass:
    """Write a classification and its reference lookups onto a hit row, returning the class written.

    The return value is read back out of the row rather than taken from ``classification``, so a
    caller cannot count one class and publish another. The annotator is required, not optional: a
    row carrying a class but no gene name would publish two different column sets.
    """
    species_index_missing = not annotator.has_index_for(row)
    hit_symbol = annotator.symbol_for(row)

    row[HIT_CLASS_COLUMN] = _persisted_class(classification, species_index_missing).value
    row[MATCHED_SYMBOL_COLUMN] = classification.matched_symbol or UNKNOWN_SYMBOL
    row[SYMBOL_LOOKUP_MISSING_COLUMN] = classification.symbol_lookup_missing
    row[HIT_SYMBOL_COLUMN] = hit_symbol or UNKNOWN_SYMBOL
    row[HIT_SYMBOL_MISSING_COLUMN] = hit_symbol is None
    row[SPECIES_INDEX_MISSING_COLUMN] = species_index_missing
    row[ORTHOLOG_EVIDENCE_COLUMN] = (
        classification.ortholog_evidence.value
        if classification.ortholog_evidence is not None
        else ORTHOLOG_EVIDENCE_NOT_APPLICABLE
    )
    return hit_class_of(row)


def _persisted_class(classification: HitClassification, species_index_missing: bool) -> HitClass:
    """Downgrade an ``off_target`` fallthrough to UNDETERMINED when no reference could be consulted.

    ON_TARGET, ORTHOLOG and REPEAT were each decided by positive evidence (a transcript ID, a
    resolved orthologue gene ID or matched symbol, a repeat-flagged guide) and stand whatever the
    reference inventory looks like.
    OFF_TARGET is the only verdict reached by exclusion, so it is the only one a missing index
    invalidates.
    """
    if species_index_missing and classification.hit_class is HitClass.OFF_TARGET:
        return HitClass.UNDETERMINED
    return classification.hit_class


def hit_class_of(row: Mapping[str, Any]) -> HitClass:
    """Read the persisted class off a hit row.

    Raises:
        KeyError: The row was never annotated.
        ValueError: The row carries a value outside the taxonomy.
    """
    return HitClass(str(row[HIT_CLASS_COLUMN]))


def is_annotated(row: Mapping[str, Any]) -> bool:
    """True when this row carries a usable value in **every** classification column.

    Checks the class, both symbols and the three flag columns, not ``hit_class`` alone. Presence of a
    key is not enough: a table read back from disk can carry an empty cell, and treating that as
    annotated republished the blank — for ``hit_class`` that raised ``ValueError`` on the next read
    after the file had already been overwritten, and for the others it published ``hit_symbol=''`` /
    ``hit_symbol_missing=''`` / ``species_index_missing=''``, which ``AggregatedOffTargetSchema``
    rejects on the very table the same run publishes against it. A partially annotated row is
    re-annotated by the orphan pass, which stays the single repair path.

    ``ortholog_evidence`` is checked the same way, against its own vocabulary rather than for mere
    non-emptiness: a table written before that column existed carries it blank, and treating such a
    row as annotated republished the blank, which the schema's ``isin`` then rejects on the very
    table the same run publishes.
    """
    if str(row.get(HIT_CLASS_COLUMN) or "") not in {member.value for member in HitClass}:
        return False
    if not str(row.get(MATCHED_SYMBOL_COLUMN) or "").strip():
        return False
    if not str(row.get(HIT_SYMBOL_COLUMN) or "").strip():
        return False
    if str(row.get(ORTHOLOG_EVIDENCE_COLUMN) or "") not in ORTHOLOG_EVIDENCE_CELL_VALUES:
        return False
    return all(
        _is_parseable_flag(row.get(column))
        for column in (SYMBOL_LOOKUP_MISSING_COLUMN, HIT_SYMBOL_MISSING_COLUMN, SPECIES_INDEX_MISSING_COLUMN)
    )


def accumulate_hit_class(
    row: Mapping[str, Any],
    counts: HitClassCounts,
    species_bucket: MutableMapping[str, int] | None,
    hit_species: str,
) -> HitClass:
    """Add one persisted row to the candidate-level and, optionally, per-species counters.

    Every quantity here is derived from the row's own columns, which is what keeps the per-hit
    table and the per-candidate totals from disagreeing. The counter field names are the
    ``HitClass`` values verbatim, so there is no class-to-counter mapping to get wrong.

    ``species_bucket=None`` skips the per-species tally, for a caller counting alignments once
    each rather than once per candidate that shares them.
    """
    hit_class = hit_class_of(row)
    field = hit_class.value
    setattr(counts, field, getattr(counts, field) + 1)
    if species_bucket is not None:
        species_bucket[field] = species_bucket.get(field, 0) + 1

    if _as_bool(row.get(SYMBOL_LOOKUP_MISSING_COLUMN)):
        counts.symbol_lookup_missing += 1
        if species_bucket is not None:
            species_bucket["symbol_lookup_missing"] = species_bucket.get("symbol_lookup_missing", 0) + 1

    if _as_bool(row.get(SPECIES_INDEX_MISSING_COLUMN)):
        counts.no_species_index += 1
        if species_bucket is not None:
            species_bucket[SPECIES_INDEX_MISSING_COLUMN] = species_bucket.get(SPECIES_INDEX_MISSING_COLUMN, 0) + 1

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
    """Rewrite an off-target TSV with the classification columns appended.

    Values are written back as the strings they were read as, so only the new columns change. The
    replacement is atomic and replaces the directory entry, which matters because Nextflow's
    publishDir may have left a symlink into the work directory here.

    Args:
        tsv_path: The TSV to rewrite in place.
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


_FLAG_TRUE: frozenset[str] = frozenset({"true", "1", "yes"})
_FLAG_FALSE: frozenset[str] = frozenset({"false", "0", "no"})


def _as_bool(value: Any) -> bool:
    """Interpret a persisted flag, which may be a real bool or the string a TSV round-tripped."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in _FLAG_TRUE


def _is_parseable_flag(value: Any) -> bool:
    """Whether a flag cell says something. A blank cell parses as False but means "not written"."""
    if isinstance(value, bool):
        return True
    return str(value or "").strip().lower() in (_FLAG_TRUE | _FLAG_FALSE)
