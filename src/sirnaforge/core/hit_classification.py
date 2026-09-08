"""Four-way hit classifier for multi-species off-target analysis.

This module classifies screening hits into four mutually exclusive categories:
ON_TARGET, ORTHOLOG, REPEAT, and OFF_TARGET. The classifier is pure (no I/O,
no alignment) and operates on pre-computed indices and hit metadata.

Classification precedence:
    1. ON_TARGET - hit is on the query gene in the query species
    2. ORTHOLOG - hit is on an ortholog in a different species (symbol match)
    3. REPEAT - guide overlaps a repeat element
    4. OFF_TARGET - everything else

On-target and ortholog deliberately outrank repeat: a repeat-flagged guide's
hits on its own gene or orthologs are still classified as such, preserving the
decomposition needed for accurate specificity reporting.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from sirnaforge.core.repeat_detection import normalize_guide_sequence
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcript_index import TranscriptGeneIndex


class HitClass(str, Enum):
    """Mutually exclusive hit classification categories."""

    ON_TARGET = "on_target"
    ORTHOLOG = "ortholog"
    REPEAT = "repeat"
    OFF_TARGET = "off_target"


@dataclass(frozen=True)
class ClassificationContext:
    """Immutable context for hit classification.

    Attributes:
        query_gene_ids: Version-stripped gene IDs of the query gene, uppercased for comparison.
        query_gene_symbols: Uppercased gene symbols of the query gene.
        on_target_transcript_ids: Version-stripped transcript IDs directly present in the input.
        query_species: Canonical name of the species the query gene belongs to.
        index: Multi-species transcript→gene index for symbol/gene lookups.
        repeat_flagged_guides: Normalised guide sequences flagged as repeat-overlapping.
        requested_species: Canonical names of species the user requested to screen.
            (Used downstream for conservation scoring; does not affect classification verdicts.)
    """

    query_gene_ids: frozenset[str]
    query_gene_symbols: frozenset[str]
    on_target_transcript_ids: frozenset[str]
    query_species: str
    index: TranscriptGeneIndex
    repeat_flagged_guides: frozenset[str]
    requested_species: frozenset[str]


@dataclass(frozen=True)
class HitClassification:
    """Result of classifying one hit.

    Attributes:
        hit_class: The assigned class (ON_TARGET, ORTHOLOG, REPEAT, or OFF_TARGET).
        matched_symbol: Gene symbol that matched (for ON_TARGET or ORTHOLOG), or None.
        symbol_lookup_missing: True when the ortholog check could not run because the hit
            species' index carries no symbol for that transcript. Reported so a thin
            annotation is visible rather than silently misclassified.
        species_index_missing: True when the ortholog check could not run because there is
            no index for the hit species at all. A missing reference and a thin annotation
            are different problems, so they are counted separately.
    """

    hit_class: HitClass
    matched_symbol: str | None
    symbol_lookup_missing: bool
    species_index_missing: bool = False


@dataclass
class HitClassCounts:
    """Aggregate counts of classified hits with addable arithmetic.

    Attributes:
        on_target: Count of hits classified as ON_TARGET.
        ortholog: Count of hits classified as ORTHOLOG.
        repeat: Count of hits classified as REPEAT.
        off_target: Count of hits classified as OFF_TARGET.
        symbol_lookup_missing: Count of hits where ortholog check failed due to missing symbol.
        no_species_index: Count of hits where the species has no index at all.
        ortholog_species: Set of canonical species names with at least one ortholog hit.
    """

    on_target: int = 0
    ortholog: int = 0
    repeat: int = 0
    off_target: int = 0
    symbol_lookup_missing: int = 0
    no_species_index: int = 0
    ortholog_species: frozenset[str] = field(default_factory=frozenset)

    def add(self, other: HitClassCounts) -> HitClassCounts:
        """Aggregate two count sets, returning a new instance.

        Args:
            other: Another HitClassCounts to add.

        Returns:
            New HitClassCounts with summed counts and unioned ortholog species.
        """
        return HitClassCounts(
            on_target=self.on_target + other.on_target,
            ortholog=self.ortholog + other.ortholog,
            repeat=self.repeat + other.repeat,
            off_target=self.off_target + other.off_target,
            symbol_lookup_missing=self.symbol_lookup_missing + other.symbol_lookup_missing,
            no_species_index=self.no_species_index + other.no_species_index,
            ortholog_species=frozenset(self.ortholog_species | other.ortholog_species),
        )

    def __add__(self, other: HitClassCounts) -> HitClassCounts:
        """Support + operator for aggregation."""
        return self.add(other)


def classify_hit(
    hit: Mapping[str, Any],
    guide_sequence: str,
    context: ClassificationContext,
) -> HitClassification:
    """Classify one screening hit into the four-way taxonomy.

    This function is PURE: given the same inputs, it returns the same output.
    It performs no I/O, no logging, no mutation, and no alignment. This design
    makes the classification testable without running BWA.

    Classification applies to hits at all mismatch levels. The caller is
    responsible for passing only the hits it wants classified (typically those
    below a mismatch threshold). The redefined off_target_count includes hits
    at all mismatch levels that are not on-target, ortholog, or repeat-mediated.

    Args:
        hit: Screening hit dict with keys 'rname' (transcript ID), 'species',
            'nm' (mismatch count), 'seed_mismatches'.
        guide_sequence: The guide sequence that produced this hit.
        context: Classification context with query gene info and indices.

    Returns:
        HitClassification with the assigned class, matched symbol (if any),
        and shortfall flags.

    Note:
        Hit species strings arriving from screening may be any alias; they are
        normalized to canonical form before comparison. An empty or None species
        label means "query species" (preserving the current workflow behavior
        where untagged hits are treated as on-target candidates).

        Gene IDs are species-specific and cannot match across species. Symbol
        comparison is case-insensitive (mouse 'Tp53' matches human 'TP53') to
        handle ortholog naming conventions.
    """
    # Extract and normalize hit metadata
    hit_transcript_raw = str(hit.get("rname", ""))
    if not hit_transcript_raw:
        # No transcript ID means we cannot classify it; default to OFF_TARGET
        return HitClassification(
            hit_class=HitClass.OFF_TARGET,
            matched_symbol=None,
            symbol_lookup_missing=False,
        )

    # Strip version suffixes for comparison (e.g., "ENST00000123456.9" -> "ENST00000123456")
    hit_transcript_id = _strip_version(hit_transcript_raw)

    hit_species_raw = hit.get("species")
    if hit_species_raw is None or str(hit_species_raw).strip() == "":
        # Empty/None species label → treated as query species
        hit_species = context.query_species
    else:
        hit_species = normalize_species_name(str(hit_species_raw))

    normalized_guide = normalize_guide_sequence(guide_sequence)

    # PRECEDENCE ORDER (tests must pin this exact order):

    # 1. ON_TARGET: hit is in the query species and matches the query gene
    if hit_species == context.query_species:
        on_target = _classify_query_species(hit_transcript_id, context)
        if on_target is not None:
            return on_target
        symbol_lookup_missing = False
        species_index_missing = False
    else:
        # 2. ORTHOLOG: different species, symbol matches (case-insensitively)
        ortholog, symbol_lookup_missing, species_index_missing = _classify_other_species(
            hit_transcript_id, hit_species, context
        )
        if ortholog is not None:
            return ortholog

    # 3. REPEAT: guide is flagged as overlapping a repeat element. Falls through from the
    # shortfall branches above deliberately -- a guide we could not check for orthology is
    # still a repeat-overlapping guide, and returning early there would hide that.
    hit_class = HitClass.REPEAT if normalized_guide in context.repeat_flagged_guides else HitClass.OFF_TARGET

    # 4. OFF_TARGET: everything else
    return HitClassification(
        hit_class=hit_class,
        matched_symbol=None,
        symbol_lookup_missing=symbol_lookup_missing,
        species_index_missing=species_index_missing,
    )


def _classify_query_species(hit_transcript_id: str, context: ClassificationContext) -> HitClassification | None:
    """Return an ON_TARGET verdict for a query-species hit, or None if it is not on-target."""
    if hit_transcript_id in context.on_target_transcript_ids:
        # Matched by transcript ID, so no symbol was consulted.
        return HitClassification(hit_class=HitClass.ON_TARGET, matched_symbol=None, symbol_lookup_missing=False)

    species_index = context.index.for_species(context.query_species)
    if species_index is None:
        return None

    hit_symbol = species_index.symbol_for(hit_transcript_id)
    hit_gene_id = species_index.gene_id_for(hit_transcript_id)
    if (hit_gene_id and hit_gene_id in context.query_gene_ids) or (
        hit_symbol and hit_symbol in context.query_gene_symbols
    ):
        return HitClassification(hit_class=HitClass.ON_TARGET, matched_symbol=hit_symbol, symbol_lookup_missing=False)
    return None


def _classify_other_species(
    hit_transcript_id: str, hit_species: str, context: ClassificationContext
) -> tuple[HitClassification | None, bool, bool]:
    """Return an ORTHOLOG verdict for a non-query-species hit, plus the two shortfall flags.

    Gene IDs are species-specific and cannot match across species, so orthology is decided on
    the gene symbol alone. A symbol we do not have is never guessed at.
    """
    species_index = context.index.for_species(hit_species)
    if species_index is None:
        return None, False, True

    hit_symbol = species_index.symbol_for(hit_transcript_id)
    if hit_symbol is None:
        return None, True, False

    if hit_symbol in context.query_gene_symbols:
        verdict = HitClassification(
            hit_class=HitClass.ORTHOLOG,
            matched_symbol=hit_symbol,
            symbol_lookup_missing=False,
        )
        return verdict, False, False
    return None, False, False


def _strip_version(identifier: str) -> str:
    """Strip Ensembl version suffixes (e.g., .9) for comparison.

    Args:
        identifier: Ensembl ID with optional version suffix.

    Returns:
        ID without version suffix.

    Note:
        This duplicates the function in transcript_index.py to avoid circular
        imports and keep this module dependency-light.
    """
    return re.sub(r"\.\d+$", "", identifier.strip())
