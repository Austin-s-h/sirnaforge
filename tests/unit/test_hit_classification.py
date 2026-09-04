"""Unit tests for four-way hit classification.

Tests the pure classifier logic without requiring BWA, filesystem I/O, or
real alignment. All tests use synthetic indices and TP53 as the documented
example gene.
"""

from typing import Any

import pytest

from sirnaforge.core import hit_classification
from sirnaforge.core.hit_classification import (
    ClassificationContext,
    HitClass,
    HitClassCounts,
    classify_hit,
)
from sirnaforge.core.repeat_detection import normalize_guide_sequence
from sirnaforge.data.transcript_index import (
    SpeciesTranscriptIndex,
    TranscriptGeneIndex,
    TranscriptRecord,
)


def _build_synthetic_index() -> TranscriptGeneIndex:
    """Build a minimal multi-species index for testing.

    Human TP53:
        - ENST00000269305 (canonical): gene ENSG00000141510, symbol TP53
        - ENST00000445888 (isoform): gene ENSG00000141510, symbol TP53

    Mouse Tp53 (ortholog, uppercased to TP53 on ingest):
        - ENSMUST00000108658: gene ENSMUSG00000059552, symbol TP53

    Rat (no TP53 annotation):
        - ENSRNOT00000012345: gene ENSRNOG00000054321, symbol Gapdh
        - ENSRNOT00000067890: gene ENSRNOG00000011111, (NO SYMBOL)
    """
    index = TranscriptGeneIndex()

    # Human index
    human_index = SpeciesTranscriptIndex(species="human")
    human_index.add(
        TranscriptRecord(
            transcript_id="ENST00000269305",
            gene_id="ENSG00000141510",
            gene_symbol="TP53",
            biotype="protein_coding",
        )
    )
    human_index.add(
        TranscriptRecord(
            transcript_id="ENST00000445888",
            gene_id="ENSG00000141510",
            gene_symbol="TP53",
            biotype="protein_coding",
        )
    )
    index._indices["human"] = human_index

    # Mouse index with ortholog symbol (uppercased on ingest, like human)
    mouse_index = SpeciesTranscriptIndex(species="mouse")
    mouse_index.add(
        TranscriptRecord(
            transcript_id="ENSMUST00000108658",
            gene_id="ENSMUSG00000059552",
            gene_symbol="TP53",  # Same symbol as human after uppercasing
            biotype="protein_coding",
        )
    )
    index._indices["mouse"] = mouse_index

    # Rat index with a transcript missing a symbol
    rat_index = SpeciesTranscriptIndex(species="rat")
    rat_index.add(
        TranscriptRecord(
            transcript_id="ENSRNOT00000012345",
            gene_id="ENSRNOG00000054321",
            gene_symbol="GAPDH",
            biotype="protein_coding",
        )
    )
    rat_index.add(
        TranscriptRecord(
            transcript_id="ENSRNOT00000067890",
            gene_id="ENSRNOG00000011111",
            gene_symbol=None,  # No symbol annotation
            biotype="protein_coding",
        )
    )
    index._indices["rat"] = rat_index

    return index


def _make_context(
    query_gene_ids: frozenset[str] | None = None,
    query_gene_symbols: frozenset[str] | None = None,
    on_target_transcript_ids: frozenset[str] | None = None,
    repeat_guides: frozenset[str] | None = None,
) -> ClassificationContext:
    """Build a minimal context for TP53 testing."""
    return ClassificationContext(
        query_gene_ids=frozenset(["ENSG00000141510"]) if query_gene_ids is None else query_gene_ids,
        query_gene_symbols=frozenset(["TP53"]) if query_gene_symbols is None else query_gene_symbols,
        on_target_transcript_ids=frozenset(["ENST00000269305"])
        if on_target_transcript_ids is None
        else on_target_transcript_ids,
        query_species="human",
        index=_build_synthetic_index(),
        repeat_flagged_guides=frozenset() if repeat_guides is None else repeat_guides,
        requested_species=frozenset(["human", "mouse", "rat"]),
    )


def _make_hit(
    rname: str,
    species: str | None = None,
    nm: int = 0,
    seed_mismatches: int = 0,
) -> dict[str, Any]:
    """Build a minimal hit dict."""
    hit: dict[str, Any] = {
        "rname": rname,
        "nm": nm,
        "seed_mismatches": seed_mismatches,
    }
    if species is not None:
        hit["species"] = species
    return hit


@pytest.mark.unit
def test_on_target_by_transcript_id():
    """ON_TARGET when hit transcript ID is in on_target_transcript_ids."""
    context = _make_context()
    hit = _make_hit(rname="ENST00000269305", species="human")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ON_TARGET
    assert result.matched_symbol is None  # Matched by transcript ID, not symbol
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_on_target_by_gene_id():
    """ON_TARGET when hit's gene ID matches query gene IDs via index."""
    context = _make_context(on_target_transcript_ids=frozenset())  # No transcript IDs
    hit = _make_hit(rname="ENST00000269305", species="human")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ON_TARGET
    assert result.matched_symbol == "TP53"
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_on_target_by_gene_symbol():
    """ON_TARGET when hit's symbol matches query symbols via index."""
    context = _make_context(
        query_gene_ids=frozenset(),  # No gene IDs
        on_target_transcript_ids=frozenset(),  # No transcript IDs
    )
    hit = _make_hit(rname="ENST00000269305", species="human")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ON_TARGET
    assert result.matched_symbol == "TP53"
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_on_target_sibling_isoform():
    """ON_TARGET for a sibling isoform via shared gene ID.

    The isoform ENST00000445888 is NOT in on_target_transcript_ids (simulating
    a gene search that returned only the canonical transcript), but it shares
    gene ID ENSG00000141510 with the canonical, so the index recognizes it as
    on-target.
    """
    context = _make_context(on_target_transcript_ids=frozenset(["ENST00000269305"]))
    hit = _make_hit(rname="ENST00000445888", species="human")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ON_TARGET
    assert result.matched_symbol == "TP53"
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_ortholog_symbol_match():
    """ORTHOLOG when hit is in a different species with matching symbol."""
    context = _make_context()
    hit = _make_hit(rname="ENSMUST00000108658", species="mouse")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ORTHOLOG
    assert result.matched_symbol == "TP53"
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_ortholog_case_insensitive():
    """ORTHOLOG when symbol differs ONLY in capitalisation.

    This test verifies that symbols are case-normalized during index building.
    Since transcript_index.py uppercases symbols on ingest, comparing mouse "Tp53"
    (uppercased to "TP53") against human "TP53" succeeds.
    """
    # The default index already has mouse with "TP53" (uppercased on ingest)
    context = _make_context()
    hit = _make_hit(rname="ENSMUST00000108658", species="mouse")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ORTHOLOG
    assert result.matched_symbol == "TP53"
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_ortholog_no_symbol_annotation():
    """OFF_TARGET with symbol_lookup_missing=True when ortholog check cannot run.

    Rat transcript ENSRNOT00000067890 has no symbol annotation. We cannot
    determine if it's an ortholog, so it stays OFF_TARGET and increments the
    shortfall counter.
    """
    context = _make_context()
    hit = _make_hit(rname="ENSRNOT00000067890", species="rat")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.OFF_TARGET
    assert result.matched_symbol is None
    assert result.symbol_lookup_missing is True


@pytest.mark.unit
def test_no_species_index():
    """OFF_TARGET when the species has no index at all.

    Rhesus is requested but not indexed. This is distinct from "indexed but
    symbol missing" and is tracked separately in HitClassCounts.no_species_index.
    """
    context = _make_context()
    hit = _make_hit(rname="ENSMMUT00000012345", species="macaque")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.OFF_TARGET
    assert result.matched_symbol is None
    assert result.symbol_lookup_missing is False  # Not symbol-missing; it's no-index


@pytest.mark.unit
def test_repeat_flagged():
    """REPEAT when guide is in repeat_flagged_guides."""
    repeat_guide = "AAAAAAAAAAAAAAAAAAAA"
    context = _make_context(repeat_guides=frozenset([repeat_guide.upper()]))
    hit = _make_hit(rname="ENSRNOT00000012345", species="rat")
    result = classify_hit(hit, guide_sequence=repeat_guide, context=context)

    assert result.hit_class == HitClass.REPEAT
    assert result.matched_symbol is None
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_precedence_on_target_beats_repeat():
    """ON_TARGET outranks REPEAT: a repeat-flagged guide's hit on its own gene is ON_TARGET."""
    repeat_guide = "ACGTACGTACGTACGTACGTA"
    context = _make_context(repeat_guides=frozenset([repeat_guide.upper()]))
    hit = _make_hit(rname="ENST00000269305", species="human")
    result = classify_hit(hit, guide_sequence=repeat_guide, context=context)

    assert result.hit_class == HitClass.ON_TARGET
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_precedence_ortholog_beats_repeat():
    """ORTHOLOG outranks REPEAT: a repeat-flagged guide's hit on its ortholog is ORTHOLOG."""
    repeat_guide = "ACGTACGTACGTACGTACGTA"
    context = _make_context(repeat_guides=frozenset([repeat_guide.upper()]))
    hit = _make_hit(rname="ENSMUST00000108658", species="mouse")
    result = classify_hit(hit, guide_sequence=repeat_guide, context=context)

    assert result.hit_class == HitClass.ORTHOLOG
    assert result.matched_symbol == "TP53"
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_off_target_default():
    """OFF_TARGET for a hit that doesn't match any higher-precedence category."""
    context = _make_context()
    hit = _make_hit(rname="ENSRNOT00000012345", species="rat")  # Gapdh, not TP53
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.OFF_TARGET
    assert result.matched_symbol is None
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_empty_species_label_treated_as_query_species():
    """A hit with empty/None species is treated as the query species (ON_TARGET if matched)."""
    context = _make_context()

    # Test None species
    hit_none = _make_hit(rname="ENST00000269305", species=None)
    result_none = classify_hit(hit_none, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)
    assert result_none.hit_class == HitClass.ON_TARGET

    # Test empty string species
    hit_empty = _make_hit(rname="ENST00000269305", species="")
    result_empty = classify_hit(hit_empty, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)
    assert result_empty.hit_class == HitClass.ON_TARGET


@pytest.mark.unit
def test_version_suffixed_transcript_id():
    """Version suffixes (.9) are stripped before comparison."""
    context = _make_context()
    hit = _make_hit(rname="ENST00000269305.14", species="human")  # Version suffix
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.ON_TARGET
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_guide_normalisation():
    """Guide sequences are normalised (uppercase, U->T) before comparison."""
    repeat_guide_rna = "uuuuuuuuuuuuuuuuuuuu"  # RNA with lowercase u
    repeat_guide_dna = "TTTTTTTTTTTTTTTTTTTT"  # DNA uppercase
    context = _make_context(repeat_guides=frozenset([repeat_guide_dna]))

    hit = _make_hit(rname="ENSRNOT00000012345", species="rat")
    result = classify_hit(hit, guide_sequence=repeat_guide_rna, context=context)

    assert result.hit_class == HitClass.REPEAT


@pytest.mark.unit
def test_hit_class_counts_aggregation():
    """HitClassCounts aggregates correctly with add() and + operator."""
    counts1 = HitClassCounts(
        on_target=5,
        ortholog=3,
        repeat=2,
        off_target=10,
        symbol_lookup_missing=1,
        no_species_index=0,
        ortholog_species=frozenset(["mouse"]),
    )

    counts2 = HitClassCounts(
        on_target=2,
        ortholog=1,
        repeat=0,
        off_target=7,
        symbol_lookup_missing=0,
        no_species_index=1,
        ortholog_species=frozenset(["rat"]),
    )

    # Test .add()
    result_add = counts1.add(counts2)
    assert result_add.on_target == 7
    assert result_add.ortholog == 4
    assert result_add.repeat == 2
    assert result_add.off_target == 17
    assert result_add.symbol_lookup_missing == 1
    assert result_add.no_species_index == 1
    assert result_add.ortholog_species == frozenset(["mouse", "rat"])

    # Test + operator
    result_plus = counts1 + counts2
    assert result_plus.on_target == 7
    assert result_plus.ortholog == 4
    assert result_plus.repeat == 2
    assert result_plus.off_target == 17
    assert result_plus.symbol_lookup_missing == 1
    assert result_plus.no_species_index == 1
    assert result_plus.ortholog_species == frozenset(["mouse", "rat"])


@pytest.mark.unit
def test_classification_multiple_hits_aggregate():
    """Aggregate classifications across multiple hits to produce summary counts."""
    context = _make_context()

    # Simulate classifying several hits
    hits_and_guides = [
        (_make_hit("ENST00000269305", "human"), "ACGTACGTACGTACGTACGTA"),  # ON_TARGET
        (_make_hit("ENSMUST00000108658", "mouse"), "ACGTACGTACGTACGTACGTA"),  # ORTHOLOG
        (_make_hit("ENSRNOT00000012345", "rat"), "ACGTACGTACGTACGTACGTA"),  # OFF_TARGET
        (_make_hit("ENSRNOT00000067890", "rat"), "GCTAGCTAGCTAGCTAGCTAG"),  # OFF_TARGET + symbol_lookup_missing
    ]

    aggregate = HitClassCounts()
    for hit, guide in hits_and_guides:
        classification = classify_hit(hit, guide, context)

        if classification.hit_class == HitClass.ON_TARGET:
            aggregate.on_target += 1
        elif classification.hit_class == HitClass.ORTHOLOG:
            aggregate.ortholog += 1
            # Track ortholog species
            hit_species = hit.get("species")
            if hit_species:
                aggregate.ortholog_species = frozenset(aggregate.ortholog_species | {hit_species})
        elif classification.hit_class == HitClass.REPEAT:
            aggregate.repeat += 1
        elif classification.hit_class == HitClass.OFF_TARGET:
            aggregate.off_target += 1

        if classification.symbol_lookup_missing:
            aggregate.symbol_lookup_missing += 1

    assert aggregate.on_target == 1
    assert aggregate.ortholog == 1
    assert aggregate.repeat == 0
    assert aggregate.off_target == 2
    assert aggregate.symbol_lookup_missing == 1
    assert aggregate.ortholog_species == frozenset(["mouse"])


@pytest.mark.unit
def test_no_transcript_id():
    """A hit with no transcript ID defaults to OFF_TARGET."""
    context = _make_context()
    hit = _make_hit(rname="", species="human")
    result = classify_hit(hit, guide_sequence="ACGTACGTACGTACGTACGTA", context=context)

    assert result.hit_class == HitClass.OFF_TARGET
    assert result.matched_symbol is None
    assert result.symbol_lookup_missing is False


@pytest.mark.unit
def test_species_index_missing_is_reported_separately_from_missing_symbol():
    """A missing reference and a thin annotation are different problems, reported separately."""
    context = _make_context()

    # Macaque was never indexed at all.
    no_index = classify_hit(_make_hit("ENSMMUT00000099999", species="rhesus"), "ACGTACGTACGTACGTACGTA", context)
    assert no_index.hit_class is HitClass.OFF_TARGET
    assert no_index.species_index_missing is True
    assert no_index.symbol_lookup_missing is False

    # Rat is indexed, but this transcript carries no gene_symbol.
    no_symbol = classify_hit(_make_hit("ENSRNOT00000067890", species="rat"), "ACGTACGTACGTACGTACGTA", context)
    assert no_symbol.hit_class is HitClass.OFF_TARGET
    assert no_symbol.symbol_lookup_missing is True
    assert no_symbol.species_index_missing is False


@pytest.mark.unit
def test_repeat_verdict_survives_an_unresolvable_ortholog_check():
    """A guide we could not check for orthology is still a repeat-overlapping guide.

    Regression: both shortfall branches used to return OFF_TARGET before the repeat check ran,
    so a repeat-flagged guide hitting an unindexed species was silently downgraded.
    """
    guide = "ACGTACGTACGTACGTACGTA"
    context = _make_context(repeat_guides=frozenset([guide]))

    # Species has no index at all.
    no_index = classify_hit(_make_hit("ENSMMUT00000099999", species="rhesus"), guide, context)
    assert no_index.hit_class is HitClass.REPEAT
    assert no_index.species_index_missing is True

    # Species is indexed but the transcript has no symbol.
    no_symbol = classify_hit(_make_hit("ENSRNOT00000067890", species="rat"), guide, context)
    assert no_symbol.hit_class is HitClass.REPEAT
    assert no_symbol.symbol_lookup_missing is True

    # And the shortfall must not be invented where the check did run cleanly.
    resolved = classify_hit(_make_hit("ENSRNOT00000012345", species="rat"), guide, context)
    assert resolved.hit_class is HitClass.REPEAT
    assert resolved.symbol_lookup_missing is False
    assert resolved.species_index_missing is False


@pytest.mark.unit
def test_classifier_shares_one_normaliser_with_the_repeat_detector():
    """The REPEAT class only fires if both modules normalise guides identically."""
    assert hit_classification.normalize_guide_sequence is normalize_guide_sequence

    # An RNA guide flagged by the detector (which stores the DNA form) still classifies as REPEAT.
    context = _make_context(repeat_guides=frozenset([normalize_guide_sequence("acguacguacguacguacgua")]))
    result = classify_hit(_make_hit("ENSRNOT00000012345", species="rat"), "ACGUACGUACGUACGUACGUA", context)
    assert result.hit_class is HitClass.REPEAT
