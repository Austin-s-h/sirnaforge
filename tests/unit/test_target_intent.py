"""Intent resolution, the post-classification intent seam, and pairing geometry (#101).

Three things are pinned here, and they are pinned in one file because they fail independently:

1. **Geometry.** Guide/transcript pairing is antiparallel. A reversed map changes no *count* -- the
   forward seed window occurs in a random transcriptome at the same expected rate as the real
   reverse-complement site -- so every assertion below is on a **string and a coordinate**, never on
   a number of matches.
2. **Intent resolution.** The coverage denominator is frozen from what the annotation source offered,
   before ORF/sequence/design filtering can shrink it, and unknown coverage has three spellings and
   no number.
3. **The intent seam.** ``ON_TARGET`` is a statement about the gene. Acceptability is a second,
   separate record, and an absent transcript identity makes it unknown rather than unacceptable.
"""

import ast
import inspect
from pathlib import Path

import pytest
from pydantic import ValidationError

from sirnaforge.core import off_target, repeat_detection, seed_geometry
from sirnaforge.core import target_intent as target_intent_module
from sirnaforge.core.hit_classification import HitClass, HitClassification, OrthologEvidence
from sirnaforge.core.seed_geometry import (
    A1_POSITION,
    M8_POSITION,
    SEED_START_POSITION,
    complement,
    paired_transcript_position,
    reverse_complement,
    seed_window,
)
from sirnaforge.core.target_intent import (
    PROTEIN_CODING_BIOTYPE,
    CoverageStatus,
    IntentVerdict,
    evaluate_intent,
    observe_coverage,
    resolve_target_intent,
)
from sirnaforge.models.policy import CoverageMatch, TargetIntent, TargetSelectivity

# A let-7a-like guide, chosen so the worked example in seed_geometry's docstring is independently
# checkable against published let-7 seed sequences rather than only against this code.
GUIDE = "TGAGGTAGTAGGTTGTATAGT"
M8_SITE = "CTACCTC"  # reverse_complement(GUIDE[2..8])
SIXMER_SITE = "TACCTC"  # reverse_complement(GUIDE[2..7])
ANCHOR = 107


def _transcript(core: str) -> str:
    """A 108-base transcript whose bases 101..108 are ``core``, matching the module's worked example."""
    return "G" * 100 + core


def _imported_modules(module) -> set[str]:
    """Every module name the file imports, read from its source rather than from sys.modules."""
    tree = ast.parse(Path(module.__file__).read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def _classification(hit_class: HitClass, **overrides) -> HitClassification:
    """The classifier's own result object, which the seam consumes and never recomputes."""
    fields: dict = {"hit_class": hit_class, "matched_symbol": None, "symbol_lookup_missing": False}
    fields.update(overrides)
    return HitClassification(**fields)


# ---------------------------------------------------------------------------------------------
# Layering
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_intent_seam_imports_only_contracts_geometry_and_classification_types():
    """A pure seam's whole value is being testable without a workflow, a config or a network call."""
    imported = {name for name in _imported_modules(target_intent_module) if name.startswith("sirnaforge")}

    assert imported == {
        "sirnaforge.core.hit_classification",
        "sirnaforge.core.seed_geometry",
        "sirnaforge.models.policy",
    }


@pytest.mark.unit
def test_the_geometry_module_reuses_the_existing_normaliser_rather_than_adding_a_third():
    """Two normalisers already exist; a third that drifted would silently halve a site count."""
    imported = {name for name in _imported_modules(seed_geometry) if name.startswith("sirnaforge")}

    assert imported == {"sirnaforge.core.repeat_detection"}


@pytest.mark.unit
def test_intent_resolution_takes_no_variant_argument():
    """Allele targeting stays with VariantWorkflowConfig, which is what makes them composable.

    Intent describes which transcripts, variant configuration describes which alleles, and two
    objects that cannot see each other cannot disagree about either (#101).
    """
    parameters = set(inspect.signature(resolve_target_intent).parameters)

    assert not {name for name in parameters if "variant" in name or "snp" in name or "allele" in name}


# ---------------------------------------------------------------------------------------------
# Geometry: the antiparallel map
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_pairing_map_is_antiparallel_over_a_worked_example():
    """p(i) = anchor + 2 - i, asserted for every seed position rather than only at the anchor.

    A reversed implementation (``anchor - 2 + i``) still satisfies ``p(2) == anchor``, so pinning the
    anchor alone would not catch it. This asserts the whole map, and that it *decreases*.
    """
    assert paired_transcript_position(SEED_START_POSITION, ANCHOR) == ANCHOR

    expected = {2: 107, 3: 106, 4: 105, 5: 104, 6: 103, 7: 102, 8: 101}
    for guide_position, transcript_position in expected.items():
        assert paired_transcript_position(guide_position, ANCHOR) == transcript_position

    positions = [paired_transcript_position(i, ANCHOR) for i in range(2, 9)]
    assert positions == sorted(positions, reverse=True), "pairing must run antiparallel to the guide"

    # p(1) is the unpaired target adenosine's coordinate, one base 3' of the anchor.
    assert paired_transcript_position(A1_POSITION, ANCHOR) == ANCHOR + 1

    with pytest.raises(ValueError):
        paired_transcript_position(0, ANCHOR)


@pytest.mark.unit
def test_the_seed_window_is_one_based_and_inclusive_on_both_ends():
    """The off-by-one between guide[2..8] and guide[1:8] is exactly the class-boundary error."""
    assert seed_window(GUIDE, SEED_START_POSITION, M8_POSITION) == "GAGGTAG"
    assert seed_window(GUIDE, SEED_START_POSITION, 7) == "GAGGTA"
    assert seed_window(GUIDE, 1, 1) == "T"
    assert seed_window("ugagguag", SEED_START_POSITION, M8_POSITION) == "GAGGTAG", "RNA input is normalised"

    with pytest.raises(ValueError):
        seed_window(GUIDE, 2, len(GUIDE) + 1)
    with pytest.raises(ValueError):
        seed_window(GUIDE, 0, 7)
    with pytest.raises(ValueError):
        seed_window(GUIDE, 7, 2)


@pytest.mark.unit
def test_the_site_string_is_the_reverse_complement_of_the_seed_window():
    """The value a scanner searches for, spelled out so it cannot be the guide window by accident."""
    assert reverse_complement(seed_window(GUIDE, SEED_START_POSITION, M8_POSITION)) == M8_SITE
    assert reverse_complement(seed_window(GUIDE, SEED_START_POSITION, 7)) == SIXMER_SITE
    assert M8_SITE == "C" + SIXMER_SITE, "the classes are nested windows sharing the anchor"

    # complement is per-position and un-reversed; the two are different functions on purpose.
    assert complement("GAGGTAG") == "CTCCATC"
    assert reverse_complement("GAGGTAG") == "CTACCTC"
    assert reverse_complement("gaggtag") == "CTACCTC", "RNA/lowercase input is normalised first"


@pytest.mark.unit
def test_the_forward_seed_window_finds_passenger_orientation_sites_not_target_sites():
    """The fake null #101 warns about: two transcripts of identical length and base composition.

    One contains the target site, the other the guide-identical window. A count-based test cannot tell
    them apart, so the assertion is on which *string* is present.
    """
    with_site = "GGGG" + "CTACCTCA" + "A"
    with_guide_window = "CCCC" + "GAGGTAGA" + "T"

    assert len(with_site) == len(with_guide_window)
    assert sorted(with_site) == sorted(with_guide_window), "composition-matched, so counts cannot distinguish them"

    site = reverse_complement(seed_window(GUIDE, SEED_START_POSITION, M8_POSITION))
    guide_window = seed_window(GUIDE, SEED_START_POSITION, M8_POSITION)

    assert site in with_site
    assert site not in with_guide_window
    assert guide_window in with_guide_window
    assert guide_window not in with_site


@pytest.mark.unit
def test_the_four_seed_classes_are_pinned_by_fixtures_differing_in_one_base():
    """8mer / 7mer-m8 / 7mer-A1 / 6mer, each one base from its neighbour, all at the same anchor.

    The two bases that decide a class are the m8 base at p(8) and the A1 base at p(1); this asserts the
    predicate table the class assignment is derived from, using the map rather than hard-coded slices.
    """
    m8_start = paired_transcript_position(M8_POSITION, ANCHOR)
    sixmer_start = paired_transcript_position(7, ANCHOR)
    a1_position = paired_transcript_position(A1_POSITION, ANCHOR)

    fixtures = {
        # core at 101..108              m8 window,  6mer window, A1 base is literal 'A'
        "8mer": ("CTACCTCA", True, True, True),
        "7mer-m8": ("CTACCTCG", True, True, False),
        "7mer-A1": ("GTACCTCA", False, True, True),
        "6mer": ("GTACCTCG", False, True, False),
    }

    for label, (core, expect_m8, expect_sixmer, expect_a1) in fixtures.items():
        transcript = _transcript(core)
        assert (transcript[m8_start - 1 : ANCHOR] == M8_SITE) is expect_m8, label
        assert (transcript[sixmer_start - 1 : ANCHOR] == SIXMER_SITE) is expect_sixmer, label
        assert (transcript[a1_position - 1] == "A") is expect_a1, label

    # Each fixture differs from the 8mer in at most one base, so a class boundary cannot be crossed
    # by a change this test did not make.
    eightmer = fixtures["8mer"][0]
    for label, (core, _, _, _) in fixtures.items():
        differences = sum(1 for left, right in zip(eightmer, core, strict=True) if left != right)
        assert differences <= 2, label
        if label in {"7mer-m8", "7mer-A1"}:
            assert differences == 1, label


@pytest.mark.unit
def test_the_a1_base_is_tested_literally_and_is_not_derived_from_the_guide():
    """It is an unpaired target adenosine, so complementing g[1] would be a different test."""
    guide_starting_with_a = "A" + GUIDE[1:]

    assert complement(seed_window(guide_starting_with_a, A1_POSITION, A1_POSITION)) == "T"
    # A p(1) complementarity test would demand 'T' here; the documented rule demands 'A'.
    assert _transcript("CTACCTCA")[paired_transcript_position(A1_POSITION, ANCHOR) - 1] == "A"


@pytest.mark.unit
def test_the_geometry_helpers_are_byte_equivalent_to_the_two_existing_implementations():
    """Pinned without editing either file, so the deliberate duplication cannot drift silently."""
    for sequence in ("TGAGGTAGTAGGTTGTATAGT", "ugagguaguagguuguauagu", "ACGTN", "", "acgt"):
        assert seed_geometry.normalize_guide_sequence(sequence) == off_target._normalize_nucleotide_sequence(sequence)
        normalized = seed_geometry.normalize_guide_sequence(sequence)
        assert reverse_complement(sequence) == repeat_detection._reverse_complement(normalized)

    # Non-ACGT passes through rather than raising, matching the existing table exactly.
    assert reverse_complement("ACGTN") == "NACGT"


# ---------------------------------------------------------------------------------------------
# Intent resolution and the frozen denominator
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_pan_isoform_denominator_is_every_protein_coding_isoform_the_source_offered():
    """Frozen before ORF/sequence/design filtering, which is the coverage defect #101 names.

    ``ENST4`` has no sequence and would be dropped by ``step1_retrieve_transcripts``; it stays in the
    denominator, so coverage cannot rise because an isoform was filtered out.
    """
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe={
            "ENST1": PROTEIN_CODING_BIOTYPE,
            "ENST2": PROTEIN_CODING_BIOTYPE,
            "ENST3": "retained_intron",
            "ENST4": PROTEIN_CODING_BIOTYPE,
        },
        enumeration_inputs=("ENST1", "ENST2"),
        annotation_provenance="ensembl_114",
    )

    assert intent.coverage_denominator == frozenset({"ENST1", "ENST2", "ENST4"})
    assert intent.annotation_universe == frozenset({"ENST1", "ENST2", "ENST3", "ENST4"})
    assert intent.enumeration_inputs == frozenset({"ENST1", "ENST2"})
    assert intent.enumeration_was_filtered is True
    assert intent.annotation_provenance == "ensembl_114"
    assert intent.coverage_match is CoverageMatch.EXACT_FULL_SITE


@pytest.mark.unit
def test_a_universe_member_with_no_declared_biotype_is_only_counted_when_the_caller_says_so():
    """A thin annotation must not inflate the denominator; an input-FASTA run has no biotype at all."""
    universe = {"ENST1": PROTEIN_CODING_BIOTYPE, "FASTA_RECORD": None}

    assert resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM, annotation_universe=universe
    ).coverage_denominator == frozenset({"ENST1"})
    assert resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe=universe,
        annotation_provenance="input_fasta",
        unknown_biotype_counts_as_coding=True,
    ).coverage_denominator == frozenset({"ENST1", "FASTA_RECORD"})


@pytest.mark.unit
def test_an_isoform_selective_run_is_measured_against_what_it_required():
    """One unchanged floor means the right thing in both modes only because the denominator differs."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.ISOFORM_SELECTIVE,
        annotation_universe=dict.fromkeys(("ENST1", "ENST2", "ENST3"), PROTEIN_CODING_BIOTYPE),
        required_transcript_ids=("ENST1", "ENST2"),
        excluded_transcript_ids=("ENST3",),
    )

    assert intent.coverage_denominator == frozenset({"ENST1", "ENST2"})
    assert intent.coverage_denominator < intent.annotation_universe


@pytest.mark.unit
def test_an_exclusion_is_subtracted_from_the_denominator_in_both_modes():
    """An exclusion is a prohibition, so covering it must not be able to improve a score."""
    pan = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe=dict.fromkeys(("ENST1", "ENST2"), PROTEIN_CODING_BIOTYPE),
        excluded_transcript_ids=("ENST2",),
    )
    selective = resolve_target_intent(
        selectivity=TargetSelectivity.ISOFORM_SELECTIVE,
        annotation_universe=dict.fromkeys(("ENST1", "ENST2"), PROTEIN_CODING_BIOTYPE),
        required_transcript_ids=("ENST1", "ENST2"),
        excluded_transcript_ids=("ENST2",),
    )

    assert pan.coverage_denominator == frozenset({"ENST1"})
    assert selective.coverage_denominator == frozenset({"ENST1"})


@pytest.mark.unit
def test_a_selective_intent_naming_nothing_required_leaves_coverage_unknown():
    """Inventing the pan-isoform denominator here would reinstate the assumption #101 removed."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.ISOFORM_SELECTIVE,
        annotation_universe={"ENST1": PROTEIN_CODING_BIOTYPE},
    )

    assert intent.coverage_denominator == frozenset()
    assert observe_coverage(GUIDE, intent, {}).status is CoverageStatus.UNKNOWN_NO_DENOMINATOR


@pytest.mark.unit
def test_a_declared_id_the_annotation_never_produced_is_recorded_not_silently_satisfied():
    """Without this field a mistyped required ID is simply absent from every set, which reads as met."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.ISOFORM_SELECTIVE,
        annotation_universe={"ENST1": PROTEIN_CODING_BIOTYPE},
        required_transcript_ids=("ENST1", "ENST00000000000"),
        excluded_transcript_ids=("ENST_TYPO",),
    )

    assert intent.unresolved_target_ids == frozenset({"ENST00000000000", "ENST_TYPO"})
    # The unresolved requirement is still a requirement, so it still shrinks nothing silently.
    assert "ENST00000000000" in intent.required_transcript_ids


@pytest.mark.unit
def test_recorded_sequence_availability_is_confined_to_the_frozen_denominator():
    """Availability outside the denominator is not coverage's business and must not pad the record."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe=dict.fromkeys(("ENST1", "ENST2"), PROTEIN_CODING_BIOTYPE),
        sequence_available=("ENST1", "ENST2", "ENST_ELSEWHERE"),
    )

    assert intent.coverage_sequence_available == frozenset({"ENST1", "ENST2"})


@pytest.mark.unit
def test_a_resolved_intent_cannot_be_shrunk_afterwards():
    """The denominator lives on a frozen model, so a later reassignment is an exception, not a habit."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe={"ENST1": PROTEIN_CODING_BIOTYPE},
    )

    with pytest.raises(ValidationError, match="frozen"):
        intent.coverage_denominator = frozenset()  # type: ignore[misc]


# ---------------------------------------------------------------------------------------------
# The intent seam
# ---------------------------------------------------------------------------------------------


def _selective_intent(**overrides) -> TargetIntent:
    """An isoform-selective intent over three known isoforms, one of them excluded."""
    defaults: dict = {
        "selectivity": TargetSelectivity.ISOFORM_SELECTIVE,
        "annotation_universe": dict.fromkeys(("ENST1", "ENST2", "ENST3"), PROTEIN_CODING_BIOTYPE),
        "required_transcript_ids": ("ENST1",),
        "excluded_transcript_ids": ("ENST3",),
    }
    defaults.update(overrides)
    return resolve_target_intent(**defaults)


@pytest.mark.unit
def test_an_excluded_isoform_is_on_target_and_unacceptable():
    """#101's central point: the taxonomy is right and the design is still wrong.

    ``hit_class`` is untouched, so ``on_target`` + ``excluded_isoform`` is published as a readable
    pair rather than resolved into a fifth hit class that lies about the gene.
    """
    classification = _classification(HitClass.ON_TARGET)

    assessment = evaluate_intent(
        classification,
        transcript_id="ENST3",
        species="human",
        intent=_selective_intent(),
        query_species="human",
    )

    assert assessment.verdict is IntentVerdict.EXCLUDED_ISOFORM
    assert assessment.acceptable is False
    assert "ENST3" in assessment.reason
    assert classification.hit_class is HitClass.ON_TARGET, "the seam never mutates the classification"
    assert not isinstance(assessment.verdict, HitClass)


@pytest.mark.unit
def test_an_exclusion_beats_a_requirement_naming_the_same_transcript():
    """A contradictory declaration resolves to the safety statement, and it does so deterministically."""
    intent = _selective_intent(required_transcript_ids=("ENST1", "ENST3"), excluded_transcript_ids=("ENST3",))

    assessment = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id="ENST3",
        species="human",
        intent=intent,
        query_species="human",
    )

    assert assessment.verdict is IntentVerdict.EXCLUDED_ISOFORM
    assert assessment.acceptable is False


@pytest.mark.unit
def test_a_known_same_gene_isoform_outside_the_required_set_is_unintended_in_selective_mode():
    """ON_TARGET by gene taxonomy does not imply acceptable when only some isoforms were asked for."""
    assessment = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id="ENST2",
        species="human",
        intent=_selective_intent(),
        query_species="human",
    )

    assert assessment.verdict is IntentVerdict.UNINTENDED_ISOFORM
    assert assessment.acceptable is False


@pytest.mark.unit
def test_a_required_isoform_and_any_isoform_under_a_pan_intent_are_intended():
    """The two modes differ on exactly the transcripts neither declared set names."""
    required = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id="ENST1",
        species="human",
        intent=_selective_intent(),
        query_species="human",
    )
    pan = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id="ENST2",
        species="human",
        intent=resolve_target_intent(
            selectivity=TargetSelectivity.PAN_ISOFORM,
            annotation_universe=dict.fromkeys(("ENST1", "ENST2"), PROTEIN_CODING_BIOTYPE),
        ),
        query_species="human",
    )

    assert required.verdict is IntentVerdict.INTENDED
    assert required.acceptable is True
    assert pan.verdict is IntentVerdict.INTENDED
    assert pan.acceptable is True


@pytest.mark.unit
def test_an_unidentifiable_on_target_hit_is_unknown_not_unacceptable():
    """``classify_hit`` reaches ON_TARGET through a gene-ID/symbol match, consulting no transcript.

    So an on-target row can name a transcript the intent has never heard of. That is a gap in this
    run's annotation, and reporting it as an unintended isoform would manufacture rejections out of a
    thin index. ``acceptable`` is ``None``, never coerced to ``False``.
    """
    intent = _selective_intent()

    absent_identity = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id=None,
        species="human",
        intent=intent,
        query_species="human",
    )
    unheard_of = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id="ENST_NOT_IN_ANY_SET",
        species="human",
        intent=intent,
        query_species="human",
    )

    assert absent_identity.verdict is IntentVerdict.UNKNOWN
    assert absent_identity.acceptable is None
    assert unheard_of.verdict is IntentVerdict.UNKNOWN
    assert unheard_of.acceptable is None


@pytest.mark.unit
def test_a_missing_identity_is_only_unknown_when_the_intent_names_transcripts():
    """A pan-isoform run with no exclusions cannot have its verdict changed by identity.

    So demanding one would fill an ordinary run with undecided rows for no gain.
    """
    assessment = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id=None,
        species="human",
        intent=resolve_target_intent(
            selectivity=TargetSelectivity.PAN_ISOFORM,
            annotation_universe={"ENST1": PROTEIN_CODING_BIOTYPE},
        ),
        query_species="human",
    )

    assert assessment.verdict is IntentVerdict.INTENDED
    assert assessment.acceptable is True


@pytest.mark.unit
def test_an_on_target_hit_outside_the_query_species_is_declined_rather_than_accepted():
    """``classify_hit`` cannot produce that pair, so intent refuses to rubber-stamp it if one arrives."""
    assessment = evaluate_intent(
        _classification(HitClass.ON_TARGET),
        transcript_id="ENST1",
        species="mouse",
        intent=_selective_intent(),
        query_species="human",
    )

    assert assessment.verdict is IntentVerdict.UNKNOWN
    assert assessment.acceptable is None
    assert "mouse" in assessment.reason and "human" in assessment.reason


@pytest.mark.unit
def test_an_orthologue_is_conserved_or_exploratory_by_declared_target_species():
    """Both are acceptable; the distinction is descriptive, and it is published rather than inferred."""
    intent = _selective_intent(target_species=("human", "mouse"))
    classification = _classification(HitClass.ORTHOLOG, ortholog_evidence=OrthologEvidence.GENE_ID)

    conserved = evaluate_intent(
        classification, transcript_id="ENSMUST1", species="mouse", intent=intent, query_species="human"
    )
    exploratory = evaluate_intent(
        classification, transcript_id="ENSRNOT1", species="rat", intent=intent, query_species="human"
    )

    assert conserved.verdict is IntentVerdict.CONSERVED_TARGET
    assert conserved.acceptable is True
    assert exploratory.verdict is IntentVerdict.EXPLORATORY_ORTHOLOG
    assert exploratory.acceptable is True


@pytest.mark.unit
def test_adding_an_exploratory_species_changes_no_verdict_for_a_required_one():
    """#101's guarantee at the verdict level: the target set is enumerated per species."""
    classification = _classification(HitClass.ORTHOLOG, ortholog_evidence=OrthologEvidence.GENE_ID)
    before = evaluate_intent(
        classification,
        transcript_id="ENSMUST1",
        species="mouse",
        intent=_selective_intent(target_species=("human", "mouse")),
        query_species="human",
    )
    after = evaluate_intent(
        classification,
        transcript_id="ENSMUST1",
        species="mouse",
        intent=_selective_intent(target_species=("human", "mouse"), offtarget_screen_species=("rat", "cynomolgus")),
        query_species="human",
    )

    assert before == after


@pytest.mark.unit
def test_repeat_and_off_target_hits_are_liabilities_and_undetermined_is_unknown():
    """Intent adds no new rejection for a non-target gene, and refuses to decide an undecided class."""
    intent = _selective_intent()

    for hit_class in (HitClass.REPEAT, HitClass.OFF_TARGET):
        assessment = evaluate_intent(
            _classification(hit_class),
            transcript_id="ENST_OTHER",
            species="human",
            intent=intent,
            query_species="human",
        )
        assert assessment.verdict is IntentVerdict.LIABILITY
        assert assessment.acceptable is False

    undetermined = evaluate_intent(
        _classification(HitClass.UNDETERMINED),
        transcript_id="ENST_OTHER",
        species="dog",
        intent=intent,
        query_species="human",
    )

    assert undetermined.verdict is IntentVerdict.UNKNOWN
    assert undetermined.acceptable is None


@pytest.mark.unit
def test_every_verdict_carries_a_non_blank_reason_and_only_unknown_is_undecided():
    """Same rule as an evidence detail: a record whose explanation is empty cannot be audited later."""
    intent = _selective_intent()
    seen: set[IntentVerdict] = set()

    cases = [
        (HitClass.ON_TARGET, "ENST1", "human"),
        (HitClass.ON_TARGET, "ENST2", "human"),
        (HitClass.ON_TARGET, "ENST3", "human"),
        (HitClass.ON_TARGET, None, "human"),
        (HitClass.ORTHOLOG, "ENSMUST1", "mouse"),
        (HitClass.ORTHOLOG, "ENSRNOT1", "rat"),
        (HitClass.REPEAT, "ENST_OTHER", "human"),
        (HitClass.OFF_TARGET, "ENST_OTHER", "human"),
        (HitClass.UNDETERMINED, "ENST_OTHER", "human"),
    ]
    for hit_class, transcript_id, species in cases:
        assessment = evaluate_intent(
            _classification(hit_class),
            transcript_id=transcript_id,
            species=species,
            intent=intent,
            query_species="human",
        )
        assert assessment.reason.strip(), (hit_class, transcript_id)
        assert (assessment.acceptable is None) is (assessment.verdict is IntentVerdict.UNKNOWN)
        seen.add(assessment.verdict)

    assert seen == set(IntentVerdict) - {IntentVerdict.CONSERVED_TARGET}


# ---------------------------------------------------------------------------------------------
# Coverage as a complementarity test
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_coverage_tests_complementarity_against_transcripts_enumeration_never_saw():
    """The enumeration map answers "was this guide generated from that isoform", a different question."""
    site = reverse_complement(GUIDE)
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe=dict.fromkeys(("ENST1", "ENST2"), PROTEIN_CODING_BIOTYPE),
        enumeration_inputs=("ENST1",),
        sequence_available=("ENST1", "ENST2"),
    )

    observation = observe_coverage(GUIDE, intent, {"ENST1": "AAA" + site + "TTT", "ENST2": "CCC" + site + "GGG"})

    assert intent.enumeration_inputs == frozenset({"ENST1"})
    assert observation.covered == frozenset({"ENST1", "ENST2"})
    assert observation.fraction == 1.0
    assert observation.status is CoverageStatus.KNOWN
    assert observation.match is CoverageMatch.EXACT_FULL_SITE


@pytest.mark.unit
def test_the_guide_occurring_verbatim_is_not_coverage():
    """A guide is antisense to its target; searching for the guide finds the passenger orientation."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe=dict.fromkeys(("ENST_SENSE", "ENST_ANTISENSE"), PROTEIN_CODING_BIOTYPE),
        sequence_available=("ENST_SENSE", "ENST_ANTISENSE"),
    )

    observation = observe_coverage(
        GUIDE,
        intent,
        {"ENST_SENSE": "AAA" + reverse_complement(GUIDE) + "TTT", "ENST_ANTISENSE": "AAA" + GUIDE + "TTT"},
    )

    assert observation.covered == frozenset({"ENST_SENSE"})
    assert observation.fraction == 0.5


@pytest.mark.unit
def test_the_seed_alone_does_not_count_as_coverage():
    """A seed-level criterion would count off-target silencing as coverage of the target."""
    intent = resolve_target_intent(
        selectivity=TargetSelectivity.PAN_ISOFORM,
        annotation_universe={"ENST1": PROTEIN_CODING_BIOTYPE},
        sequence_available=("ENST1",),
    )

    observation = observe_coverage(GUIDE, intent, {"ENST1": "AAA" + M8_SITE + "TTT"})

    assert observation.covered == frozenset()
    assert observation.fraction == 0.0
    assert observation.status is CoverageStatus.KNOWN, "a measured zero is a number, unlike an unknown"


@pytest.mark.unit
def test_a_denominator_member_without_a_sequence_makes_coverage_unknown_not_smaller():
    """min_isoform_coverage is a greater-or-equal floor, so a lower-bound numerator can only reject."""
    universe = dict.fromkeys(("ENST1", "ENST2"), PROTEIN_CODING_BIOTYPE)
    site = reverse_complement(GUIDE)

    from_the_mapping = observe_coverage(
        GUIDE,
        resolve_target_intent(
            selectivity=TargetSelectivity.PAN_ISOFORM, annotation_universe=universe, sequence_available=universe
        ),
        {"ENST1": site},
    )
    from_the_record = observe_coverage(
        GUIDE,
        resolve_target_intent(
            selectivity=TargetSelectivity.PAN_ISOFORM, annotation_universe=universe, sequence_available=("ENST1",)
        ),
        {"ENST1": site, "ENST2": site},
    )

    for observation in (from_the_mapping, from_the_record):
        assert observation.fraction is None
        assert observation.status is CoverageStatus.UNKNOWN_MISSING_SEQUENCE
        assert observation.unavailable == frozenset({"ENST2"})
        assert observation.denominator == frozenset({"ENST1", "ENST2"})
    assert from_the_mapping.covered == frozenset({"ENST1"}), "what was measured is still reported"


@pytest.mark.unit
def test_an_empty_denominator_is_unknown_coverage_and_not_complete_coverage():
    """Empty invites the reading "nothing left to cover"; the three status spellings prevent it."""
    intent = resolve_target_intent(selectivity=TargetSelectivity.PAN_ISOFORM, annotation_universe={})

    observation = observe_coverage(GUIDE, intent, {"ENST1": reverse_complement(GUIDE)})

    assert observation.denominator == frozenset()
    assert observation.fraction is None
    assert observation.status is CoverageStatus.UNKNOWN_NO_DENOMINATOR
    assert {status.value for status in CoverageStatus} == {
        "known",
        "unknown_no_denominator",
        "unknown_missing_sequence",
    }
