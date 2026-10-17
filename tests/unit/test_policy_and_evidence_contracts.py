"""The shared 0.7.1 contracts: run policy vocabulary and screening plan/evidence types.

These modules are types, not behaviour. Four downstream packages (#99, #100, #101, #103) fork off
them, so what is pinned here is the vocabulary and the deliberate absence of resolution: no
defaults that make a policy decision, no producer, no consumer.
"""

import ast
import dataclasses
import inspect
from pathlib import Path

import pytest
from pydantic import BaseModel, ValidationError

from sirnaforge.core.hit_classification import HitClass, HitClassCounts
from sirnaforge.models import evidence as evidence_module
from sirnaforge.models import policy as policy_module
from sirnaforge.models.evidence import (
    EVIDENCE_SCHEMA_VERSION,
    EvidenceStatus,
    HitCountCell,
    HitCountMatrix,
    ObservedCount,
    ObservedCounts,
    ScreeningEvidence,
    ScreeningEvidenceEntry,
    ScreeningPlan,
    ScreeningPlanEntry,
)
from sirnaforge.models.policy import (
    ChannelRequirement,
    EvidenceRequirements,
    FilterAction,
    FilterComparator,
    FilterDescriptor,
    FilterEvaluation,
    FilterScope,
    FilterVerdict,
    OrthologyAssertion,
    OrthologyEvidenceSource,
    Requiredness,
    RunMode,
    ScreeningChannel,
    TargetIntent,
    TargetSelectivity,
    UnknownEvidenceAction,
)

DIGEST = "sha256:abc"


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


def _models_defined_in(module) -> list[type[BaseModel]]:
    """Pydantic models declared in this module, so a new one cannot skip the config checks."""
    return [
        obj
        for _, obj in inspect.getmembers(module, inspect.isclass)
        if issubclass(obj, BaseModel) and obj is not BaseModel and obj.__module__ == module.__name__
    ]


def _intent(**overrides) -> TargetIntent:
    """A pan-isoform intent; selectivity has no default, so every construction states one."""
    return TargetIntent(selectivity=TargetSelectivity.PAN_ISOFORM, **overrides)


@pytest.mark.unit
def test_contract_modules_do_not_import_workflow():
    """These are data contracts: a dependency on the orchestrator would invert the layering."""
    for module in (policy_module, evidence_module):
        assert not any("workflow" in name for name in _imported_modules(module))


@pytest.mark.unit
def test_evidence_module_depends_only_on_the_policy_vocabulary():
    """Channel names and filter scope are shared, so evidence imports policy and never the reverse."""
    assert "sirnaforge.models.policy" in _imported_modules(evidence_module)
    assert not any(name.startswith("sirnaforge.models.evidence") for name in _imported_modules(policy_module))


@pytest.mark.unit
def test_contract_models_forbid_unknown_keys_and_are_frozen():
    """A misnamed key must not become an unset field: #99 reads field-set provenance off these models."""
    for module in (policy_module, evidence_module):
        for model in _models_defined_in(module):
            assert model.model_config.get("extra") == "forbid", f"{model.__name__} accepts unknown keys"
            assert model.model_config.get("frozen") is True, f"{model.__name__} is mutable"


@pytest.mark.unit
def test_a_misnamed_key_is_rejected_rather_than_silently_dropped():
    """The concrete failure this guards: a typo that reads as an omitted option rather than an error."""
    with pytest.raises(ValidationError):
        FilterScope(specie=frozenset({"human"}))  # type: ignore[call-arg]
    with pytest.raises(ValidationError):
        ObservedCount(values=3)  # type: ignore[call-arg]


@pytest.mark.unit
def test_run_modes_are_the_three_declared_states():
    """design-only, exploratory and qualified, spelled once."""
    assert {mode.value for mode in RunMode} == {"design_only", "exploratory", "qualified"}


@pytest.mark.unit
def test_action_and_evaluation_are_separate_vocabularies():
    """A filter's action is what it does; its evaluation is what it found."""
    assert {action.value for action in FilterAction} == {"off", "warn", "fail"}
    assert {result.value for result in FilterEvaluation} == {"pass", "fail", "unknown", "not_evaluated"}
    # "fail" is the only word the two share, and it means different things in each.
    assert FilterAction.FAIL.value == FilterEvaluation.FAIL.value
    assert FilterAction.FAIL is not FilterEvaluation.FAIL


@pytest.mark.unit
def test_unknown_and_not_evaluated_are_separate_first_class_evaluations():
    """Undecided and never-ran are different claims, and neither may borrow PASS or FAIL."""
    decided = (FilterEvaluation.PASS, FilterEvaluation.FAIL)
    assert FilterEvaluation.UNKNOWN not in decided
    assert FilterEvaluation.NOT_EVALUATED not in decided
    assert FilterEvaluation.NOT_EVALUATED is not FilterEvaluation.UNKNOWN


@pytest.mark.unit
def test_unknown_evidence_action_cannot_be_off():
    """OFF means never evaluated, so it is not an available response to an evaluated UNKNOWN."""
    assert {action.value for action in UnknownEvidenceAction} == {"warn", "fail"}
    assert {action.value for action in UnknownEvidenceAction} < {action.value for action in FilterAction}
    with pytest.raises(ValidationError):
        EvidenceRequirements(unknown_evidence_action="off")  # type: ignore[arg-type]


@pytest.mark.unit
def test_comparators_cover_both_directions_and_apply_themselves():
    """Direction used to live only in a max_/min_ prefix; here it is a value with one implementation."""
    assert {comparator.value for comparator in FilterComparator} == {"le", "lt", "ge", "gt"}
    assert FilterComparator.LE.passes(15, 15) and not FilterComparator.LE.passes(16, 15)
    assert FilterComparator.LT.passes(14, 15) and not FilterComparator.LT.passes(15, 15)
    assert FilterComparator.GE.passes(0.65, 0.65) and not FilterComparator.GE.passes(0.64, 0.65)
    assert FilterComparator.GT.passes(1, 0) and not FilterComparator.GT.passes(0, 0)


@pytest.mark.unit
def test_a_boolean_fail_if_any_flag_is_expressible_as_a_ceiling_of_zero():
    """fail_on_high_risk_mirna has no threshold field today; as data it is "at most zero"."""
    descriptor = FilterDescriptor(
        filter_id="fail_on_high_risk_mirna",
        column="mirna_hits_high_risk",
        comparator=FilterComparator.LE,
        threshold=0,
        scope=FilterScope(species=frozenset({"human"})),
        action=FilterAction.FAIL,
    )

    assert descriptor.comparator.passes(0, descriptor.threshold)
    assert not descriptor.comparator.passes(1, descriptor.threshold)


@pytest.mark.unit
def test_filter_scope_expresses_species_mismatch_ceiling_and_hit_class_as_data():
    """D3's "human hits at nm <= 2, genuine off-targets only" with nothing encoded in a field name."""
    scope = FilterScope(
        species=frozenset({"human"}),
        max_mismatches=2,
        hit_classes=frozenset({HitClass.OFF_TARGET.value}),
    )

    assert scope.species == frozenset({"human"})
    assert scope.max_mismatches == 2
    assert HitClass.ON_TARGET.value not in scope.hit_classes
    assert HitClass.REPEAT.value not in scope.hit_classes


@pytest.mark.unit
def test_filter_scope_accepts_every_hit_class_value_and_rejects_a_blank_one():
    """The strings are HitClass values; a blank member would narrow a scope to nothing real."""
    every_class = FilterScope(hit_classes=frozenset(member.value for member in HitClass))
    assert len(every_class.hit_classes) == len(list(HitClass))

    with pytest.raises(ValidationError):
        FilterScope(hit_classes=frozenset({""}))
    with pytest.raises(ValidationError):
        FilterScope(species=frozenset({" "}))


@pytest.mark.unit
def test_an_unrestricted_scope_is_the_default_and_must_still_be_stated():
    """A gate with no scope restriction says so explicitly; scope has no default on the descriptor."""
    assert FilterScope() == FilterScope(species=frozenset(), max_mismatches=None, hit_classes=frozenset())
    with pytest.raises(ValidationError):
        FilterDescriptor(
            filter_id="gc_content_max",
            column="gc_content",
            comparator=FilterComparator.LE,
            threshold=52.0,
            action=FilterAction.FAIL,
        )  # type: ignore[call-arg]


@pytest.mark.unit
def test_a_verdict_carries_the_descriptor_it_was_decided_under():
    """#103 re-thresholds in the browser, so a verdict read against another threshold is a wrong report."""
    verdict = FilterVerdict(
        descriptor=FilterDescriptor(
            filter_id="max_off_target_count",
            column="off_target_count",
            comparator=FilterComparator.LE,
            threshold=15,
            scope=FilterScope(
                species=frozenset({"human"}), max_mismatches=2, hit_classes=frozenset({HitClass.OFF_TARGET.value})
            ),
            action=FilterAction.FAIL,
        ),
        observed=40,
        evaluation=FilterEvaluation.FAIL,
    )

    assert verdict.descriptor.threshold == 15
    assert verdict.evaluation is FilterEvaluation.FAIL
    # The five keys #103 generates its gate descriptors from, straight off the model.
    assert set(verdict.descriptor.model_dump()) >= {"column", "comparator", "threshold", "action", "scope"}


@pytest.mark.unit
def test_a_decided_verdict_must_carry_both_numbers_it_was_decided_from():
    """PASS or FAIL with no observed value is the blank-reads-as-a-verdict failure, spelled in a model."""
    descriptor = FilterDescriptor(
        filter_id="max_transcriptome_hits_0mm",
        column="transcriptome_hits_0mm",
        comparator=FilterComparator.LE,
        threshold=1,
        scope=FilterScope(species=frozenset({"human"})),
        action=FilterAction.FAIL,
    )

    with pytest.raises(ValidationError):
        FilterVerdict(descriptor=descriptor, observed=None, evaluation=FilterEvaluation.PASS)

    undecided = FilterVerdict(descriptor=descriptor, observed=None, evaluation=FilterEvaluation.UNKNOWN)
    assert undecided.observed is None


@pytest.mark.unit
def test_an_undeclared_threshold_cannot_produce_a_verdict():
    """max_transcriptome_seed_perfect defaults to None; that is unevaluable, not a pass."""
    descriptor = FilterDescriptor(
        filter_id="max_transcriptome_seed_perfect",
        column="transcriptome_hits_seed_0mm",
        comparator=FilterComparator.LE,
        threshold=None,
        scope=FilterScope(),
        action=FilterAction.WARN,
    )

    with pytest.raises(ValidationError):
        FilterVerdict(descriptor=descriptor, observed=3, evaluation=FilterEvaluation.PASS)
    assert FilterVerdict(descriptor=descriptor, observed=3, evaluation=FilterEvaluation.NOT_EVALUATED).observed == 3


@pytest.mark.unit
def test_a_descriptor_that_cannot_reproduce_its_own_verdict_is_rejected():
    """#103 re-applies the descriptor client-side, so a descriptor that disagrees is a wrong report."""
    descriptor = FilterDescriptor(
        filter_id="max_off_target_count",
        column="off_target_count",
        comparator=FilterComparator.LE,
        threshold=15,
        scope=FilterScope(species=frozenset({"human"}), max_mismatches=2),
        action=FilterAction.FAIL,
    )

    assert FilterVerdict(descriptor=descriptor, observed=3, evaluation=FilterEvaluation.PASS).observed == 3
    with pytest.raises(ValidationError, match="3 le 15"):
        FilterVerdict(descriptor=descriptor, observed=3, evaluation=FilterEvaluation.FAIL)
    with pytest.raises(ValidationError):
        FilterVerdict(descriptor=descriptor, observed=40, evaluation=FilterEvaluation.PASS)


@pytest.mark.unit
def test_an_off_filter_reports_not_evaluated_and_nothing_else():
    """D9's asymmetry gate switched off: a word in the column, no verdict, and the value still visible."""
    off_descriptor = FilterDescriptor(
        filter_id="min_asymmetry_score",
        column="asymmetry_score",
        comparator=FilterComparator.GE,
        threshold=0.65,
        scope=FilterScope(),
        action=FilterAction.OFF,
    )

    retained = FilterVerdict(descriptor=off_descriptor, observed=0.41, evaluation=FilterEvaluation.NOT_EVALUATED)
    assert retained.evaluation.value == "not_evaluated"
    assert retained.observed == 0.41

    for evaluation in (FilterEvaluation.PASS, FilterEvaluation.FAIL, FilterEvaluation.UNKNOWN):
        with pytest.raises(ValidationError):
            FilterVerdict(descriptor=off_descriptor, observed=0.41, evaluation=evaluation)


@pytest.mark.unit
def test_a_verdict_does_not_decide_what_an_unknown_costs():
    """That needs EvidenceRequirements, so no property here folds action and evaluation into a rejection."""
    assert not hasattr(FilterVerdict, "rejects")
    assert not hasattr(FilterVerdict, "is_rejected")


@pytest.mark.unit
def test_target_intent_forces_the_caller_to_declare_selectivity():
    """Pan-isoform and isoform-selective disagree on the transcripts neither set names (#101 AC 1)."""
    with pytest.raises(ValidationError):
        TargetIntent()  # type: ignore[call-arg]

    assert {value.value for value in TargetSelectivity} == {"pan_isoform", "isoform_selective"}
    assert _intent().selectivity is TargetSelectivity.PAN_ISOFORM
    assert (
        TargetIntent(selectivity=TargetSelectivity.ISOFORM_SELECTIVE).selectivity is TargetSelectivity.ISOFORM_SELECTIVE
    )


@pytest.mark.unit
def test_target_intent_defaults_to_human_targets_and_no_screening_species():
    """The documented default is human alone as a target, with screening left undeclared."""
    intent = _intent()

    assert intent.target_species == frozenset({"human"})
    assert intent.offtarget_screen_species == frozenset()
    assert intent.required_transcript_ids == frozenset()
    assert intent.excluded_transcript_ids == frozenset()


@pytest.mark.unit
def test_target_and_screening_species_are_independent_sets():
    """A species may be a target, a screening subject, both, or neither."""
    intent = _intent(
        target_species=frozenset({"human", "cynomolgus"}),
        offtarget_screen_species=frozenset({"human", "mouse", "rat"}),
    )

    assert intent.target_species - intent.offtarget_screen_species == {"cynomolgus"}
    assert intent.offtarget_screen_species - intent.target_species == {"mouse", "rat"}
    assert intent.target_species & intent.offtarget_screen_species == {"human"}


@pytest.mark.unit
def test_coverage_denominator_is_separate_from_what_enumeration_saw():
    """Coverage must not silently shrink to the isoforms that survived design filtering."""
    intent = _intent(
        annotation_universe=frozenset({"ENST1", "ENST2", "ENST3", "ENST4"}),
        coverage_denominator=frozenset({"ENST1", "ENST2", "ENST3"}),
        enumeration_inputs=frozenset({"ENST1", "ENST2"}),
        annotation_provenance="ensembl_114",
    )

    assert intent.coverage_denominator < intent.annotation_universe
    assert intent.enumeration_was_filtered is True
    assert (
        _intent(
            coverage_denominator=frozenset({"ENST1"}), enumeration_inputs=frozenset({"ENST1"})
        ).enumeration_was_filtered
        is False
    )


@pytest.mark.unit
def test_an_unfrozen_denominator_makes_filtering_unknown_not_false():
    """With no denominator there is nothing to prove enumeration was complete against."""
    assert _intent().enumeration_was_filtered is None
    assert _intent(coverage_denominator=frozenset({"ENST1"})).enumeration_was_filtered is None


@pytest.mark.unit
def test_orthology_provenance_separates_an_explicit_mapping_from_symbol_equality():
    """Symbol matching is a heuristic; required conservation must be able to say which it rests on."""
    mapped = OrthologyAssertion(species="mouse", gene_id="ENSMUSG1", source=OrthologyEvidenceSource.EXPLICIT_MAPPING)
    guessed = OrthologyAssertion(
        species="mouse", gene_symbol="Tp53", source=OrthologyEvidenceSource.SYMBOL_EQUALITY, confidence=0.5
    )

    assert mapped.is_validated_orthology
    assert not guessed.is_validated_orthology
    with pytest.raises(ValidationError):
        OrthologyAssertion(species="mouse", source=OrthologyEvidenceSource.SYMBOL_EQUALITY)

    intent = _intent(target_species=frozenset({"human", "mouse"}), orthology_evidence=(mapped, guessed))
    assert len(intent.orthology_evidence) == 2


@pytest.mark.unit
def test_target_intent_is_immutable():
    """Intent is resolved once by its owner and read everywhere else."""
    intent = _intent()
    with pytest.raises(ValidationError):
        intent.target_species = frozenset({"mouse"})


@pytest.mark.unit
def test_evidence_requirements_forces_the_caller_to_choose_the_unknown_action():
    """What an UNKNOWN does to a candidate is a policy decision, so there is no default."""
    with pytest.raises(ValidationError):
        EvidenceRequirements()  # type: ignore[call-arg]

    requirements = EvidenceRequirements(unknown_evidence_action=UnknownEvidenceAction.WARN)

    assert requirements.unknown_evidence_action is UnknownEvidenceAction.WARN
    assert requirements.channel_requirements == ()


@pytest.mark.unit
def test_requiredness_is_declared_per_channel_and_species():
    """A channel required in one species and exploratory in another is one run, not two axes."""
    requirements = EvidenceRequirements(
        channel_requirements=(
            ChannelRequirement(
                channel=ScreeningChannel.TRANSCRIPTOME, species="human", requiredness=Requiredness.REQUIRED
            ),
            ChannelRequirement(
                channel=ScreeningChannel.TRANSCRIPTOME, species="mouse", requiredness=Requiredness.EXPLORATORY
            ),
            ChannelRequirement(
                channel=ScreeningChannel.MIRNA_SEED, species="human", requiredness=Requiredness.EXPLORATORY
            ),
        ),
        unknown_evidence_action=UnknownEvidenceAction.FAIL,
    )

    assert requirements.requiredness_of(ScreeningChannel.TRANSCRIPTOME, "human") is Requiredness.REQUIRED
    assert requirements.requiredness_of(ScreeningChannel.TRANSCRIPTOME, "mouse") is Requiredness.EXPLORATORY
    assert requirements.required_pairs == frozenset({("transcriptome", "human")})
    # An undeclared pair is undeclared, not exploratory: what it means is the caller's decision.
    assert requirements.requiredness_of(ScreeningChannel.TRANSCRIPT_SEED, "human") is None


@pytest.mark.unit
def test_adding_an_exploratory_species_cannot_change_a_required_pair():
    """#101's guarantee, held by enumerating pairs rather than intersecting two species sets."""
    human_required = ChannelRequirement(
        channel=ScreeningChannel.TRANSCRIPTOME, species="human", requiredness=Requiredness.REQUIRED
    )
    before = EvidenceRequirements(
        channel_requirements=(human_required,), unknown_evidence_action=UnknownEvidenceAction.FAIL
    )
    after = EvidenceRequirements(
        channel_requirements=(
            human_required,
            ChannelRequirement(
                channel=ScreeningChannel.TRANSCRIPTOME, species="rat", requiredness=Requiredness.EXPLORATORY
            ),
        ),
        unknown_evidence_action=UnknownEvidenceAction.FAIL,
    )

    assert before.required_pairs == after.required_pairs


@pytest.mark.unit
def test_one_requiredness_position_per_pair():
    """Two positions on one pair would restore the ambiguity this type exists to remove."""
    with pytest.raises(ValidationError):
        EvidenceRequirements(
            channel_requirements=(
                ChannelRequirement(
                    channel=ScreeningChannel.MIRNA_SEED, species="human", requiredness=Requiredness.REQUIRED
                ),
                ChannelRequirement(
                    channel=ScreeningChannel.MIRNA_SEED, species="human", requiredness=Requiredness.EXPLORATORY
                ),
            ),
            unknown_evidence_action=UnknownEvidenceAction.WARN,
        )


@pytest.mark.unit
def test_plan_and_evidence_entries_carry_no_requiredness_of_their_own():
    """EvidenceRequirements is the single authority; a restatement here could contradict it invisibly."""
    assert "required" not in ScreeningPlanEntry.model_fields
    assert "required" not in ScreeningEvidenceEntry.model_fields


@pytest.mark.unit
def test_screening_channels_name_the_three_liability_sources():
    """Transcript-seed sites are a separate channel from known-miRNA resemblance."""
    assert {channel.value for channel in ScreeningChannel} == {
        "transcriptome",
        "mirna_seed",
        "transcript_seed",
    }


@pytest.mark.unit
def test_evidence_status_distinguishes_all_four_outcomes():
    """Never-asked and asked-and-failed must not collapse into one shape."""
    assert {status.value for status in EvidenceStatus} == {
        "complete",
        "failed",
        "not_requested",
        "censored",
    }


@pytest.mark.unit
def test_observed_counts_start_unknown_not_zero_at_every_unit():
    """An unobserved count is None; a fabricated zero is the failure mode this prevents."""
    counts = ObservedCounts()

    for unit in (counts.sites, counts.distinct_transcripts, counts.distinct_genes, counts.unresolved_gene_sites):
        assert unit.value is None
        assert unit.is_lower_bound is False
        assert unit.truncated is False
        assert unit.cap is None


@pytest.mark.unit
def test_site_transcript_and_gene_counts_have_independent_cap_semantics():
    """A per-query alignment cap and a gene-level cap truncate different things (#101)."""
    counts = ObservedCounts(
        sites=ObservedCount(value=500, is_lower_bound=True, cap=500, truncated=True),
        distinct_transcripts=ObservedCount(value=112),
        distinct_genes=ObservedCount(value=61, is_lower_bound=True),
        unresolved_gene_sites=ObservedCount(value=9),
    )

    assert counts.sites.truncated and counts.sites.cap == 500
    assert counts.distinct_transcripts.truncated is False
    # Genes are a lower bound because nine sites had no resolvable gene, not because of a cap.
    assert counts.distinct_genes.is_lower_bound and counts.distinct_genes.cap is None
    assert counts.unresolved_gene_sites.value == 9


@pytest.mark.unit
def test_the_hit_count_matrix_declares_the_scope_it_is_complete_over():
    """Inside the declared scope an absent cell is a real zero; outside it the matrix says nothing."""
    matrix = HitCountMatrix(
        scope=FilterScope(species=frozenset({"human", "mouse"}), max_mismatches=3),
        cells=(
            HitCountCell(species="human", mismatches=0, hit_class=HitClass.ON_TARGET.value, sites=34),
            HitCountCell(species="human", mismatches=2, hit_class=HitClass.OFF_TARGET.value, sites=18),
            HitCountCell(species="mouse", mismatches=2, hit_class=HitClass.ORTHOLOG.value, sites=4),
        ),
    )

    assert matrix.scope.max_mismatches == 3
    assert sum(cell.sites for cell in matrix.cells if cell.hit_class == HitClass.OFF_TARGET.value) == 18
    assert HitCountMatrix.model_validate_json(matrix.model_dump_json()) == matrix


@pytest.mark.unit
def test_the_hit_count_matrix_rejects_a_duplicated_cell():
    """Two cells at one coordinate would make a total depend on iteration order."""
    cell = HitCountCell(species="human", mismatches=1, hit_class=HitClass.OFF_TARGET.value, sites=3)
    with pytest.raises(ValidationError):
        HitCountMatrix(scope=FilterScope(), cells=(cell, cell))


@pytest.mark.unit
def test_plan_and_evidence_entries_share_one_identity_including_the_guide_set():
    """Reconciliation is the consumer's job, but both sides must key on the same tuple."""
    planned = ScreeningPlanEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        reference_id="ensembl_114_cdna",
        guide_set_digest=DIGEST,
        search_settings={"bwa_k": 12, "max_hits": 500},
    )
    observed = ScreeningEvidenceEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        reference_id="ensembl_114_cdna",
        guide_set_digest=DIGEST,
        status=EvidenceStatus.FAILED,
        detail="bwa-mem2 index build ran out of memory",
    )

    assert planned.key == observed.key
    assert planned.search_settings["bwa_k"] == 12


@pytest.mark.unit
def test_two_guide_sets_against_one_reference_are_different_evidence():
    """Without the digest in the key, one screen's counts can be attributed to the other's guides."""
    common = {"channel": ScreeningChannel.TRANSCRIPTOME, "species": "human", "reference_id": "ensembl_114_cdna"}
    first = ScreeningPlanEntry(guide_set_digest="sha256:aaa", **common)
    second = ScreeningPlanEntry(guide_set_digest="sha256:bbb", **common)

    assert first.key != second.key
    with pytest.raises(ValidationError):
        ScreeningPlanEntry(**common)  # type: ignore[arg-type]


@pytest.mark.unit
def test_a_submitted_guide_set_differing_from_the_plan_is_representable():
    """The point of carrying both digests is that the mismatch can be detected rather than assumed away."""
    entry = ScreeningEvidenceEntry(
        channel=ScreeningChannel.MIRNA_SEED,
        species="human",
        guide_set_digest=DIGEST,
        submitted_guide_digest="sha256:different",
        status=EvidenceStatus.COMPLETE,
    )

    assert entry.submitted_guide_digest != entry.guide_set_digest


@pytest.mark.unit
def test_plan_and_evidence_carry_a_schema_version_and_round_trip():
    """The payload is versioned so a report and its run cannot be silently mismatched."""
    plan = ScreeningPlan(
        entries=(
            ScreeningPlanEntry(channel=ScreeningChannel.TRANSCRIPTOME, species="human", guide_set_digest=DIGEST),
            ScreeningPlanEntry(channel=ScreeningChannel.MIRNA_SEED, species="human", guide_set_digest=DIGEST),
        )
    )
    evidence = ScreeningEvidence(
        entries=(
            ScreeningEvidenceEntry(
                channel=ScreeningChannel.MIRNA_SEED,
                species="human",
                guide_set_digest=DIGEST,
                status=EvidenceStatus.COMPLETE,
                counts=ObservedCounts(sites=ObservedCount(value=0)),
            ),
        )
    )

    assert plan.schema_version == EVIDENCE_SCHEMA_VERSION
    assert evidence.schema_version == EVIDENCE_SCHEMA_VERSION
    assert ScreeningPlan.model_validate_json(plan.model_dump_json()) == plan
    assert ScreeningEvidence.model_validate_json(evidence.model_dump_json()) == evidence
    # A completed channel that found nothing is a real zero, distinct from an unobserved None.
    assert evidence.entries[0].counts.sites.value == 0


@pytest.mark.unit
def test_a_truncated_count_must_declare_itself_a_lower_bound():
    """Discarded members are uncounted members, so truncated without is_lower_bound is a total that isn't."""
    with pytest.raises(ValidationError):
        ObservedCount(value=500, cap=500, truncated=True)

    honest = ObservedCount(value=500, cap=500, truncated=True, is_lower_bound=True)
    assert honest.is_lower_bound and honest.truncated


@pytest.mark.unit
def test_a_count_cannot_exceed_the_cap_it_declares():
    """A value above its own cap means the declared cap was not the cap in force."""
    with pytest.raises(ValidationError):
        ObservedCount(value=501, cap=500)

    assert ObservedCount(value=500, cap=500).value == 500


@pytest.mark.unit
def test_failed_and_censored_evidence_must_say_why():
    """A rejection with no reason stops reading as a rejection; both statuses exist to be acted on."""
    common = {
        "channel": ScreeningChannel.TRANSCRIPTOME,
        "species": "mouse",
        "guide_set_digest": DIGEST,
    }
    for status in (EvidenceStatus.FAILED, EvidenceStatus.CENSORED):
        with pytest.raises(ValidationError):
            ScreeningEvidenceEntry(status=status, **common)  # type: ignore[arg-type]
        with pytest.raises(ValidationError):
            ScreeningEvidenceEntry(status=status, detail="   ", **common)  # type: ignore[arg-type]

    assert ScreeningEvidenceEntry(status=EvidenceStatus.FAILED, detail="index build OOM", **common).detail  # type: ignore[arg-type]
    # COMPLETE and NOT_REQUESTED need no reason: nothing went wrong to explain.
    assert ScreeningEvidenceEntry(status=EvidenceStatus.COMPLETE, **common).detail is None  # type: ignore[arg-type]


@pytest.mark.unit
def test_never_requested_evidence_cannot_carry_an_observed_count():
    """A zero on a search nobody ran is the fabricated zero the whole counts model exists to prevent."""
    common = {
        "channel": ScreeningChannel.TRANSCRIPTOME,
        "species": "rat",
        "guide_set_digest": DIGEST,
        "status": EvidenceStatus.NOT_REQUESTED,
    }
    with pytest.raises(ValidationError):
        ScreeningEvidenceEntry(counts=ObservedCounts(sites=ObservedCount(value=0)), **common)  # type: ignore[arg-type]
    with pytest.raises(ValidationError):
        ScreeningEvidenceEntry(counts=ObservedCounts(distinct_genes=ObservedCount(value=3)), **common)  # type: ignore[arg-type]

    unrequested = ScreeningEvidenceEntry(**common)  # type: ignore[arg-type]
    assert unrequested.counts.observed_units == ()


@pytest.mark.unit
def test_an_empty_plan_is_representable():
    """A design-only run plans nothing, and that is not an error to construct."""
    assert ScreeningPlan().entries == ()
    assert ScreeningEvidence().entries == ()


@pytest.mark.unit
def test_hit_class_counts_field_names_match_the_hit_class_values():
    """accumulate_hit_class does setattr(counts, hit_class.value); mypy cannot see that coupling.

    #101 is chartered to add HitClass members. Without this assertion the first new member fails as
    an AttributeError inside #100's aggregation loop instead of here.
    """
    counter_fields = {field.name for field in dataclasses.fields(HitClassCounts)}
    missing = {member.value for member in HitClass} - counter_fields

    assert not missing, f"HitClassCounts has no counter field for HitClass member(s): {sorted(missing)}"
