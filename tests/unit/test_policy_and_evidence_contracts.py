"""The shared 0.7.1 contracts: run policy vocabulary and screening plan/evidence types.

These modules are types, not behaviour. Four downstream packages (#99, #100, #101, #103) fork off
them, so what is pinned here is the vocabulary and the deliberate absence of resolution: no
defaults that make a policy decision, no producer, no consumer.
"""

import ast
from pathlib import Path

import pytest
from pydantic import ValidationError

from sirnaforge.models import evidence as evidence_module
from sirnaforge.models import policy as policy_module
from sirnaforge.models.evidence import (
    EVIDENCE_SCHEMA_VERSION,
    EvidenceStatus,
    ObservedCounts,
    ScreeningEvidence,
    ScreeningEvidenceEntry,
    ScreeningPlan,
    ScreeningPlanEntry,
)
from sirnaforge.models.policy import (
    EvidenceRequirements,
    FilterAction,
    FilterEvaluation,
    RunMode,
    ScreeningChannel,
    TargetIntent,
)


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


@pytest.mark.unit
def test_contract_modules_do_not_import_workflow():
    """These are data contracts: a dependency on the orchestrator would invert the layering."""
    for module in (policy_module, evidence_module):
        assert not any("workflow" in name for name in _imported_modules(module))


@pytest.mark.unit
def test_evidence_module_depends_only_on_the_policy_vocabulary():
    """Channel names are shared, so evidence imports policy and never the reverse."""
    assert "sirnaforge.models.policy" in _imported_modules(evidence_module)
    assert not any(name.startswith("sirnaforge.models.evidence") for name in _imported_modules(policy_module))


@pytest.mark.unit
def test_run_modes_are_the_three_declared_states():
    """design-only, exploratory and qualified, spelled once."""
    assert {mode.value for mode in RunMode} == {"design_only", "exploratory", "qualified"}


@pytest.mark.unit
def test_action_and_evaluation_are_separate_vocabularies():
    """A filter's action is what it does; its evaluation is what it found."""
    assert {action.value for action in FilterAction} == {"off", "warn", "fail"}
    assert {result.value for result in FilterEvaluation} == {"pass", "fail", "unknown"}
    # "fail" is the only word the two share, and it means different things in each.
    assert FilterAction.FAIL.value == FilterEvaluation.FAIL.value
    assert FilterAction.FAIL is not FilterEvaluation.FAIL


@pytest.mark.unit
def test_unknown_is_a_first_class_evaluation():
    """An undecided filter must be expressible without borrowing pass or fail."""
    assert FilterEvaluation.UNKNOWN not in (FilterEvaluation.PASS, FilterEvaluation.FAIL)


@pytest.mark.unit
def test_target_intent_defaults_to_human_targets_and_no_screening_species():
    """The documented default is human alone as a target, with screening left undeclared."""
    intent = TargetIntent()

    assert intent.target_species == frozenset({"human"})
    assert intent.offtarget_screen_species == frozenset()
    assert intent.required_transcript_ids == frozenset()
    assert intent.excluded_transcript_ids == frozenset()


@pytest.mark.unit
def test_target_and_screening_species_are_independent_sets():
    """A species may be a target, a screening subject, both, or neither."""
    intent = TargetIntent(
        target_species=frozenset({"human", "cynomolgus"}),
        offtarget_screen_species=frozenset({"human", "mouse", "rat"}),
    )

    assert intent.target_species - intent.offtarget_screen_species == {"cynomolgus"}
    assert intent.offtarget_screen_species - intent.target_species == {"mouse", "rat"}
    assert intent.target_species & intent.offtarget_screen_species == {"human"}


@pytest.mark.unit
def test_target_intent_is_immutable():
    """Intent is resolved once by its owner and read everywhere else."""
    intent = TargetIntent()
    with pytest.raises(ValidationError):
        intent.target_species = frozenset({"mouse"})


@pytest.mark.unit
def test_evidence_requirements_forces_the_caller_to_choose_the_unknown_action():
    """What an UNKNOWN does to a candidate is a policy decision, so there is no default."""
    with pytest.raises(ValidationError):
        EvidenceRequirements()

    requirements = EvidenceRequirements(
        required_channels=frozenset({ScreeningChannel.TRANSCRIPTOME}),
        exploratory_channels=frozenset({ScreeningChannel.MIRNA_SEED}),
        required_species=frozenset({"human"}),
        unknown_evidence_action=FilterAction.WARN,
    )

    assert requirements.unknown_evidence_action is FilterAction.WARN
    assert requirements.exploratory_species == frozenset()


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
def test_observed_counts_start_unknown_not_zero():
    """An unobserved count is None; a fabricated zero is the failure mode this prevents."""
    counts = ObservedCounts()

    assert counts.hits_retained is None
    assert counts.is_lower_bound is False
    assert counts.truncated is False
    assert counts.hit_cap is None


@pytest.mark.unit
def test_observed_counts_can_declare_a_lower_bound():
    """A capped search reports what it kept and says the total may be larger."""
    counts = ObservedCounts(hits_retained=500, is_lower_bound=True, hit_cap=500, truncated=True)

    assert counts.hits_retained == 500
    assert counts.is_lower_bound
    assert counts.truncated


@pytest.mark.unit
def test_plan_and_evidence_entries_share_one_identity():
    """Reconciliation is the consumer's job, but both sides must key on the same tuple."""
    planned = ScreeningPlanEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        reference_id="ensembl_114_cdna",
        required=True,
        guide_set_digest="sha256:abc",
        search_settings={"bwa_k": 12, "max_hits": 500},
    )
    observed = ScreeningEvidenceEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        reference_id="ensembl_114_cdna",
        required=True,
        status=EvidenceStatus.FAILED,
        detail="bwa-mem2 index build ran out of memory",
    )

    assert planned.key == observed.key
    assert planned.search_settings["bwa_k"] == 12


@pytest.mark.unit
def test_plan_and_evidence_carry_a_schema_version_and_round_trip():
    """The payload is versioned so a report and its run cannot be silently mismatched."""
    plan = ScreeningPlan(
        entries=(
            ScreeningPlanEntry(channel=ScreeningChannel.TRANSCRIPTOME, species="human", required=True),
            ScreeningPlanEntry(channel=ScreeningChannel.MIRNA_SEED, species="human"),
        )
    )
    evidence = ScreeningEvidence(
        entries=(
            ScreeningEvidenceEntry(
                channel=ScreeningChannel.MIRNA_SEED,
                species="human",
                status=EvidenceStatus.COMPLETE,
                counts=ObservedCounts(hits_retained=0),
            ),
        )
    )

    assert plan.schema_version == EVIDENCE_SCHEMA_VERSION
    assert evidence.schema_version == EVIDENCE_SCHEMA_VERSION
    assert ScreeningPlan.model_validate_json(plan.model_dump_json()) == plan
    assert ScreeningEvidence.model_validate_json(evidence.model_dump_json()) == evidence
    # A completed channel that found nothing is a real zero, distinct from an unobserved None.
    assert evidence.entries[0].counts.hits_retained == 0


@pytest.mark.unit
def test_an_empty_plan_is_representable():
    """A design-only run plans nothing, and that is not an error to construct."""
    assert ScreeningPlan().entries == ()
    assert ScreeningEvidence().entries == ()
