"""The run resolves an intent, freezes a denominator, and keeps three channels apart (#101, S4).

This is the wiring slice's test file. The pure seams it drives are tested on their own terms in
``tests/unit/test_target_intent.py`` and ``tests/unit/test_transcript_seed_channel.py``; what is
asserted here is that ``workflow.py`` calls them at the right moment, with the right inputs, and
publishes the answers where a consumer can check them.

Four things #101 says must not happen, each pinned below:

1. **The coverage denominator must not shrink.** It was assigned from ``step1``'s *filtered* return
   value, so an isoform dropped for want of a sequence or by ORF validation left the denominator and
   coverage rose for a guide that covered nothing more. It is now a field on one frozen instance,
   reached through a read-only property, so a later assignment is an ``AttributeError``.
2. **Coverage must not be the enumeration map.** ``_store_guide_to_transcripts`` records which
   transcripts a guide was *enumerated on*, which differs from "which transcripts does this guide
   silence" on exactly the transcripts enumeration never saw.
3. **Missing annotation must read as unknown, never as a number.** An unavailable sequence makes
   coverage undecidable, and an undecidable gate withholds rather than rejects.
4. **A seed site must not be counted as an alignment.** The seed rows arrive on their own key with a
   positive channel discriminator; nothing about them reaches ``off_target_count``.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

import asyncio
import inspect
from pathlib import Path
from typing import Any

import pytest
from rich.progress import Progress

from sirnaforge.cli import workflow as workflow_command
from sirnaforge.config.run_policy import FILTER_SPECS, EntryPoint, resolve_run_policy
from sirnaforge.core.hit_classification import HitClass, HitClassCounts, HitClassification
from sirnaforge.core.screening_evidence import collect_evidence
from sirnaforge.core.seed_geometry import reverse_complement
from sirnaforge.core.target_intent import CoverageStatus, IntentVerdict, resolve_target_intent
from sirnaforge.core.transcript_seed import SeedScanScope, SiteRegion, TranscriptSeedScanResult
from sirnaforge.data.base import DatabaseType, GeneInfo, TranscriptInfo
from sirnaforge.data.gene_search import GeneSearchResult
from sirnaforge.models.evidence import EvidenceStatus, ObservedCounts
from sirnaforge.models.policy import (
    FilterEvaluation,
    FilterStage,
    ScreeningChannel,
    TargetSelectivity,
)
from sirnaforge.models.sirna import OffTargetFilterCriteria, SiRNACandidate
from sirnaforge.workflow import (
    POST_SCREEN_FILTER_CHANNELS,
    SiRNAWorkflow,
    WorkflowConfig,
)

#: A 21-mer guide, and the site a transcript must carry to be covered by it: the guide's full-length
#: reverse complement, because a guide is antisense to its target. Searching for the guide itself
#: would find passenger-orientation sites at the same rate and mean nothing.
GUIDE = "ACGTACGTACGTACGTACGTA"
SITE = reverse_complement(GUIDE)

REQUIRED = "ENST00000000001"
SIBLING = "ENST00000000002"
NON_CODING = "ENST00000000003"


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    stated: dict[str, object] | None = None,
    selectivity: TargetSelectivity | None = None,
    required: list[str] | None = None,
    excluded: list[str] | None = None,
    transcript_seed_scope: str | None = None,
) -> SiRNAWorkflow:
    """A workflow whose policy is resolved through the one resolver, as every real entry point does."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        stated=stated,
        query_species="human",
        screen_species=["human"],
    )
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        resolved_policy=policy,
        screen_species=["human"],
        query_species="human",
        selectivity=selectivity,
        required_transcripts=required,
        excluded_transcripts=excluded,
        transcript_seed_scope=transcript_seed_scope,
    )
    return SiRNAWorkflow(config)


def _candidate(candidate_id: str = "probe") -> SiRNACandidate:
    """A minimal candidate carrying enough component scores to be re-scored after screening."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id=REQUIRED,
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=reverse_complement(GUIDE),
        length=len(GUIDE),
        gc_content=47.6,
        asymmetry_score=0.5,
        design_score=50.0,
        component_scores={"target_accessibility": 0.5, "asymmetry": 0.5, "gc_content": 0.5},
        screen_query_id=candidate_id,
    )


def _transcript(transcript_id: str, *, biotype: str | None, sequence: str | None) -> TranscriptInfo:
    return TranscriptInfo(
        transcript_id=transcript_id,
        transcript_type=biotype,
        gene_id="ENSG00000141510",
        gene_name="TP53",
        sequence=sequence,
        length=len(sequence) if sequence else None,
        database=DatabaseType.ENSEMBL,
    )


def _padded(site: str | None = None) -> str:
    """A synthetic transcript of fixed composition, optionally carrying one guide site."""
    filler = "AAACCCGGGTTT" * 4
    return f"{filler}{site or ''}{filler}"


# ──────────────────────────────────────────────────────
#  The frozen denominator
# ──────────────────────────────────────────────────────


@pytest.mark.unit
def test_the_coverage_denominator_is_frozen_before_orf_and_design_filtering(tmp_path: Path) -> None:
    """The universe intent is built from is the RETRIEVED set, not the set step1 returns.

    ``step1_retrieve_transcripts`` returns only transcripts that are protein-coding *and* carry a
    sequence. Building intent from its return value is what let an isoform leave the denominator
    silently: a protein-coding sibling with no fetched sequence simply stopped being something the
    guide had to cover, and coverage rose without the guide changing. So the universe has to come
    from the recorded pre-filter set, and it has to still contain the non-coding isoform -- which the
    biotype rule then keeps out of the *denominator* on the record rather than by never having seen it.
    """
    workflow = _workflow(tmp_path, "frozen")
    workflow._retrieved_annotation_universe = {
        REQUIRED: "protein_coding",
        SIBLING: "protein_coding",
        NON_CODING: "retained_intron",
    }
    workflow._target_sequences = {REQUIRED: _padded(SITE), SIBLING: _padded(SITE)}

    # What step1 actually returns: the one transcript that is coding AND has a sequence.
    survivors = [_transcript(REQUIRED, biotype="protein_coding", sequence=_padded(SITE))]
    intent = workflow._freeze_target_intent(survivors)

    assert NON_CODING in intent.annotation_universe
    assert intent.annotation_universe == frozenset({REQUIRED, SIBLING, NON_CODING})
    # The sequence-less protein-coding sibling stays in the denominator; the non-coding isoform is
    # excluded by biotype, which is a decision on the record rather than an omission.
    assert intent.coverage_denominator == frozenset({REQUIRED, SIBLING})
    assert workflow._protein_coding_transcript_ids == frozenset({REQUIRED, SIBLING})
    assert workflow._protein_coding_transcript_count == 2


@pytest.mark.unit
def test_step1_records_the_unfiltered_universe_without_changing_its_return_value(tmp_path: Path) -> None:
    """Step 1 keeps its signature and its filter, and records the wider set as a side effect.

    The signature is load-bearing: a great many tests patch this method, so widening its return type
    would have been a repo-wide change. The recording is what makes the freeze possible at all --
    ``run_complete_workflow`` never sees the unfiltered set otherwise.
    """
    workflow = _workflow(tmp_path, "step1")
    retrieved = [
        _transcript(REQUIRED, biotype="protein_coding", sequence=_padded(SITE)),
        _transcript(SIBLING, biotype="protein_coding", sequence=None),
        _transcript(NON_CODING, biotype="retained_intron", sequence=_padded()),
    ]

    class _Searcher:
        def query_species(self, _database: object) -> str:
            return "human"

        async def search_gene(self, *_args: object, **_kwargs: object) -> Any:
            return GeneSearchResult(
                query="TP53",
                database=DatabaseType.ENSEMBL,
                gene_info=GeneInfo(gene_id="ENSG00000141510", gene_name="TP53", database=DatabaseType.ENSEMBL),
                transcripts=retrieved,
            )

    workflow.gene_searcher = _Searcher()  # type: ignore[assignment]
    workflow._annotation_client = None

    with Progress() as progress:
        returned = asyncio.run(workflow.step1_retrieve_transcripts(progress))

    assert [t.transcript_id for t in returned] == [REQUIRED]
    assert set(workflow._retrieved_annotation_universe) == {REQUIRED, SIBLING, NON_CODING}
    assert workflow._retrieved_annotation_universe[NON_CODING] == "retained_intron"
    # A transcript with no sequence is in the universe but not in the sequence map, which is what
    # makes its coverage unknown rather than zero.
    assert set(workflow._target_sequences) == {REQUIRED, NON_CODING}


@pytest.mark.unit
def test_nothing_downstream_can_shrink_the_frozen_denominator(tmp_path: Path) -> None:
    """The denominator is a read-only property, and freezing twice is refused.

    This is the guarantee, not a convention: ``_protein_coding_transcript_ids`` used to be a mutable
    set attribute that any later step could reassign, and the fix is that assigning it now raises.
    Freezing twice is refused for the same reason -- two resolutions would produce two denominators,
    and the point of freezing is that there is exactly one.
    """
    workflow = _workflow(tmp_path, "readonly")
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    workflow._freeze_target_intent([])

    with pytest.raises(AttributeError):
        workflow._protein_coding_transcript_ids = {REQUIRED}  # type: ignore[misc]
    with pytest.raises(AttributeError):
        workflow._protein_coding_transcript_count = 1  # type: ignore[misc]
    with pytest.raises(RuntimeError, match="already frozen"):
        workflow._freeze_target_intent([])

    assert workflow._protein_coding_transcript_ids == frozenset({REQUIRED, SIBLING})


@pytest.mark.unit
def test_declaring_required_transcripts_makes_the_run_isoform_selective(tmp_path: Path) -> None:
    """Naming targets IS selectivity, so it is inferred rather than defaulted away.

    Left to default to pan-isoform, a run that named two required transcripts would have had its
    coverage measured against every isoform of the gene -- and a perfectly selective guide would
    then score as poor coverage and be rejected by the very floor meant to protect it.
    """
    workflow = _workflow(tmp_path, "selective", required=[f"{REQUIRED}.9"])
    assert workflow.config.selectivity is TargetSelectivity.ISOFORM_SELECTIVE
    # Version-stripped through the same function the universe is keyed with, or the declared id would
    # match nothing and be recorded as unresolved.
    assert workflow.config.required_transcripts == (REQUIRED,)

    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    intent = workflow._freeze_target_intent([])
    assert intent.coverage_denominator == frozenset({REQUIRED})
    assert intent.unresolved_target_ids == frozenset()


@pytest.mark.unit
def test_a_declared_transcript_the_annotation_never_offered_is_recorded_as_unresolved(tmp_path: Path) -> None:
    """A typo'd id constrains nothing, and absence from every set is indistinguishable from compliance."""
    workflow = _workflow(tmp_path, "unresolved", required=[REQUIRED], excluded=["ENST00000000000"])
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding"}
    intent = workflow._freeze_target_intent([])

    assert intent.unresolved_target_ids == frozenset({"ENST00000000000"})
    assert workflow._target_intent_payload()["unresolved_target_ids"] == ["ENST00000000000"]


# ──────────────────────────────────────────────────────
#  Coverage: the numerator, and the two ways it is unknown
# ──────────────────────────────────────────────────────


@pytest.mark.unit
def test_coverage_counts_a_required_transcript_enumeration_never_saw(tmp_path: Path) -> None:
    """The numerator is complementarity against sequence, not the enumeration map.

    ``_guide_to_transcripts`` records that this guide was enumerated on ``REQUIRED`` alone. The
    sibling carries the same site verbatim, so the guide covers it -- and the old numerator, an
    intersection with the enumeration map, would have reported half the coverage that exists. The two
    answers differ on exactly the transcripts enumeration never saw, which are the ones an
    isoform-selective design cares most about.
    """
    workflow = _workflow(tmp_path, "numerator")
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    workflow._target_sequences = {REQUIRED: _padded(SITE), SIBLING: _padded(SITE)}
    workflow._freeze_target_intent([])
    workflow._guide_to_transcripts = {GUIDE: frozenset({REQUIRED})}

    observation = workflow._observe_candidate_coverage(_candidate())
    assert observation is not None
    assert observation.covered == frozenset({REQUIRED, SIBLING})
    assert observation.fraction == pytest.approx(1.0)
    assert observation.status is CoverageStatus.KNOWN

    candidate = _candidate()
    workflow._score_candidate_post_screen(candidate, _empty_hit_counts(), frozenset())
    assert candidate.isoform_coverage == pytest.approx(1.0)
    assert candidate.intent_coverage_status == CoverageStatus.KNOWN.value


@pytest.mark.unit
def test_a_guide_present_in_only_one_denominator_member_reports_a_partial_fraction(tmp_path: Path) -> None:
    """The complementarity test is a real test, so a guide that covers one of two reports 0.5."""
    workflow = _workflow(tmp_path, "partial")
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    workflow._target_sequences = {REQUIRED: _padded(SITE), SIBLING: _padded()}
    workflow._freeze_target_intent([])

    observation = workflow._observe_candidate_coverage(_candidate())
    assert observation is not None
    assert observation.covered == frozenset({REQUIRED})
    assert observation.fraction == pytest.approx(0.5)


@pytest.mark.unit
def test_missing_sequence_makes_coverage_unknown_not_zero(tmp_path: Path) -> None:
    """An unavailable sequence withholds the fraction, and the gate reports UNKNOWN, never FAIL.

    ``min_isoform_coverage`` is a greater-or-equal floor, so a numerator missing an untested
    transcript is a lower bound -- and a low lower bound can only wrongly REJECT. Publishing 0.0 (or
    even the honest 0.5 over a denominator of 2) would turn "we could not fetch that sequence" into
    evidence against the guide. Unknown never rejects; it withholds from a qualified shortlist, which
    is a statement about the run rather than about the candidate.
    """
    workflow = _workflow(tmp_path, "unknown", stated={"min_isoform_coverage": 0.75})
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    # SIBLING is in the denominator and this run holds no sequence for it.
    workflow._target_sequences = {REQUIRED: _padded(SITE)}
    workflow._freeze_target_intent([])

    observation = workflow._observe_candidate_coverage(_candidate())
    assert observation is not None
    assert observation.unavailable == frozenset({SIBLING})
    assert observation.fraction is None
    assert observation.status is CoverageStatus.UNKNOWN_MISSING_SEQUENCE

    candidate = _candidate()
    workflow._score_candidate_post_screen(candidate, _empty_hit_counts(), frozenset())

    # The gate first, because the gate's verdict is the consequence that matters: a fabricated
    # fraction of 0.5 against a floor of 0.75 would be a decided FAIL, and it would reject.
    rejected = workflow._apply_isoform_coverage_gate(candidate)
    assert candidate.filter_verdicts["min_isoform_coverage"] == FilterEvaluation.UNKNOWN.value
    assert rejected is False
    assert candidate.filter_observed["min_isoform_coverage"] is None

    assert candidate.isoform_coverage is None
    assert candidate.intent_coverage_status == CoverageStatus.UNKNOWN_MISSING_SEQUENCE.value

    payload = workflow._target_intent_payload()
    assert payload["coverage_status"] == CoverageStatus.UNKNOWN_MISSING_SEQUENCE.value
    assert payload["coverage_missing_sequence"] == [SIBLING]
    # No fraction anywhere in the published block: the three unknown spellings are words.
    assert "isoform_coverage" not in payload


@pytest.mark.unit
def test_an_empty_denominator_is_unknown_rather_than_complete_coverage(tmp_path: Path) -> None:
    """An isoform-selective intent that required nothing names nothing to cover.

    The reading an empty set invites is "covered everything", which is exactly backwards. It is
    reported as its own unknown spelling so the fix -- declare a target set -- is readable from the
    summary.
    """
    workflow = _workflow(tmp_path, "empty", selectivity=TargetSelectivity.ISOFORM_SELECTIVE)
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding"}
    workflow._target_sequences = {REQUIRED: _padded(SITE)}
    workflow._freeze_target_intent([])

    observation = workflow._observe_candidate_coverage(_candidate())
    assert observation is not None
    assert observation.fraction is None
    assert observation.status is CoverageStatus.UNKNOWN_NO_DENOMINATOR
    assert workflow._target_intent_payload()["coverage_status"] == CoverageStatus.UNKNOWN_NO_DENOMINATOR.value


@pytest.mark.unit
def test_the_published_block_derives_enumeration_filtering_from_what_was_observed(tmp_path: Path) -> None:
    """``enumeration_was_filtered`` reads the recorded map, so a later stage's loss cannot be hidden.

    Intent is frozen before ORF validation, so the enumeration inputs it holds are what design was
    *about* to be offered. ORF validation and a failed design batch both remove transcripts after
    that, so a published flag derived from the frozen field alone could claim an enumeration was
    complete when it was not.
    """
    workflow = _workflow(tmp_path, "enumeration")
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    workflow._target_sequences = {REQUIRED: _padded(SITE), SIBLING: _padded(SITE)}
    workflow._freeze_target_intent(
        [
            _transcript(REQUIRED, biotype="protein_coding", sequence=_padded(SITE)),
            _transcript(SIBLING, biotype="protein_coding", sequence=_padded(SITE)),
        ]
    )
    assert workflow._target_intent is not None
    assert workflow._target_intent.enumeration_was_filtered is False

    # Design only ever enumerated on one of the two, which the frozen field cannot know.
    workflow._guide_to_transcripts = {GUIDE: frozenset({REQUIRED})}
    payload = workflow._target_intent_payload()
    assert payload["enumeration_transcripts_observed"] == [REQUIRED]
    assert payload["enumeration_was_filtered"] is True


# ──────────────────────────────────────────────────────
#  Intent evaluation after classification
# ──────────────────────────────────────────────────────


def _empty_hit_counts() -> Any:
    return HitClassCounts()


def _on_target(transcript_id: str) -> tuple[dict[str, Any], HitClassification]:
    """One on-target alignment row and the classifier verdict that goes with it."""
    row: dict[str, Any] = {"qname": "probe", "rname": transcript_id, "species": "human", "nm": 0}
    return row, HitClassification(hit_class=HitClass.ON_TARGET, matched_symbol="TP53", symbol_lookup_missing=False)


@pytest.mark.unit
def test_an_excluded_isoform_hit_stays_on_target_and_is_written_beside_the_class(tmp_path: Path) -> None:
    """``hit_class`` is untouched and ``intent_verdict`` is written next to it.

    On-target is a statement about the gene and it is true: this hit IS on the target gene. It is also
    forbidden, because the caller excluded that isoform. Folding the second into the first would have
    needed a fifth hit class that lies about the gene, so both are published and the pair
    ``on_target`` + ``excluded_isoform`` reads without contradiction.
    """
    workflow = _workflow(tmp_path, "excluded", excluded=[SIBLING])
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    workflow._freeze_target_intent([])

    row, classification = _on_target(f"{SIBLING}.4")
    verdict = workflow._record_intent_verdict(row, classification, "human", "human")

    assert verdict is IntentVerdict.EXCLUDED_ISOFORM
    assert row["intent_verdict"] == "excluded_isoform"
    assert row["intent_reason"]
    assert "hit_class" not in row  # written by annotate_hit_row, and not by this method


@pytest.mark.unit
def test_an_unknown_intent_verdict_feeds_neither_counter(tmp_path: Path) -> None:
    """A transcript the intent has never heard of is a gap, not a finding.

    ``classify_hit`` reaches ON_TARGET through a gene-id or symbol match and consults no transcript
    identity at all, so an on-target row can legitimately name a transcript no declared set contains.
    Counting that as an unintended isoform would manufacture rejections out of a thin index, which is
    the same defect as a missing UTR annotation reading as a clean UTR screen.
    """
    workflow = _workflow(tmp_path, "unknownintent", required=[REQUIRED])
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding"}
    workflow._freeze_target_intent([])

    row, classification = _on_target("ENST00000999999")
    verdict = workflow._record_intent_verdict(row, classification, "human", "human")
    assert verdict is IntentVerdict.UNKNOWN
    assert row["intent_verdict"] == "unknown"

    known_sibling_row, sibling_classification = _on_target(SIBLING)
    workflow._retrieved_annotation_universe[SIBLING] = "protein_coding"
    # The intent is already frozen, so the sibling is still outside its universe: same verdict.
    assert workflow._record_intent_verdict(known_sibling_row, sibling_classification, "human", "human") is (
        IntentVerdict.UNKNOWN
    )


@pytest.mark.unit
def test_no_intent_resolved_writes_no_verdict_at_all(tmp_path: Path) -> None:
    """A direct caller that froze no intent has no verdict to report, so the column is left alone."""
    workflow = _workflow(tmp_path, "nointent")
    row, classification = _on_target(REQUIRED)
    assert workflow._record_intent_verdict(row, classification, "human", "human") is None
    assert "intent_verdict" not in row


@pytest.mark.unit
def test_an_excluded_isoform_hit_rejects_through_the_declared_gate(tmp_path: Path) -> None:
    """The rejection is expressed by a gate reading a counter, never by rewriting the taxonomy.

    ``max_excluded_isoform_hits`` ships fail at 0 because it restates something the caller declared,
    and #101's rule is that a required limit uses fail. The counter is what the gate reads, so the
    hit table and the candidate row agree by construction.
    """
    workflow = _workflow(tmp_path, "gate", excluded=[SIBLING])
    workflow._retrieved_annotation_universe = {REQUIRED: "protein_coding", SIBLING: "protein_coding"}
    workflow._freeze_target_intent([])
    workflow._screened_species_scope = frozenset({"human"})

    candidate = _candidate()
    should_fail, status = workflow._check_offtarget_filters(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        OffTargetFilterCriteria(),
        candidate,
        excluded_isoform_hits=1,
        complete_pairs=frozenset({(ScreeningChannel.TRANSCRIPTOME.value, "human")}),
    )
    assert should_fail is True
    assert status is SiRNACandidate.FilterStatus.EXCLUDED_ISOFORM
    assert candidate.filter_verdicts["max_excluded_isoform_hits"] == FilterEvaluation.FAIL.value
    # And with nothing excluded the same gate passes trivially, so a default run rejects nothing new.
    clean = _candidate("clean")
    should_fail, _ = workflow._check_offtarget_filters(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        OffTargetFilterCriteria(),
        clean,
        excluded_isoform_hits=0,
        complete_pairs=frozenset({(ScreeningChannel.TRANSCRIPTOME.value, "human")}),
    )
    assert should_fail is False
    assert clean.filter_verdicts["max_excluded_isoform_hits"] == FilterEvaluation.PASS.value


@pytest.mark.unit
def test_the_unintended_isoform_ceiling_is_off_until_a_caller_states_one(tmp_path: Path) -> None:
    """Same-gene-but-not-required is a consequence of asking for selectivity, not a prohibition.

    The excluded set is where a prohibition is stated. So this ceiling ships with no threshold and
    resolves to ``off``; naming one opts in, and ``_filter_action_for`` promotes it to fail.
    """
    workflow = _workflow(tmp_path, "unintended")
    workflow._screened_species_scope = frozenset({"human"})
    complete = frozenset({(ScreeningChannel.TRANSCRIPTOME.value, "human")})

    reported = _candidate("reported")
    workflow._check_offtarget_filters(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        OffTargetFilterCriteria(),
        reported,
        unintended_isoform_hits=3,
        complete_pairs=complete,
    )
    assert reported.filter_verdicts["max_unintended_isoform_hits"] == FilterEvaluation.NOT_EVALUATED.value
    assert reported.filter_observed["max_unintended_isoform_hits"] == 3

    gated = _candidate("gated")
    should_fail, status = workflow._check_offtarget_filters(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        OffTargetFilterCriteria(max_unintended_isoform_hits=1),
        gated,
        unintended_isoform_hits=3,
        complete_pairs=complete,
    )
    assert should_fail is True
    assert status is SiRNACandidate.FilterStatus.UNINTENDED_ISOFORM


@pytest.mark.unit
def test_every_declared_post_screen_gate_names_its_channel_and_the_intent_gates_are_query_scoped() -> None:
    """The wiring trap, pinned: an unscoped intent gate would report UNKNOWN on a partial run.

    ``_gate_evidence_pairs`` falls back to every screened species for a gate with an unrestricted
    scope, so an intent gate that did not declare ``scope_query_species`` would demand transcriptome
    evidence for every screened species and go undecidable whenever a secondary species failed --
    while counting only query-species rows, since ``evaluate_intent`` returns UNKNOWN for an
    on-target row in any other species.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human")
    post_screen = {
        resolved.filter_id for resolved in policy.filters if resolved.descriptor.stage is FilterStage.POST_SCREEN
    }
    # Only min_isoform_coverage reads annotation rather than a screening channel.
    assert post_screen - set(POST_SCREEN_FILTER_CHANNELS) == {"min_isoform_coverage"}

    intent_gates = ("max_excluded_isoform_hits", "max_unintended_isoform_hits")
    seed_gates = ("max_transcript_seed_sites", "max_transcript_seed_transcripts", "max_transcript_seed_genes")
    for filter_id in intent_gates:
        assert POST_SCREEN_FILTER_CHANNELS[filter_id] == frozenset({ScreeningChannel.TRANSCRIPTOME})
    for filter_id in seed_gates:
        assert POST_SCREEN_FILTER_CHANNELS[filter_id] == frozenset({ScreeningChannel.TRANSCRIPT_SEED})

    by_id = {spec.filter_id: spec for spec in FILTER_SPECS}
    for filter_id in intent_gates + seed_gates:
        assert by_id[filter_id].scope_query_species is True, filter_id


# ──────────────────────────────────────────────────────
#  The transcript-seed channel: a third channel, kept apart
# ──────────────────────────────────────────────────────


def _write_seed_tsv(path: Path, *, guide_id: str, transcript_id: str, gene_id: str | None) -> None:
    """One published seed-site row, in ``TranscriptSeedSiteSchema``'s shape."""
    header = (
        "guide_id\tguide_sequence\tqueried_strand\tspecies\ttranscript_id\ttranscript_version\tgene_id\t"
        "site_start\tsite_end\tanchor_position\tsite_class\tsite_strand\tregion\tcoordinate_system\t"
        "annotation_provenance"
    )
    row = (
        f"{guide_id}\t{GUIDE}\tguide\thuman\t{transcript_id}\t3\t{gene_id or ''}\t"
        "101\t107\t107\t7mer-m8\ttranscript_sense\tfull_cdna\ttranscript_cdna_1based\tensembl_cdna:test"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{header}\n{row}\n")


@pytest.mark.unit
def test_a_seed_site_row_is_not_counted_as_an_alignment(tmp_path: Path) -> None:
    """Seed rows land on their own key, with a positive channel discriminator.

    A seed site has no ``nm``, no ``cigar`` and no ``mapq``. Ingested through ``_ingest_row`` it would
    arrive as an alignment whose missing ``nm`` defaults to 0 -- a perfect-match liability -- and
    would inflate ``off_target_count``, the mismatch strata and the four-way class tally, all of which
    read ``results[qname]["hits"]`` positionally. So the parser puts them in a SIBLING key and stamps
    ``channel`` on each row, rather than leaving a third channel's rows to be recognised by the
    absence of ``mirna_id``/``database`` the way miRNA rows already are.
    """
    output_dir = tmp_path / "results"
    aggregated = output_dir / "aggregated"
    aggregated.mkdir(parents=True)
    (aggregated / "combined_offtargets.tsv").write_text(
        f"qname\trname\tspecies\tnm\tseed_mismatches\nprobe\t{SIBLING}\thuman\t0\t0\n"
    )
    _write_seed_tsv(
        output_dir / "human_transcript_seed_sites.tsv", guide_id="probe", transcript_id=SIBLING, gene_id="ENSG1"
    )

    workflow = _workflow(tmp_path, "parse", transcript_seed_scope="full_cdna")
    parsed = asyncio.run(workflow._parse_nextflow_results(output_dir))

    # Exactly one alignment, and it is the alignment row.
    assert parsed["results"]["probe"]["off_target_count"] == 1
    assert len(parsed["results"]["probe"]["hits"]) == 1
    assert all("site_class" not in hit for hit in parsed["results"]["probe"]["hits"])

    seed_rows = parsed["transcript_seed"]["probe"]
    assert len(seed_rows) == 1
    assert seed_rows[0]["channel"] == ScreeningChannel.TRANSCRIPT_SEED.value
    assert seed_rows[0]["site_class"] == "7mer-m8"
    assert parsed["transcript_seed_files"]


@pytest.mark.unit
def test_the_three_seed_units_are_counted_and_reported_separately(tmp_path: Path) -> None:
    """Sites, distinct transcripts and distinct genes are three answers, and an unresolved gene is kept.

    One transcript carries many sites and one gene many transcripts, so no two of the three are
    derivable from each other. A site whose transcript resolved to no gene is counted in
    ``transcript_seed_unresolved_gene_sites`` rather than dropped, which is what makes the gene count
    a declared lower bound instead of a quiet undercount.
    """
    workflow = _workflow(tmp_path, "units", transcript_seed_scope="full_cdna")
    workflow._transcript_seed_observed = True
    workflow._transcript_seed_rows = {
        "probe": [
            {"species": "human", "transcript_id": SIBLING, "gene_id": "ENSG1"},
            {"species": "human", "transcript_id": SIBLING, "gene_id": "ENSG1"},
            {"species": "human", "transcript_id": NON_CODING, "gene_id": None},
            {"species": "mouse", "transcript_id": "ENSMUST1", "gene_id": "ENSMUSG1"},
        ]
    }

    candidate = _candidate()
    counts = workflow._transcript_seed_gate_counts(candidate)

    assert counts.sites == 3  # query species only
    assert counts.distinct_transcripts == 2
    assert counts.distinct_genes == 1  # the unresolved one is not invented into a gene
    assert candidate.transcript_seed_unresolved_gene_sites == 1
    assert candidate.transcript_seed_sites_total == 4
    assert candidate.transcript_seed_sites_query == 3
    assert candidate.transcript_seed_genes_query == 1


@pytest.mark.unit
def test_the_seed_ceilings_are_unknown_rather_than_zero_when_no_scan_answered(tmp_path: Path) -> None:
    """A channel that never ran must not pass a ceiling on a zero nothing measured.

    The channel is opt-in, so the honest cell for a run that did not request it is an unobserved
    value, which ``evaluate_gate`` turns into UNKNOWN the moment a ceiling exists. A zero would be a
    claim about a scan that never happened -- the same defect as a UTR-only request with no UTR
    annotation publishing a clean UTR screen.
    """
    workflow = _workflow(tmp_path, "noscan")
    workflow._screened_species_scope = frozenset({"human"})
    candidate = _candidate()

    counts = workflow._transcript_seed_gate_counts(candidate)
    assert counts.sites is None
    assert candidate.transcript_seed_sites_query is None

    workflow._check_offtarget_filters(
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        OffTargetFilterCriteria(max_transcript_seed_sites=5),
        candidate,
        transcript_seed_sites=counts.sites,
        transcript_seed_transcripts=counts.distinct_transcripts,
        transcript_seed_genes=counts.distinct_genes,
        complete_pairs=frozenset({(ScreeningChannel.TRANSCRIPTOME.value, "human")}),
    )
    assert candidate.filter_verdicts["max_transcript_seed_sites"] == FilterEvaluation.UNKNOWN.value
    assert candidate.filter_observed["max_transcript_seed_sites"] is None


@pytest.mark.unit
def test_a_refused_scope_leaves_the_seed_columns_unobserved(tmp_path: Path) -> None:
    """A FAILED scan is "we asked and could not answer", which is not a measurement of zero.

    Nothing in this repository carries UTR intervals in transcript coordinates, so a UTR-only request
    is refused by the scanner. Its consequence here has to be an unobserved column, or a
    transcript-seed ceiling would pass against a clean UTR screen that never ran.
    """
    workflow = _workflow(tmp_path, "refused", transcript_seed_scope="utr3")
    workflow._transcript_seed_results["human"] = TranscriptSeedScanResult(
        sites=(),
        counts=ObservedCounts(),
        status=EvidenceStatus.FAILED,
        detail="3' UTR scope requested and no UTR annotation was available for human",
        scope=SeedScanScope(
            species="human",
            region=SiteRegion.UTR3,
            classes=frozenset(),
            max_sites_per_guide=None,
            reference_id=None,
        ),
        submitted_guides=1,
        processed_guides=0,
    )
    # Deliberately NOT setting _transcript_seed_observed: a FAILED scan never sets it.
    candidate = _candidate()
    assert workflow._transcript_seed_gate_counts(candidate).sites is None
    assert candidate.transcript_seed_sites_query is None


@pytest.mark.unit
def test_the_channel_is_not_requested_unless_a_scope_is_declared(tmp_path: Path) -> None:
    """Opt-in, and the plan says so: a run that did not ask records no transcript-seed unit at all.

    That is what keeps ``NOT_REQUESTED`` distinguishable from ``FAILED``, and what makes a default
    run's plan identical to the one it produced before the channel existed.
    """
    default_run = _workflow(tmp_path, "default")
    assert default_run.config.transcript_seed_scope is None
    assert default_run._requested_transcript_seed_units() == ()
    assert default_run._run_transcript_seed_scan([_candidate()]) == {"status": "not_requested"}

    opted_in = _workflow(tmp_path, "optedin", transcript_seed_scope="full_cdna")
    assert opted_in._requested_transcript_seed_units() == (("human", None),)


@pytest.mark.unit
def test_the_scan_reuses_the_materialised_reference_and_never_fetches_one(tmp_path: Path) -> None:
    """No reference means no scan and no counts -- never a download, and never a fabricated zero.

    The scan reads the cDNA copy screening already materialised, exactly as repeat detection does.
    With nothing materialised there is nothing to read, and the honest outcome is a skip whose
    columns stay unobserved.
    """
    workflow = _workflow(tmp_path, "noref", transcript_seed_scope="full_cdna")
    summary = workflow._run_transcript_seed_scan([_candidate()])
    assert summary == {"status": "skipped", "reason": "reference_unavailable", "species": "human"}
    assert workflow._transcript_seed_observed is False


@pytest.mark.unit
def test_the_scan_publishes_sites_evidence_and_a_summary(tmp_path: Path) -> None:
    """End to end over a synthetic reference: sites, a schema-shaped TSV, and a #100 envelope.

    The envelope matters as much as the sites: it is what puts the (channel, species) pair into
    ``completed_pairs``, and therefore what lets a stated ceiling reach a decided verdict instead of
    UNKNOWN.
    """
    reference = tmp_path / "human_cdna.fa"
    reference.write_text(f">{SIBLING}.3 gene:ENSG1\n{_padded(SITE)}\n")

    workflow = _workflow(tmp_path, "scan", transcript_seed_scope="full_cdna")
    workflow._species_cdna_fasta["human"] = reference
    workflow._guide_set_digest = "sha256:" + "0" * 64
    summary = workflow._run_transcript_seed_scan([_candidate()])

    assert summary["status"] == EvidenceStatus.COMPLETE.value
    assert summary["sites"] >= 1
    assert workflow._transcript_seed_observed is True

    sites_file = Path(str(summary["detail_files"]["sites"]))
    assert sites_file.exists()
    header = sites_file.read_text().splitlines()[0].split("\t")
    assert header[0] == "guide_id"
    assert "channel" not in header  # a discriminator on the row, not a property of a site

    envelopes = collect_evidence(sites_file.parent)
    assert [envelope.entry.channel for envelope in envelopes] == [ScreeningChannel.TRANSCRIPT_SEED]
    assert envelopes[0].entry.species == "human"


@pytest.mark.unit
def test_the_summary_block_carries_the_frozen_lists_and_the_match_semantics(tmp_path: Path) -> None:
    """#101 requires the numerator/denominator ids, the unit, the provenance and the match be reported.

    A free-text match description would drift from the code that implements it, so the published
    value is the enum member's own word.
    """
    workflow = _workflow(tmp_path, "summary", required=[REQUIRED], excluded=[SIBLING])
    workflow._retrieved_annotation_universe = {
        REQUIRED: "protein_coding",
        SIBLING: "protein_coding",
        NON_CODING: "retained_intron",
    }
    workflow._target_sequences = {REQUIRED: _padded(SITE)}
    workflow._freeze_target_intent([])

    payload = workflow._target_intent_payload()
    assert payload["resolved"] is True
    assert payload["selectivity"] == TargetSelectivity.ISOFORM_SELECTIVE.value
    assert payload["coverage_denominator"] == [REQUIRED]
    assert payload["excluded_transcript_ids"] == [SIBLING]
    assert payload["annotation_universe"] == sorted([REQUIRED, SIBLING, NON_CODING])
    assert payload["coverage_unit"] == "transcript"
    assert payload["coverage_match"] == "exact_full_site"
    assert payload["coverage_status"] == CoverageStatus.KNOWN.value
    assert payload["annotation_provenance"] == "ensembl_gene_query"


@pytest.mark.unit
def test_an_unresolved_intent_publishes_that_it_is_unresolved_rather_than_an_empty_intent(tmp_path: Path) -> None:
    """A run that never froze an intent says so, instead of publishing empty lists that read as facts."""
    workflow = _workflow(tmp_path, "notresolved")
    assert workflow._target_intent_payload() == {"resolved": False}


@pytest.mark.unit
def test_the_intent_takes_no_variant_argument_and_variant_resolution_is_untouched() -> None:
    """Two axes, two owners: which transcripts versus which alleles, and neither reads the other.

    #101 is explicit that target intent must compose with ``VariantWorkflowConfig`` rather than open a
    second SNP configuration system. Two objects that cannot see each other cannot disagree, so this
    pins the absence of the coupling rather than its behaviour.
    """
    parameters = set(inspect.signature(resolve_target_intent).parameters)
    assert not {name for name in parameters if "variant" in name or "snp" in name or "allele" in name}

    freeze = set(inspect.signature(SiRNAWorkflow._freeze_target_intent).parameters)
    assert freeze == {"self", "transcripts"}

    # And the variant step's own signature is unchanged by this slice.
    variant_step = inspect.signature(SiRNAWorkflow.resolve_variants_step)
    assert set(variant_step.parameters) == {"self", "progress"}


@pytest.mark.unit
def test_the_cli_exposes_the_intent_and_the_seed_controls() -> None:
    """The declarations reach the user through the same plumbing every other setting uses.

    Named options rather than a policy file only: an exclusion is a safety statement and has to be
    typeable on the command line beside the run it constrains.
    """
    parameters = inspect.signature(workflow_command).parameters
    for name in (
        "required_transcripts",
        "excluded_transcripts",
        "selectivity",
        "transcript_seed_scope",
        "max_transcript_seed_sites",
        "max_transcript_seed_transcripts",
        "max_transcript_seed_genes",
    ):
        assert name in parameters, name
