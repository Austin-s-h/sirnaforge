"""workflow.py produces #100's screening evidence, and eligibility reads it (W1).

Before this slice nothing in the tree ever emitted a ``ScreeningEvidence``. ``reference_summary``
carried only transcriptome/screening/scope/screening_plan, the plan itself was recorded *after*
``_resolve_screening_references`` had already dropped every reference it could not fetch, and
``_completed_evidence_pairs`` -- the one set the whole eligibility engine compares against -- was
hand-rolled from the ``_species_with_alignment_evidence`` species list. The consequences pinned
here, each independently:

1. A species the run asked for and never screened produced no record at all: it was absent from the
   plan (built from what survived), absent from the evidence (there was none), and its zero hits
   read as clean.
2. The sequence-only fallback returned ``{"status": "partial", "method": "basic"}`` and nothing
   else, so a run that aligned against no reference at all published the same evidence as one that
   screened cleanly: none, which read as unrestricted.
3. A reference rejected by ``_reject_reference`` appeared only inside
   ``to_metadata()['rejections']`` -- a provenance note, not a completeness fact anything joined on.

All fixtures are synthetic; TP53 is the documented public example gene.
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

from sirnaforge.core.screening_evidence import (
    EvidenceProducer,
    EvidenceSource,
    completed_pairs,
    guide_set_digest,
    write_evidence,
)
from sirnaforge.data.mirna_manager import MiRNADatabaseManager
from sirnaforge.models.evidence import EvidenceStatus, ScreeningEvidenceEntry, ScreeningPlan
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate
from sirnaforge.pipeline import NextflowConfig
from sirnaforge.utils.cli_inputs import resolve_species_inputs
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

GUIDE = "ATGCGATGCGATGCGATGCGC"

#: Offline orthologue evidence, so nothing here touches Ensembl Compara.
_ORTHOLOG_MAPPING_FIXTURE = Path(__file__).resolve().parent / "data" / "ortholog_mapping_synthetic.json"


def _candidate(candidate_id: str = "cand_1") -> SiRNACandidate:
    """A minimal candidate carrying one guide."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE[::-1],
        gc_content=57.1,
        length=len(GUIDE),
        asymmetry_score=0.7,
        composite_score=50.0,
    )


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    species: list[str] | None = None,
    transcriptome_indices: str | None = None,
    mirna_species: list[str] | None = None,
) -> SiRNAWorkflow:
    """A TP53 workflow that resolves orthologues from a file rather than over the network.

    A transcriptome reference is configured unless the test names its own: with none at all the run
    requested no transcriptome screen, and a channel nobody asked for has no plan entry by design.
    """
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        screen_species=species or ["human"],
        design_params=DesignParameters(),
        transcriptome_indices=transcriptome_indices,
        transcriptome_fasta=None if transcriptome_indices else str(_reference_fasta(tmp_path)),
        mirna_species=mirna_species,
        ortholog_mapping_file=_ORTHOLOG_MAPPING_FIXTURE,
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    return workflow


def _reference_fasta(tmp_path: Path) -> Path:
    """A cDNA file standing in for a configured screening reference; nothing here resolves it."""
    path = tmp_path / "reference_cdna.fa"
    if not path.exists():
        path.write_text(">ENST00000000001 cdna gene:ENSG00000000001 gene_symbol:TP53\n" + GUIDE * 3 + "\n")
    return path


def _guides(tmp_path: Path) -> Path:
    """The staged guide FASTA every digest in a run is keyed on."""
    path = tmp_path / "input_candidates.fasta"
    path.write_text(f">cand_1\n{GUIDE}\n")
    return path


def _legacy_results(workflow: SiRNAWorkflow, *, species_screened: list[str]) -> Path:
    """A results directory in the pre-#100 shape: an aggregate summary, and no evidence anywhere."""
    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    (aggregated / "combined_summary.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "species_screened": species_screened,
                "species_analyzed": species_screened,
                "hits_per_species": dict.fromkeys(species_screened, 0),
                "total_candidates": 1,
            }
        )
    )
    return results_dir


def _entries_by_key(payload: dict) -> dict[tuple[str, str], dict]:
    """The published evidence entries, keyed by (channel, species)."""
    return {(entry["channel"], entry["species"]): entry for entry in payload["evidence"]["entries"]}


# ---------------------------------------------------------------------------
# 1. The plan is built from what was REQUESTED, before resolution can drop it
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_plan_holds_a_species_whose_reference_never_resolves(tmp_path: Path) -> None:
    """A plan derived from the resolved set cannot represent a missing reference at all.

    ``ScreeningReferenceSet.plan()`` walks ``self.references``, so recording the plan after
    resolution silently answered a different question: not "what did this run ask for" but "what
    did it manage to fetch". The species with no usable index has to be in the plan, or nothing
    downstream can say it failed.
    """
    workflow = _workflow(tmp_path, "requested_plan", transcriptome_indices="rat:/nonexistent")
    params: dict[str, object] = {"max_hits": 100}

    workflow._record_screening_plan(_guides(tmp_path), params)
    resolved = asyncio.run(workflow._resolve_screening_references(params))

    assert resolved is False, "the only reference requested here cannot resolve"
    assert workflow._screening_references.species == (), "nothing survived resolution"
    planned = {(entry.channel.value, entry.species): entry for entry in workflow._screening_plan.entries}
    assert ("transcriptome", "rat") in planned
    assert planned[("transcriptome", "rat")].reference_id == "/nonexistent", "the request's own identity"
    assert planned[("transcriptome", "rat")].search_settings == {"max_hits": 100}
    # The same plan restated from the resolved set: what W1 replaced, and what it must not lose.
    assert workflow._screening_references.plan(guide_set_digest="x").entries == ()
    assert [entry.species for entry in workflow._screening_references.requested_plan(guide_set_digest="x").entries] == [
        "rat"
    ]


@pytest.mark.unit
def test_the_plan_is_serialized_and_threaded_to_the_pipeline(tmp_path: Path) -> None:
    """The pipeline must reconcile against the same expectation the workflow holds.

    ``--evidence_plan`` is how the plan reaches ``aggregate_results_cli``; without it the pipeline
    can only reconcile against the species that happened to build an index, which is the derivation
    #100 exists to remove.
    """
    workflow = _workflow(tmp_path, "threaded_plan", species=["human"])
    params: dict[str, object] = {}

    guides = _guides(tmp_path)
    workflow._record_screening_plan(guides, params)

    plan_file = Path(str(params["evidence_plan"]))
    assert plan_file == workflow.config.output_dir / "screening_plan.json"
    on_disk = ScreeningPlan.model_validate_json(plan_file.read_text())
    assert on_disk == workflow._screening_plan
    digest = guide_set_digest(guides)
    assert {entry.guide_set_digest for entry in on_disk.entries} == {digest}
    assert workflow._guide_set_digest == digest


@pytest.mark.unit
def test_the_plan_does_not_key_the_shared_nextflow_work_directory(tmp_path: Path) -> None:
    """Threading the plan must not cost every run its cached indices.

    The work-dir cache key is built from the pipeline parameters, and the plan's path lives inside
    this run's own output directory -- so keying on it would give each run a private work dir and
    rebuild multi-gigabyte BWA indices every time. Nextflow stages the plan as a path input, so it
    still re-runs the aggregation task when the plan's content changes.
    """
    workflow = _workflow(tmp_path, "cache_key")
    nf_config = NextflowConfig(profile="local")

    without = workflow._prepare_nextflow_cache(
        nf_config=nf_config, screen_species=["human"], additional_params={}, pipeline_revision="rev"
    )
    with_plan = workflow._prepare_nextflow_cache(
        nf_config=nf_config,
        screen_species=["human"],
        additional_params={"evidence_plan": str(workflow.config.output_dir / "screening_plan.json")},
        pipeline_revision="rev",
    )

    assert with_plan["cache_key"] == without["cache_key"]


@pytest.mark.unit
def test_an_unrequested_channel_is_not_requested_rather_than_failed(tmp_path: Path) -> None:
    """Requestedness is the plan's to state, not a producer's silence.

    ``_mirna_channel_completed`` returns False both for a channel nobody asked for and for one that
    was asked for and published nothing. The plan separates them: no entry at all versus an entry
    with no evidence.
    """
    unrequested = _workflow(tmp_path, "mirna_off", species=["human"])
    requested = _workflow(tmp_path, "mirna_on", species=["human"], mirna_species=["human"])
    for workflow in (unrequested, requested):
        workflow._record_screening_plan(_guides(tmp_path), {})

    mirna_key = (ScreeningChannel.MIRNA_SEED.value, "human")
    assert mirna_key not in {(entry.channel.value, entry.species) for entry in unrequested._screening_plan.entries}
    assert mirna_key in {(entry.channel.value, entry.species) for entry in requested._screening_plan.entries}

    for workflow, expected in ((unrequested, EvidenceStatus.NOT_REQUESTED), (requested, EvidenceStatus.FAILED)):
        results_dir = _legacy_results(workflow, species_screened=["human"])
        reconciliation = workflow._reconcile_screening_evidence(
            results_dir, screened_species=["human"], mirna_screened=False
        )
        entry = next(
            entry for entry in reconciliation.evidence.entries if (entry.channel.value, entry.species) == mirna_key
        )
        assert entry.status is expected
        assert entry.counts.observed_units == (), "an unscreened channel observed nothing, not zero"


# ---------------------------------------------------------------------------
# 2. _process_nextflow_results publishes the reconciliation and feeds eligibility from it
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_species_the_screen_never_covered_is_published_as_failed_evidence(tmp_path: Path) -> None:
    """The published record must name the unscreened species, and eligibility must read it from there.

    ``reference_summary`` had no ``screening_evidence`` key at all, and the pairs the gates compare
    against were built from the aggregate's species list -- so a requested species that produced
    nothing was simply absent from both, and its absence was indistinguishable from a clean zero.
    """
    workflow = _workflow(tmp_path, "failed_species", species=["human", "rat"])
    workflow._record_screening_plan(_guides(tmp_path), {})
    results_dir = _legacy_results(workflow, species_screened=["human"])

    outcome = asyncio.run(workflow._process_nextflow_results([_candidate()], results_dir, {"status": "completed"}))
    published = workflow._summarize_screening_references()["screening_evidence"]
    entries = _entries_by_key(published)

    assert entries[("transcriptome", "rat")]["status"] == EvidenceStatus.FAILED.value
    assert entries[("transcriptome", "rat")]["detail"], "a failure with no reason is not actionable"
    assert entries[("transcriptome", "human")]["status"] == EvidenceStatus.COMPLETE.value
    # The heuristic survives as the fallback, and says so: nothing published a per-unit envelope.
    assert published["sources"]["transcriptome|human|" + str(workflow._guide_set_digest)] == (
        EvidenceSource.LEGACY_SUMMARY.value
    )
    assert workflow._completed_evidence_pairs == frozenset({("transcriptome", "human")})
    # The whole reconciliation, not its ``evidence``: the plan has to be in scope where the join key
    # is projected down to (channel, species) and loses the guide-set digest.
    assert workflow._completed_evidence_pairs == completed_pairs(workflow._screening_evidence)
    # The aggregate itself claims nothing was missing, so the run-level word is still "completed":
    # the evidence is what disagrees with it, per unit. Reconciling the two words is W4's (#100).
    assert outcome["status"] == "completed"


@pytest.mark.unit
def test_the_evidence_the_pipeline_published_wins_over_the_summary_heuristic(tmp_path: Path) -> None:
    """A per-unit envelope is first-hand; the species list in a summary is an inference.

    The envelope is what carries a *failed* unit, so a summary that lists a species as screened must
    not be able to overrule it -- that is exactly the "optional output presence is not completion
    evidence" rule, applied to the summary itself.
    """
    from sirnaforge.core.screening_evidence import EvidenceProducer, write_evidence  # noqa: PLC0415
    from sirnaforge.models.evidence import ScreeningEvidenceEntry  # noqa: PLC0415

    workflow = _workflow(tmp_path, "envelope_wins", species=["human"])
    workflow._record_screening_plan(_guides(tmp_path), {})
    results_dir = _legacy_results(workflow, species_screened=["human"])
    write_evidence(
        results_dir / "transcriptome",
        producer=EvidenceProducer.OFFTARGET_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.TRANSCRIPTOME,
            species="human",
            guide_set_digest=str(workflow._guide_set_digest),
            status=EvidenceStatus.FAILED,
            detail="no usable BWA-MEM2 index at the prefix handed to this task",
        ),
    )

    asyncio.run(workflow._process_nextflow_results([_candidate()], results_dir, {"status": "completed"}))

    entries = _entries_by_key(workflow._summarize_screening_references()["screening_evidence"])
    assert entries[("transcriptome", "human")]["status"] == EvidenceStatus.FAILED.value
    assert workflow._completed_evidence_pairs == frozenset()


# ---------------------------------------------------------------------------
# 3. The paths that align nothing say so, per planned unit
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_basic_fallback_reports_every_planned_unit_as_failed(tmp_path: Path) -> None:
    """The sequence-only fallback published no evidence at all, which read as no restriction.

    It aligns against no reference, so every unit the run planned failed -- and that has to be
    stated per unit, because "partial" is a run-level word that no gate and no eligibility rule
    reads.
    """
    workflow = _workflow(tmp_path, "basic_fallback", species=["human", "rat"], mirna_species=["human"])
    workflow._record_screening_plan(_guides(tmp_path), {})

    outcome = asyncio.run(workflow._basic_offtarget_analysis([_candidate()]))

    entries = _entries_by_key(outcome["screening_evidence"])
    planned = {(entry.channel.value, entry.species) for entry in workflow._screening_plan.entries}
    assert planned == {("transcriptome", "human"), ("transcriptome", "rat"), ("mirna_seed", "human")}
    for key in planned:
        assert entries[key]["status"] == EvidenceStatus.FAILED.value, key
        assert entries[key]["detail"]
    assert workflow._completed_evidence_pairs == frozenset()


@pytest.mark.unit
def test_a_rejected_reference_reconciles_as_failed_under_its_own_plan_key(tmp_path: Path) -> None:
    """A rejection was provenance only; it has to be evidence, keyed identically to the plan entry.

    ``_reject_reference`` records the species in ``_species_screening_shortfalls`` and in
    ``to_metadata()['rejections']``, neither of which anything joins on. The reconciled entry must
    carry the plan entry's own key -- channel, species, reference id and guide-set digest -- or the
    failure describes a unit nobody planned.
    """
    workflow = _workflow(tmp_path, "rejected_reference", species=["rat"], transcriptome_indices="rat:/nonexistent")
    params: dict[str, object] = {}
    workflow._record_screening_plan(_guides(tmp_path), params)
    asyncio.run(workflow._resolve_screening_references(params))

    assert "rat" in workflow._species_screening_shortfalls
    results_dir = workflow.config.output_dir / "off_target" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    reconciliation = workflow._reconcile_screening_evidence(results_dir, screened_species=[], mirna_screened=False)

    plan_entry = next(entry for entry in workflow._screening_plan.entries if entry.species == "rat")
    evidence_entry = next(entry for entry in reconciliation.evidence.entries if entry.species == "rat")
    assert evidence_entry.key == plan_entry.key
    assert evidence_entry.status is EvidenceStatus.FAILED
    assert evidence_entry.detail


# ---------------------------------------------------------------------------
# 4. The miRNA channel is planned and emitted in ONE species vocabulary
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_mirna_channel_is_planned_and_emitted_in_one_species_vocabulary(tmp_path: Path) -> None:
    """The default path published miRNA evidence under a key nothing had planned.

    ``resolve_species_inputs`` answers with miRNA *database codes* -- ``--species human,rhesus`` comes
    back as ``['hsa', 'mml']`` -- and those codes reached ``nextflow_config['mirna_species']`` raw, so
    ``run_mirna_seed_analysis`` wrote ``mirna_seed_hsa_evidence.json`` while the plan and the declared
    ``EvidenceRequirements`` were keyed on ``human``. Nothing joined: every mirna_seed unit reconciled
    FAILED on an ordinary ``sirnaforge workflow TP53 --species human`` run, and the envelope that did
    exist landed in ``unplanned``.

    The two vocabularies are pinned to each other rather than each to a literal, so neither side can
    drift again without this failing.
    """
    resolved = resolve_species_inputs(species="human,rhesus", mirna_db="mirgenedb", mirna_species=None)
    assert resolved.mirna_species == ["hsa", "mml"], "the resolver's own answer is database codes"

    workflow = _workflow(
        tmp_path,
        "mirna_vocabulary",
        species=resolved.screen_species,
        mirna_species=resolved.mirna_species,
    )
    workflow._record_screening_plan(_guides(tmp_path), {})

    # The parameter the .nf module splits into run_mirna_seed_analysis' own per-species loop.
    parameter_species = str(workflow.config.nextflow_config["mirna_species"]).split(",")
    planned_mirna = [
        entry.species for entry in workflow._screening_plan.entries if entry.channel is ScreeningChannel.MIRNA_SEED
    ]
    assert planned_mirna == ["human", "macaque"], "the plan speaks canonical species"
    assert parameter_species == planned_mirna, "one vocabulary, or the join key cannot match"

    # Verified, not assumed: a canonical name and its code resolve to the same source, so speaking the
    # plan's vocabulary costs the screen nothing.
    for canonical, code in zip(planned_mirna, resolved.mirna_species, strict=True):
        by_name = MiRNADatabaseManager.get_source_configuration("mirgenedb", canonical)
        by_code = MiRNADatabaseManager.get_source_configuration("mirgenedb", code)
        assert by_name is not None
        assert by_name == by_code
        assert by_name.cache_key() == by_code.cache_key(), "the same cached database, under either name"

    # And what the emitter writes now reconciles against the plan instead of arriving unplanned.
    results_dir = workflow.config.output_dir / "off_target" / "results"
    for species in parameter_species:
        write_evidence(
            results_dir / "mirna",
            producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
            entry=ScreeningEvidenceEntry(
                channel=ScreeningChannel.MIRNA_SEED,
                species=species,
                guide_set_digest=str(workflow._guide_set_digest),
                status=EvidenceStatus.COMPLETE,
            ),
        )
    reconciliation = workflow._reconcile_screening_evidence(results_dir, screened_species=[], mirna_screened=True)

    mirna_status = {
        entry.species: entry.status
        for entry in reconciliation.evidence.entries
        if entry.channel is ScreeningChannel.MIRNA_SEED
    }
    assert mirna_status == dict.fromkeys(planned_mirna, EvidenceStatus.COMPLETE)
    assert reconciliation.unplanned == (), "no envelope describes a unit nobody planned"
    assert {(ScreeningChannel.MIRNA_SEED.value, species) for species in planned_mirna} <= completed_pairs(
        reconciliation
    ), "and eligibility can join on it"
