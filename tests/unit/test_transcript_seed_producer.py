"""Deterministic tests for the Nextflow-called transcript-seed producer."""

import json
from types import SimpleNamespace

from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.core.screening_evidence import collect_evidence
from sirnaforge.core.target_intent import TargetIntentRequest
from sirnaforge.core.transcript_seed import run_transcript_seed_analysis
from sirnaforge.models.evidence import EvidenceStatus
from sirnaforge.models.policy import Requiredness, ScreeningChannel, TargetSelectivity
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.models.transcript_seed import (
    TRANSCRIPT_SEED_SITE_IDENTITY,
    TRANSCRIPT_SEED_SUMMARY_SCHEMA_VERSION,
    TranscriptSequenceScope,
)
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig


def test_producer_scans_full_reference_and_publishes_counts_and_reference_checksum(tmp_path) -> None:
    """A full cDNA scan emits real site observations and a complete evidence envelope."""
    guides = tmp_path / "guides.fa"
    guides.write_text(">guide-1\nAACCTGTAAG\n")
    reference = tmp_path / "human.fa"
    reference.write_text(">tx-1\nCCCTACAGGTACCC\n>tx-2\nCCCCCCCCCCCC\n")

    result = run_transcript_seed_analysis(guides, reference, species="human", output_dir=tmp_path / "out")

    assert result["scope"] == TranscriptSequenceScope.FULL_CDNA.value
    assert result["reference_sha256"]
    assert result["site_count"] == 1
    assert result["transcript_count"] == 1
    assert result["transcripts_scanned"] == 2
    assert result["per_guide"]["guide-1"]["site_count"] == 1
    assert result["per_guide"]["guide-1"]["gene_count_is_lower_bound"] is True
    assert result["per_guide"]["guide-1"]["unresolved_gene_site_count"] == 1
    envelopes = collect_evidence(tmp_path / "out")
    assert len(envelopes) == 1
    assert envelopes[0].entry.channel is ScreeningChannel.TRANSCRIPT_SEED
    assert envelopes[0].entry.status is EvidenceStatus.COMPLETE
    assert envelopes[0].entry.counts.sites.value == 1
    assert envelopes[0].entry.counts.distinct_genes.value == 0
    assert envelopes[0].entry.counts.distinct_genes.is_lower_bound is True
    assert envelopes[0].entry.counts.unresolved_gene_sites.value == 1
    summary = json.loads((tmp_path / "out" / "transcript_seed_human_summary.json").read_text())
    assert summary["reference_sha256"] == result["reference_sha256"]
    assert summary["schema_version"] == TRANSCRIPT_SEED_SUMMARY_SCHEMA_VERSION
    assert summary["site_identity"] == TRANSCRIPT_SEED_SITE_IDENTITY


def test_producer_resolves_seed_site_genes_from_ensembl_cdna_headers(tmp_path) -> None:
    """The standard workflow FASTA carries gene identity even without a separate UTR manifest."""
    guides = tmp_path / "guides.fa"
    guides.write_text(">guide-1\nAACCTGTAAG\n")
    reference = tmp_path / "human.fa"
    reference.write_text(
        ">ENST00000000001.7 cdna chromosome:GRCh38:1:1:14:1 gene:ENSG00000000001.4 gene_biotype:protein_coding\n"
        "CCCTACAGGTACCC\n"
    )

    result = run_transcript_seed_analysis(guides, reference, species="human", output_dir=tmp_path / "out")

    site = result["sites"][0]
    assert site["gene_id"] == "ENSG00000000001"
    assert result["gene_count"] == 1
    assert result["gene_count_is_lower_bound"] is False
    assert result["per_guide"]["guide-1"]["gene_count"] == 1
    assert result["per_guide"]["guide-1"]["gene_count_is_lower_bound"] is False
    envelope = collect_evidence(tmp_path / "out")[0]
    assert envelope.entry.counts.distinct_genes.value == 1
    assert envelope.entry.counts.distinct_genes.is_lower_bound is False


def test_per_species_seed_limits_fail_on_lower_bound_exceedance_but_withhold_unproven_pass(tmp_path) -> None:
    """Seed site, transcript, and gene ceilings use per-species counts and lower-bound logic."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        stated={
            "species_liability_limits": {
                "human": {
                    "transcript_seed_sites": 10,
                    "transcript_seed_transcripts": 5,
                    "transcript_seed_genes": 3,
                    "action": "fail",
                }
            }
        },
    )
    workflow = SiRNAWorkflow(
        WorkflowConfig(
            output_dir=tmp_path / "workflow",
            gene_query="TP53",
            resolved_policy=policy,
            transcript_seed_species=["human"],
            transcript_seed_required_species=["human"],
        )
    )
    candidate = SiRNACandidate(
        id="candidate",
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence="AACCTGTAAGAACCTGTAAGA",
        passenger_sequence="TTTTTTTTTTTTTTTTTTTTT",
        gc_content=50.0,
        length=21,
        asymmetry_score=0.7,
        screen_query_id="guide",
    )
    result_dir = tmp_path / "results" / "transcript_seed"
    result_dir.mkdir(parents=True)
    summary_path = result_dir / "transcript_seed_human_summary.json"

    def publish(gene_count: int) -> None:
        summary_path.write_text(
            json.dumps(
                {
                    "species": "human",
                    "schema_version": TRANSCRIPT_SEED_SUMMARY_SCHEMA_VERSION,
                    "site_identity": TRANSCRIPT_SEED_SITE_IDENTITY,
                    "per_guide": {
                        "guide": {
                            "status": "complete",
                            "site_count": 12,
                            "transcript_count": 6,
                            "gene_count": gene_count,
                            "gene_count_is_lower_bound": True,
                        }
                    },
                }
            )
        )

    workflow._completed_evidence_pairs = frozenset({("transcript_seed", "human")})
    workflow._candidate_id_to_representative = {}
    publish(0)
    workflow._collect_transcript_seed_observations(result_dir.parent, [candidate])
    assert candidate.filter_verdicts["species_liability_human_transcript_seed_sites"] == "fail"
    assert candidate.filter_verdicts["species_liability_human_transcript_seed_transcripts"] == "fail"
    assert candidate.filter_verdicts["species_liability_human_transcript_seed_genes"] == "unknown"
    assert candidate.filter_observed["species_liability_human_transcript_seed_genes"] is None

    over_limit = candidate.model_copy(deep=True)
    over_limit.filter_verdicts.clear()
    over_limit.filter_observed.clear()
    publish(4)
    workflow._collect_transcript_seed_observations(result_dir.parent, [over_limit])
    assert over_limit.filter_verdicts["species_liability_human_transcript_seed_genes"] == "fail"


def test_producer_marks_missing_three_prime_utr_annotation_unavailable(tmp_path) -> None:
    """A plain cDNA FASTA cannot publish an empty 3′UTR result as complete."""
    guides = tmp_path / "guides.fa"
    guides.write_text(">guide-1\nAACCTGTAAG\n")
    reference = tmp_path / "human.fa"
    reference.write_text(">tx-1\nCCCTACAGGTACCC\n")

    result = run_transcript_seed_analysis(
        guides,
        reference,
        species="human",
        output_dir=tmp_path / "out",
        scope=TranscriptSequenceScope.THREE_PRIME_UTR,
    )

    assert result["site_count"] is None
    assert result["transcripts_unavailable"] == ["tx-1"]
    assert result["per_guide"]["guide-1"]["status"] == "unknown"
    envelope = collect_evidence(tmp_path / "out")[0]
    assert envelope.entry.status is EvidenceStatus.FAILED
    assert envelope.entry.counts.sites.value is None


def test_workflow_maps_only_reconciled_guide_observations_and_writes_provenance(tmp_path) -> None:
    """Candidate output only inherits a complete observation for its submitted guide and species."""
    guide = "ATGCGATGCGATGCGATGCGC"
    candidate = SiRNACandidate(
        id="candidate-1",
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=guide,
        passenger_sequence=guide[::-1],
        gc_content=57.1,
        length=len(guide),
        asymmetry_score=0.7,
        composite_score=50.0,
        screen_query_id="submitted-guide",
    )
    result_dir = tmp_path / "results" / "transcript_seed"
    result_dir.mkdir(parents=True)
    mouse_summary = {
        "species": "mouse",
        "schema_version": TRANSCRIPT_SEED_SUMMARY_SCHEMA_VERSION,
        "site_identity": TRANSCRIPT_SEED_SITE_IDENTITY,
        "scope": "full_cdna",
        "reference_sha256": "abc123",
        "transcripts_scanned": 900,
        "per_guide": {"submitted-guide": {"status": "complete", "site_count": 4}},
    }
    mouse_summary["schema_version"] = "1"
    mouse_summary.pop("site_identity")
    mouse_summary_path = result_dir / "transcript_seed_mouse_summary.json"
    mouse_summary_path.write_text(json.dumps(mouse_summary))
    workflow = object.__new__(SiRNAWorkflow)
    workflow.config = SimpleNamespace(
        transcript_seed_species=("mouse", "human"),
        transcript_seed_required_species=("mouse",),
        transcript_seed_scope=TranscriptSequenceScope.FULL_CDNA,
        max_transcript_seed_sites=5,
        output_dir=tmp_path / "workflow",
    )
    workflow._completed_evidence_pairs = frozenset({("transcript_seed", "mouse")})
    workflow._candidate_id_to_representative = {"candidate-1": "submitted-guide"}

    legacy = workflow._collect_transcript_seed_observations(tmp_path / "results", [candidate])
    legacy_observation = legacy["per_candidate"]["candidate-1"]
    assert legacy_observation["required_species_site_count"] is None
    assert legacy_observation["by_species"]["mouse"] == {"status": "unknown", "site_count": None}

    mouse_summary["schema_version"] = TRANSCRIPT_SEED_SUMMARY_SCHEMA_VERSION
    mouse_summary["site_identity"] = TRANSCRIPT_SEED_SITE_IDENTITY
    mouse_summary_path.write_text(json.dumps(mouse_summary))
    candidate.filter_verdicts.clear()
    candidate.filter_observed.clear()
    result = workflow._collect_transcript_seed_observations(tmp_path / "results", [candidate])

    assert result["status"] == "partial"
    observation = result["per_candidate"]["candidate-1"]
    assert observation["required_species_site_count"] == 4
    assert observation["by_species"]["mouse"]["site_count"] == 4
    assert observation["by_species"]["human"] == {"status": "unknown", "site_count": None}
    assert candidate.filter_verdicts["max_transcript_seed_sites"] == "pass"
    assert json.loads((tmp_path / "workflow/transcript_seed/guide_observations.json").read_text()) == result

    # Exploratory species are still reported when complete, but their much larger count must not
    # enter the required-species ceiling or change candidate eligibility.
    (result_dir / "transcript_seed_human_summary.json").write_text(
        json.dumps(
            {
                "species": "human",
                "schema_version": TRANSCRIPT_SEED_SUMMARY_SCHEMA_VERSION,
                "site_identity": TRANSCRIPT_SEED_SITE_IDENTITY,
                "scope": "full_cdna",
                "reference_sha256": "def456",
                "transcripts_scanned": 1200,
                "per_guide": {
                    "submitted-guide": {
                        "status": "complete",
                        "site_count": 400,
                        "transcript_count": 200,
                        "gene_count": 60,
                    }
                },
            }
        )
    )
    workflow._completed_evidence_pairs = frozenset({("transcript_seed", "mouse"), ("transcript_seed", "human")})
    exploratory = candidate.model_copy(deep=True)
    exploratory.filter_verdicts.clear()
    exploratory.filter_observed.clear()
    complete = workflow._collect_transcript_seed_observations(tmp_path / "results", [exploratory])
    exploratory_observation = complete["per_candidate"]["candidate-1"]
    assert complete["status"] == "complete"
    assert exploratory_observation["required_species_site_count"] == 4
    assert exploratory_observation["by_species"]["human"]["site_count"] == 400
    assert exploratory.filter_verdicts["max_transcript_seed_sites"] == "pass"

    withheld = candidate.model_copy(deep=True)
    withheld.filter_verdicts.clear()
    withheld.filter_observed.clear()
    workflow._completed_evidence_pairs = frozenset()
    missing_required = workflow._collect_transcript_seed_observations(tmp_path / "results", [withheld])
    assert missing_required["per_candidate"]["candidate-1"]["required_species_site_count"] is None
    assert withheld.filter_verdicts["max_transcript_seed_sites"] == "unknown"


def test_workflow_config_forwards_seed_channel_independently_and_freezes_target_request(tmp_path) -> None:
    """Seed species do not inherit target species, and requiredness reaches the evidence policy."""
    intent = TargetIntentRequest(selectivity=TargetSelectivity.PAN_ISOFORM, target_species=frozenset({"human"}))
    config = WorkflowConfig(
        output_dir=tmp_path / "workflow",
        gene_query="TP53",
        screen_species=["human"],
        mirna_species=["human"],
        transcriptome_fasta=str(tmp_path / "human.fa"),
        target_intent_request=intent,
        transcript_seed_species=["mouse", "human"],
        transcript_seed_required_species=["mouse"],
        transcript_seed_scope=TranscriptSequenceScope.FULL_CDNA,
        max_transcript_seed_sites=10,
    )

    requirements = {
        entry.species: entry.requiredness
        for entry in config.resolved_policy.evidence_requirements.channel_requirements
        if entry.channel is ScreeningChannel.TRANSCRIPT_SEED
    }
    assert config.target_intent_request is intent
    assert requirements == {"mouse": Requiredness.REQUIRED, "human": Requiredness.EXPLORATORY}
    assert config.nextflow_config["transcript_seed_species"] == "mouse,human"
    assert config.nextflow_config["transcript_seed_scope"] == "full_cdna"
    assert config.nextflow_config["transcript_seed_site_identity"] == TRANSCRIPT_SEED_SITE_IDENTITY
    assert config.max_transcript_seed_sites == 10
