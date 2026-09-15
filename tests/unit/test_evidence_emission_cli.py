"""``pipeline/nextflow_cli.py``: the per-unit evidence emitter and the aggregate reconciler (#100).

``offtarget_analysis_cli`` already owns the failed-vs-completed decision for one species -- these
tests pin that it writes a #100 evidence envelope on *both* the missing-index branch and the
aligned branch, not only the one that reaches ``run_bwa_alignment_analysis``. ``aggregate_results_cli``
and ``aggregate_mirna_results_cli`` used to treat "nothing staged to aggregate" as a silent
``final_summary.txt``-only exit; these tests pin that ``evidence.json`` is now written
unconditionally, on every return path, so a pipeline that aborted before a single analysis file
existed still leaves a machine-readable failure record behind instead of vanishing.
"""

import json
from pathlib import Path

import pytest

from sirnaforge.core.off_target import ExecutionOutcome
from sirnaforge.core.screening_evidence import (
    EvidenceProducer,
    EvidenceSource,
    build_plan,
    parse_reconciliation_payload,
    read_evidence,
    write_evidence,
)
from sirnaforge.models.evidence import EvidenceStatus, ObservedCounts, ScreeningEvidenceEntry
from sirnaforge.models.off_target import OffTargetHit
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.pipeline.nextflow_cli import (
    aggregate_mirna_results_cli,
    aggregate_results_cli,
    offtarget_analysis_cli,
)

GUIDE = "ATGCGATGCGATGCGATGCGC"


def _candidates_fasta(tmp_path: Path) -> Path:
    """A real candidates FASTA, so a missing index is the only thing wrong with the call."""
    path = tmp_path / "candidates.fasta"
    path.write_text(f">cand_1\n{GUIDE}\n")
    return path


# ---------------------------------------------------------------------------
# offtarget_analysis_cli -- the missing-index branch writes real evidence
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_missing_index_writes_a_failed_transcriptome_evidence_envelope(tmp_path):
    """A named index prefix that never resolved must reconcile as failed, not as unobserved."""
    staged = tmp_path / "staged"
    result = offtarget_analysis_cli(
        species="human",
        index_prefix="/nonexistent/prefix",
        candidates_file=str(_candidates_fasta(tmp_path)),
        output_dir=str(staged),
    )

    assert result["status"] == "failed"
    envelope = read_evidence(staged / "transcriptome_human_evidence.json")
    assert envelope is not None
    assert envelope.producer is EvidenceProducer.OFFTARGET_ANALYSIS
    assert envelope.entry.channel is ScreeningChannel.TRANSCRIPTOME
    assert envelope.entry.species == "human"
    assert envelope.entry.status is EvidenceStatus.FAILED
    assert envelope.entry.detail
    assert "index" in envelope.entry.detail


@pytest.mark.unit
def test_a_missing_index_records_the_submitted_guide_digest_and_count(tmp_path):
    """The submitted digest is derivable even on the failure branch.

    No alignment ran, but the guide set that was going to be screened is fully known from the
    staged FASTA.
    """
    staged = tmp_path / "staged"
    candidates = _candidates_fasta(tmp_path)
    offtarget_analysis_cli(
        species="human",
        index_prefix="/nonexistent/prefix",
        candidates_file=str(candidates),
        output_dir=str(staged),
    )

    envelope = read_evidence(staged / "transcriptome_human_evidence.json")
    assert envelope is not None
    assert envelope.entry.guide_set_digest == envelope.entry.submitted_guide_digest
    assert envelope.entry.submitted_guides == 1
    assert envelope.entry.processed_guides is None
    assert envelope.entry.counts == ObservedCounts()


# ---------------------------------------------------------------------------
# offtarget_analysis_cli -- the aligned branch delegates to run_bwa_alignment_analysis
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_completed_alignment_also_writes_a_transcriptome_evidence_envelope(tmp_path, monkeypatch):
    """A clean zero-hit screen reconciles as complete, not as absent.

    The aligned branch opts ``run_bwa_alignment_analysis`` into evidence via its own
    ``evidence_dir`` seam.
    """
    monkeypatch.setattr(
        "sirnaforge.core.off_target.BwaAnalyzer.analyze_sequences_with_outcome",
        lambda _self, sequences: (
            [],
            ExecutionOutcome(
                completed=True,
                submitted=len(sequences),
                processed=len(sequences),
                retained_hits=0,
                pre_cap_hits=0,
                cap=None,
                truncated=False,
            ),
        ),
    )
    monkeypatch.setattr("sirnaforge.pipeline.nextflow_cli.validate_index_files", lambda *_a, **_k: True)
    staged = tmp_path / "staged"

    result = offtarget_analysis_cli(
        species="human",
        index_prefix=str(tmp_path / "index" / "human"),
        candidates_file=str(_candidates_fasta(tmp_path)),
        output_dir=str(staged),
    )

    assert result["status"] == "completed"
    envelope = read_evidence(staged / "transcriptome_human_evidence.json")
    assert envelope is not None
    assert envelope.entry.status is EvidenceStatus.COMPLETE
    assert envelope.entry.counts.sites.value == 0
    assert envelope.entry.counts.sites.is_lower_bound is False


# ---------------------------------------------------------------------------
# aggregate_results_cli -- the zero-files branch (the contract's headline red)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_zero_files_branch_still_writes_a_reconciliation_file(tmp_path, monkeypatch):
    """#100 requires ``evidence.json`` too, not only ``final_summary.txt``.

    Every requested species reconciles as failed even though nothing was ever staged.
    """
    monkeypatch.chdir(tmp_path)
    output_dir = tmp_path / "aggregated"

    result = aggregate_results_cli(
        transcriptome_species="human,rat",
        output_dir=str(output_dir),
        analysis_files=[],
        summary_files=[],
    )

    assert result["status"] == "empty"
    reconciliation_file = output_dir / "evidence.json"
    assert reconciliation_file.exists()

    payload = json.loads(reconciliation_file.read_text())
    reconciliation = parse_reconciliation_payload(payload)
    assert reconciliation is not None
    statuses = {entry.species: entry.status for entry in reconciliation.evidence.entries}
    assert statuses == {"human": EvidenceStatus.FAILED, "rat": EvidenceStatus.FAILED}
    for entry in reconciliation.evidence.entries:
        assert entry.detail


@pytest.mark.unit
def test_the_zero_files_branch_return_dict_carries_the_same_evidence(tmp_path, monkeypatch):
    """The in-memory return value must agree with what was written to disk."""
    monkeypatch.chdir(tmp_path)
    output_dir = tmp_path / "aggregated"

    result = aggregate_results_cli(
        transcriptome_species="human,rat",
        output_dir=str(output_dir),
        analysis_files=[],
        summary_files=[],
    )

    assert result["reconciliation_file"] == str(output_dir / "evidence.json")
    entries = result["evidence"]["entries"]
    assert {entry["species"] for entry in entries} == {"human", "rat"}
    assert all(entry["status"] == "failed" for entry in entries)


@pytest.mark.unit
def test_a_real_evidence_plan_is_preferred_over_the_local_fallback(tmp_path, monkeypatch):
    """A real reference identity survives to the output when a real plan is threaded through.

    ``workflow.py``'s plan -- not a locally-derived stand-in -- is what gets reconciled.
    """
    monkeypatch.chdir(tmp_path)
    plan = build_plan(
        guide_set_digest="abc123abc123abcd",
        transcriptome=[("human", "ensembl_human_cdna")],
        mirna_species=(),
    )
    plan_path = tmp_path / "screening_plan.json"
    plan_path.write_text(plan.model_dump_json())
    output_dir = tmp_path / "aggregated"

    result = aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(output_dir),
        analysis_files=[],
        summary_files=[],
        evidence_plan=str(plan_path),
    )

    (entry,) = result["evidence"]["entries"]
    assert entry["reference_id"] == "ensembl_human_cdna"
    assert entry["guide_set_digest"] == "abc123abc123abcd"


# ---------------------------------------------------------------------------
# aggregate_results_cli -- the aggregated branch also reconciles unconditionally
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_aggregated_branch_also_writes_a_reconciliation_file(tmp_path, monkeypatch):
    """Even when real files were staged and aggregated, evidence.json is not optional."""
    monkeypatch.chdir(tmp_path)
    flat = tmp_path / "staged_flat"
    flat.mkdir()
    analysis = flat / "human_analysis.tsv"
    analysis.write_text(OffTargetHit.tsv_header() + "\n")
    output_dir = tmp_path / "aggregated"

    result = aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(output_dir),
        analysis_files=[str(analysis)],
        summary_files=[],
    )

    assert result["status"] == "completed"
    reconciliation_file = output_dir / "evidence.json"
    assert reconciliation_file.exists()
    assert result["reconciliation_file"] == str(reconciliation_file)


@pytest.mark.unit
def test_an_envelope_staged_alongside_analysis_files_is_picked_up_as_complete(tmp_path, monkeypatch):
    """A real per-species envelope staged by OFFTARGET_ANALYSIS reconciles as complete evidence.

    Not as the synthesized-failed fallback the local plan would otherwise produce.
    """
    monkeypatch.chdir(tmp_path)
    flat = tmp_path / "staged_flat"
    flat.mkdir()
    analysis = flat / "human_analysis.tsv"
    analysis.write_text(OffTargetHit.tsv_header() + "\n")
    write_evidence(
        flat,
        producer=EvidenceProducer.OFFTARGET_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.TRANSCRIPTOME,
            species="human",
            guide_set_digest="abc123abc123abcd",
            status=EvidenceStatus.COMPLETE,
        ),
    )
    output_dir = tmp_path / "aggregated"

    result = aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(output_dir),
        analysis_files=[str(analysis)],
        summary_files=[],
    )

    (entry,) = result["evidence"]["entries"]
    assert entry["status"] == "complete"
    payload = json.loads((output_dir / "evidence.json").read_text())
    reconciliation = parse_reconciliation_payload(payload)
    assert reconciliation is not None
    assert reconciliation.sources[("transcriptome", "human", "abc123abc123abcd")] is EvidenceSource.ENVELOPE


# ---------------------------------------------------------------------------
# aggregate_mirna_results_cli -- the same additive treatment
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_aggregate_mirna_results_cli_also_writes_a_reconciliation_file(tmp_path):
    """A caller of this entry point gets the same evidence contract as the transcriptome channel.

    ``aggregate_results.nf`` never actually calls this function, but it should not be left with a
    silently weaker evidence contract than ``aggregate_results_cli``.
    """
    results_dir = tmp_path / "mirna_results"
    results_dir.mkdir()
    output_dir = tmp_path / "aggregated_mirna"

    result = aggregate_mirna_results_cli(
        mirna_db="mirgenedb",
        mirna_species="human,mouse",
        results_dir=str(results_dir),
        output_dir=str(output_dir),
    )

    assert result["status"] == "completed"
    reconciliation_file = output_dir / "evidence.json"
    assert reconciliation_file.exists()
    assert result["reconciliation_file"] == str(reconciliation_file)
    statuses = {entry["species"]: entry["status"] for entry in result["evidence"]["entries"]}
    assert statuses == {"human": "failed", "mouse": "failed"}


# ---------------------------------------------------------------------------
# Compatibility: the six/seven pre-existing return keys must survive untouched.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_aggregate_results_cli_keeps_its_six_existing_return_keys(tmp_path, monkeypatch):
    """``evidence``/``reconciliation_file`` are additive: none of the six original keys move."""
    monkeypatch.chdir(tmp_path)
    result = aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(tmp_path / "aggregated"),
        analysis_files=[],
        summary_files=[],
    )
    for key in ("status", "analysis_files_processed", "summary_files_processed", "species", "output_dir"):
        assert key in result


@pytest.mark.unit
def test_aggregate_mirna_results_cli_keeps_its_six_existing_return_keys(tmp_path):
    """Same additive guarantee as ``aggregate_results_cli``, for the miRNA entry point."""
    results_dir = tmp_path / "mirna_results"
    results_dir.mkdir()
    result = aggregate_mirna_results_cli(
        mirna_db="mirgenedb",
        mirna_species="human",
        results_dir=str(results_dir),
        output_dir=str(tmp_path / "aggregated_mirna"),
    )
    for key in ("status", "mirna_database", "species", "total_hits", "candidates_analyzed", "output_dir"):
        assert key in result
