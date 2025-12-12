"""Unit tests for running the workflow with an input FASTA file."""

import json
from pathlib import Path

import pytest

from sirnaforge.data.base import FastaUtils
from sirnaforge.models.modifications import StrandMetadata
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.utils.control_candidates import DIRTY_CONTROL_LABEL
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig, run_sirna_workflow


def _write_test_fasta(path: Path) -> None:
    sequences = [("trans1", "ATG" + "A" * 300 + "TAA"), ("trans2", "ATG" + "C" * 300 + "TGA")]
    FastaUtils.save_sequences_fasta(sequences, path)


def _load_patisiran_metadata() -> tuple[StrandMetadata, StrandMetadata]:
    examples_dir = Path(__file__).parent.parent.parent / "examples" / "modification_patterns"
    onpattro_file = examples_dir / "fda_approved_onpattro.json"

    with onpattro_file.open() as fh:
        data = json.load(fh)

    sequences = data.get("sequences", {})
    guide_meta = StrandMetadata.model_validate(sequences["patisiran_ttr_guide"])
    passenger_meta = StrandMetadata.model_validate(sequences["patisiran_ttr_passenger"])
    return guide_meta, passenger_meta


def _make_candidate(candidate_id: str, guide_seq: str, passes_filters: object = True) -> SiRNACandidate:
    """Create a minimal SiRNACandidate for tests."""
    passenger_seq = guide_seq.translate(str.maketrans("ATCG", "TAGC"))
    gc_content = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq) * 100

    return SiRNACandidate(
        id=candidate_id,
        transcript_id="trans1",
        position=1,
        guide_sequence=guide_seq,
        passenger_sequence=passenger_seq,
        gc_content=gc_content,
        length=len(guide_seq),
        asymmetry_score=0.5,
        component_scores={"empirical": 1.0},
        composite_score=90.0,
        passes_filters=passes_filters,
    )


@pytest.mark.asyncio
@pytest.mark.unit
@pytest.mark.ci
async def test_workflow_runs_from_fasta(tmp_path, monkeypatch):
    """Test complete workflow execution from FASTA input file."""
    # Prepare input FASTA
    fasta = tmp_path / "test_input.fasta"
    _write_test_fasta(fasta)

    # Mock ORFAnalyzer.analyze_transcript to return a minimal object with has_valid_orf
    class DummyAnalysis:
        def __init__(self, tid):
            self.transcript_id = tid
            self.sequence_length = 306
            self.gc_content = 50.0
            self.orfs = []
            self.longest_orf = None
            self.has_valid_orf = True

    async def fake_analyze_transcript(self, transcript):
        return DummyAnalysis(transcript.transcript_id)

    monkeypatch.setattr("sirnaforge.data.orf_analysis.ORFAnalyzer.analyze_transcript", fake_analyze_transcript)

    # Mock SiRNADesigner.design_from_file to return a simple DesignResult using embedded Patisiran metadata
    def fake_design_from_file(self, path):
        guide_meta, passenger_meta = _load_patisiran_metadata()
        guide_seq = guide_meta.sequence.upper()
        passenger_seq = passenger_meta.sequence.upper()

        gc_content = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq) * 100

        cand = SiRNACandidate(
            id="c1",
            transcript_id="trans1",
            position=1,
            guide_sequence=guide_seq,
            passenger_sequence=passenger_seq,
            gc_content=gc_content,
            length=len(guide_seq),
            asymmetry_score=0.5,
            duplex_stability=-20.0,
            paired_fraction=0.1,
            component_scores={"empirical": 1.0},
            composite_score=100.0,
            guide_metadata=guide_meta,
            passenger_metadata=passenger_meta,
        )
        return DesignResult(
            input_file=str(path),
            parameters=DesignParameters(),
            candidates=[cand],
            top_candidates=[cand],
            total_sequences=2,
            total_candidates=1,
            filtered_candidates=1,
            processing_time=0.1,
            tool_versions={"test": "1.0"},
        )

    monkeypatch.setattr("sirnaforge.core.design.SiRNADesigner.design_from_file", fake_design_from_file)

    async def fake_offtarget(self, design_results):  # noqa: ARG001
        return {"status": "skipped", "reason": "unit_test"}

    monkeypatch.setattr("sirnaforge.workflow.SiRNAWorkflow.step5_offtarget_analysis", fake_offtarget)

    # Run workflow
    results = await run_sirna_workflow(gene_query="test", output_dir=str(tmp_path / "out"), input_fasta=str(fasta))

    # Assertions on results structure
    assert "transcript_summary" in results
    assert results["transcript_summary"]["total_transcripts"] == 2
    assert results["design_summary"]["total_candidates"] == 1


@pytest.mark.asyncio
@pytest.mark.unit
@pytest.mark.ci
async def test_empty_fasta_triggers_design_error(tmp_path, monkeypatch):
    """Test that empty FASTA files trigger appropriate design errors."""
    # Create empty FASTA
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    # Mock SiRNADesigner.design_from_file to raise when given empty input
    def fake_design_from_file(self, path):
        raise ValueError("No input sequences for design")

    monkeypatch.setattr("sirnaforge.core.design.SiRNADesigner.design_from_file", fake_design_from_file)

    with pytest.raises(ValueError):
        await run_sirna_workflow(gene_query="empty", output_dir=str(tmp_path / "out_empty"), input_fasta=str(fasta))


@pytest.mark.asyncio
@pytest.mark.unit
@pytest.mark.ci
async def test_invalid_sequence_in_fasta_raises(tmp_path):
    """Test that FASTA files with invalid sequences raise validation errors."""
    # Prepare FASTA with invalid characters
    fasta = tmp_path / "bad.fasta"
    FastaUtils.save_sequences_fasta([("bad1", "ATGNNNXYZTAA")], fasta)

    # Running the workflow should raise a validation error when constructing TranscriptInfo
    with pytest.raises(ValueError):
        await run_sirna_workflow(gene_query="bad", output_dir=str(tmp_path / "out_bad"), input_fasta=str(fasta))


@pytest.mark.asyncio
@pytest.mark.unit
@pytest.mark.ci
async def test_single_sequence_workflow_success(tmp_path, monkeypatch):
    """Test successful workflow execution with single FASTA sequence."""
    # Single valid FASTA sequence
    fasta = tmp_path / "single.fasta"
    FastaUtils.save_sequences_fasta([("single", "ATG" + "A" * 300 + "TAA")], fasta)

    # Minimal analysis stub
    class DummyAnalysis:
        def __init__(self, tid):
            self.transcript_id = tid
            self.sequence_length = 306
            self.gc_content = 50.0
            self.orfs = []
            self.longest_orf = None
            self.has_valid_orf = True

    async def fake_analyze_transcript(self, transcript):
        return DummyAnalysis(transcript.transcript_id)

    monkeypatch.setattr("sirnaforge.data.orf_analysis.ORFAnalyzer.analyze_transcript", fake_analyze_transcript)

    # Fake designer returns empty-but-valid DesignResult
    def fake_design_from_file(self, path):
        guide_meta, passenger_meta = _load_patisiran_metadata()
        guide_seq = guide_meta.sequence.upper()
        passenger_seq = passenger_meta.sequence.upper()

        gc_content = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq) * 100

        cand = SiRNACandidate(
            id="c1",
            transcript_id="single",
            position=1,
            guide_sequence=guide_seq,
            passenger_sequence=passenger_seq,
            gc_content=gc_content,
            length=len(guide_seq),
            asymmetry_score=0.5,
            duplex_stability=-20.0,
            paired_fraction=0.1,
            component_scores={"empirical": 1.0},
            composite_score=100.0,
            guide_metadata=guide_meta,
            passenger_metadata=passenger_meta,
        )

        return DesignResult(
            input_file=str(path),
            parameters=DesignParameters(),
            candidates=[cand],
            top_candidates=[cand],
            total_sequences=1,
            total_candidates=1,
            filtered_candidates=1,
            processing_time=0.1,
            tool_versions={"test": "1.0"},
        )

    monkeypatch.setattr("sirnaforge.core.design.SiRNADesigner.design_from_file", fake_design_from_file)

    async def fake_offtarget(self, design_results):  # noqa: ARG001
        return {"status": "skipped", "reason": "unit_test"}

    monkeypatch.setattr("sirnaforge.workflow.SiRNAWorkflow.step5_offtarget_analysis", fake_offtarget)

    results = await run_sirna_workflow(
        gene_query="single", output_dir=str(tmp_path / "out_single"), input_fasta=str(fasta)
    )

    assert results["transcript_summary"]["total_transcripts"] == 1
    assert results["design_summary"]["total_candidates"] == 1


@pytest.mark.unit
def test_offtarget_selection_includes_dirty_controls(tmp_path):
    """Dirty control candidates must appear in off-target FASTA selection."""
    params = DesignParameters(top_n=1)
    config = WorkflowConfig(output_dir=tmp_path / "out_dirty", gene_query="tp53", design_params=params)
    workflow = SiRNAWorkflow(config)

    primary = _make_candidate("c_primary", "ATGCGATGCGATGCGATGCGC")
    dirty = _make_candidate(
        "c_dirty",
        "TTGCGATGCGATGCGATGCGA",
        passes_filters=SiRNACandidate.FilterStatus.DIRTY_CONTROL,
    )
    dirty.quality_issues = [DIRTY_CONTROL_LABEL]

    design_result = DesignResult(
        input_file="input",
        parameters=params,
        candidates=[primary, dirty],
        top_candidates=[primary, dirty],
        total_sequences=1,
        total_candidates=2,
        filtered_candidates=1,
        processing_time=0.1,
        tool_versions={},
    )

    selected = workflow._select_candidates_for_offtarget(design_result)

    assert primary in selected
    assert dirty in selected
    assert len(selected) == 2


def test_load_offtarget_aggregates_reads_summaries(tmp_path):
    """Aggregated Nextflow summaries should be surfaced for downstream schemas."""
    config = WorkflowConfig(output_dir=tmp_path / "out_agg", gene_query="tp53")
    workflow = SiRNAWorkflow(config)

    results_dir = config.output_dir / "off_target" / "results"
    aggregated_dir = results_dir / "aggregated"
    aggregated_dir.mkdir(parents=True, exist_ok=True)

    transcriptome_summary = {"total_results": 5, "human_hits": 2, "other_species_hits": 3}
    mirna_summary = {"total_mirna_hits": 4, "human_hits": 1, "other_species_hits": 3}

    (aggregated_dir / "combined_summary.json").write_text(json.dumps(transcriptome_summary))
    (aggregated_dir / "combined_mirna_summary.json").write_text(json.dumps(mirna_summary))

    aggregates = workflow._load_offtarget_aggregates(results_dir)

    assert aggregates["transcriptome"]["human_hits"] == 2
    assert aggregates["mirna"]["other_species_hits"] == 3


@pytest.mark.asyncio
@pytest.mark.unit
async def test_process_nextflow_results_includes_aggregated(monkeypatch, tmp_path):
    """Workflow summaries must include aggregated human vs other breakdowns."""
    config = WorkflowConfig(output_dir=tmp_path / "out_process", gene_query="tp53")
    workflow = SiRNAWorkflow(config)

    candidate = _make_candidate("cand1", "ATGCGATGCGATGCGATGCGC")

    async def fake_parse(_self, output_dir):  # noqa: ARG001
        return {
            "status": "completed",
            "method": "nextflow",
            "output_dir": str(output_dir),
            "results": {candidate.id: {"off_target_count": 0, "off_target_score": 0.0, "hits": []}},
        }

    monkeypatch.setattr(SiRNAWorkflow, "_parse_nextflow_results", fake_parse)

    sentinel = {"transcriptome": {"human_hits": 1, "other_species_hits": 0}}
    monkeypatch.setattr(SiRNAWorkflow, "_load_offtarget_aggregates", lambda *_: sentinel)

    output_dir = config.output_dir / "off_target" / "results"
    output = await workflow._process_nextflow_results([candidate], output_dir, {"status": "completed"})

    assert output["aggregated"] == sentinel


@pytest.mark.asyncio
@pytest.mark.unit
async def test_process_nextflow_results_flags_missing_species(monkeypatch, tmp_path):
    """Missing transcriptome species should downgrade workflow status to partial."""
    config = WorkflowConfig(output_dir=tmp_path / "out_process_missing", gene_query="tp53")
    workflow = SiRNAWorkflow(config)

    candidate = _make_candidate("cand_missing", "ATGCGATGCGATGCGATGCGC")

    async def fake_parse(_self, output_dir):  # noqa: ARG001
        return {
            "status": "completed",
            "method": "nextflow",
            "output_dir": str(output_dir),
            "results": {candidate.id: {"off_target_count": 0, "off_target_score": 0.0, "hits": []}},
        }

    monkeypatch.setattr(SiRNAWorkflow, "_parse_nextflow_results", fake_parse)

    sentinel = {
        "transcriptome": {
            "human_hits": 0,
            "other_species_hits": 0,
            "missing_species": ["human"],
            "species_analyzed": ["human"],
            "hits_per_species": {"human": 0},
        }
    }
    monkeypatch.setattr(SiRNAWorkflow, "_load_offtarget_aggregates", lambda *_: sentinel)

    output_dir = config.output_dir / "off_target" / "results"
    result = await workflow._process_nextflow_results([candidate], output_dir, {"status": "completed"})

    assert result["status"] == "partial"
    assert sentinel["transcriptome"] in result["aggregated"].values()
    assert any("No transcriptome alignment files" in warning for warning in result.get("warnings", []))
