"""Unit tests for running the workflow with an input FASTA file."""

import json
from pathlib import Path

import pytest

from sirnaforge.data.base import FastaUtils
from sirnaforge.models.modifications import StrandMetadata
from sirnaforge.models.sirna import DesignParameters, DesignResult, SiRNACandidate
from sirnaforge.workflow import run_sirna_workflow


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


@pytest.mark.asyncio
@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
async def test_workflow_runs_from_fasta(tmp_path, monkeypatch):
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
@pytest.mark.local_python
@pytest.mark.ci
async def test_empty_fasta_triggers_design_error(tmp_path, monkeypatch):
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
@pytest.mark.local_python
@pytest.mark.ci
async def test_invalid_sequence_in_fasta_raises(tmp_path):
    # Prepare FASTA with invalid characters
    fasta = tmp_path / "bad.fasta"
    FastaUtils.save_sequences_fasta([("bad1", "ATGNNNXYZTAA")], fasta)

    # Running the workflow should raise a validation error when constructing TranscriptInfo
    with pytest.raises(ValueError):
        await run_sirna_workflow(gene_query="bad", output_dir=str(tmp_path / "out_bad"), input_fasta=str(fasta))


@pytest.mark.asyncio
@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
async def test_single_sequence_workflow_success(tmp_path, monkeypatch):
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
