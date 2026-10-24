"""Pinned references pass through configuration, task acquisition and report evidence."""

import hashlib
import json
from pathlib import Path

import pytest

from sirnaforge.config.mirna_reference import read_reference_evidence, stage_reference_manifest
from sirnaforge.core.off_target import run_mirna_seed_analysis
from sirnaforge.core.screening_evidence import read_evidence
from sirnaforge.data.base import DatabaseType, TranscriptInfo
from sirnaforge.data.mirna_manager import MiRNADatabaseManager
from sirnaforge.models.evidence import EvidenceStatus
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

pytestmark = pytest.mark.unit


def reference_manifest(tmp_path: Path, sequence: str = "A" * 21) -> Path:
    """Build a checksummed synthetic mature reference."""
    fasta = tmp_path / "mature.fa"
    fasta.write_text(f">hsa-test\n{sequence}\n")
    manifest = tmp_path / "pins.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "references": [
                    {
                        "source": "mirbase",
                        "species": "human",
                        "release": "test-release",
                        "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
                        "path": fasta.name,
                    }
                ],
            }
        )
    )
    return manifest


def test_config_stages_local_references_for_isolated_task(tmp_path):
    """Staged FASTAs remain available after original local input removal."""
    manifest = reference_manifest(tmp_path)
    config = WorkflowConfig(gene_query="TP53", output_dir=str(tmp_path / "run"), mirna_reference_manifest=manifest)
    staged = Path(config.nextflow_config["mirna_reference_manifest"])
    entry = json.loads(staged.read_text())["references"][0]
    manifest.unlink()
    (tmp_path / "mature.fa").unlink()
    assert Path(entry["path"]).read_text().startswith(">hsa-test")
    assert entry["release"] == "test-release"
    assert hashlib.sha256(Path(entry["path"]).read_bytes()).hexdigest() == entry["sha256"]


@pytest.mark.parametrize("positive", [False, True])
def test_pinned_task_publishes_actual_reference_and_completed_membership(tmp_path, monkeypatch, positive):
    """The actual scanner attests reference provenance and submitted sequence membership."""
    monkeypatch.setenv("SIRNAFORGE_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(MiRNADatabaseManager, "_download_file", lambda *_a, **_kw: pytest.fail("network forbidden"))
    staged = stage_reference_manifest(reference_manifest(tmp_path), tmp_path / "staged")
    candidates = tmp_path / "guides.fa"
    candidates.write_text(">guide\n" + ("A" * 21 if positive else "C" * 21) + "\n")
    output = tmp_path / "run" / "off_target" / "results" / "mirna"
    run_mirna_seed_analysis(
        candidates,
        "batch",
        "mirbase",
        ["human"],
        output,
        backend="exhaustive_python",
        evidence_dir=output,
        reference_manifest=staged,
    )
    contract = read_reference_evidence(tmp_path / "run", source="mirbase", species=["human"], skipped=False)
    evidence = contract["mirna"]
    assert evidence["screening"]["status"] == "complete"
    assert evidence["screening"]["guide_sequences"] == ["A" * 21 if positive else "C" * 21]
    assert evidence["sources"][0]["references"]["release"] == "test-release"
    summary = json.loads((output / "batch_mirna_summary.json").read_text())
    assert bool(summary["total_hits"]) is positive


def test_unavailable_source_emits_failed_evidence_not_zero_completion(tmp_path, monkeypatch):
    """Failure emits unknown counts and cannot claim a completed zero-hit screen."""
    monkeypatch.setenv("SIRNAFORGE_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(MiRNADatabaseManager, "_download_file", lambda *_a, **_kw: None)
    monkeypatch.setattr("sirnaforge.data.mirna_manager.time.sleep", lambda *_a: None)
    candidates = tmp_path / "guides.fa"
    candidates.write_text(">guide\n" + "C" * 21 + "\n")
    output = tmp_path / "result"
    run_mirna_seed_analysis(
        candidates, "batch", "mirbase", ["human"], output, backend="exhaustive_python", evidence_dir=output
    )
    contract = json.loads((output / "mirna_reference_evidence.json").read_text())["mirna"]
    assert contract["acquisition"]["status"] == "failed"
    assert contract["screening"]["status"] == "failed"
    assert read_evidence(output / "mirna_seed_human_evidence.json").entry.status is EvidenceStatus.FAILED


def test_skipped_run_does_not_reuse_stale_task_completion(tmp_path):
    """Explicit skips override artifacts from an earlier run in the same directory."""
    artifact = tmp_path / "off_target" / "results" / "mirna" / "mirna_reference_evidence.json"
    artifact.parent.mkdir(parents=True)
    artifact.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "mirna": {
                    "source": "mirbase",
                    "requested_species": ["human"],
                    "screening": {"status": "complete"},
                },
            }
        )
    )
    contract = read_reference_evidence(tmp_path, source="mirbase", species=["human"], skipped=True)
    assert contract["mirna"]["screening"]["status"] == "skipped"


@pytest.mark.parametrize("authority", [None, "Ensembl REST lookup"])
def test_target_identity_requires_real_canonical_authority(tmp_path, authority):
    """A provider's heuristic boolean cannot establish canonical identity without provenance."""
    workflow = SiRNAWorkflow(WorkflowConfig(output_dir=tmp_path / "run", gene_query="TP53"))
    transcript = TranscriptInfo(
        transcript_id="ENST1",
        gene_id="ENSG1",
        gene_name="TP53",
        database=DatabaseType.ENSEMBL,
        is_canonical=True,
        canonical_source=authority,
        canonical_release="test-release" if authority else None,
    )
    workflow._record_target_identity([transcript])
    identity = workflow._target_transcript_metadata[0]
    assert identity["is_canonical"] is (True if authority else None)
    assert identity["canonical_source"] == authority


@pytest.mark.parametrize("allowed,digest", [(False, "current"), (True, "other")])
def test_producer_evidence_requires_current_success_and_guide_digest(tmp_path, allowed, digest):
    """A stale artifact cannot attest a failed invocation or a different submitted guide set."""
    artifact = tmp_path / "off_target" / "results" / "mirna" / "mirna_reference_evidence.json"
    artifact.parent.mkdir(parents=True)
    artifact.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "mirna": {
                    "source": "mirbase",
                    "requested_species": ["human"],
                    "screening": {"status": "complete", "guide_set_digest": "current"},
                },
            }
        )
    )
    contract = read_reference_evidence(
        tmp_path, source="mirbase", species=["human"], skipped=False, allow_artifact=allowed, expected_digest=digest
    )
    assert contract["mirna"]["screening"]["status"] == "unknown"
