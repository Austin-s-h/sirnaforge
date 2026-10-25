"""Fresh pinned references reach the installed CLI, Nextflow producer, and report."""

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
from release_artifacts import assert_report_rendered

pytestmark = [pytest.mark.integration, pytest.mark.runs_in_container]


@pytest.mark.parametrize("valid", [True, False])
def test_fresh_pinned_mirbase_workflow(tmp_path: Path, valid: bool) -> None:
    """A verified local pin completes offline; a bad digest remains failed evidence."""
    reference = tmp_path / "mature.fa"
    reference.write_text(">hsa-synthetic-fixture\nACGUACGUACGUACGUACGUA\n")
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    pin = tmp_path / "references.json"
    pin.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "references": [
                    {
                        "source": "mirbase",
                        "species": "human",
                        "release": "synthetic-offline-fixture",
                        "sha256": digest if valid else "0" * 64,
                        "path": reference.name,
                    }
                ],
            }
        )
    )
    design_input = tmp_path / "transcript.fa"
    design_input.write_text(">SYNTHETIC_TRANSCRIPT\n" + "ACGT" * 40 + "\n")
    output = tmp_path / "run"
    result = subprocess.run(
        [
            "sirnaforge",
            "workflow",
            "PINNED_FIXTURE",
            "--input-fasta",
            str(design_input),
            "--output-dir",
            str(output),
            "--run-mode",
            "exploratory",
            "--species",
            "human",
            "--mirna-species",
            "human",
            "--mirna-db",
            "mirbase",
            "--mirna-reference-manifest",
            str(pin),
            "--top-n",
            "5",
        ],
        env={
            **os.environ,
            "SIRNAFORGE_CACHE_DIR": str(tmp_path / "fresh-cache"),
            "SIRNAFORGE_ANNOTATION_PROVIDER": "none",
        },
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    evidence_path = output / "off_target/results/mirna/mirna_reference_evidence.json"
    assert evidence_path.is_file(), result.stdout + result.stderr
    evidence = json.loads(evidence_path.read_text())["mirna"]
    expected = "complete" if valid else "failed"
    assert evidence["acquisition"]["status"] == expected
    assert evidence["screening"]["status"] == expected
    if valid:
        metadata = evidence["sources"][0]["references"]
        assert metadata["release"] == "synthetic-offline-fixture"
        assert metadata["raw_sha256"] == digest
        assert metadata["acquisition_mode"] == "pinned_import"
        assert evidence["screening"]["submitted_guides"] > 0
    else:
        assert evidence["acquisition"]["reason"]
    assert_report_rendered(output)
