"""Nextflow Docker integration tests.

Tests that verify Nextflow workflows work correctly in Docker environments
with proper resource constraints and Docker-in-Docker functionality.
"""

import json
import re
import subprocess
import tempfile
from pathlib import Path

import pytest

import sirnaforge.pipeline.nextflow.config
from sirnaforge.models.evidence import EVIDENCE_SCHEMA_VERSION, EvidenceStatus
from sirnaforge.pipeline.nextflow.runner import NextflowRunner


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_nextflow_available():
    """Test that Nextflow is available in the Docker container."""
    try:
        result = subprocess.run(["nextflow", "-version"], capture_output=True, text=True, timeout=30, check=True)
        assert "nextflow" in result.stdout.lower()

        # Check version compatibility (should be 25.x or higher)
        lines = result.stdout.strip().split("\n")
        for line in lines:
            if "version" in line.lower():
                # Extract version number
                version_match = re.search(r"(\d+)\.(\d+)\.(\d+)", line)
                if version_match:
                    major = int(version_match.group(1))
                    assert major >= 25, f"Nextflow version too old: {line}"
                break

    except (FileNotFoundError, subprocess.CalledProcessError):
        pytest.skip("Nextflow not available - run this test in Docker container with Nextflow")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_nextflow_local_profile():
    """Test that Nextflow can use local profiles correctly."""
    try:
        # Test that local profile is available
        result = subprocess.run(
            ["nextflow", "config", "-profile", "local", "-show-profiles"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,  # May fail if no config file present, that's ok
        )

        # If nextflow config works, it should not crash with import errors
        if result.returncode != 0 and any(term in result.stderr.lower() for term in ["import", "module", "class"]):
            pytest.fail(f"Nextflow config failed with import error: {result.stderr}")

    except FileNotFoundError:
        pytest.skip("Nextflow not available")
    except subprocess.TimeoutExpired:
        pytest.skip("Nextflow config command timed out")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_sirnaforge_nextflow_workflow_syntax():
    """Test that siRNAforge Nextflow workflow has valid syntax."""
    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        try:
            # Use the runner to get the workflow file
            runner = NextflowRunner()
            main_nf = runner.get_main_workflow()

            if not main_nf.exists():
                pytest.skip(f"Nextflow workflow not found at {main_nf}")

            # Test syntax with dry-run
            result = subprocess.run(
                [
                    "nextflow",
                    "run",
                    str(main_nf),
                    "--input",
                    "/dev/null",  # Invalid input to trigger early exit
                    "--dry-run",
                ],
                capture_output=True,
                text=True,
                cwd=work_dir,
                timeout=60,
                check=False,
            )

            # Should fail due to invalid input, not syntax errors
            if "Channel.fromList" in result.stderr:
                pytest.fail("Nextflow syntax error - Channel.fromList compatibility issue")

            # Check for other syntax errors
            error_indicators = [
                "compilation failed",
                "syntax error",
                "unexpected token",
                "groovy.lang.MissingMethodException",
            ]

            for indicator in error_indicators:
                if indicator.lower() in result.stderr.lower():
                    pytest.fail(f"Nextflow syntax error detected: {result.stderr}")

        except (ImportError, AttributeError, TypeError):
            pytest.skip("siRNAforge workflows module not available for import")
        except FileNotFoundError:
            pytest.skip("Nextflow not available")
        except subprocess.TimeoutExpired:
            pytest.skip("Nextflow syntax check timed out")


@pytest.mark.integration
@pytest.mark.runs_in_container
@pytest.mark.requires_nextflow
def test_stub_run_emits_evidence_and_cannot_look_like_a_screen():
    """A real `-stub-run` materialises #100's evidence files, and they admit nothing was aligned.

    This replaces `test_sirnaforge_nextflow_minimal_execution`, which ran
    `sirnaforge workflow TEST --input <fasta> ...`. There is no `--input` option on that command --
    there never has been -- so every run of it got a Typer usage error, and the test passed anyway:
    its only failure paths were substring matches on stderr, and "no such option" matched none of
    them. It asserted nothing about Nextflow while claiming to be the minimal-execution test. The
    embedded-Nextflow path it meant to cover is exercised with assertions by
    `test_workflow_modes.py::test_minimal_toy_workflow`.

    What is covered here that nothing else covers: the `.nf` side of #100 is otherwise asserted by
    regex over the workflow text (`tests/unit/test_evidence_plan_threading.py`, whose docstring says
    no test there invokes a real nextflow binary). Regex cannot see a renamed emit glob or a stub
    block that stops satisfying its declared outputs -- an execution can. Measured at ~4s: no
    aligner, no reference data, no network.

    The second half of the name is the contract that matters. A stub run must be *legible* as a
    non-result: `producer: stub`, `source: synthesized`, and a FAILED unit whose detail says no
    aligner executed. If a stub ever emitted `complete`, every consumer that trusts the evidence
    contract would read a screen that never happened.
    """
    runner = NextflowRunner()
    main_workflow = runner.get_main_workflow()
    assert main_workflow.is_file(), f"embedded workflow missing from the image: {main_workflow}"

    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)
        candidates = work_dir / "candidates.fasta"
        candidates.write_text(">candidate_1\nAUGAAAGUGAACUACAACUGU\n>candidate_2\nAUGCCAGUGAACUACAACUGU\n")
        outdir = work_dir / "out"

        result = subprocess.run(
            [
                "nextflow",
                "run",
                str(main_workflow),
                "-stub-run",
                "--input",
                str(candidates),
                "--outdir",
                str(outdir),
                "--mirna_species",
                "human",
            ],
            capture_output=True,
            text=True,
            cwd=work_dir,
            timeout=300,
            check=False,
        )
        assert result.returncode == 0, (
            f"stub run failed:\nSTDOUT: {result.stdout[-2000:]}\nSTDERR: {result.stderr[-2000:]}"
        )

        # The aggregate emits the plan/evidence envelope even with no plan supplied. Its shape is the
        # thing every 0.7.1 reader joins on.
        envelope = json.loads((outdir / "aggregated" / "evidence.json").read_text())
        assert envelope["schema_version"] == EVIDENCE_SCHEMA_VERSION
        assert set(envelope) >= {"plan", "evidence", "sources", "unplanned"}, sorted(envelope)

        # The per-unit evidence the stub block declares as an output must actually materialise.
        stub_files = sorted((outdir / "mirna").glob("mirna_seed_*_evidence.json"))
        assert stub_files, f"stub emitted no per-unit evidence: {sorted((outdir / 'mirna').glob('*'))}"

        stub_evidence = json.loads(stub_files[0].read_text())
        assert stub_evidence["producer"] == "stub"
        assert stub_evidence["source"] == "synthesized"
        entry = stub_evidence["entry"]
        assert entry["status"] == EvidenceStatus.FAILED.value, (
            f"a stub run must not publish a completed unit: {entry['status']}"
        )
        assert entry["detail"], "a stub unit must say why it is not a result"


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_nextflow_config_generation():
    """Test that siRNAforge can generate valid Nextflow configuration."""
    try:
        # Test NextflowConfig class functionality
        config = sirnaforge.pipeline.nextflow.config.NextflowConfig()

        # Test environment detection
        env_info = config.get_environment_info()
        assert hasattr(env_info, "running_in_docker")
        assert hasattr(env_info, "docker_available")
        assert hasattr(env_info, "recommended_profile")

        # Test profile selection
        profile = config.get_execution_profile()
        assert profile in ["docker", "local", "test", "conda", "singularity"]

        # Test argument generation
        with tempfile.NamedTemporaryFile(suffix=".fasta") as tf:
            args = config.get_nextflow_args(
                input_file=Path(tf.name),
                output_dir=Path("/tmp/test"),
                screen_species=["human"],
                include_test_profile=True,
            )

            assert "--input" in args
            assert tf.name in args
            assert "--outdir" in args

        print(f"✓ NextflowConfig working correctly, profile: {profile}")

    except ImportError as e:
        pytest.skip(f"siRNAforge NextflowConfig not available: {e}")
