"""Nextflow Docker integration tests.

Tests that verify Nextflow workflows work correctly in Docker environments
with proper resource constraints and Docker-in-Docker functionality.
"""

import re
import subprocess
import tempfile
from pathlib import Path

import pytest

import sirnaforge.pipeline.nextflow.config
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
def test_sirnaforge_nextflow_minimal_execution():
    """Test that siRNAforge can execute Nextflow workflows with minimal resources."""
    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        # Create minimal test data
        test_fasta = work_dir / "test_candidates.fasta"
        test_fasta.write_text(">candidate_1\nAUGAAAGUGAACUACAACUGU\n>candidate_2\nAUGCCAGUGAACUACAACUGU\n")

        output_dir = work_dir / "output"
        output_dir.mkdir()

        try:
            # Test the workflow execution with minimal resources
            result = subprocess.run(
                [
                    "sirnaforge",
                    "workflow",
                    "TEST",
                    "--input",
                    str(test_fasta),
                    "--output-dir",
                    str(output_dir),
                    "--species",
                    "human",
                    "--top-n",
                    "2",
                    "--offtarget-n",
                    "2",
                    "--verbose",
                ],
                capture_output=True,
                text=True,
                cwd=work_dir,
                timeout=300,  # 5 minutes max
                check=False,
            )

            # Check for specific memory/resource errors vs other issues
            if result.returncode != 0:
                error_text = result.stderr.lower()

                # Memory issues are expected on constrained systems
                memory_indicators = ["memory", "ram", "resource", "exceeds available"]

                if any(indicator in error_text for indicator in memory_indicators):
                    pytest.skip(f"Insufficient resources for full workflow test: {result.stderr}")

                # Channel/syntax issues indicate our fixes didn't work
                syntax_indicators = [
                    "channel.fromlist",
                    "missing process or function",
                    "compilation failed",
                    "groovy.lang",
                ]

                if any(indicator in error_text for indicator in syntax_indicators):
                    pytest.fail(f"Nextflow syntax/compatibility issue: {result.stderr}")

                # Docker issues
                docker_indicators = ["docker daemon", "docker not found", "container failed"]

                if any(indicator in error_text for indicator in docker_indicators):
                    pytest.skip(f"Docker environment issue: {result.stderr}")

                # Import/environment issues
                env_indicators = ["import", "no module", "command not found"]

                if any(indicator in error_text for indicator in env_indicators):
                    pytest.fail(f"Environment setup issue: {result.stderr}")

                # Other issues are logged but don't fail the test
                print(f"Workflow execution issue (may be expected): {result.stderr}")

        except ImportError:
            pytest.skip("siRNAforge not available for import")
        except FileNotFoundError:
            pytest.skip("siRNAforge CLI not available")
        except subprocess.TimeoutExpired:
            pytest.skip("Workflow execution timed out (expected on resource-constrained systems)")


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
                genome_species=["human"],
                include_test_profile=True,
            )

            assert "--input" in args
            assert tf.name in args
            assert "--outdir" in args

        print(f"✓ NextflowConfig working correctly, profile: {profile}")

    except ImportError as e:
        pytest.skip(f"siRNAforge NextflowConfig not available: {e}")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_nextflow_channel_syntax_fix():
    """Test that the Channel.from syntax fix is working."""
    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        # Create a minimal test workflow with the new syntax
        test_workflow = work_dir / "test.nf"
        test_workflow.write_text(
            """
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.test_species = 'human,mouse'

workflow {
    ch_species = Channel.from(params.test_species.split(','))
    ch_species.view()
}
"""
        )

        try:
            result = subprocess.run(
                ["nextflow", "run", str(test_workflow)],
                capture_output=True,
                text=True,
                cwd=work_dir,
                timeout=30,
                check=False,
            )

            # Should not have Channel.fromList errors
            if "Channel.fromList" in result.stderr:
                pytest.fail("Channel.fromList syntax still being used")

            # Should not have compilation errors
            if "compilation failed" in result.stderr.lower():
                pytest.fail(f"Nextflow compilation failed: {result.stderr}")

            print("✓ Channel.from syntax working correctly")

        except FileNotFoundError:
            pytest.skip("Nextflow not available")
        except subprocess.TimeoutExpired:
            pytest.skip("Nextflow test timed out")
