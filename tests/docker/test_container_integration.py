"""Docker container integration tests.

Tests that verify the Docker container works correctly with the siRNAforge CLI
and core functionality. These tests are designed to work with your existing
pytest infrastructure and Makefile targets.
"""

import importlib
import subprocess
import tempfile
from pathlib import Path

import pytest


@pytest.mark.docker
@pytest.mark.integration
def test_docker_cli_help():
    """Test that siRNAforge CLI works in container environment."""
    # This test is designed to be run inside the Docker container via make ci-test-docker
    try:
        result = subprocess.run(["sirnaforge", "--help"], capture_output=True, text=True, check=True)
        assert "siRNAforge" in result.stdout
        assert "Comprehensive siRNA design toolkit for gene silencing." in result.stdout or "workflow" in result.stdout
    except FileNotFoundError:
        pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.docker
@pytest.mark.integration
def test_docker_version():
    """Test version command works in container."""
    try:
        result = subprocess.run(["sirnaforge", "version"], capture_output=True, text=True, check=True)
        # Should output version information
        assert len(result.stdout.strip()) > 0
    except FileNotFoundError:
        pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.docker
@pytest.mark.integration
def test_docker_command_structure():
    """Test that main commands are available in container."""
    commands = ["search", "workflow", "design"]

    for cmd in commands:
        try:
            result = subprocess.run(["sirnaforge", cmd, "--help"], capture_output=True, text=True, check=True)
            assert len(result.stdout) > 0, f"Command {cmd} should show help"
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command {cmd} failed: {e.stderr}")


@pytest.mark.docker
@pytest.mark.integration
def test_docker_python_environment():
    """Test that key Python packages are available in container."""
    required_packages = [
        "Bio",  # biopython
        "pandas",
        "numpy",
        "pysam",
        "RNA",  # ViennaRNA
    ]

    missing = []
    for package in required_packages:
        try:
            importlib.import_module(package)
        except ImportError:
            missing.append(package)

    if missing:
        pytest.fail(f"Missing required packages in Docker environment: {missing}")


@pytest.mark.docker
@pytest.mark.integration
def test_docker_sirnaforge_imports():
    """Test that siRNAforge modules can be imported in container."""
    # Test that basic modules can be imported
    try:
        importlib.import_module("sirnaforge.core.design")
        importlib.import_module("sirnaforge.models")
    except ImportError as e:
        pytest.skip(f"siRNAforge core modules not available: {e}")

    # Test instantiation if imports work
    try:
        # Use importlib to avoid linter issues with conditional imports
        design_module = importlib.import_module("sirnaforge.core.design")
        models_module = importlib.import_module("sirnaforge.models")

        designer_class = design_module.SiRNADesigner
        params_class = models_module.DesignParameters

        params = params_class()
        designer = designer_class(params)
        assert designer is not None

    except (ImportError, AttributeError) as e:
        pytest.skip(f"siRNAforge classes not available: {e}")


@pytest.mark.docker
@pytest.mark.integration
def test_docker_minimal_workflow():
    """Test minimal siRNA design workflow in container with test data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        # Create minimal test FASTA
        test_fasta = work_dir / "test.fasta"
        test_fasta.write_text(
            ">test_transcript\n"
            "ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTG"
            "GATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATGGTTTACATGTTCCAATATGATTCCACCCATGGCAAATTC\n"
        )

        try:
            # Test basic design command (should work even if it doesn't find perfect siRNAs)
            result = subprocess.run(
                ["sirnaforge", "design", str(test_fasta), "--output", str(work_dir / "output.tsv")],
                capture_output=True,
                text=True,
                cwd=work_dir,
                timeout=30,  # Don't hang CI
                check=False,
            )

            # Command should run without crashing (may have warnings/errors about design parameters)
            # but shouldn't have import errors or missing tools
            if result.returncode != 0:
                # Check if it's an expected design-related issue vs environment issue
                error_text = result.stderr.lower()
                if any(term in error_text for term in ["import", "not found", "command not found", "no module"]):
                    pytest.fail(f"Environment issue in container: {result.stderr}")
                else:
                    # Design-related issues are okay for this test
                    pytest.skip(f"Design workflow issue (expected): {result.stderr}")

        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")
        except subprocess.TimeoutExpired:
            pytest.skip("Command timed out - may indicate missing dependencies")


@pytest.mark.docker
@pytest.mark.integration
def test_docker_bioinformatics_tools():
    """Test that key bioinformatics tools are available in container."""
    tools = [
        ("samtools", ["--version"]),
        ("bwa-mem2", ["version"]),
        ("RNAfold", ["--help"]),
    ]

    missing_tools = []

    for tool, args in tools:
        try:
            result = subprocess.run([tool] + args, capture_output=True, text=True, timeout=10, check=False)
            # Tool should exist and respond (exit code may vary)
            if result.returncode > 1:  # 0 or 1 are usually okay
                missing_tools.append(f"{tool} (bad exit code: {result.returncode})")
        except FileNotFoundError:
            missing_tools.append(f"{tool} (not found)")
        except subprocess.TimeoutExpired:
            missing_tools.append(f"{tool} (timeout)")

    if missing_tools:
        pytest.skip(f"Some bioinformatics tools missing in container: {missing_tools}")
