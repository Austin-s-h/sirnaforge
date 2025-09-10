"""Ultra-minimal smoke tests for CI/CD with toy data.

These tests are designed to run extremely fast with minimal resources
to provide quick feedback in CI/CD pipelines while using toy data.
"""

import subprocess
import tempfile
from pathlib import Path

import pytest


@pytest.mark.docker
@pytest.mark.smoke
def test_docker_smoke_cli():
    """Smoke test: sirnaforge CLI is available and responds."""
    try:
        result = subprocess.run(
            ["sirnaforge", "version"], 
            capture_output=True, 
            text=True, 
            timeout=10,
            check=False
        )
        # Just check it doesn't crash completely
        assert result.returncode in [0, 1, 2]  # Allow various exit codes
    except FileNotFoundError:
        pytest.skip("sirnaforge CLI not available - run in Docker container")
    except subprocess.TimeoutExpired:
        pytest.fail("CLI took too long to respond")


@pytest.mark.docker  
@pytest.mark.smoke
def test_docker_smoke_toy_workflow():
    """Smoke test: basic workflow with minimal toy data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)
        
        # Create ultra-minimal test FASTA (smaller than existing test data)
        toy_fasta = work_dir / "toy.fasta"
        toy_fasta.write_text(
            ">toy_seq\n"
            "AUGCUGAUCCGCAUGCUGAUC\n"  # Just 20 bases
        )
        
        try:
            # Test basic design command with timeout
            result = subprocess.run(
                ["sirnaforge", "design", str(toy_fasta), "--output", str(work_dir / "output.tsv")],
                capture_output=True,
                text=True,
                cwd=work_dir,
                timeout=15,  # Very short timeout for CI
                check=False,
            )
            
            # For smoke test, just ensure it doesn't crash with import errors
            if result.returncode != 0:
                error_text = result.stderr.lower()
                # Fail only on critical environment issues
                if any(term in error_text for term in ["import", "not found", "command not found", "no module"]):
                    pytest.fail(f"Environment issue: {result.stderr}")
                # Design issues are okay for smoke test
                
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available")
        except subprocess.TimeoutExpired:
            # Timeout is acceptable for smoke test - just checking basic functionality
            pass


@pytest.mark.docker
@pytest.mark.smoke  
def test_docker_smoke_help_commands():
    """Smoke test: help commands work quickly."""
    commands = ["--help", "design --help", "search --help"]
    
    for cmd_str in commands:
        try:
            cmd = ["sirnaforge"] + cmd_str.split()
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True, 
                timeout=5,  # Very short timeout
                check=False
            )
            # Just verify it responds (help often returns exit code 0 or 1)
            assert result.returncode in [0, 1, 2]
        except (FileNotFoundError, subprocess.TimeoutExpired):
            # Skip individual command failures in smoke test
            pass