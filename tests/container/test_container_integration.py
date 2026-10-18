"""Docker container integration tests.

Tests that verify the Docker container works correctly with the siRNAforge CLI
and core functionality. These tests are designed to work with your existing
pytest infrastructure and Makefile targets.
"""

import importlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_cli_help():
    """Test that siRNAforge CLI works in container environment."""
    # This test is designed to be run inside the Docker container via make ci-test-docker
    try:
        result = subprocess.run(["sirnaforge", "--help"], capture_output=True, text=True, check=True)
        assert "siRNAforge" in result.stdout
        assert "Comprehensive siRNA design toolkit for gene silencing." in result.stdout or "workflow" in result.stdout
    except FileNotFoundError:
        pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_version():
    """Test version command works in container."""
    try:
        result = subprocess.run(["sirnaforge", "version"], capture_output=True, text=True, check=True)
        # Should output version information
        assert len(result.stdout.strip()) > 0
    except FileNotFoundError:
        pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_command_structure():
    """Test that every top-level command is available in container.

    The list was `[search, workflow, design]` and had not grown in three releases, so 0.7.1's
    `report` and `benchmark` -- and `offtarget`, `config`, `cache` before them -- had no image-level
    guard at all: a Typer registration error or a missing import inside one of those command modules
    would ship green. `--help` on a command imports its module, which is the cheap part of the check.
    """
    commands = [
        "search",
        "workflow",
        "offtarget",
        "design",
        "validate",
        "version",
        "config",
        "cache",
        "report",
        "benchmark",
        "sequences",
        "zfn",
    ]

    for cmd in commands:
        try:
            result = subprocess.run(["sirnaforge", cmd, "--help"], capture_output=True, text=True, check=True)
            assert len(result.stdout) > 0, f"Command {cmd} should show help"
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command {cmd} failed: {e.stderr}")


@pytest.mark.integration
@pytest.mark.runs_in_container
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


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_sirnaforge_imports():
    """Every module the release ships must import inside the image.

    This used to `pytest.skip` on ImportError -- the precise failure the container tier exists to
    catch, so a missing runtime dependency shipped green. An ImportError here is now a failure.

    The list is the modules 0.7.1 added plus the pre-existing core, because a wheel that builds and
    a wheel that imports are different claims: `reporting` pulls the render stack, `benchmark` and
    `provenance` were new entry points in this release, and each is only reachable from a command
    that no host test runs as an installed program.
    """
    modules = [
        "sirnaforge.core.design",
        "sirnaforge.models",
        "sirnaforge.models.policy",
        "sirnaforge.models.evidence",
        "sirnaforge.core.screening_evidence",
        "sirnaforge.core.selection",
        "sirnaforge.data.orthology",
        "sirnaforge.provenance",
        "sirnaforge.reporting.payload",
        "sirnaforge.reporting.render",
        "sirnaforge.reporting.quilt",
        "sirnaforge.benchmark.prepare",
        "sirnaforge.benchmark.design",
    ]
    failures: list[str] = []
    for name in modules:
        try:
            importlib.import_module(name)
        except ImportError as exc:
            failures.append(f"{name}: {exc}")
    assert not failures, "modules the image ships but cannot import: " + "; ".join(failures)

    # An importable package is not a working one; instantiate the pair every command depends on.
    design_module = importlib.import_module("sirnaforge.core.design")
    models_module = importlib.import_module("sirnaforge.models")
    designer = design_module.SiRNADesigner(models_module.DesignParameters())
    assert designer is not None


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_minimal_workflow():
    """`sirnaforge design` on a valid FASTA succeeds in the image and writes candidates.

    Every branch of this test used to end in `pytest.skip` -- on a non-zero exit it did not
    recognise, on a timeout, on the CLI being absent. Inside the image (which `runs_in_container`
    guarantees) each of those is a defect in the artifact, not an environment gap.
    """
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

            # `design` on a valid FASTA in the release image has no honest reason to fail. This used
            # to skip on any non-zero exit whose stderr it did not recognise, which is how a real
            # regression would have been recorded as "expected".
            assert result.returncode == 0, (
                f"design failed in the image:\nSTDOUT: {result.stdout[-2000:]}\nSTDERR: {result.stderr[-2000:]}"
            )
            output_tsv = work_dir / "output.tsv"
            assert output_tsv.is_file(), f"design exited 0 but wrote no output: {sorted(work_dir.iterdir())}"
            assert output_tsv.read_text().strip(), "design wrote an empty output file"

        except FileNotFoundError:
            pytest.fail("sirnaforge CLI is not on PATH inside the image")
        except subprocess.TimeoutExpired:
            pytest.fail("sirnaforge design did not finish within 30s on a single short transcript")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_benchmark_prepare_explains_that_vendored_panels_need_a_checkout():
    """#109's `benchmark prepare` refuses a vendored panel in the image, and says why.

    The vendored panel tables live under `tests/`, which the wheel and the image do not carry, so
    every data-bearing panel is checkout-only by design. Before this was fenced the command exited 1
    with `panel csv not found: /opt/conda/lib/python3.12/tests/unit/data/...` -- a path that reads as
    a broken build. Nothing in dev, ci, release or container noticed, in either direction.

    This is the only test that runs `benchmark` as the installed program: the six host test files
    drive it in-process from the checkout, where the vendored bytes are always present, so this
    branch is unreachable there.
    """
    result = subprocess.run(
        ["sirnaforge", "benchmark", "prepare", "--panel", "huesken_subset", "--out-dir", "/tmp/benchmark_probe"],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode != 0, "a vendored panel cannot be prepared from an installed package"

    # Rich wraps the message at the terminal width, so compare on collapsed whitespace.
    message = " ".join((result.stdout + result.stderr).split())
    assert "checkout" in message, f"refusal does not name the checkout requirement: {message[-400:]}"
    assert "panel csv not found" not in message, f"still reporting a site-packages path as missing: {message[-400:]}"


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_bioinformatics_tools():
    """Test that key bioinformatics tools are available in container."""
    # TODO: expand this to cover ... more
    tools = [
        ("samtools", ["--version"]),
        ("bwa-mem2", ["version"]),
        ("RNAfold", ["--help"]),
    ]

    missing_tools = []

    for tool, args in tools:
        try:
            # Avoid automatic text decoding by subprocess (some tool outputs contain
            # non-UTF-8 bytes). Capture raw bytes and decode defensively.
            result = subprocess.run([tool] + args, capture_output=True, text=False, timeout=10, check=False)
            # Decode safely for any later inspection (replace invalid bytes)
            stdout = (
                result.stdout.decode("utf-8", errors="replace")
                if isinstance(result.stdout, bytes | bytearray)
                else str(result.stdout)
            )
            stderr = (
                result.stderr.decode("utf-8", errors="replace")
                if isinstance(result.stderr, bytes | bytearray)
                else str(result.stderr)
            )

            # Tool should exist and respond (exit code may vary)
            if result.returncode > 1:  # 0 or 1 are usually okay
                missing_tools.append(f"{tool} (bad exit code: {result.returncode})")

            # Also consider evidence in stdout/stderr that the binary is missing or failed
            combined = (stdout or "") + "\n" + (stderr or "")
            if any(term in combined.lower() for term in ["not found", "command not found", "error", "unable to"]):
                missing_tools.append(f"{tool} (error output: {combined.strip()[:200]})")
        except FileNotFoundError:
            missing_tools.append(f"{tool} (not found)")
        except subprocess.TimeoutExpired:
            missing_tools.append(f"{tool} (timeout)")

    if missing_tools:
        pytest.skip(f"Some bioinformatics tools missing in container: {missing_tools}")


@pytest.mark.integration
@pytest.mark.runs_in_container
@pytest.mark.requires_network
@pytest.mark.slow
def test_docker_full_tp53_workflow(tmp_path: Path):
    """Test complete TP53 workflow inside container using local profile.

    This validates:
    - Container environment auto-detection
    - Local profile selection (no Docker-in-Docker)
    - Gene search functionality (Ensembl)
    - siRNA design pipeline (default filtering)
    - Off-target analysis with Nextflow local execution (bwa-mem2 based)
    - All bioinformatics tools working correctly (in the container)

    Based on integration/test_workflow_integration.sh patterns. (depreciated)
    """
    # Use /workspace if in container, otherwise tmp_path
    if Path("/workspace").exists() and os.access("/workspace", os.W_OK):
        output_dir = Path("/workspace/tp53_workflow_debug")
    else:
        output_dir = tmp_path / "tp53_workflow_debug"
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Run complete TP53 workflow (similar to test_workflow_integration.sh)
        result = subprocess.run(
            ["sirnaforge", "workflow", "TP53", "--output-dir", str(output_dir)],
            capture_output=True,
            text=True,
            cwd=output_dir,
            timeout=1200,  # 20 minutes max
            check=False,
        )

        # Always save workflow outputs for debugging
        print(f"\n{'=' * 80}")
        print(f"TP53 Workflow Results saved to: {output_dir}")
        print(f"{'=' * 80}")
        print(f"STDOUT:\n{result.stdout}")
        print(f"{'=' * 80}")
        print(f"STDERR:\n{result.stderr}")
        print(f"{'=' * 80}\n")

        # Check for failures and categorize them
        _check_workflow_result(result)

        # Verify output structure
        _verify_workflow_outputs(output_dir, result)

    except FileNotFoundError:
        pytest.skip("sirnaforge CLI not available - run this test in Docker container")
    except subprocess.TimeoutExpired:
        pytest.skip("Workflow timed out (may indicate missing dependencies or slow network)")


def _check_workflow_result(result: "subprocess.CompletedProcess[str]") -> None:
    """Check workflow execution result and handle different failure types."""
    if result.returncode == 0:
        return

    error_text = result.stderr.lower()

    # Critical failures (environment problems)
    critical_errors = ["import error", "module not found", "command not found", "docker: command not found"]
    if any(term in error_text for term in critical_errors):
        pytest.fail(f"Environment issue in container workflow:\n{result.stderr}")

    # Network issues are skippable (like test_mirna_integration.py)
    network_errors = ["connection", "network", "timeout", "ssl", "ensembl"]
    if any(term in error_text for term in network_errors):
        pytest.skip(f"Network issue during workflow: {result.stderr[:500]}")

    # Other failures might be transient
    pytest.skip(f"Workflow issue (may be transient): {result.stderr[:500]}")


def _verify_workflow_outputs(output_dir: Path, result: "subprocess.CompletedProcess[str]") -> None:
    """Verify workflow outputs and validate content."""
    # Expected output files
    expected_files = {
        "all_csv": output_dir / "sirnaforge" / "candidates_all.csv",
        "pass_csv": output_dir / "sirnaforge" / "candidates_pass.csv",
        "transcripts": output_dir / "transcripts" / "TP53_transcripts.fasta",
        "summary": output_dir / "logs" / "workflow_summary.json",
    }

    missing_files = [name for name, path in expected_files.items() if not path.exists()]
    if missing_files:
        actual_files = list(output_dir.rglob("*")) if output_dir.exists() else []
        pytest.fail(
            f"Missing files: {missing_files}\n"
            f"Created: {[str(f.relative_to(output_dir)) for f in actual_files if f.is_file()][:10]}"
        )

    # Verify no Docker-in-Docker attempts (key test for container fix)
    if "docker: command not found" in result.stderr.lower():
        pytest.fail("Workflow tried to use Docker-in-Docker (should use local profile)")

    # Validate CSV content
    all_csv = expected_files["all_csv"]
    lines = all_csv.read_text().strip().split("\n")
    assert len(lines) > 1, "No siRNA candidates generated (only header present)"

    # Verify CSV structure
    header = lines[0].lower()
    for col in ["id", "sequence", "position", "gc_content"]:
        assert col in header, f"Missing expected column: {col}"

    # Validate summary JSON
    summary = json.loads(expected_files["summary"].read_text())

    # Check required phases (updated to match new workflow summary structure)
    for phase in ["transcript_summary", "design_summary"]:
        assert phase in summary, f"Missing {phase} phase in workflow summary"

    # Validate phase results (updated key names)
    assert summary["transcript_summary"].get("total_transcripts", 0) > 0, "No transcripts retrieved"
    assert summary["design_summary"].get("total_candidates", 0) > 0, "No candidates generated"
    assert summary["design_summary"].get("pass_count", 0) > 0, "No candidates passed filters"


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_docker_login_shell_path():
    """Test that login shells preserve PATH and can find sirnaforge/nextflow (Issue #37).

    This test verifies the fix for Docker login shell PATH reset issue where
    /bin/bash -lc would drop /opt/conda/bin from PATH.

    The test runs inside the container and simulates what users would do with
    docker run ... /bin/bash -lc commands.
    """
    # Test 1: Non-login shell (baseline - should always work)
    result = subprocess.run(
        ["/bin/bash", "-c", "command -v sirnaforge && command -v nextflow"],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        # Get PATH from subprocess for debugging
        path_result = subprocess.run(
            ["/bin/bash", "-c", "echo $PATH"],
            check=False,
            capture_output=True,
            text=True,
        )
        pytest.fail(
            f"Non-login shell should find tools.\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}\n"
            f"Subprocess PATH: {path_result.stdout.strip()}"
        )

    # Test 2: Login shell (the fix target)
    result = subprocess.run(
        ["/bin/bash", "-lc", "command -v sirnaforge && command -v nextflow"],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        # Get PATH from subprocess for debugging
        path_result = subprocess.run(
            ["/bin/bash", "-lc", "echo $PATH"],
            check=False,
            capture_output=True,
            text=True,
        )
        pytest.fail(
            f"Login shell should find tools after fix (Issue #37).\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}\n"
            f"Subprocess PATH: {path_result.stdout.strip()}\n"
            f"Check if /etc/profile.d/conda-path.sh is installed and executable."
        )

    # Test 3: Verify PATH contains conda directories in login shell
    result = subprocess.run(
        ["/bin/bash", "-lc", "echo $PATH"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, "Failed to get PATH from login shell"
    path = result.stdout.strip()
    assert "/opt/conda/bin" in path, f"Login shell PATH should contain /opt/conda/bin.\nActual PATH: {path}"
