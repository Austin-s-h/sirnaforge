"""Enhanced Docker test profiles and configurations.

This module provides different testing profiles for various Docker environments
and use cases, making it easier to test siRNAforge in different scenarios.
"""

import os
import subprocess
from dataclasses import dataclass
from typing import Optional

import pytest


@dataclass
class DockerTestProfile:
    """Configuration for different Docker testing scenarios."""

    name: str
    description: str
    required_tools: list[str]
    test_markers: list[str]
    environment_vars: dict[str, str]
    timeout_multiplier: float = 1.0

    def should_skip(self) -> Optional[str]:
        """Check if this profile should be skipped based on available tools."""
        missing_tools = []
        for tool in self.required_tools:
            try:
                result = subprocess.run(
                    [tool, "--version"] if tool != "nextflow" else [tool, "-version"],
                    capture_output=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=5,
                    check=False,  # We handle the return code manually
                )
                if result.returncode != 0:
                    missing_tools.append(tool)
            except (FileNotFoundError, subprocess.TimeoutExpired):
                missing_tools.append(tool)

        if missing_tools:
            return f"Missing required tools: {', '.join(missing_tools)}"
        return None

    def get_pytest_args(self) -> list[str]:
        """Get pytest arguments for this profile."""
        args = ["-v"]
        if self.test_markers:
            args.extend(["-m", " or ".join(self.test_markers)])
        return args


# Define comprehensive test profiles
TEST_PROFILES = {
    "minimal": DockerTestProfile(
        name="minimal",
        description="Basic functionality tests - no external tools required",
        required_tools=[],
        test_markers=["not docker and not integration and not slow"],
        environment_vars={},
    ),
    "bioinformatics": DockerTestProfile(
        name="bioinformatics",
        description="Tests requiring bioinformatics tools (samtools, bwa-mem2, RNAfold)",
        required_tools=["samtools", "bwa-mem2", "RNAfold"],
        test_markers=["docker"],
        environment_vars={},
    ),
    "workflow": DockerTestProfile(
        name="workflow",
        description="Complete workflow tests including Nextflow pipelines",
        required_tools=["nextflow", "docker"],
        test_markers=["integration", "pipeline"],
        environment_vars={"NXF_WORK": "/tmp/nextflow_work", "NXF_TEMP": "/tmp/nextflow_temp"},
        timeout_multiplier=2.0,
    ),
    "performance": DockerTestProfile(
        name="performance",
        description="Performance and stress tests",
        required_tools=["samtools", "bwa-mem2"],
        test_markers=["performance", "slow"],
        environment_vars={"SIRNAFORGE_PERF_MODE": "1"},
        timeout_multiplier=3.0,
    ),
    "comprehensive": DockerTestProfile(
        name="comprehensive",
        description="All tests including slow and integration tests",
        required_tools=["samtools", "bwa-mem2", "RNAfold", "nextflow"],
        test_markers=[],  # Run all tests
        environment_vars={"SIRNAFORGE_TEST_ALL": "1", "NXF_WORK": "/tmp/nextflow_work"},
        timeout_multiplier=5.0,
    ),
}


class DockerTestRunner:
    """Runner for executing tests in different Docker profiles."""

    def __init__(self, profile: DockerTestProfile):
        """Initialize Docker test runner with profile configuration."""
        self.profile = profile

    def run_tests(self) -> subprocess.CompletedProcess:
        """Run tests for this profile."""
        skip_reason = self.profile.should_skip()
        if skip_reason:
            pytest.skip(f"Profile '{self.profile.name}' skipped: {skip_reason}")

        # Set environment variables
        env = dict(os.environ)
        env.update(self.profile.environment_vars)

        # Build pytest command
        cmd = ["uv", "run", "pytest"] + self.profile.get_pytest_args()

        # Apply timeout multiplier
        timeout = 300 * self.profile.timeout_multiplier  # 5 minutes base

        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
            check=False,  # We handle the return code in the calling code
        )


# Example test configurations for different scenarios
EXAMPLE_CONFIGURATIONS = {
    "gene_targeting": {
        "description": "Test siRNA design for specific gene targeting",
        "input_files": ["examples/sample_transcripts.fasta"],
        "parameters": {"gene": "TP53", "top_n": 10, "genome_species": ["human"]},
    },
    "multi_species": {
        "description": "Test off-target analysis across multiple species",
        "input_files": ["examples/sample_transcripts.fasta"],
        "parameters": {"genome_species": ["human", "mouse", "rat"], "top_n_offtarget": 20},
    },
    "large_dataset": {
        "description": "Performance test with large transcript datasets",
        "input_files": ["examples/large_transcripts.fasta"],
        "parameters": {"top_n": 100, "batch_size": 50},
    },
    "custom_scoring": {
        "description": "Test with custom scoring weights",
        "input_files": ["examples/sample_transcripts.fasta"],
        "parameters": {
            "scoring_weights": {"off_target": 0.4, "accessibility": 0.3, "asymmetry": 0.2, "gc_content": 0.1}
        },
    },
}


def get_available_profiles() -> dict[str, DockerTestProfile]:
    """Get profiles that can run in the current environment."""
    available = {}
    for name, profile in TEST_PROFILES.items():
        if profile.should_skip() is None:
            available[name] = profile
    return available


def run_profile_tests(profile_name: str) -> subprocess.CompletedProcess:
    """Run tests for a specific profile."""
    if profile_name not in TEST_PROFILES:
        raise ValueError(f"Unknown profile: {profile_name}")

    profile = TEST_PROFILES[profile_name]
    runner = DockerTestRunner(profile)
    return runner.run_tests()


# Convenience functions for common test scenarios
def run_minimal_functionality():
    """Test basic siRNAforge functionality."""
    return run_profile_tests("minimal")


def run_bioinformatics_pipeline():
    """Test bioinformatics tool integration."""
    return run_profile_tests("bioinformatics")


def run_complete_workflow():
    """Test complete siRNA design workflow."""
    return run_profile_tests("workflow")


def run_performance_baseline():
    """Test performance characteristics."""
    return run_profile_tests("performance")
