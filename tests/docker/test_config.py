"""Test configuration profiles for different resource scenarios.

This module provides different test configurations optimized for various
environments and resource constraints.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class TestProfile:
    """Configuration for different test scenarios."""

    name: str
    description: str
    pytest_args: list[str]
    docker_cpus: float
    docker_memory: str
    docker_memory_swap: str
    max_fail: int
    timeout: Optional[int] = None


# Predefined test profiles
TEST_PROFILES = {
    "lightweight": TestProfile(
        name="lightweight",
        description="Fast, minimal resource tests for development",
        pytest_args=["-q", "-n", "1", "-m", "lightweight or docker", "--maxfail=3"],
        docker_cpus=1.0,
        docker_memory="1g",
        docker_memory_swap="2g",
        max_fail=3,
        timeout=300,
    ),
    "development": TestProfile(
        name="development",
        description="Balanced tests for development workflow",
        pytest_args=["-v", "-n", "1", "--maxfail=5"],
        docker_cpus=2.0,
        docker_memory="4g",
        docker_memory_swap="6g",
        max_fail=5,
        timeout=600,
    ),
    "ci": TestProfile(
        name="ci",
        description="Comprehensive tests for CI/CD pipelines",
        pytest_args=["-v", "-n", "2", "--cov=sirnaforge", "--cov-report=term-missing"],
        docker_cpus=4.0,
        docker_memory="8g",
        docker_memory_swap="12g",
        max_fail=10,
        timeout=1800,
    ),
    "fast": TestProfile(
        name="fast",
        description="Quick tests excluding slow operations",
        pytest_args=["-q", "-n", "1", "-m", "not slow", "--maxfail=3"],
        docker_cpus=1.0,
        docker_memory="2g",
        docker_memory_swap="3g",
        max_fail=3,
        timeout=300,
    ),
    "minimal": TestProfile(
        name="minimal",
        description="Minimal smoke tests only",
        pytest_args=["-q", "-n", "1", "-k", "test_docker_cli_help or test_docker_version", "--maxfail=1"],
        docker_cpus=0.5,
        docker_memory="512m",
        docker_memory_swap="1g",
        max_fail=1,
        timeout=120,
    ),
}


def get_profile(name: str) -> TestProfile:
    """Get a test profile by name."""
    if name not in TEST_PROFILES:
        available = list(TEST_PROFILES.keys())
        raise ValueError(f"Unknown profile '{name}'. Available: {available}")
    return TEST_PROFILES[name]


def list_profiles() -> dict[str, str]:
    """List all available test profiles with descriptions."""
    return {name: profile.description for name, profile in TEST_PROFILES.items()}


def get_docker_run_command(profile: TestProfile, image_name: str = "sirnaforge:latest") -> list[str]:
    """Generate Docker run command for a profile."""
    return [
        "docker",
        "run",
        "--rm",
        "--cpus",
        str(profile.docker_cpus),
        "--memory",
        profile.docker_memory,
        "--memory-swap",
        profile.docker_memory_swap,
        "-v",
        "$(pwd):/workspace",
        "-w",
        "/workspace",
        image_name,
        "uv",
        "run",
        "pytest",
    ] + profile.pytest_args


# Quick access functions
def lightweight_profile() -> TestProfile:
    """Get lightweight profile for development."""
    return get_profile("lightweight")


def development_profile() -> TestProfile:
    """Get development profile."""
    return get_profile("development")


def ci_profile() -> TestProfile:
    """Get CI profile."""
    return get_profile("ci")
