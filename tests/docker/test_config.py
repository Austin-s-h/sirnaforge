"""Test configuration profiles and Nextflow evaluation helpers.

Provides two Docker testing profiles (smoke + production) and exposes the
Nextflow configurations that the workflow can auto-detect so tests/scripts can
validate their availability before attempting container launches.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from sirnaforge.pipeline.nextflow.config import EnvironmentInfo, NextflowConfig


def _get_available_memory_gb() -> float:
    """Detect available system memory in GB.

    Checks cgroup limits (for containers) first, then falls back to /proc/meminfo.
    """
    # Check cgroup v2 memory limit
    cgroup_v2_path = Path("/sys/fs/cgroup/memory.max")
    if cgroup_v2_path.exists():
        try:
            with cgroup_v2_path.open() as f:
                limit = f.read().strip()
                if limit != "max":
                    return int(limit) / (1024**3)
        except (OSError, ValueError):
            pass

    # Check cgroup v1 memory limit
    cgroup_v1_path = Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")
    if cgroup_v1_path.exists():
        try:
            with cgroup_v1_path.open() as f:
                limit = int(f.read().strip())
                # Ignore unreasonably large limits (indicates no limit)
                if limit < (1 << 60):
                    return limit / (1024**3)
        except (OSError, ValueError):
            pass

    # Fall back to total system memory
    meminfo_path = Path("/proc/meminfo")
    try:
        with meminfo_path.open() as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return kb / (1024**2)
    except (OSError, ValueError):
        pass

    # Default fallback
    return 8.0


def _calculate_safe_memory_gb(total_gb: float, buffer_gb: float = 0.5, min_gb: int = 1) -> int:
    """Calculate safe memory limit with buffer."""
    safe_gb = math.floor(total_gb - buffer_gb)
    return max(safe_gb, min_gb)


ConfigFactory = Callable[[], NextflowConfig]


@dataclass
class NextflowConfigEvaluation:
    """Summary of a Nextflow configuration detected from the workflow."""

    name: str
    requested_profile: str
    recommended_profile: str
    docker_available: bool
    running_in_docker: bool
    max_cpus: int
    max_memory: str
    max_time: str
    env_summary: str
    docker_image: str | None = None
    environment: EnvironmentInfo | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize evaluation for logging or JSON reporting."""
        return {
            "name": self.name,
            "requested_profile": self.requested_profile,
            "recommended_profile": self.recommended_profile,
            "docker_available": self.docker_available,
            "running_in_docker": self.running_in_docker,
            "max_cpus": self.max_cpus,
            "max_memory": self.max_memory,
            "max_time": self.max_time,
            "docker_image": self.docker_image,
            "env_summary": self.env_summary,
        }


def _default_nextflow_factories() -> dict[str, ConfigFactory]:
    """Factories mirroring the workflow's available Nextflow configurations."""
    return {
        "auto": NextflowConfig.auto_configure,
        "testing": NextflowConfig.for_testing,
        "production": NextflowConfig.for_production,
    }


def _evaluate_factories(factories: Mapping[str, ConfigFactory]) -> dict[str, NextflowConfigEvaluation]:
    evaluations: dict[str, NextflowConfigEvaluation] = {}
    for name, factory in factories.items():
        config = factory()
        env_info = config.get_environment_info()
        evaluations[name] = NextflowConfigEvaluation(
            name=name,
            requested_profile=config.profile,
            recommended_profile=env_info.recommended_profile,
            docker_available=env_info.docker_available,
            running_in_docker=env_info.running_in_docker,
            max_cpus=config.max_cpus,
            max_memory=config.max_memory,
            max_time=config.max_time,
            docker_image=config.docker_image if env_info.recommended_profile == "docker" else None,
            env_summary=env_info.get_execution_summary(),
            environment=env_info,
        )
    return evaluations


@functools.lru_cache(maxsize=1)
def evaluate_nextflow_configurations() -> dict[str, NextflowConfigEvaluation]:
    """Evaluate workflow-backed Nextflow configurations once per session."""
    return _evaluate_factories(_default_nextflow_factories())


def get_nextflow_configuration(name: str) -> NextflowConfigEvaluation:
    """Fetch a specific evaluated Nextflow configuration by name."""
    configs = evaluate_nextflow_configurations()
    if name not in configs:
        available = ", ".join(sorted(configs.keys()))
        raise ValueError(f"Unknown Nextflow configuration '{name}'. Available: {available}")
    return configs[name]


@dataclass
class TestProfile:
    """Configuration for Docker-based test execution."""

    name: str
    description: str
    pytest_args: list[str]
    docker_cpus: float
    docker_memory: str
    docker_memory_swap: str
    max_fail: int
    timeout: int | None = None

    def available_nextflow_profiles(self) -> dict[str, NextflowConfigEvaluation]:
        """Expose workflow-driven Nextflow configuration evaluations."""
        return evaluate_nextflow_configurations()


# Prevent pytest from treating the dataclass as a test container
TestProfile.__test__ = False


# Predefined test profiles
TEST_PROFILES = {
    "smoke": TestProfile(
        name="smoke",
        description="Quick smoke tests with fixed 2GB RAM limit",
        pytest_args=["-q", "-n", "1", "-m", "smoke", "--maxfail=1"],
        docker_cpus=1.0,
        docker_memory="2g",
        docker_memory_swap="2g",
        max_fail=1,
        timeout=120,
    ),
    "production": TestProfile(
        name="production",
        description="Full test suite with auto-detected RAM limit",
        pytest_args=["-v", "--maxfail=10"],
        docker_cpus=0,  # 0 means no limit (use all available)
        docker_memory=f"{_calculate_safe_memory_gb(_get_available_memory_gb())}g",
        docker_memory_swap=f"{_calculate_safe_memory_gb(_get_available_memory_gb())}g",
        max_fail=10,
        timeout=3600,
    ),
}


def get_profile(name: str) -> TestProfile:
    """Get a test profile by name.

    Args:
        name: Profile name ('smoke' or 'production')

    Returns:
        TestProfile instance

    Raises:
        ValueError: If profile name is unknown
    """
    if name not in TEST_PROFILES:
        available = list(TEST_PROFILES.keys())
        raise ValueError(f"Unknown profile '{name}'. Available: {available}")
    profile = TEST_PROFILES[name]
    # Ensure Nextflow configurations are evaluated whenever a profile is requested
    profile.available_nextflow_profiles()
    return profile


def get_docker_run_command(profile: TestProfile, image_name: str = "sirnaforge:latest") -> list[str]:
    """Generate Docker run command for a profile.

    Args:
        profile: TestProfile instance
        image_name: Docker image name/tag

    Returns:
        Command list suitable for subprocess.run()
    """
    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        "$(pwd):/workspace",
        "-w",
        "/workspace",
    ]

    # Add CPU limit if specified
    if profile.docker_cpus > 0:
        cmd.extend(["--cpus", str(profile.docker_cpus)])

    # Add memory limits
    cmd.extend(
        [
            "--memory",
            profile.docker_memory,
            "--memory-swap",
            profile.docker_memory_swap,
        ]
    )

    # Add image and pytest command
    cmd.extend(
        [
            image_name,
            "uv",
            "run",
            "pytest",
        ]
        + profile.pytest_args
    )

    return cmd


# Convenience accessors
def smoke_profile() -> TestProfile:
    """Get smoke test profile (2GB RAM)."""
    return get_profile("smoke")


def production_profile() -> TestProfile:
    """Get production test profile (auto-detected RAM)."""
    return get_profile("production")
