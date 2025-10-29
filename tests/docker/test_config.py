"""Test configuration profiles for Docker-based testing.

Provides two simple profiles:
- smoke: Quick validation with fixed 2GB RAM limit
- production: Full tests with auto-detected RAM based on system capacity
"""

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


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
    timeout: Optional[int] = None


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
    return TEST_PROFILES[name]


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
