"""Resource defaults must fit small CI runners and quota-limited containers."""

from pathlib import Path

import pytest

from sirnaforge.pipeline.nextflow.config import NextflowConfig

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("host_cpus", "affinity", "quota_v2", "quota_v1", "expected"),
    [
        (32, 32, "max 100000", None, 16),
        (4, 4, "max 100000", None, 4),
        (32, 4, "max 100000", None, 4),
        (32, 32, "400000 100000", None, 4),
        (32, 32, "250000 100000", None, 2),
        (32, 32, "50000 100000", None, 1),
        (32, 32, None, (200000, 100000), 2),
        (32, 32, "invalid", (-1, 100000), 16),
        (None, 1, None, None, 1),
    ],
)
def test_default_cpu_limit_respects_available_resources(monkeypatch, host_cpus, affinity, quota_v2, quota_v1, expected):
    """Affinity and quotas constrain defaults in both standard and production constructors."""
    monkeypatch.setattr("os.cpu_count", lambda: host_cpus)
    monkeypatch.setattr("os.sched_getaffinity", lambda _pid: set(range(affinity)), raising=False)
    files = {"/sys/fs/cgroup/cpu.max": quota_v2}
    if quota_v1:
        files.update(
            {
                "/sys/fs/cgroup/cpu/cpu.cfs_quota_us": str(quota_v1[0]),
                "/sys/fs/cgroup/cpu/cpu.cfs_period_us": str(quota_v1[1]),
            }
        )

    def read(path, *args, **kwargs):
        value = files.get(str(path))
        if value is None:
            raise FileNotFoundError(path)
        return value

    monkeypatch.setattr(Path, "read_text", read)
    assert NextflowConfig().max_cpus == expected
    assert NextflowConfig.for_production().max_cpus == expected


def test_cpu_detection_falls_back_when_affinity_is_unavailable(monkeypatch):
    """Non-Linux systems still use their detected processor count."""
    monkeypatch.setattr("os.cpu_count", lambda: 4)
    monkeypatch.delattr("os.sched_getaffinity", raising=False)
    monkeypatch.setattr(Path, "read_text", lambda *_args, **_kwargs: "max 100000")
    assert NextflowConfig().max_cpus == 4


def test_explicit_production_limits_are_preserved():
    """Explicit settings override defaults without duplicate keyword arguments."""
    config = NextflowConfig.for_production(max_cpus=2, max_memory="4.GB", profile="local")
    assert (config.max_cpus, config.max_memory, config.profile) == (2, "4.GB", "local")


def test_invalid_cpu_limit_is_rejected():
    """Reject a CPU ceiling that cannot schedule any task."""
    with pytest.raises(ValueError, match="max_cpus must be at least 1"):
        NextflowConfig(max_cpus=0)
