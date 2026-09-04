"""Shared pytest fixtures for all test modules."""

import os
import shutil
import socket
import ssl
import subprocess
from functools import lru_cache
from pathlib import Path

import pytest

from sirnaforge.core.off_target import build_bwa_index

_NEXTFLOW_PROBE_TIMEOUT = 30.0

# Set by the Makefile's `docker run` invocations. Preferred over sniffing the
# filesystem because only Docker (/.dockerenv) and Podman (/run/.containerenv)
# leave a marker file -- containerd, CRI-O and Kubernetes leave none, and a
# missed detection would skip the entire container tier while still exiting 0.
_CONTAINER_ENV_FLAG = "SIRNAFORGE_IN_CONTAINER"
_CONTAINER_MARKER_FILES = (Path("/.dockerenv"), Path("/run/.containerenv"))

# Every `requires_network` test in this suite ultimately talks to Ensembl REST.
_NETWORK_PROBE_HOST = "rest.ensembl.org"
_NETWORK_PROBE_PORT = 443
_NETWORK_PROBE_TIMEOUT = 5.0


def _running_in_container() -> bool:
    """Whether this interpreter is running inside the sirnaforge container image."""
    if os.environ.get(_CONTAINER_ENV_FLAG, "").strip().lower() in {"1", "true", "yes"}:
        return True
    return any(marker.exists() for marker in _CONTAINER_MARKER_FILES)


@lru_cache(maxsize=1)
def _nextflow_runnable() -> bool:
    """Whether Nextflow can actually execute, not merely whether the launcher exists.

    ``shutil.which`` is not enough: the launcher is a small script that fetches the
    framework JAR into NXF_HOME on first use, so on an image without that JAR baked
    in -- or on a network that cannot reach nextflow.io -- ``which`` succeeds while
    every real invocation exits non-zero. Mirrors the product's own check in
    ``NextflowRunner.validate_installation`` so tests skip exactly when the pipeline
    would report ``nextflow_unavailable``.
    """
    launcher = shutil.which("nextflow")
    if launcher is None:
        return False
    try:
        return (
            subprocess.run(  # noqa: S603 - fixed argv, resolved executable
                [launcher, "-version"],
                capture_output=True,
                timeout=_NEXTFLOW_PROBE_TIMEOUT,
                check=False,
            ).returncode
            == 0
        )
    except (OSError, subprocess.SubprocessError):
        return False


@lru_cache(maxsize=1)
def _network_available() -> bool:
    """Whether a *verified* TLS handshake with Ensembl REST is possible.

    Deliberately performs the handshake with ``ssl.create_default_context()`` --
    the same trust store aiohttp uses -- rather than a bare TCP connect. A plain
    socket connect succeeds behind a TLS-intercepting corporate proxy even though
    every subsequent HTTPS request fails certificate verification, so a TCP-only
    probe would report "online" and let the tests fail confusingly.

    Only transport reachability is probed; no sirnaforge code path is exercised,
    so a genuine client regression still fails the test rather than skipping it.

    Context construction is inside the ``try`` on purpose: it reads SSL_CERT_FILE,
    so a stale value (the proxy-CA workaround in docs/developer/testing_guide.md
    points at /tmp, which gets reaped) must skip rather than error in setup.
    """
    try:
        context = ssl.create_default_context()
        with (
            socket.create_connection(
                (_NETWORK_PROBE_HOST, _NETWORK_PROBE_PORT), timeout=_NETWORK_PROBE_TIMEOUT
            ) as raw_sock,
            context.wrap_socket(raw_sock, server_hostname=_NETWORK_PROBE_HOST),
        ):
            return True
    except OSError:  # covers socket.timeout, gaierror and ssl.SSLError
        return False


@pytest.fixture
def toy_genome_path():
    """Path to toy transcriptome database for testing off-target analysis.

    This provides a small transcriptome database suitable for testing
    BWA-based off-target analysis against cDNA sequences (resource-intensive mode).
    Note: This is transcriptome (cDNA) alignment, not genomic DNA.
    """
    test_data_dir = Path(__file__).parent / "unit" / "data"
    genome_path = test_data_dir / "toy_transcriptome_db.fasta"

    if not genome_path.exists():
        pytest.skip(f"Toy genome database not found: {genome_path}")

    return genome_path


@pytest.fixture(scope="session")
def toy_genome_index_prefix(tmp_path_factory):
    """Build a BWA index for the toy transcriptome database once per test session."""
    test_data_dir = Path(__file__).parent / "unit" / "data"
    genome_fasta = test_data_dir / "toy_transcriptome_db.fasta"

    if not genome_fasta.exists():
        pytest.skip(f"Toy genome database not found: {genome_fasta}")

    index_dir = tmp_path_factory.mktemp("toy_transcriptome_index")
    index_prefix = Path(index_dir) / "toy_transcriptome"

    if not all(index_prefix.with_suffix(suffix).exists() for suffix in (".amb", ".ann", ".bwt.2bit.64", ".pac", ".sa")):
        try:
            build_bwa_index(genome_fasta, index_prefix)
        except Exception as exc:  # pragma: no cover - requires bwa tooling
            pytest.skip(f"Unable to build toy genome index: {exc}")

    return index_prefix


@pytest.fixture
def realistic_transcripts_fasta():
    """Path to the packaged realistic transcript FASTA used across tests."""
    data_path = Path(__file__).parent / "unit" / "data" / "realistic_transcripts.fasta"

    if not data_path.exists():
        pytest.skip(f"Realistic transcript FASTA not found: {data_path}")

    return data_path


@pytest.fixture
def genome_config_for_nextflow(toy_genome_path):
    """Nextflow-compatible transcriptome configuration for tests.

    Enables the optional resource-intensive transcriptome off-target analysis
    (BWA alignment against cDNA sequences). Without this, only lightweight
    miRNA seed match analysis runs.

    Returns a dict with genome parameters that can be passed to Nextflow workflows.
    Note: 'genome' in parameter names refers to the Nextflow convention, not genomic DNA.
    """
    return {
        "--genome_species": "test_species",
        "--genome_fastas": f"test_species:{toy_genome_path}",
    }


@pytest.fixture
def mirna_only_config_for_nextflow():
    """Nextflow configuration for lightweight miRNA-only analysis.

    This skips the resource-intensive genome/transcriptome analysis (8-60GB RAM)
    and only runs lightweight miRNA seed match comparison (<1GB RAM).

    Returns a dict with minimal parameters for miRNA-only mode.
    """
    return {}


@pytest.fixture
def nextflow_test_work_dir(tmp_path):
    """Persistent work directory for Nextflow tests.

    Returns a path that exists outside Docker container's /tmp to ensure
    artifacts persist for inspection on host.
    """
    work_dir = tmp_path / "nextflow_work"
    work_dir.mkdir(exist_ok=True, parents=True)
    return work_dir


def pytest_runtest_setup(item):
    """Enforce the requirement markers so unmet requirements skip, not fail.

    `runs_in_container`, `requires_nextflow` and `requires_network` describe an
    environment a test cannot supply for itself. Without this gate the tier targets
    stay green only because they filter the markers out (`make test-ci`,
    `make test-release-host`), while `make test` runs them anywhere and reports
    environment gaps as assertion failures.

    Each probe asks whether the dependency actually *works*, not whether it is
    nominally present -- a resolvable hostname or an executable on PATH both lie.
    """
    marker_names = {mark.name for mark in item.iter_markers()}

    if "runs_in_container" in marker_names and not _running_in_container():
        pytest.skip("runs_in_container test: only valid inside the sirnaforge image - use 'make docker-test'")

    if "requires_nextflow" in marker_names and not _nextflow_runnable():
        pytest.skip(
            "requires_nextflow test: `nextflow -version` does not succeed here (launcher present is not enough)"
        )

    if "requires_network" in marker_names and not _network_available():
        pytest.skip(
            f"requires_network test: no verified TLS route to {_NETWORK_PROBE_HOST} "
            "(offline, or a TLS-intercepting proxy whose root CA is missing from Python's "
            "trust store - see docs/developer/testing_guide.md)"
        )


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config, items):
    """Auto-assign tier markers based on test type.

    Tier hierarchy:
    - dev: Fast unit tests for development iteration
    - ci: Smoke tests and container validation for CI/CD
    - release: Integration tests and heavy workloads for release validation
    """
    for item in items:
        marker_names = {mark.name for mark in item.iter_markers()}

        # Smoke tests → CI tier
        if "smoke" in marker_names:
            if "ci" not in marker_names:
                item.add_marker(pytest.mark.ci)
            continue

        # Container tests → CI tier (validate Docker image)
        if "runs_in_container" in marker_names:
            if "ci" not in marker_names:
                item.add_marker(pytest.mark.ci)
            # Container integration tests also get release tier
            if "integration" in marker_names and "release" not in marker_names:
                item.add_marker(pytest.mark.release)
            continue

        # Heavy workloads → Release tier
        release_triggers = {
            "integration",
            "requires_docker",
            "requires_nextflow",
            "requires_tools",
            "slow",
        }
        if marker_names & release_triggers:
            if "release" not in marker_names:
                item.add_marker(pytest.mark.release)
            continue

        # Unit tests → Dev tier
        if "unit" in marker_names:
            if "dev" not in marker_names:
                item.add_marker(pytest.mark.dev)
            continue

        # Default: untagged tests → Dev tier
        if "dev" not in marker_names and "release" not in marker_names and "ci" not in marker_names:
            item.add_marker(pytest.mark.dev)
