"""Shared pytest fixtures for all test modules."""

from pathlib import Path

import pytest

from sirnaforge.core.off_target import build_bwa_index


@pytest.fixture
def toy_genome_path():
    """Path to toy genome database for testing off-target analysis.

    This provides a small transcriptome database suitable for testing
    genome/transcriptome off-target analysis (resource-intensive mode).
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
def genome_config_for_nextflow(toy_genome_path):
    """Nextflow-compatible genome configuration for tests.

    Enables the optional resource-intensive genome/transcriptome off-target analysis.
    Without this, only lightweight miRNA seed match analysis runs.

    Returns a dict with genome parameters that can be passed to Nextflow workflows.
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
