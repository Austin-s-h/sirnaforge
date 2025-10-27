"""Shared pytest fixtures for all test modules."""

from pathlib import Path

import pytest


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
    return {
        # No genome_fastas or genome_indices provided
        # This triggers miRNA-only analysis mode
    }


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
    """Normalize environment-tier markers for consistent workflows."""
    release_aliases = {"integration", "docker_integration", "pipeline", "slow", "nextflow"}

    for item in items:
        marker_names = {mark.name for mark in item.iter_markers()}

        # Ensure smoke implies ci
        if "smoke" in marker_names and "ci" not in marker_names:
            item.add_marker("ci")

        # Promote heavy workloads to release tier automatically
        if marker_names & release_aliases or "release" in marker_names:
            item.add_marker("release")
            continue

        # Skip automatically tagging explicit CI tests as dev
        if "ci" in marker_names:
            continue

        # Default bucket for normal local development runs
        if "dev" not in marker_names:
            item.add_marker("dev")
