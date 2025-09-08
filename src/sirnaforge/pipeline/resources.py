"""
Resource Management for siRNAforge Pipeline

This module provides utilities for managing test data and pipeline resources.
"""

import importlib.resources as pkg_resources
from pathlib import Path

from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


def get_resource_path(resource_name: str) -> Path:
    """
    Get path to a pipeline resource file.

    Args:
        resource_name: Name of the resource file

    Returns:
        Path to the resource file

    Raises:
        FileNotFoundError: If resource is not found
    """
    try:
        # Try to get from package resources first
        resources = pkg_resources.files("sirnaforge.pipeline.resources")
        resource_path = resources / resource_name
        try:
            # Check if the resource exists using the traversable path
            if hasattr(resource_path, "exists") and resource_path.exists():
                return Path(str(resource_path))
        except (AttributeError, OSError):
            # Fallback for different importlib versions
            pass
    except (ImportError, AttributeError):
        pass

    # Fallback to file-based lookup
    package_dir = Path(__file__).parent
    resource_path = package_dir / "resources" / resource_name

    if not resource_path.exists():
        raise FileNotFoundError(f"Resource not found: {resource_name}")

    return resource_path


def get_test_data_path(filename: str) -> Path:
    """
    Get path to test data file.

    Args:
        filename: Name of the test data file

    Returns:
        Path to the test data file

    Raises:
        FileNotFoundError: If test data file is not found
    """
    return get_resource_path(f"test_data/{filename}")


def get_workflow_path(workflow_name: str = "main.nf") -> Path:
    """
    Get path to embedded Nextflow workflow.

    Args:
        workflow_name: Name of the workflow file (default: main.nf)

    Returns:
        Path to the workflow file

    Raises:
        FileNotFoundError: If workflow is not found
    """
    package_dir = Path(__file__).parent
    workflow_path = package_dir / "nextflow" / "workflows" / workflow_name

    if not workflow_path.exists():
        raise FileNotFoundError(f"Workflow not found: {workflow_name}")

    return workflow_path


def list_test_data() -> list[str]:
    """
    List available test data files.

    Returns:
        List of test data filenames
    """
    try:
        test_data_dir = get_resource_path("test_data")
        return [f.name for f in test_data_dir.glob("*") if f.is_file()]
    except FileNotFoundError:
        logger.warning("Test data directory not found")
        return []


def validate_test_data() -> dict[str, bool]:
    """
    Validate that required test data files are present.

    Returns:
        Dictionary mapping filenames to availability status
    """
    required_files = [
        "test_candidates.fasta",
        "test_transcriptome.fasta",
        "test_mirna_seeds.fasta",
        "genomes.yaml",
    ]

    status = {}
    for filename in required_files:
        try:
            get_test_data_path(filename)
            status[filename] = True
        except FileNotFoundError:
            status[filename] = False
            logger.warning(f"Missing test data file: {filename}")

    return status


class ResourceManager:
    """Manage pipeline resources and test data."""

    def __init__(self) -> None:
        """Initialize resource manager."""
        self.package_dir = Path(__file__).parent
        self.resources_dir = self.package_dir / "resources"
        self.test_data_dir = self.resources_dir / "test_data"

    def ensure_test_data(self) -> bool:
        """
        Ensure test data directory exists and contains required files.

        Returns:
            True if all required test data is available
        """
        if not self.test_data_dir.exists():
            logger.warning(f"Test data directory not found: {self.test_data_dir}")
            return False

        validation = validate_test_data()
        missing_files = [f for f, available in validation.items() if not available]

        if missing_files:
            logger.warning(f"Missing test data files: {missing_files}")
            return False

        logger.info("All required test data files are available")
        return True

    def get_test_config(self) -> dict[str, str]:
        """
        Get test configuration with paths to test data.

        Returns:
            Dictionary containing test data paths
        """
        try:
            return {
                "test_candidates": str(get_test_data_path("test_candidates.fasta")),
                "test_transcriptome": str(get_test_data_path("test_transcriptome.fasta")),
                "test_mirna_seeds": str(get_test_data_path("test_mirna_seeds.fasta")),
                "genomes_config": str(get_test_data_path("genomes.yaml")),
            }
        except FileNotFoundError as e:
            logger.error(f"Failed to create test config: {e}")
            return {}
