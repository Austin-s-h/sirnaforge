"""Example configurations for different siRNAforge testing scenarios.

This file provides ready-to-use configurations for common testing patterns
and scenarios in siRNAforge development and CI/CD pipelines.
"""

from pathlib import Path
from typing import Any

# Base directory for test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "examples"

# Example test configurations
EXAMPLE_CONFIGS = {
    "basic_functionality": {
        "description": "Basic siRNA design functionality test",
        "input_file": str(TEST_DATA_DIR / "sample_transcripts.fasta"),
        "parameters": {
            "top_n": 5,
            "sirna_length": 21,
            "gc_min": 30.0,
            "gc_max": 52.0,
        },
        "expected_outputs": ["sirna_candidates.tsv"],
        "timeout": 60,
    },
    "gene_targeting": {
        "description": "Gene-specific siRNA targeting test",
        "gene_query": "TP53",
        "database": "ensembl",
        "parameters": {
            "top_n_candidates": 10,
            "top_n_offtarget": 5,
            "genome_species": ["human"],
            "gc_min": 35.0,
            "gc_max": 50.0,
        },
        "expected_outputs": [
            "transcripts/TP53_transcripts.fasta",
            "sirnaforge/sirna_candidates.tsv",
            "off_target/results/",
        ],
        "timeout": 300,
    },
    "multi_species_offtarget": {
        "description": "Multi-species off-target analysis test",
        "input_file": str(TEST_DATA_DIR / "sample_transcripts.fasta"),
        "parameters": {
            "top_n_candidates": 3,
            "top_n_offtarget": 2,
            "genome_species": ["human", "mouse", "rat"],
            "gc_min": 40.0,
            "gc_max": 60.0,
        },
        "expected_outputs": [
            "sirnaforge/sirna_candidates.tsv",
            "off_target/results/",
        ],
        "timeout": 600,
    },
    "large_dataset": {
        "description": "Large dataset performance test",
        "input_file": str(TEST_DATA_DIR / "large_transcripts.fasta"),
        "parameters": {
            "top_n_candidates": 50,
            "batch_size": 20,
            "gc_min": 30.0,
            "gc_max": 70.0,
        },
        "expected_outputs": ["sirnaforge/sirna_candidates.tsv"],
        "timeout": 1200,
    },
    "custom_scoring": {
        "description": "Custom scoring weights test",
        "input_file": str(TEST_DATA_DIR / "sample_transcripts.fasta"),
        "parameters": {
            "top_n": 8,
            "scoring_weights": {
                "off_target": 0.4,
                "accessibility": 0.3,
                "asymmetry": 0.2,
                "gc_content": 0.1,
            },
        },
        "expected_outputs": ["sirnaforge/sirna_candidates.tsv"],
        "timeout": 120,
    },
    "validation_only": {
        "description": "Input validation and FASTA checking test",
        "input_file": str(TEST_DATA_DIR / "sample_transcripts.fasta"),
        "command": "validate",
        "expected_outputs": [],
        "timeout": 30,
    },
}

# CI/CD pipeline configurations
CI_CONFIGS = {
    "unit_tests": {
        "description": "Fast unit tests for CI",
        "pytest_args": ["-v", "-m", "not slow and not integration"],
        "timeout": 300,
        "required_tools": [],
    },
    "integration_tests": {
        "description": "Integration tests with bioinformatics tools",
        "pytest_args": ["-v", "-m", "integration"],
        "timeout": 900,
        "required_tools": ["samtools", "bwa-mem2", "RNAfold"],
    },
    "docker_tests": {
        "description": "Docker container integration tests",
        "pytest_args": ["-v", "-m", "docker"],
        "timeout": 600,
        "required_tools": ["docker"],
    },
    "performance_tests": {
        "description": "Performance and load tests",
        "pytest_args": ["-v", "-m", "performance or slow"],
        "timeout": 1800,
        "required_tools": ["samtools", "bwa-mem2"],
    },
    "full_suite": {
        "description": "Complete test suite",
        "pytest_args": ["-v", "--cov=sirnaforge", "--cov-report=html"],
        "timeout": 3600,
        "required_tools": ["samtools", "bwa-mem2", "RNAfold", "nextflow"],
    },
}

# Docker-specific test configurations
DOCKER_CONFIGS = {
    "container_health": {
        "description": "Basic container health checks",
        "commands": [
            "sirnaforge --help",
            "sirnaforge version",
            "samtools --version",
            "bwa-mem2 version",
            "RNAfold --version",
        ],
        "timeout": 60,
    },
    "workflow_integration": {
        "description": "Complete workflow in container",
        "commands": ["sirnaforge workflow TP53 --output-dir /tmp/test_output --top-n 3 --genome-species human"],
        "timeout": 600,
        "expected_files": [
            "/tmp/test_output/sirnaforge/sirna_candidates.tsv",
            "/tmp/test_output/workflow_summary.json",
        ],
    },
    "tool_validation": {
        "description": "Validate all bioinformatics tools",
        "tools": {
            "samtools": ["--version"],
            "bwa-mem2": ["version"],
            "RNAfold": ["--version"],
            "nextflow": ["-version"],
        },
        "timeout": 120,
    },
}


def get_example_config(name: str) -> dict[str, Any]:
    """Get a specific example configuration."""
    if name not in EXAMPLE_CONFIGS:
        raise ValueError(f"Unknown example configuration: {name}")
    return EXAMPLE_CONFIGS[name]


def get_ci_config(name: str) -> dict[str, Any]:
    """Get a specific CI configuration."""
    if name not in CI_CONFIGS:
        raise ValueError(f"Unknown CI configuration: {name}")
    return CI_CONFIGS[name]


def get_docker_config(name: str) -> dict[str, Any]:
    """Get a specific Docker configuration."""
    if name not in DOCKER_CONFIGS:
        raise ValueError(f"Unknown Docker configuration: {name}")
    return DOCKER_CONFIGS[name]


def list_available_configs() -> dict[str, list[str]]:
    """List all available configurations by category."""
    return {
        "examples": list(EXAMPLE_CONFIGS.keys()),
        "ci": list(CI_CONFIGS.keys()),
        "docker": list(DOCKER_CONFIGS.keys()),
    }


# Quick test scenarios for development
QUICK_TESTS = {
    "smoke_test": {
        "description": "Quick smoke test to verify basic functionality",
        "commands": [
            "sirnaforge --help",
            "sirnaforge version",
            "sirnaforge validate examples/sample_transcripts.fasta",
        ],
        "timeout": 30,
    },
    "design_test": {
        "description": "Quick siRNA design test",
        "commands": ["sirnaforge design examples/sample_transcripts.fasta --output /tmp/quick_test.tsv --top-n 3"],
        "timeout": 60,
        "cleanup": ["/tmp/quick_test.tsv"],
    },
}
