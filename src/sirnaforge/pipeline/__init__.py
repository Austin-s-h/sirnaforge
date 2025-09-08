"""
siRNAforge Pipeline Module

Integrated Nextflow pipeline for comprehensive off-target analysis.
This module provides a Python interface to the embedded Nextflow workflows.
"""

from .nextflow.config import NextflowConfig
from .nextflow.runner import NextflowRunner
from .resources import get_resource_path, get_test_data_path

__all__ = [
    "NextflowRunner",
    "NextflowConfig",
    "get_test_data_path",
    "get_resource_path",
]
