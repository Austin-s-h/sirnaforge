"""
Nextflow Configuration Management

This module handles configuration for Nextflow workflows, including
Docker settings, resource management, and parameter validation.
"""

import subprocess
from pathlib import Path
from typing import Any, Optional

from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class NextflowConfig:
    """Configuration manager for Nextflow workflows."""

    def __init__(
        self,
        docker_image: str = "ghcr.io/austin-s-h/sirnaforge:latest",
        profile: str = "docker",
        work_dir: Optional[Path] = None,
        max_cpus: int = 16,
        max_memory: str = "128.GB",
        max_time: str = "240.h",
        **kwargs: Any,
    ) -> None:
        """
        Initialize Nextflow configuration.

        Args:
            docker_image: Docker container image to use
            profile: Nextflow profile (docker, singularity, conda, local)
            work_dir: Working directory for Nextflow execution
            max_cpus: Maximum CPU cores
            max_memory: Maximum memory allocation
            max_time: Maximum execution time
            **kwargs: Additional configuration parameters
        """
        self.docker_image = docker_image
        self.profile = profile
        self.work_dir = work_dir or Path.cwd() / "nextflow_work"
        self.max_cpus = max_cpus
        self.max_memory = max_memory
        self.max_time = max_time
        self.extra_params = kwargs

    def get_nextflow_args(
        self,
        input_file: Path,
        output_dir: Path,
        genome_species: list[str],
        additional_params: Optional[dict[str, Any]] = None,
    ) -> list[str]:
        """
        Generate Nextflow command arguments.

        Args:
            input_file: Input FASTA file path
            output_dir: Output directory
            genome_species: List of genome species to analyze
            additional_params: Additional parameters to pass

        Returns:
            List of command arguments for Nextflow
        """
        args = [
            "--input",
            str(input_file),
            "--outdir",
            str(output_dir),
            "--genome_species",
            ",".join(genome_species),
            "-profile",
            self.profile,
            "-w",
            str(self.work_dir),
            "-resume",
        ]

        # Add Docker image if using Docker profile
        if self.profile == "docker":
            args.extend(["-with-docker", self.docker_image])

        # Add resource limits
        args.extend(
            [
                "--max_cpus",
                str(self.max_cpus),
                "--max_memory",
                f"'{self.max_memory}'",
                "--max_time",
                f"'{self.max_time}'",
            ]
        )

        # Add extra parameters from initialization
        for key, value in self.extra_params.items():
            if isinstance(value, bool):
                if value:
                    args.append(f"--{key}")
            else:
                args.extend([f"--{key}", str(value)])

        # Add additional runtime parameters
        if additional_params:
            for key, value in additional_params.items():
                if isinstance(value, bool):
                    if value:
                        args.append(f"--{key}")
                else:
                    args.extend([f"--{key}", str(value)])

        return args

    def create_config_file(self, config_path: Path) -> Path:
        """
        Create a custom Nextflow configuration file.

        Args:
            config_path: Path where to create the config file

        Returns:
            Path to the created configuration file
        """
        config_content = f"""
// Generated Nextflow configuration for siRNAforge
params {{
    max_cpus = {self.max_cpus}
    max_memory = '{self.max_memory}'
    max_time = '{self.max_time}'
}}

process {{
    container = '{self.docker_image}'
}}

{self.profile} {{
    enabled = true
}}
"""

        config_path.write_text(config_content)
        logger.info(f"Created Nextflow config at {config_path}")
        return config_path

    def validate_docker_available(self) -> bool:
        """
        Check if Docker is available for execution.

        Returns:
            True if Docker is available and accessible
        """
        try:
            subprocess.run(["docker", "version"], capture_output=True, timeout=10, check=True)
            logger.debug("Docker is available")
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            logger.warning("Docker is not available")
            return False

    def get_execution_profile(self) -> str:
        """
        Get the appropriate execution profile based on available tools.

        Returns:
            Recommended execution profile
        """
        if self.profile == "docker" and self.validate_docker_available():
            return "docker"
        if self.profile == "docker" and not self.validate_docker_available():
            logger.warning("Docker requested but not available, falling back to local")
            return "local"
        return self.profile

    @classmethod
    def for_testing(cls) -> "NextflowConfig":
        """
        Create a configuration optimized for testing.

        Returns:
            NextflowConfig instance with test-friendly settings
        """
        return cls(
            profile="test",
            max_cpus=2,
            max_memory="6.GB",
            max_time="6.h",
            max_hits=100,
        )

    @classmethod
    def for_production(cls, docker_image: Optional[str] = None) -> "NextflowConfig":
        """
        Create a configuration optimized for production use.

        Args:
            docker_image: Override default Docker image

        Returns:
            NextflowConfig instance with production settings
        """
        return cls(
            docker_image=docker_image or "ghcr.io/austin-s-h/sirnaforge:latest",
            profile="docker",
            max_cpus=16,
            max_memory="128.GB",
            max_time="240.h",
        )
