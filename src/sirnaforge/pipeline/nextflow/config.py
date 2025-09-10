"""
Nextflow Configuration Management

This module handles configuration for Nextflow workflows, including
Docker settings, resource management, and parameter validation.

Simple Usage Examples:

    # Auto-configure based on environment (easiest)
    config = NextflowConfig.auto_configure()

    # Production settings
    config = NextflowConfig.for_production()

    # Testing settings
    config = NextflowConfig.for_testing()

    # Local Docker testing (uses local image built by 'make docker')
    config = NextflowConfig.for_local_docker_testing()
"""

import os
import shutil
import subprocess  # nosec B404
from pathlib import Path
from typing import Any, Optional

from pydantic import BaseModel, Field

from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class EnvironmentInfo(BaseModel):
    """Information about the current execution environment."""

    running_in_docker: bool = Field(description="Whether the current process is running inside a Docker container")
    docker_available: bool = Field(description="Whether Docker is available and functional for Nextflow execution")
    requested_profile: str = Field(description="The originally requested execution profile")
    recommended_profile: str = Field(description="The profile recommended based on environment detection")
    docker_image: Optional[str] = Field(default=None, description="Docker image to use (if applicable)")
    profile_override_reason: Optional[str] = Field(default=None, description="Reason for profile override (if any)")

    def is_profile_overridden(self) -> bool:
        """Check if the recommended profile differs from the requested profile."""
        return self.requested_profile != self.recommended_profile

    def get_execution_summary(self) -> str:
        """Get a human-readable summary of the execution environment."""
        summary = f"Profile: {self.recommended_profile}"

        if self.is_profile_overridden():
            summary += f" (overridden from {self.requested_profile})"
            if self.profile_override_reason:
                summary += f" - {self.profile_override_reason}"

        if self.running_in_docker:
            summary += " | Running in container"

        if self.docker_available and self.recommended_profile == "docker":
            summary += f" | Using Docker image: {self.docker_image}"

        return summary


def _get_executable_path(tool_name: str) -> Optional[str]:
    """Get the full path to an executable, ensuring it exists."""
    path = shutil.which(tool_name)
    if path is None:
        logger.warning(f"Tool '{tool_name}' not found in PATH")
    return path


def _validate_command_args(cmd: list[str]) -> None:
    """Validate command arguments for subprocess execution."""
    if not cmd:
        raise ValueError("Command list cannot be empty")

    executable = cmd[0]
    if not executable:
        raise ValueError("Executable path cannot be empty")

    # Ensure we have an absolute path to the executable
    if not Path(executable).is_absolute():
        raise ValueError(f"Executable must be an absolute path: {executable}")


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
        include_test_profile: bool = False,
    ) -> list[str]:
        """
        Generate Nextflow command arguments.

        Args:
            input_file: Input FASTA file path
            output_dir: Output directory
            genome_species: List of genome species to analyze
            additional_params: Additional parameters to pass
            include_test_profile: Whether to include 'test' profile for integration testing

        Returns:
            List of command arguments for Nextflow
        """
        # Ensure all paths are absolute to avoid working directory issues
        abs_input_file = input_file.resolve()
        abs_output_dir = output_dir.resolve()
        abs_work_dir = self.work_dir.resolve()
        args = [
            "--input",
            str(abs_input_file),
            "--outdir",
            str(abs_output_dir),
            "--genome_species",
            ",".join(genome_species),
            "-profile",
            self.profile,
            "-w",
            str(abs_work_dir),
            "-resume",
        ]

        # Add test profile for local integration testing
        if include_test_profile and os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes"):
            args.extend(["-profile", "test"])

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
        Check if Docker is available for Nextflow execution.

        This checks if Docker can be used by Nextflow to run containers.
        Note: This is different from running tests inside Docker containers.

        Returns:
            True if Docker is available and accessible for Nextflow
        """
        try:
            # Get absolute path to docker executable
            docker_path = _get_executable_path("docker")
            if not docker_path:
                logger.warning("Docker executable not found in PATH")
                return False

            cmd = [docker_path, "version"]
            _validate_command_args(cmd)
            subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603

            # Additional check: try to run a simple container to ensure Docker daemon is accessible
            # This helps distinguish between Docker being installed vs Docker daemon being available
            test_cmd = [docker_path, "run", "--rm", "hello-world"]
            subprocess.run(test_cmd, capture_output=True, timeout=30, check=True)  # nosec B603

            logger.debug("Docker is available and functional for Nextflow")
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            logger.warning("Docker is not available or not functional for Nextflow")
            return False

    def is_running_in_docker(self) -> bool:
        """
        Check if we're currently running inside a Docker container.

        This is useful for determining the appropriate execution profile
        when running tests or workflows.

        Returns:
            True if running inside a Docker container
        """
        try:
            # Check for Docker-specific files
            if Path("/.dockerenv").exists():
                return True

            # Check cgroup for Docker container indicators
            cgroup_content = Path("/proc/1/cgroup").read_text()
            return "docker" in cgroup_content or "containerd" in cgroup_content
        except (FileNotFoundError, OSError):
            return False

    def _is_singularity_available(self) -> bool:
        """Check if Singularity is available."""
        try:
            singularity_path = _get_executable_path("singularity")
            if not singularity_path:
                return False
            cmd = [singularity_path, "--version"]
            _validate_command_args(cmd)
            subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            return False

    def _is_conda_available(self) -> bool:
        """Check if Conda or uv is available for environment management."""
        # First check for uv (preferred for this project)
        try:
            uv_path = _get_executable_path("uv")
            if uv_path:
                cmd = [uv_path, "--version"]
                _validate_command_args(cmd)
                subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603
                return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            pass

        # Fallback to conda
        try:
            conda_path = _get_executable_path("conda")
            if not conda_path:
                return False
            cmd = [conda_path, "--version"]
            _validate_command_args(cmd)
            subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            return False

    def get_execution_profile(self) -> str:
        """
        Get the appropriate execution profile based on available tools and environment.

        This method considers:
        1. Environment variables (SIRNAFORGE_USE_LOCAL_EXECUTION)
        2. Whether we're running inside a Docker container (for testing)
        3. Whether Docker is available for Nextflow execution
        4. Availability of Singularity or Conda as fallbacks
        5. The requested profile

        Returns:
            Recommended execution profile
        """
        recommended_profile: Optional[str] = None

        # Check for explicit environment variable to force local execution
        if os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes"):
            logger.info("SIRNAFORGE_USE_LOCAL_EXECUTION set, using local execution profile")
            recommended_profile = "local"

        # If we're running inside a Docker container (e.g., for testing),
        # we might not have access to Docker daemon, so use local execution
        if recommended_profile is None and self.is_running_in_docker():
            logger.info("Running inside Docker container, using local execution profile")
            recommended_profile = "local"

        # If Docker profile is requested, check if Docker is available
        if recommended_profile is None and self.profile == "docker":
            if self.validate_docker_available():
                logger.debug("Using Docker profile for Nextflow execution")
                recommended_profile = "docker"
            else:
                logger.warning("Docker requested but not available, falling back to alternatives")

        # Check for alternative container runtimes
        if recommended_profile is None:
            if self._is_singularity_available():
                logger.info("Using Singularity profile as Docker alternative")
                recommended_profile = "singularity"
            elif self._is_conda_available():
                logger.info("Using Conda profile as container alternative")
                recommended_profile = "conda"

        # Fallback to local if no container runtime available and a container profile was requested
        if recommended_profile is None:
            if self.profile in ["docker", "singularity", "conda"]:
                logger.warning(f"Requested profile '{self.profile}' not available, falling back to local")
                recommended_profile = "local"
            else:
                supported_profiles = ["local", "docker", "singularity", "conda", "test"]
                if self.profile not in supported_profiles:
                    raise ValueError(
                        f"Requested profile '{self.profile}' is not supported. Supported: {supported_profiles}"
                    )
                logger.debug(f"Using requested profile: {self.profile}")
                recommended_profile = self.profile

        return recommended_profile

    def get_environment_info(self) -> EnvironmentInfo:
        """
        Get information about the current execution environment.

        This provides structured information about Docker availability,
        profile selection, and environment detection.

        Returns:
            EnvironmentInfo model with environment details
        """
        running_in_docker = self.is_running_in_docker()
        docker_available = self.validate_docker_available()
        recommended_profile = self.get_execution_profile()

        # Determine override reason
        override_reason = None
        if self.profile != recommended_profile:
            if running_in_docker:
                override_reason = "Running inside container"
            elif os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION"):
                override_reason = "Environment variable SIRNAFORGE_USE_LOCAL_EXECUTION set"
            elif self.profile == "docker" and not docker_available:
                override_reason = "Docker not available"

        return EnvironmentInfo(
            running_in_docker=running_in_docker,
            docker_available=docker_available,
            requested_profile=self.profile,
            recommended_profile=recommended_profile,
            docker_image=self.docker_image if recommended_profile == "docker" else None,
            profile_override_reason=override_reason,
        )

    @classmethod
    def for_testing(cls) -> "NextflowConfig":
        """
        Create a configuration optimized for testing.

        This automatically detects if we're running in Docker and adjusts accordingly.
        Uses uv/conda for environment management when available.

        Returns:
            NextflowConfig instance with test-friendly settings
        """
        # For testing, use test profile by default, but allow environment override
        instance = cls(
            profile="test",  # Use test profile which runs locally with uv/conda
            max_cpus=2,
            max_memory="6.GB",
            max_time="6.h",
            max_hits=100,
        )

        # Auto-detect and adjust profile only if environment variable is set or in container
        detected_profile = instance.get_execution_profile()
        if detected_profile != instance.profile and (
            os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes")
            or instance.is_running_in_docker()
        ):
            instance.profile = detected_profile
        return instance

    @classmethod
    def for_production(cls, **kwargs: Any) -> "NextflowConfig":
        """
        Create a configuration optimized for production use.

        This uses Docker by default for reproducible execution with full resources.

        Args:
            **kwargs: Additional configuration parameters to override defaults

        Returns:
            NextflowConfig instance with production settings
        """
        return cls(
            docker_image="ghcr.io/austin-s-h/sirnaforge:latest",
            profile="docker",
            max_cpus=16,
            max_memory="128.GB",
            max_time="240.h",
            **kwargs,
        )

    @classmethod
    def auto_configure(cls, **kwargs: Any) -> "NextflowConfig":
        """
        Auto-configure Nextflow settings based on environment detection.

        This method automatically detects available tools and selects the best profile.

        Args:
            **kwargs: Additional configuration parameters to override defaults

        Returns:
            NextflowConfig instance with auto-detected settings
        """
        # Start with default production settings
        instance = cls.for_production(**kwargs)

        # Auto-detect the best profile
        detected_profile = instance.get_execution_profile()
        instance.profile = detected_profile

        logger.info(f"Auto-configured Nextflow with profile: {detected_profile}")
        return instance
