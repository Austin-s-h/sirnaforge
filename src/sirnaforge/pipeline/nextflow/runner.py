"""
Nextflow Pipeline Runner

This module provides a Python interface to execute Nextflow workflows
for siRNA off-target analysis with proper Docker integration.

Simple Usage Examples:

    # Auto-configured runner (easiest)
    runner = NextflowRunner.create()
    result = runner.run_sync(
        input_file=Path("candidates.fasta"),
        output_dir=Path("results")
    )

    # Local Docker testing (uses local image built by 'make docker')
    runner = NextflowRunner.for_local_docker_testing()
    result = runner.run_sync(...)

    # Async version
    result = await runner.run(
        input_file=Path("candidates.fasta"),
        output_dir=Path("results")
    )
"""

import asyncio
import shutil
import subprocess  # nosec B404
from pathlib import Path
from typing import Any, Optional

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from sirnaforge.utils.logging_utils import get_logger

from .config import NextflowConfig

logger = get_logger(__name__)


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


console = Console()


class NextflowRunner:
    """Execute Nextflow workflows from Python with proper error handling."""

    def __init__(self, config: Optional[NextflowConfig] = None):
        """
        Initialize Nextflow runner.

        Args:
            config: NextflowConfig instance, creates auto-configured if None
        """
        self.config = config or NextflowConfig.auto_configure()
        self.workflow_dir = self._get_workflow_dir()

    def _get_workflow_dir(self) -> Path:
        """Get the directory containing embedded Nextflow workflows."""
        # Get the package directory and find the workflows
        package_dir = Path(__file__).parent
        workflow_dir = package_dir / "workflows"

        if not workflow_dir.exists():
            raise FileNotFoundError(f"Nextflow workflows not found at {workflow_dir}")

        return workflow_dir

    def get_main_workflow(self) -> Path:
        """Get path to the main Nextflow workflow."""
        main_nf = self.workflow_dir / "main.nf"
        if not main_nf.exists():
            raise FileNotFoundError(f"Main workflow not found at {main_nf}")
        return main_nf

    async def run(
        self, input_file: Path, output_dir: Path, genome_species: Optional[list[str]] = None, **kwargs: Any
    ) -> dict[str, Any]:
        """
        Simple method to run Nextflow workflow with auto-validation and defaults.

        Args:
            input_file: Path to input FASTA file
            output_dir: Output directory for results
            genome_species: List of genome species (defaults to ["human", "rat", "rhesus"])
            **kwargs: Additional parameters passed to run_offtarget_analysis

        Returns:
            Dictionary containing execution results and metadata

        Raises:
            NextflowExecutionError: If workflow execution fails
        """
        # Auto-validate installation
        validation = self.validate_installation()
        if not validation.get("nextflow", False):
            raise NextflowExecutionError("Nextflow is not available. Please install Nextflow.")
        if not validation.get("workflow_files", False):
            raise NextflowExecutionError("Nextflow workflow files not found.")

        # Set defaults
        genome_species = genome_species or ["human", "rat", "rhesus"]

        # Run the analysis
        return await self.run_offtarget_analysis(
            input_file=input_file, output_dir=output_dir, genome_species=genome_species, **kwargs
        )

    def run_sync(
        self, input_file: Path, output_dir: Path, genome_species: Optional[list[str]] = None, **kwargs: Any
    ) -> dict[str, Any]:
        """
        Synchronous version of run() for simpler usage without async/await.

        Args:
            input_file: Path to input FASTA file
            output_dir: Output directory for results
            genome_species: List of genome species (defaults to ["human", "rat", "rhesus"])
            **kwargs: Additional parameters passed to run_offtarget_analysis

        Returns:
            Dictionary containing execution results and metadata
        """
        # Auto-validate installation
        validation = self.validate_installation()
        if not validation.get("nextflow", False):
            raise NextflowExecutionError("Nextflow is not available. Please install Nextflow.")
        if not validation.get("workflow_files", False):
            raise NextflowExecutionError("Nextflow workflow files not found.")

        # Set defaults
        genome_species = genome_species or ["human", "rat", "rhesus"]

        # Run synchronously
        return asyncio.run(
            self.run_offtarget_analysis(
                input_file=input_file, output_dir=output_dir, genome_species=genome_species, **kwargs
            )
        )

    async def run_offtarget_analysis(
        self,
        input_file: Path,
        output_dir: Path,
        genome_species: list[str],
        additional_params: Optional[dict[str, Any]] = None,
        show_progress: bool = True,
    ) -> dict[str, Any]:
        """
        Run the off-target analysis Nextflow workflow.

        Args:
            input_file: Path to siRNA candidates FASTA file
            output_dir: Output directory for results
            genome_species: List of genome species to analyze against
            additional_params: Additional parameters for the workflow
            show_progress: Whether to show progress indicators

        Returns:
            Dictionary containing execution results and metadata

        Raises:
            NextflowExecutionError: If workflow execution fails
        """
        # Validate inputs
        abs_input_file = input_file.resolve()
        if not abs_input_file.exists():
            raise FileNotFoundError(f"Input file not found: {abs_input_file}")

        output_dir.mkdir(parents=True, exist_ok=True)
        abs_output_dir = output_dir.resolve()

        # Get workflow path
        workflow_path = self.get_main_workflow()

        # Prepare command arguments with absolute paths
        args = self.config.get_nextflow_args(
            input_file=abs_input_file,
            output_dir=abs_output_dir,
            genome_species=genome_species,
            additional_params=additional_params,
        )

        # Build full command
        cmd = ["nextflow", "run", str(workflow_path)] + args

        logger.info(f"Executing Nextflow workflow from {Path.cwd()}: {' '.join(cmd)}")
        logger.debug(f"Input file exists: {abs_input_file.exists()} at {abs_input_file}")
        logger.debug(f"Output directory: {abs_output_dir}")
        logger.debug(f"Workflow path: {workflow_path}")

        if show_progress:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                task = progress.add_task("Running Nextflow off-target analysis...", total=None)

                try:
                    result = await self._run_subprocess(cmd)
                    progress.update(task, description="✅ Nextflow execution completed")

                except subprocess.CalledProcessError as e:
                    progress.update(task, description="❌ Nextflow execution failed")
                    raise NextflowExecutionError(f"Nextflow failed: {e}") from e
        else:
            try:
                result = await self._run_subprocess(cmd)
            except subprocess.CalledProcessError as e:
                raise NextflowExecutionError(f"Nextflow failed: {e}") from e

        # Process results
        return self._process_results(output_dir, result)

    async def _run_subprocess(self, cmd: list[str]) -> Any:
        """
        Run subprocess asynchronously with proper logging.

        Args:
            cmd: Command to execute as list of strings

        Returns:
            AsyncResult instance with subprocess-like interface

        Raises:
            subprocess.CalledProcessError: If command fails
        """
        logger.debug(f"Running command: {' '.join(cmd)}")

        # Run in thread pool to avoid blocking
        # Execute from project root directory to ensure relative paths work correctly
        project_root = Path.cwd()
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=project_root,
        )

        stdout, stderr = await process.communicate()

        # Log outputs
        if stdout:
            logger.debug(f"Nextflow stdout: {stdout.decode()}")
        if stderr:
            logger.debug(f"Nextflow stderr: {stderr.decode()}")

        # Determine if this is a real error or just a version warning
        final_return_code = process.returncode if process.returncode is not None else 0

        if process.returncode != 0:
            error_msg = stderr.decode() if stderr else "Unknown error"
            # Check if this is just a version warning (not a real error)
            if self._is_nextflow_version_warning_only(error_msg):
                logger.warning(f"Nextflow version warning: {error_msg.strip()}")
                # Treat version warnings as success
                final_return_code = 0
            else:
                logger.error(f"Nextflow failed with code {process.returncode}: {error_msg}")
                # This is a real error - raise exception
                return_code = process.returncode if process.returncode is not None else 1
                raise subprocess.CalledProcessError(return_code, cmd, output=stdout, stderr=stderr)

        # Create CompletedProcess-like object
        class AsyncResult:
            def __init__(self, returncode: int, stdout: bytes, stderr: bytes) -> None:
                self.returncode = returncode
                self.stdout = stdout
                self.stderr = stderr

        return AsyncResult(final_return_code, stdout, stderr)

    def _is_nextflow_version_warning_only(self, error_msg: str) -> bool:
        """
        Check if the error message is only a Nextflow version warning.

        Args:
            error_msg: The error message from stderr

        Returns:
            True if this is only a version warning, False if it's a real error
        """
        lines = error_msg.strip().split("\n")
        # Filter out empty lines
        non_empty_lines = [line.strip() for line in lines if line.strip()]

        if not non_empty_lines:
            return False

        # Check if all non-empty lines are version warnings
        version_warning_patterns = [
            "is available - Please consider updating your version",
            "Nextflow",  # Simple version info lines
        ]

        for line in non_empty_lines:
            is_version_related = any(pattern in line for pattern in version_warning_patterns)
            if not is_version_related:
                return False

        return True

    def _process_results(self, output_dir: Path, result: Any) -> dict[str, Any]:
        """
        Process Nextflow execution results and extract metadata.

        Args:
            output_dir: Output directory where results were written
            result: Subprocess result object

        Returns:
            Dictionary containing execution metadata and result paths
        """
        # Find output files
        output_files = {
            "combined_analyses": list(output_dir.glob("**/combined_*_analysis.tsv")),
            "combined_summary": list(output_dir.glob("**/combined_summary.json")),
            "html_report": list(output_dir.glob("**/analysis_report.html")),
            "validation_report": list(output_dir.glob("**/validation_report.txt")),
            "individual_results": list(output_dir.glob("**/individual_results/")),
        }

        # Count results
        result_counts = {
            "total_files": sum(len(files) for files in output_files.values()),
            "analysis_files": len(output_files["combined_analyses"]),
            "summary_files": len(output_files["combined_summary"]),
            "report_files": len(output_files["html_report"]),
        }

        return {
            "status": "completed",
            "output_dir": str(output_dir),
            "output_files": {k: [str(f) for f in v] for k, v in output_files.items()},
            "result_counts": result_counts,
            "returncode": result.returncode,
            "stdout": result.stdout.decode() if hasattr(result, "stdout") else "",
            "stderr": result.stderr.decode() if hasattr(result, "stderr") else "",
        }

    def validate_installation(self) -> dict[str, bool]:
        """
        Validate that Nextflow and required tools are available.

        Returns:
            Dictionary of tool availability status
        """
        tools = {}

        # Check Nextflow
        try:
            # Get absolute path to nextflow executable
            nextflow_path = _get_executable_path("nextflow")
            if not nextflow_path:
                tools["nextflow"] = False
                logger.warning("Nextflow not available - not found in PATH")
            else:
                cmd = [nextflow_path, "-version"]
                _validate_command_args(cmd)
                result = subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603
                tools["nextflow"] = True
                logger.debug(f"Nextflow version: {result.stdout.decode()}")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            tools["nextflow"] = False
            logger.warning("Nextflow not available")

        # Check Docker if using Docker profile
        tools["docker"] = self.config.validate_docker_available()

        # Check uv/conda for environment management
        try:
            uv_path = _get_executable_path("uv")
            if uv_path:
                cmd = [uv_path, "--version"]
                _validate_command_args(cmd)
                result = subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603
                tools["uv"] = True
                logger.debug(f"uv version: {result.stdout.decode()}")
            else:
                tools["uv"] = False
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            tools["uv"] = False

        # Check conda as fallback
        try:
            conda_path = _get_executable_path("conda")
            if conda_path:
                cmd = [conda_path, "--version"]
                _validate_command_args(cmd)
                result = subprocess.run(cmd, capture_output=True, timeout=10, check=True)  # nosec B603
                tools["conda"] = True
                logger.debug(f"Conda version: {result.stdout.decode()}")
            else:
                tools["conda"] = False
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            tools["conda"] = False

        # Check workflow files
        tools["workflow_files"] = self.workflow_dir.exists() and (self.workflow_dir / "main.nf").exists()

        return tools

    @classmethod
    def for_testing(cls) -> "NextflowRunner":
        """
        Create a runner configured for testing.

        This uses the test configuration which automatically detects
        if we're running in Docker and adjusts accordingly.

        Returns:
            NextflowRunner configured for testing
        """
        config = NextflowConfig.for_testing()
        return cls(config)

    @classmethod
    def create(cls, **config_kwargs: Any) -> "NextflowRunner":
        """
        Create a NextflowRunner with auto-configured settings.

        Args:
            **config_kwargs: Configuration parameters for NextflowConfig

        Returns:
            NextflowRunner with auto-configured NextflowConfig
        """
        config = NextflowConfig.auto_configure(**config_kwargs)
        return cls(config)


class NextflowExecutionError(Exception):
    """Exception raised when Nextflow execution fails."""

    def __init__(self, message: str, stdout: str = "", stderr: str = ""):
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr
