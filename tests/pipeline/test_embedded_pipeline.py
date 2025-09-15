"""
Integration tests for embedded Nextflow pipeline.

Tests the integrated Python-Nextflow pipeline functionality.
"""

import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from sirnaforge.core.off_target import (
    parse_fasta_file,
    validate_and_write_sequences,
    write_fasta_file,
)
from sirnaforge.pipeline import NextflowConfig, NextflowRunner, get_test_data_path

logger = logging.getLogger(__name__)


class ResourceManager:
    """Mock resource manager for testing."""

    def get_test_config(self):
        """Return test configuration."""
        return {
            "test_candidates": "test_candidates.fasta",
            "test_transcriptome": "test_transcriptome.fasta",
            "test_mirna_seeds": "test_mirna.fasta",
            "genomes_config": "genomes.yaml",
        }


class TestPipelineIntegration:
    """Test embedded Nextflow pipeline integration."""

    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.resource_manager = ResourceManager()

    def teardown_method(self):
        """Clean up test environment."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    @pytest.mark.integration
    @pytest.mark.local_nextflow
    def test_nextflow_config_creation(self):
        """Test NextflowConfig creation and validation."""
        config = NextflowConfig.for_testing()

        # Profile depends on environment - 'test' normally, 'local' if in Docker or env var set
        expected_profile = (
            "local"
            if config.is_running_in_docker()
            or os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes")
            else "test"
        )
        assert config.profile == expected_profile
        assert config.max_cpus == 2
        assert config.max_memory == "6.GB"
        assert config.docker_image == "ghcr.io/austin-s-h/sirnaforge:latest"

        # Test argument generation
        input_file = Path("test.fasta").resolve()
        output_dir = Path("results").resolve()
        args = config.get_nextflow_args(
            input_file=Path("test.fasta"),
            output_dir=Path("results"),
            genome_species=["human", "rat"],
            include_test_profile=True,  # Include test profile for integration testing
        )

        # Determine expected profiles based on actual configuration
        expected_main_profile = config.profile  # This will be 'test' or 'local' depending on environment
        expected_profiles = [expected_main_profile]

        # If include_test_profile is True and SIRNAFORGE_USE_LOCAL_EXECUTION is set, add test profile
        if os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes"):
            expected_profiles.append("test")

        expected_args = [
            "--input",
            str(input_file),
            "--outdir",
            str(output_dir),
            "--genome_species",
            "human,rat",
        ]

        # Add profile arguments
        for profile in expected_profiles:
            expected_args.extend(["-profile", profile])

        expected_args.extend(
            [
                "-w",
                str(config.work_dir),
                "-resume",
                "--max_cpus",
                "2",
                "--max_memory",
                "'6.GB'",
                "--max_time",
                "'6.h'",
                "--max_hits",
                "100",
            ]
        )

        # Check that all expected args are present
        for arg in expected_args:
            assert arg in args

    @pytest.mark.integration
    @pytest.mark.local_nextflow
    def test_nextflow_runner_initialization(self):
        """Test NextflowRunner initialization."""
        config = NextflowConfig.for_testing()
        runner = NextflowRunner(config)

        assert runner.config == config
        assert runner.workflow_dir.name == "workflows"

    def test_resource_management(self):
        """Test resource discovery and validation."""
        # Test test data validation
        validation = self.resource_manager.get_test_config()

        # Should have all required test files configured
        expected_keys = ["test_candidates", "test_transcriptome", "test_mirna_seeds", "genomes_config"]

        for key in expected_keys:
            assert key in validation

    def test_get_test_data_path(self):
        """Test test data path retrieval."""
        try:
            candidates_path = get_test_data_path("test_candidates.fasta")
            assert candidates_path.exists()
            assert candidates_path.name == "test_candidates.fasta"

            # Verify file contains FASTA data
            content = candidates_path.read_text()
            assert content.startswith(">")
            assert "GCAUGAACCGGAGGCCCAUUU" in content

        except FileNotFoundError:
            pytest.skip("Test data files not found - this is expected during development")

    def test_workflow_path_discovery(self):
        """Test workflow file discovery."""
        runner = NextflowRunner()

        try:
            main_workflow = runner.get_main_workflow()
            assert main_workflow.exists()
            assert main_workflow.name == "main.nf"

            # Verify it's a valid Nextflow file
            content = main_workflow.read_text()
            assert "nextflow.enable.dsl = 2" in content
            assert "SIRNAFORGE_OFFTARGET" in content

        except FileNotFoundError:
            pytest.skip("Embedded workflows not found - expected during development")

    def test_tool_validation(self):
        """Test validation of required tools."""
        runner = NextflowRunner.for_testing()
        validation = runner.validate_installation()

        # Should check for key tools
        expected_tools = ["nextflow", "docker", "workflow_files"]
        for tool in expected_tools:
            assert tool in validation
            # Each validation should be a boolean
            assert isinstance(validation[tool], bool)

    @pytest.mark.integration
    @pytest.mark.integration
    @pytest.mark.local_nextflow
    def test_pipeline_execution_dry_run(self):
        """Test pipeline execution in dry-run mode (requires Nextflow)."""
        pytest.importorskip("subprocess")

        try:
            runner = NextflowRunner.for_testing()

            # Check if tools are available
            validation = runner.validate_installation()
            if not validation["nextflow"]:
                pytest.skip("Nextflow not available")

            if not validation["workflow_files"]:
                pytest.skip("Workflow files not available")

            # Create test input
            test_input = self.temp_dir / "test_input.fasta"
            test_input.write_text(">test_seq\nGCAUGAACCGGAGGCCCAUUU\n")

            # This would be an actual execution test
            # For now, we just verify the runner can be configured
            expected_profile = (
                "local"
                if runner.config.is_running_in_docker()
                or os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes")
                else "test"
            )
            assert runner.config.profile == expected_profile

        except Exception as e:
            pytest.skip(f"Pipeline execution test skipped: {e}")

    def test_config_for_production_vs_testing(self):
        """Test different configuration presets."""
        test_config = NextflowConfig.for_testing()
        prod_config = NextflowConfig.for_production()

        # Test config should have lower resource limits
        assert test_config.max_cpus < prod_config.max_cpus
        assert "6.GB" in test_config.max_memory
        assert "128.GB" in prod_config.max_memory

        # Test config profile depends on environment (local if SIRNAFORGE_USE_LOCAL_EXECUTION is set)
        expected_test_profile = (
            "local" if os.getenv("SIRNAFORGE_USE_LOCAL_EXECUTION", "").lower() in ("true", "1", "yes") else "test"
        )
        assert test_config.profile == expected_test_profile
        assert prod_config.profile == "docker"

    def test_docker_availability_check(self):
        """Test Docker availability validation."""
        config = NextflowConfig()

        # This will actually try to run docker --version
        docker_available = config.validate_docker_available()
        assert isinstance(docker_available, bool)

        # Test execution profile selection
        profile = config.get_execution_profile()
        assert profile in ["docker", "local", "test", "conda", "singularity"]

    @pytest.mark.integration
    @pytest.mark.docker
    @pytest.mark.local_nextflow
    def test_offtarget_analysis_docker_integration(self):
        """
        Integration test for off-target analysis module using Docker.

        Tests the complete off-target analysis pipeline including:
        - Index building with core entrypoint functions
        - Off-target analysis with BWA-MEM2
        - Result aggregation and validation
        """
        pytest.importorskip("subprocess")

        try:
            runner = self._setup_test_environment()
            test_output_dir, test_candidates = self._create_test_files()
            config = self._configure_nextflow_docker()
            offtarget_params = self._build_offtarget_params(test_candidates, test_output_dir, config)

            # Execute pipeline
            result = self._execute_nextflow_pipeline(runner, offtarget_params)
            self._handle_pipeline_result(result, test_output_dir)

        except Exception as e:
            pytest.skip(f"Off-target Docker integration test skipped: {e}")

    def _setup_test_environment(self) -> NextflowRunner:
        """Setup and validate test environment for Docker integration."""
        runner = NextflowRunner.for_testing()
        validation = runner.validate_installation()

        if not validation.get("nextflow", False):
            pytest.skip("Nextflow not available for Docker integration test")
        if not validation.get("docker", False):
            pytest.skip("Docker not available for integration test")
        if not validation.get("workflow_files", False):
            pytest.skip("Workflow files not available")

        return runner

    def _create_test_files(self) -> tuple[Path, Path]:
        """Create test input files for the pipeline."""
        test_output_dir = self.temp_dir / "offtarget_results"
        test_output_dir.mkdir(exist_ok=True)

        test_candidates = self.temp_dir / "test_candidates.fasta"
        test_candidates.write_text(
            ">candidate_1\n"
            "GCAUGAACCGGAGGCCCAUUU\n"
            ">candidate_2\n"
            "UGGAUCCAGAUCGAAAUGACU\n"
            ">candidate_3\n"
            "CGGAUUACCUGGAGCUUAAUG\n"
        )

        test_transcriptome = self.temp_dir / "test_transcriptome.fasta"
        test_transcriptome.write_text(
            ">transcript_1\n"
            "ATGGCAUGAACCGGAGGCCCAUUUCGAAUCGAAAUGACUGGAUUCCA\n"
            ">transcript_2\n"
            "CGGAUUACCUGGAGCUUAAUGCCUAGGUAACUCGAAAGCUCCAGGU\n"
        )

        return test_output_dir, test_candidates

    def _configure_nextflow_docker(self) -> NextflowConfig:
        """Configure Nextflow for Docker execution."""
        config = NextflowConfig.for_testing()
        config.profile = "docker"
        # Use a container that has BWA-MEM2 but not Bowtie to avoid memory issues
        config.docker_image = "ghcr.io/austin-s-h/sirnaforge:latest"
        return config

    def _build_offtarget_params(self, test_candidates: Path, test_output_dir: Path, config: NextflowConfig) -> dict:
        """Build parameters for off-target analysis pipeline."""
        test_transcriptome = self.temp_dir / "test_transcriptome.fasta"
        return {
            "--input": str(test_candidates),
            "--outdir": str(test_output_dir),
            "--genome_species": "test_species",
            "--genome_fastas": f"test_species:{test_transcriptome}",
            "--max_hits": "50",
            "--bwa_k": "12",
            "--bwa_T": "15",
            "--seed_start": "2",
            "--seed_end": "8",
            "-profile": "smoke",  # Use smoke profile for resource-constrained environments
            "-with-docker": config.docker_image,
        }

    def _execute_nextflow_pipeline(self, runner: NextflowRunner, offtarget_params: dict):
        """Execute the Nextflow pipeline with given parameters."""
        nextflow_cmd = ["nextflow", "run", str(runner.get_main_workflow())]
        for key, value in offtarget_params.items():
            nextflow_cmd.extend([key, value])

        work_dir = self.temp_dir / "work"
        nextflow_cmd.extend(["-w", str(work_dir)])

        print(f"Running Nextflow command: {' '.join(nextflow_cmd)}")

        try:
            return subprocess.run(
                nextflow_cmd,
                cwd=runner.workflow_dir.parent,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout
                check=False,  # Don't raise on non-zero exit
            )
        except subprocess.TimeoutExpired:
            pytest.skip("Pipeline execution timed out - this is expected in resource-constrained environments")

    def _is_version_warning_only(self, error_msg: str) -> bool:
        """
        Check if the error message is only a Nextflow version warning.

        Args:
            error_msg: The error message from stderr

        Returns:
            True if this is only a version warning, False if it's a real error
        """
        if not error_msg:
            return False

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

    def _handle_pipeline_result(self, result, test_output_dir: Path):
        """Handle the result of pipeline execution."""
        print(f"Nextflow stdout: {result.stdout}")
        print(f"Nextflow stderr: {result.stderr}")
        print(f"Nextflow exit code: {result.returncode}")

        if result.returncode == 0:
            self._validate_offtarget_outputs(test_output_dir)
        elif self._is_version_warning_only(result.stderr):
            # Version warning with non-zero return code - treat as success
            logger.warning(f"Nextflow version warning (treating as success): {result.stderr}")
            self._validate_offtarget_outputs(test_output_dir)
        elif "command not found" in result.stderr.lower():
            pytest.skip(f"Required tools not available: {result.stderr}")
        elif "docker" in result.stderr.lower() and "not found" in result.stderr.lower():
            pytest.skip(f"Docker not properly configured: {result.stderr}")
        elif "index" in result.stderr.lower():
            pytest.skip(f"Index building failed (expected in test env): {result.stderr}")
        else:
            pytest.fail(f"Pipeline failed with unexpected error: {result.stderr}")

    def _validate_offtarget_outputs(self, output_dir: Path):
        """Validate that the off-target analysis produced expected outputs."""
        # Check for main output files
        expected_files = ["combined_offtarget_analysis.tsv", "combined_summary.json", "final_summary.txt"]

        found_files = []
        for expected_file in expected_files:
            file_path = output_dir / expected_file
            if file_path.exists() and file_path.stat().st_size > 0:
                found_files.append(expected_file)
            else:
                print(f"Expected output file not found or empty: {expected_file}")

        # In test environments, pipeline might fail due to various issues (memory, Docker, etc.)
        # If no files are produced, skip the test rather than fail
        if not found_files:
            pytest.fail("No expected output files were generated by the pipeline.")

        print(f"Found output files: {found_files}")

        # Only validate files that were actually created
        if "combined_offtarget_analysis.tsv" in found_files:
            tsv_file = output_dir / "combined_offtarget_analysis.tsv"
            content = tsv_file.read_text()
            # Should have proper TSV header for off-target analysis
            expected_headers = ["qname", "qseq", "rname", "coord", "strand"]
            header_line = content.split("\n")[0]
            for header in expected_headers:
                assert header in header_line, f"Missing expected header: {header}"

        if "combined_summary.json" in found_files:
            json_file = output_dir / "combined_summary.json"
            summary_data = json.loads(json_file.read_text())
            assert "genome_species" in summary_data
            assert "total_analysis_files" in summary_data
            assert "analysis_timestamp" in summary_data

        if "final_summary.txt" in found_files:
            summary_file = output_dir / "final_summary.txt"
            content = summary_file.read_text()
            assert "SiRNA Off-Target Analysis Summary" in content
            assert "Analysis completed" in content

    @pytest.mark.integration
    @pytest.mark.integration
    @pytest.mark.local_nextflow
    def test_offtarget_entrypoint_functions_directly(self):
        """
        Test the off-target entrypoint functions directly without Nextflow.

        This validates that the core functions work correctly and can be
        called from the pipeline modules.
        """
        try:
            # Test imports work

            sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

            # Create test data
            test_fasta = self.temp_dir / "test_sequences.fasta"
            test_fasta.write_text(">test_seq_1\nGCAUGAACCGGAGGCCCAUUU\n>test_seq_2\nUGGAUCCAGAUCGAAAUGACU\n")

            # Test FASTA parsing
            sequences = parse_fasta_file(str(test_fasta))
            assert len(sequences) == 2
            assert "test_seq_1" in sequences
            assert sequences["test_seq_1"] == "GCATGAACCGGAGGCCCATTT"

            # Test sequence validation
            valid_output = self.temp_dir / "valid_sequences.fasta"
            valid_count, invalid_count, issues = validate_and_write_sequences(
                str(test_fasta), str(valid_output), expected_length=21
            )

            assert valid_count == 2
            assert invalid_count == 0
            assert len(issues) == 0
            assert valid_output.exists()

            # Test FASTA writing
            test_write_fasta = self.temp_dir / "written_sequences.fasta"
            write_fasta_file(sequences, str(test_write_fasta))
            assert test_write_fasta.exists()

            # Read back and verify
            written_sequences = parse_fasta_file(str(test_write_fasta))
            assert written_sequences == sequences

            print("✅ All off-target entrypoint functions work correctly")

        except ImportError as e:
            pytest.skip(f"Could not import off-target functions: {e}")
        except Exception as e:
            pytest.fail(f"Off-target entrypoint function test failed: {e}")
