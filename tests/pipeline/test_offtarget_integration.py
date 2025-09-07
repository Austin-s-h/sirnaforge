"""Integration test for off-target analysis wrapper."""

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import pytest


@pytest.mark.integration
@pytest.mark.docker
def test_offtarget_wrapper_integration():
    """Test off-target wrapper with toy genome and queries in Docker environment."""

    # Check required tools
    try:
        subprocess.run(["bwa-mem2", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pytest.skip("bwa-mem2 not available - this test requires Docker environment")

    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        # Create toy genome
        genome_file = work_dir / "genome.fa"
        genome_file.write_text(
            ">chr1\n"
            "ATGCGTACGTTAGCTAGCTAGCTAGCTGACTGACTGATCGATGCTAGCTAGCTGATCGATGCTAGC\n"
            ">chr2\n"
            "GCTAGCTAGCTAGGCTAGCTAGCTAACGATCGATCGATCGATCGTACGATCGATCGTAGCTAGCTA\n"
        )

        # Create query sequences
        query_file = work_dir / "queries.fasta"
        query_file.write_text(">q1\nGCTAGCTAGCTAGGCTAGCTA\n>q2\nTACGTTAGCTAGCTAGCTAGC\n")

        # Build BWA-MEM2 index
        index_cmd = ["bwa-mem2", "index", str(genome_file)]
        result = subprocess.run(index_cmd, check=False, capture_output=True, text=True, cwd=work_dir)
        assert result.returncode == 0, f"BWA-MEM2 index failed: {result.stderr}"

        # Verify index files were created
        expected_index_files = [
            "genome.fa.0123",
            "genome.fa.amb",
            "genome.fa.ann",
            "genome.fa.bwt.2bit.64",
            "genome.fa.pac",
        ]
        for idx_file in expected_index_files:
            assert (work_dir / idx_file).exists(), f"Missing index file: {idx_file}"

        # Run off-target analysis
        results_dir = work_dir / "results"
        results_dir.mkdir()
        out_prefix = results_dir / "offtarget"

        # Find the offtarget_wrapper.py script
        wrapper_script = Path(__file__).parent.parent.parent / "nextflow_pipeline" / "modules" / "offtarget_wrapper.py"
        assert wrapper_script.exists(), f"offtarget_wrapper.py not found at {wrapper_script}"

        analysis_cmd = [
            "python3",
            str(wrapper_script),
            "--mode",
            "bwa-mem2-only",
            "--index-prefix",
            str(genome_file),
            "--queries",
            str(query_file),
            "--out-prefix",
            str(out_prefix),
            "--max-hits",
            "100",
            "--verbose",
        ]

        result = subprocess.run(analysis_cmd, check=False, capture_output=True, text=True, cwd=work_dir)
        assert result.returncode == 0, f"Off-target analysis failed: {result.stderr}\nStdout: {result.stdout}"

        # Verify output files exist
        tsv_file = Path(f"{out_prefix}.tsv")
        json_file = Path(f"{out_prefix}.json")

        assert tsv_file.exists(), "TSV output file not created"
        assert json_file.exists(), "JSON output file not created"

        # Validate TSV output
        tsv_content = tsv_file.read_text()
        lines = tsv_content.strip().split("\n")
        assert len(lines) >= 1, "TSV file should have at least header line"

        # Check header
        expected_headers = [
            "qname",
            "qseq",
            "rname",
            "coord",
            "strand",
            "cigar",
            "mapq",
            "AS",
            "NM",
            "mismatch_count",
            "mismatch_positions",
            "offtarget_score",
            "source",
        ]
        header = lines[0].split("\t")
        assert header == expected_headers, f"TSV header mismatch: {header}"

        # Validate JSON output
        json_content = json_file.read_text()
        results: list[dict[str, Any]] = json.loads(json_content)
        assert isinstance(results, list), "JSON should contain a list of results"

        if results:  # If we have alignments
            # Check structure of first result
            first_result = results[0]
            required_keys = ["qname", "rname", "coord", "strand", "offtarget_score", "source"]
            for key in required_keys:
                assert key in first_result, f"Missing key '{key}' in JSON result"

            # Verify source is correct
            assert first_result["source"] == "bwa-mem2", "Source should be bwa-mem2"

            # Verify off-target score is a number
            assert isinstance(first_result["offtarget_score"], (int, float)), "offtarget_score should be numeric"


@pytest.mark.integration
@pytest.mark.docker
@pytest.mark.slow
def test_offtarget_wrapper_combined_mode():
    """Test off-target wrapper in combined mode (requires both bowtie and bwa-mem2)."""

    # Check for both tools
    missing_tools = []
    for tool in ["bowtie", "bwa-mem2"]:
        try:
            subprocess.run([tool, "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing_tools.append(tool)

    if missing_tools:
        pytest.skip(f"Missing tools: {missing_tools} - this test requires full Docker environment")

    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)

        # Create slightly larger genome for combined test
        genome_file = work_dir / "genome.fa"
        genome_file.write_text(
            ">chr1\n"
            "ATGCGTACGTTAGCTAGCTAGCTAGCTGACTGACTGATCGATGCTAGCTAGCTGATCGATGCTAGC" * 3 + "\n"
            ">chr2\n"
            "GCTAGCTAGCTAGGCTAGCTAGCTAACGATCGATCGATCGATCGTACGATCGATCGTAGCTAGCTA" * 3 + "\n"
        )

        query_file = work_dir / "queries.fasta"
        query_file.write_text(">q1\nGCTAGCTAGCTAGGCTAGCTA\n>q2\nTACGTTAGCTAGCTAGCTAGC\n>q3\nATGCGTACGTTAGCTAGCTA\n")

        # Build both BWA-MEM2 and Bowtie indices
        bwa_result = subprocess.run(
            ["bwa-mem2", "index", str(genome_file)], check=False, capture_output=True, text=True, cwd=work_dir
        )
        assert bwa_result.returncode == 0, f"BWA-MEM2 index failed: {bwa_result.stderr}"

        bowtie_result = subprocess.run(
            ["bowtie-build", str(genome_file), str(genome_file)], check=False, capture_output=True, text=True, cwd=work_dir
        )
        assert bowtie_result.returncode == 0, f"Bowtie index failed: {bowtie_result.stderr}"

        # Run combined analysis
        results_dir = work_dir / "results"
        results_dir.mkdir()
        out_prefix = results_dir / "combined_offtarget"

        wrapper_script = Path(__file__).parent.parent.parent / "nextflow_pipeline" / "modules" / "offtarget_wrapper.py"

        analysis_cmd = [
            "python3",
            str(wrapper_script),
            "--mode",
            "combined",
            "--index-prefix",
            str(genome_file),
            "--queries",
            str(query_file),
            "--out-prefix",
            str(out_prefix),
            "--max-hits",
            "50",
        ]

        result = subprocess.run(analysis_cmd, check=False, capture_output=True, text=True, cwd=work_dir)
        assert result.returncode == 0, f"Combined analysis failed: {result.stderr}\nStdout: {result.stdout}"

        # Verify outputs
        tsv_file = Path(f"{out_prefix}.tsv")
        json_file = Path(f"{out_prefix}.json")

        assert tsv_file.exists(), "Combined TSV output not created"
        assert json_file.exists(), "Combined JSON output not created"

        # Check that we have results from both tools
        json_content = json.loads(json_file.read_text())
        if json_content:  # If alignments found
            sources = {result["source"] for result in json_content}
            # We should have results from at least one tool, possibly both
            assert sources.intersection({"bowtie", "bwa-mem2"}), f"No expected sources found: {sources}"
