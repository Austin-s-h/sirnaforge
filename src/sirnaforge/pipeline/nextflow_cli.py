"""
Command-line interface functions for Nextflow module integration.

These functions serve as entry points for Nextflow modules, providing
clean separation between workflow orchestration and business logic.
All complex logic is properly encapsulated in the core modules.
"""

import json
import shutil
import sys
from pathlib import Path
from typing import Any

from sirnaforge.core.off_target import (
    aggregate_mirna_results,
    aggregate_offtarget_results,
    build_bwa_index,
    parse_fasta_file,
    run_bwa_alignment_analysis,
    run_mirna_seed_analysis,
    validate_and_write_sequences,
)
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


def split_candidates_cli(input_fasta: str, output_dir: str = ".") -> dict[str, Any]:
    """
    Split multi-FASTA into individual candidate files with manifest.

    Args:
        input_fasta: Path to input FASTA file
        output_dir: Directory to write output files

    Returns:
        Dictionary with candidate count and manifest path
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Parse candidates
    sequences = parse_fasta_file(input_fasta)

    # Split into individual files
    candidate_files = []
    for i, (seq_id, sequence) in enumerate(sequences.items()):
        filename = f"candidate_{i + 1:04d}.fasta"
        filepath = output_path / filename
        with filepath.open("w") as f:
            f.write(f">{seq_id}\n{sequence}\n")
        candidate_files.append((seq_id, filename))

    # Write manifest
    manifest_path = output_path / "candidate_manifest.txt"
    with manifest_path.open("w") as f:
        f.write("sequence_id\tfilename\n")
        for seq_id, filename in candidate_files:
            f.write(f"{seq_id}\t{filename}\n")

    logger.info(f"Split {len(sequences)} candidates into individual files")

    return {
        "candidate_count": len(sequences),
        "manifest": str(manifest_path),
        "output_files": [str(output_path / fname) for _, fname in candidate_files],
    }


def prepare_candidates_cli(
    input_fasta: str, output_fasta: str = "validated_candidates.fasta", expected_length: int = 21
) -> dict[str, Any]:
    """
    Validate and prepare siRNA candidates for analysis.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to write validated sequences
        expected_length: Expected siRNA length

    Returns:
        Dictionary with validation statistics
    """
    # Validate and write sequences (returns: valid_count, invalid_count, issues)
    valid, invalid, errors = validate_and_write_sequences(
        input_file=input_fasta, output_file=output_fasta, expected_length=expected_length
    )
    total = valid + invalid

    # Write validation report
    report_path = Path("validation_report.txt")
    with report_path.open("w") as f:
        f.write(f"Total candidates: {total}\n")
        f.write(f"Valid candidates: {valid}\n")
        f.write(f"Invalid candidates: {invalid}\n")
        if errors:
            f.write("\nErrors:\n")
            for error in errors:
                f.write(f"  {error}\n")

    logger.info(f"Validated {valid} out of {total} candidates")

    return {
        "total": total,
        "valid": valid,
        "invalid": invalid,
        "output_file": output_fasta,
        "report": str(report_path),
    }


def build_bwa_index_cli(fasta_file: str, species: str, output_dir: str = ".") -> dict[str, Any]:
    """
    Build BWA-MEM2 index for genome/transcriptome.

    Args:
        fasta_file: Path to input FASTA file
        species: Species identifier
        output_dir: Directory to write index files

    Returns:
        Dictionary with index prefix path
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Build index with species-specific prefix
    index_prefix = output_path / f"{species}_index"
    result_prefix = build_bwa_index(fasta_file=fasta_file, index_prefix=str(index_prefix))

    logger.info(f"Built BWA index for {species}: {result_prefix}")

    return {
        "species": species,
        "index_prefix": str(result_prefix),
        "index_files": list(output_path.glob(f"{species}_index*")),
    }


def run_offtarget_analysis_cli(
    candidate_fasta: str,
    candidate_id: str,
    species: str,
    index_prefix: str,
    output_dir: str = ".",
    max_hits: int = 10000,
    bwa_k: int = 12,
    bwa_T: int = 15,
    seed_start: int = 2,
    seed_end: int = 8,
) -> dict[str, Any]:
    """
    Run off-target analysis for a single candidate against a genome.

    Args:
        candidate_fasta: Path to candidate FASTA file
        candidate_id: Candidate identifier
        species: Species identifier
        index_prefix: BWA index prefix path
        output_dir: Directory to write results
        max_hits: Maximum hits to report
        bwa_k: BWA seed length
        bwa_T: BWA minimum score threshold
        seed_start: Seed region start (1-based)
        seed_end: Seed region end (1-based)

    Returns:
        Dictionary with analysis results
    """
    logger.info(f"Running off-target analysis for {candidate_id} against {species}")

    # Run analysis (this creates properly named output files)
    output_path = run_bwa_alignment_analysis(
        candidates_file=candidate_fasta,
        index_prefix=index_prefix,
        species=species,
        output_dir=output_dir,
        max_hits=max_hits,
        bwa_k=bwa_k,
        bwa_T=bwa_T,
        seed_start=seed_start,
        seed_end=seed_end,
    )

    # Find the generated files
    analysis_file = list(Path(output_dir).glob(f"*{species}_analysis.tsv"))[0]
    summary_file = list(Path(output_dir).glob(f"*{species}_summary.json"))[0]

    return {
        "candidate_id": candidate_id,
        "species": species,
        "analysis_file": str(analysis_file),
        "summary_file": str(summary_file),
        "output_dir": str(output_path),
    }


def run_mirna_seed_analysis_cli(
    candidate_fasta: str,
    candidate_id: str,
    mirna_db: str = "mirgenedb",
    mirna_species: str = "human",
    output_dir: str = ".",
) -> dict[str, Any]:
    """
    Run miRNA seed match analysis for a single candidate.

    Args:
        candidate_fasta: Path to candidate FASTA file
        candidate_id: Candidate identifier
        mirna_db: miRNA database to use (mirgenedb, mirbase, etc.)
        mirna_species: Comma-separated list of species
        output_dir: Directory to write results

    Returns:
        Dictionary with analysis results
    """
    logger.info(f"Running miRNA seed match analysis for {candidate_id}")
    logger.info(f"Database: {mirna_db}, Species: {mirna_species}")

    # Parse species list
    species_list = [s.strip() for s in mirna_species.split(",") if s.strip()]

    if not species_list:
        raise ValueError("No species provided for miRNA analysis")

    # Run the analysis using the core function
    output_path = run_mirna_seed_analysis(
        candidates_file=candidate_fasta,
        candidate_id=candidate_id,
        mirna_db=mirna_db,
        mirna_species=species_list,
        output_dir=output_dir,
    )

    # Return metadata about the results
    analysis_file = output_path / f"{candidate_id}_mirna_analysis.tsv"
    summary_file = output_path / f"{candidate_id}_mirna_summary.json"

    # Load summary to get hit count
    hit_count = 0
    if summary_file.exists():
        with summary_file.open() as f:
            summary_data = json.load(f)
            hit_count = summary_data.get("total_hits", 0)

    logger.info(f"miRNA seed analysis completed for {candidate_id}: {hit_count} hits")

    return {
        "candidate_id": candidate_id,
        "total_hits": hit_count,
        "analysis_file": str(analysis_file),
        "summary_file": str(summary_file),
        "status": "completed",
    }


def aggregate_results_cli(genome_species: str, output_dir: str = ".") -> dict[str, Any]:
    """
    Aggregate off-target analysis results from multiple candidates and genomes.

    Args:
        genome_species: Comma-separated list of species
        output_dir: Directory to write aggregated results

    Returns:
        Dictionary with aggregation statistics
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Collect all analysis and summary files from current directory (Nextflow stages them)
    current_dir = Path()
    analysis_files = list(current_dir.glob("*_analysis.tsv"))
    summary_files = list(current_dir.glob("*_summary.json"))

    logger.info(f"Found {len(analysis_files)} analysis files and {len(summary_files)} summary files")

    if analysis_files or summary_files:
        # Create a temporary results directory structure for aggregation
        results_dir = Path("temp_results")
        results_dir.mkdir(exist_ok=True)

        # Organize files by species (extract from filename)
        species_list = [s.strip() for s in genome_species.split(",") if s.strip()]
        for species in species_list:
            species_dir = results_dir / species
            species_dir.mkdir(exist_ok=True)

            # Copy relevant files to species directory
            for f in analysis_files:
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

            for f in summary_files:
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

        # Run aggregation using core function
        result_path = aggregate_offtarget_results(
            results_dir=str(results_dir), output_dir=output_dir, genome_species=genome_species
        )

        logger.info(f"Aggregation completed: {result_path}")

        return {
            "status": "completed",
            "analysis_files_processed": len(analysis_files),
            "summary_files_processed": len(summary_files),
            "species": species_list,
            "output_dir": str(result_path),
        }

    # Create empty final summary
    final_summary = output_path / "final_summary.txt"
    with final_summary.open("w") as f:  # type: ignore[assignment]
        f.write("No analysis results found to aggregate\n")  # type: ignore[attr-defined]

    logger.warning("No files to aggregate")

    return {
        "status": "empty",
        "analysis_files_processed": 0,
        "summary_files_processed": 0,
        "species": [],
        "output_dir": str(output_path),
    }


def aggregate_mirna_results_cli(
    mirna_db: str,
    mirna_species: str,
    results_dir: str = ".",
    output_dir: str = ".",
) -> dict[str, Any]:
    """
    Aggregate miRNA seed analysis results from multiple candidates.

    Args:
        mirna_db: miRNA database name used for analysis
        mirna_species: Comma-separated list of species analyzed
        results_dir: Directory containing individual miRNA results
        output_dir: Directory to write aggregated results

    Returns:
        Dictionary with aggregation statistics
    """
    logger.info(f"Aggregating miRNA results from {results_dir}")

    # Run aggregation using core function
    result_path = aggregate_mirna_results(
        results_dir=results_dir, output_dir=output_dir, mirna_db=mirna_db, mirna_species=mirna_species
    )

    # Load summary to get statistics
    summary_file = result_path / "combined_mirna_summary.json"
    stats = {}
    if summary_file.exists():
        with summary_file.open() as f:
            stats = json.load(f)

    logger.info(f"miRNA aggregation completed: {result_path}")

    return {
        "status": "completed",
        "mirna_database": mirna_db,
        "species": mirna_species.split(","),
        "total_hits": stats.get("total_mirna_hits", 0),
        "candidates_analyzed": stats.get("total_candidates", 0),
        "output_dir": str(result_path),
    }


def resolve_genome_indices_cli(
    genome_species: str,
    genome_indices_override: str = "",
    output_file: str = "resolved_indices.json",
) -> dict[str, Any]:
    """
    Resolve genome indices from configuration or auto-discovery.

    Args:
        genome_species: Comma-separated list of species
        genome_indices_override: Override indices (format: species:path,species:path)
        output_file: Path to write resolved indices JSON

    Returns:
        Dictionary of species to index paths
    """
    # This is a placeholder - implement actual resolution logic
    # For now, just create an empty resolution
    indices = {}

    species_list = [s.strip() for s in genome_species.split(",") if s.strip()]

    # Parse override if provided
    if genome_indices_override:
        for entry in genome_indices_override.split(","):
            if ":" in entry:
                species, index_path = entry.split(":", 1)
                indices[species.strip()] = index_path.strip()

    # Write resolved indices to JSON
    output_path = Path(output_file)
    with output_path.open("w") as f:
        json.dump(indices, f, indent=2)

    logger.info(f"Resolved indices for species: {list(indices.keys())}")
    for species, index_path in indices.items():
        logger.info(f"  {species}: {index_path}")

    return {
        "species": species_list,
        "resolved_indices": indices,
        "output_file": str(output_path),
    }


# Main entry point for command-line execution
def main() -> None:
    """Main entry point for CLI commands."""
    if len(sys.argv) < 2:
        print("Usage: python -m sirnaforge.pipeline.nextflow_cli <command> [args...]")
        print("\nAvailable commands:")
        print("  split_candidates <input_fasta>")
        print("  prepare_candidates <input_fasta> <output_fasta> <expected_length>")
        print("  build_bwa_index <fasta_file> <species>")
        print("  run_offtarget_analysis <candidate_fasta> <candidate_id> <species> <index_prefix>")
        print("  aggregate_results <genome_species>")
        sys.exit(1)

    command = sys.argv[1]

    if command == "split_candidates":
        result = split_candidates_cli(sys.argv[2])
        print(json.dumps(result, indent=2))

    elif command == "prepare_candidates":
        result = prepare_candidates_cli(sys.argv[2], sys.argv[3], int(sys.argv[4]))
        print(json.dumps(result, indent=2))

    elif command == "build_bwa_index":
        result = build_bwa_index_cli(sys.argv[2], sys.argv[3])
        print(json.dumps(result, default=str, indent=2))

    elif command == "run_offtarget_analysis":
        result = run_offtarget_analysis_cli(
            candidate_fasta=sys.argv[2],
            candidate_id=sys.argv[3],
            species=sys.argv[4],
            index_prefix=sys.argv[5],
        )
        print(json.dumps(result, indent=2))

    elif command == "aggregate_results":
        result = aggregate_results_cli(genome_species=sys.argv[2], output_dir=".")
        print(json.dumps(result, indent=2))

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
