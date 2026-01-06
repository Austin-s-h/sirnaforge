"""Command-line entry points used by embedded Nextflow modules."""

import json
import shutil
from pathlib import Path
from typing import Any

from sirnaforge.config import DEFAULT_MIRNA_CANONICAL_SPECIES
from sirnaforge.core.off_target import aggregate_mirna_results, aggregate_offtarget_results, build_bwa_index
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

DEFAULT_MIRNA_SPECIES_ARGUMENT = ",".join(DEFAULT_MIRNA_CANONICAL_SPECIES)


def build_bwa_index_cli(fasta_file: str, species: str, output_dir: str = ".") -> dict[str, Any]:
    """Build BWA-MEM2 index for genome/transcriptome.

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


def aggregate_results_cli(  # noqa: PLR0912
    genome_species: str,
    output_dir: str = ".",
    mirna_db: str | None = None,
    mirna_species: str | None = None,
) -> dict[str, Any]:
    """Aggregate off-target analysis results from multiple candidates and genomes.

    Args:
        genome_species: Comma-separated list of species
        output_dir: Directory to write aggregated results
        mirna_db: The database that provided the reference
        mirna_species: The species code for the matching miRNA

    Returns:
        Dictionary with aggregation statistics
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Collect all analysis and summary files from current directory (Nextflow stages them)
    current_dir = Path()
    analysis_files = list(current_dir.glob("*_analysis.tsv")) + list(current_dir.glob("mirna_analysis.tsv"))
    summary_files = list(current_dir.glob("*_summary.json")) + list(current_dir.glob("mirna_summary.json"))

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
                name_lower = f.name.lower()
                if name_lower.startswith("mirna"):
                    continue
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

            for f in summary_files:
                name_lower = f.name.lower()
                if name_lower.startswith("mirna"):
                    continue
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

        # Also persist miRNA batch results (if present) under dedicated directory
        mirna_results_dir = results_dir / "mirna"
        mirna_results_dir.mkdir(exist_ok=True)
        mirna_analysis_files: list[Path] = []
        for f in analysis_files:
            if "mirna" not in f.name.lower():
                continue
            dest_name = f.name
            if not dest_name.endswith("_mirna_analysis.tsv"):
                dest_name = "batch_mirna_analysis.tsv"
            dest_path = mirna_results_dir / dest_name
            shutil.copy(f, dest_path)
            mirna_analysis_files.append(dest_path)

        for f in summary_files:
            if "mirna" not in f.name.lower():
                continue
            dest_name = f.name
            if not dest_name.endswith("_mirna_summary.json"):
                dest_name = "batch_mirna_summary.json"
            dest_path = mirna_results_dir / dest_name
            shutil.copy(f, dest_path)

        # Run aggregation using core function
        result_path = aggregate_offtarget_results(
            results_dir=str(results_dir), output_dir=output_dir, genome_species=genome_species
        )

        logger.info(f"Aggregation completed: {result_path}")

        mirna_summary: dict[str, Any] | None = None
        if mirna_analysis_files and (mirna_db or mirna_species):
            resolved_mirna_db = mirna_db or "mirgenedb"
            resolved_mirna_species = mirna_species or DEFAULT_MIRNA_SPECIES_ARGUMENT
            mirna_output = aggregate_mirna_results(
                results_dir=str(mirna_results_dir),
                output_dir=output_dir,
                mirna_db=resolved_mirna_db,
                mirna_species=resolved_mirna_species,
            )
            mirna_summary = {
                "mirna_db": resolved_mirna_db,
                "mirna_species": resolved_mirna_species,
                "output_dir": str(mirna_output),
                "analysis_files_processed": len(mirna_analysis_files),
            }

        return {
            "status": "completed",
            "analysis_files_processed": len(analysis_files),
            "summary_files_processed": len(summary_files),
            "species": species_list,
            "output_dir": str(result_path),
            "mirna": mirna_summary,
        }

    # Create empty final summary
    final_summary = output_path / "final_summary.txt"
    with final_summary.open("w") as handle:
        handle.write("No analysis results found to aggregate\n")

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
    """Aggregate miRNA seed analysis results from multiple candidates.

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
