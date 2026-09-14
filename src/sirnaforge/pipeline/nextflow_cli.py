"""Command-line entry points used by embedded Nextflow modules."""

import json
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from sirnaforge.config import DEFAULT_MIRNA_CANONICAL_SPECIES
from sirnaforge.core.off_target import (
    aggregate_mirna_results,
    aggregate_offtarget_results,
    build_bwa_index,
    parse_fasta_file,
    run_bwa_alignment_analysis,
    validate_index_files,
)
from sirnaforge.core.screening_evidence import (
    EvidenceProducer,
    Reconciliation,
    build_plan,
    collect_evidence,
    guide_set_digest,
    reconcile,
    write_evidence,
    write_reconciliation,
)
from sirnaforge.models.evidence import EvidenceStatus, ScreeningEvidenceEntry, ScreeningPlan
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

DEFAULT_MIRNA_SPECIES_ARGUMENT = ",".join(DEFAULT_MIRNA_CANONICAL_SPECIES)

#: Guide-set digest stamped on a locally-derived fallback plan, used only when no ``evidence_plan``
#: file was threaded through a call. #100's real digests are 16-character sha256 prefixes, so this
#: sentinel never collides with one: a fallback plan entry can therefore only ever explain an
#: envelope that was NOT observed (reconciling to failed), never falsely match one that was.
_NO_PLAN_DIGEST_PLACEHOLDER = "no-evidence-plan-supplied"


def _reconcile_species(
    *,
    evidence_plan: str | None,
    transcriptome_species: Sequence[str] = (),
    mirna_species: Sequence[str] = (),
    search_root: str | Path,
) -> Reconciliation:
    """Reconcile the species this call expected against whatever evidence was actually staged.

    Prefers the real plan threaded in via ``evidence_plan`` (#100's expected-plan contract, built by
    ``workflow.py`` before any reference can be dropped). With none supplied -- a bare Nextflow run,
    or a pipeline invocation that aborted before a plan file existed -- falls back to a plan built
    from this call's own species arguments, so a species with no evidence still produces a
    machine-readable failure record instead of vanishing with only a ``final_summary.txt`` behind.
    The fallback plan's digest is borrowed from whatever evidence was actually observed (a run
    screens one guide set at a time, by construction) so a real envelope still joins against it;
    with nothing observed at all there is no digest to borrow, and the placeholder is used purely
    to explain the absence -- it can never falsely match a real envelope written later.
    """
    observed = collect_evidence(search_root)
    if evidence_plan is not None:
        plan = ScreeningPlan.model_validate_json(Path(evidence_plan).read_text())
    else:
        digest = observed[0].entry.guide_set_digest if observed else _NO_PLAN_DIGEST_PLACEHOLDER
        plan = build_plan(
            guide_set_digest=digest,
            transcriptome=[(species, None) for species in transcriptome_species],
            mirna_species=mirna_species,
        )
    return reconcile(plan, observed)


def build_bwa_index_cli(fasta_file: str, species: str, output_dir: str = ".") -> dict[str, Any]:
    """Build a BWA-MEM2 index for one screening reference.

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


def offtarget_analysis_cli(
    species: str,
    index_prefix: str,
    candidates_file: str,
    output_dir: str = ".",
    max_hits: int | None = None,
    bwa_k: int = 12,
    bwa_T: int = 15,
    seed_start: int = 2,
    seed_end: int = 8,
) -> dict[str, Any]:
    """Align one species' reference for OFFTARGET_ANALYSIS, or record that it could not be aligned.

    A named index prefix that does not resolve to a usable BWA-MEM2 index used to reach the aligner,
    which logged a failure, returned no hits, and published a header-only table -- a completed screen
    with zero hits, indistinguishable from a clean one. The species is now published as an EMPTY
    analysis file plus a failed summary, which ``aggregate_offtarget_results`` reports as a per-species
    rejection rather than as a clean result. Either way, this function owns the failed-vs-completed
    decision, so it also writes the #100 ``transcriptome_<species>_evidence.json`` envelope on both
    branches -- directly here on the missing-index branch, and via ``run_bwa_alignment_analysis``'s
    ``evidence_dir`` on the aligned branch.

    Args:
        species: Species identifier this reference belongs to
        index_prefix: BWA-MEM2 index prefix as resolved for this species
        candidates_file: FASTA of candidate guides to screen
        output_dir: Directory to write ``<species>_analysis.tsv`` and ``<species>_summary.json``
        max_hits: Maximum hits per candidate (``None`` = exhaustive)
        bwa_k: BWA seed length
        bwa_T: BWA minimum score threshold
        seed_start: Seed region start (1-based)
        seed_end: Seed region end (1-based)

    Returns:
        Dictionary with the per-species status and the published file paths
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    analysis_file = output_path / f"{species}_analysis.tsv"
    summary_file = output_path / f"{species}_summary.json"

    if not validate_index_files(index_prefix, "bwa-mem2"):
        error = (
            f"No usable BWA-MEM2 index at prefix '{index_prefix}': the index files are missing, empty, "
            "or the path points at something that is not an index (a FASTA left behind by a failed "
            "index build looks exactly like this). No alignment was attempted."
        )
        logger.error(f"{species}: {error}")
        # Deliberately empty, not header-only: a header-only table is a screen that found nothing.
        analysis_file.write_bytes(b"")
        with summary_file.open("w") as handle:
            json.dump({"species": species, "status": "failed", "error": error}, handle, indent=2)
        # This branch owns the failed-vs-completed decision, so it is the one that must write the
        # #100 evidence envelope: a species whose index never resolved must never reconcile as
        # completed just because the aligner was never reached.
        digest = guide_set_digest(candidates_file)
        try:
            submitted_guides: int | None = len(parse_fasta_file(candidates_file))
        except (OSError, ValueError):
            submitted_guides = None
        write_evidence(
            output_dir,
            producer=EvidenceProducer.OFFTARGET_ANALYSIS,
            entry=ScreeningEvidenceEntry(
                channel=ScreeningChannel.TRANSCRIPTOME,
                species=species,
                guide_set_digest=digest,
                status=EvidenceStatus.FAILED,
                submitted_guide_digest=digest,
                submitted_guides=submitted_guides,
                detail=error,
            ),
        )
        return {"species": species, "status": "failed", "error": error, "analysis_file": str(analysis_file)}

    run_bwa_alignment_analysis(
        candidates_file=candidates_file,
        index_prefix=index_prefix,
        species=species,
        output_dir=output_dir,
        max_hits=max_hits,
        bwa_k=bwa_k,
        bwa_T=bwa_T,
        seed_start=seed_start,
        seed_end=seed_end,
        evidence_dir=output_dir,
    )

    # run_bwa_alignment_analysis names its outputs after the candidates file; the process declares
    # them per species.
    stem = Path(candidates_file).stem
    (output_path / f"{stem}_{species}_analysis.tsv").replace(analysis_file)
    (output_path / f"{stem}_{species}_summary.json").replace(summary_file)

    return {"species": species, "status": "completed", "analysis_file": str(analysis_file)}


def aggregate_results_cli(  # noqa: PLR0912
    transcriptome_species: str,
    output_dir: str = ".",
    mirna_db: str | None = None,
    mirna_species: str | None = None,
    analysis_files: list[str] | None = None,
    summary_files: list[str] | None = None,
    *,
    evidence_plan: str | None = None,
) -> dict[str, Any]:
    """Aggregate off-target analysis results from multiple candidates and references.

    Args:
        transcriptome_species: Comma-separated list of species screened
        output_dir: Directory to write aggregated results
        mirna_db: The database that provided the reference
        mirna_species: The species code for the matching miRNA
        analysis_files: Optional explicit list of staged analysis files to aggregate
        summary_files: Optional explicit list of staged summary files to aggregate
        evidence_plan: Path to the serialized #100 ``ScreeningPlan`` (``workflow.py``'s
            ``screening_plan.json``), additive and keyword-only. Reconciliation against it is
            written unconditionally on every return path, including the zero-files branch, so a
            pipeline that aborted before producing a single analysis file still leaves a
            machine-readable ``evidence.json`` behind instead of only a ``final_summary.txt``.

    Returns:
        Dictionary with aggregation statistics
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    species_list = [s.strip() for s in transcriptome_species.split(",") if s.strip()]

    # Collect all analysis and summary files from current directory (Nextflow stages them)
    current_dir = Path()
    resolved_analysis_files = (
        [Path(path) for path in analysis_files]
        if analysis_files is not None
        else list(current_dir.rglob("*_analysis.tsv")) + list(current_dir.rglob("mirna_analysis.tsv"))
    )
    resolved_summary_files = (
        [Path(path) for path in summary_files]
        if summary_files is not None
        else list(current_dir.rglob("*_summary.json")) + list(current_dir.rglob("mirna_summary.json"))
    )

    logger.info(f"Found {len(resolved_analysis_files)} analysis files and {len(resolved_summary_files)} summary files")

    if resolved_analysis_files or resolved_summary_files:
        # Create a temporary results directory structure for aggregation
        results_dir = Path("temp_results")
        results_dir.mkdir(exist_ok=True)

        # Organize files by species (extract from filename)
        for species in species_list:
            species_dir = results_dir / species
            species_dir.mkdir(exist_ok=True)

            # Copy relevant files to species directory
            for f in resolved_analysis_files:
                name_lower = f.name.lower()
                if name_lower.startswith("mirna"):
                    continue
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

            for f in resolved_summary_files:
                name_lower = f.name.lower()
                if name_lower.startswith("mirna"):
                    continue
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

        # Also persist miRNA batch results (if present) under dedicated directory
        mirna_results_dir = results_dir / "mirna"
        mirna_results_dir.mkdir(exist_ok=True)
        mirna_analysis_files: list[Path] = []
        for f in resolved_analysis_files:
            if "mirna" not in f.name.lower():
                continue
            dest_name = f.name
            if not dest_name.endswith("_mirna_analysis.tsv"):
                dest_name = "batch_mirna_analysis.tsv"
            dest_path = mirna_results_dir / dest_name
            shutil.copy(f, dest_path)
            mirna_analysis_files.append(dest_path)

        for f in resolved_summary_files:
            if "mirna" not in f.name.lower():
                continue
            dest_name = f.name
            if not dest_name.endswith("_mirna_summary.json"):
                dest_name = "batch_mirna_summary.json"
            dest_path = mirna_results_dir / dest_name
            shutil.copy(f, dest_path)

        # Run aggregation using core function
        result_path = aggregate_offtarget_results(
            results_dir=str(results_dir), output_dir=output_dir, transcriptome_species=transcriptome_species
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

        reconciliation = _reconcile_species(
            evidence_plan=evidence_plan,
            transcriptome_species=species_list,
            search_root=current_dir,
        )
        reconciliation_file = write_reconciliation(output_dir, reconciliation)

        return {
            "status": "completed",
            "analysis_files_processed": len(resolved_analysis_files),
            "summary_files_processed": len(resolved_summary_files),
            "species": species_list,
            "output_dir": str(result_path),
            "mirna": mirna_summary,
            "evidence": reconciliation.evidence.model_dump(mode="json"),
            "reconciliation_file": str(reconciliation_file),
        }

    # Create empty final summary
    final_summary = output_path / "final_summary.txt"
    with final_summary.open("w") as handle:
        handle.write("No analysis results found to aggregate\n")

    logger.warning("No files to aggregate")

    # Optional output presence is not completion evidence (#100): even though nothing was staged to
    # aggregate, every species this call expected to see still gets a reconciled -- here, failed --
    # entry, so a pipeline that aborted before a single analysis file existed leaves the same
    # machine-readable record behind as one that completed and legitimately found nothing.
    reconciliation = _reconcile_species(
        evidence_plan=evidence_plan,
        transcriptome_species=species_list,
        search_root=current_dir,
    )
    reconciliation_file = write_reconciliation(output_dir, reconciliation)

    return {
        "status": "empty",
        "analysis_files_processed": 0,
        "summary_files_processed": 0,
        "species": [],
        "output_dir": str(output_path),
        "evidence": reconciliation.evidence.model_dump(mode="json"),
        "reconciliation_file": str(reconciliation_file),
    }


def aggregate_mirna_results_cli(
    mirna_db: str,
    mirna_species: str,
    results_dir: str = ".",
    output_dir: str = ".",
    *,
    evidence_plan: str | None = None,
) -> dict[str, Any]:
    """Aggregate miRNA seed analysis results from multiple candidates.

    Args:
        mirna_db: miRNA database name used for analysis
        mirna_species: Comma-separated list of species analyzed
        results_dir: Directory containing individual miRNA results
        output_dir: Directory to write aggregated results
        evidence_plan: Path to the serialized #100 ``ScreeningPlan``, additive and keyword-only.
            ``aggregate_results.nf`` does not call this entry point, but it gets the same
            unconditional reconciliation treatment as ``aggregate_results_cli`` so a caller of it
            is not left with a weaker evidence contract than the transcriptome channel.

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

    reconciliation = _reconcile_species(
        evidence_plan=evidence_plan,
        mirna_species=[s.strip() for s in mirna_species.split(",") if s.strip()],
        search_root=results_dir,
    )
    reconciliation_file = write_reconciliation(output_dir, reconciliation)

    return {
        "status": "completed",
        "mirna_database": mirna_db,
        "species": mirna_species.split(","),
        "total_hits": stats.get("total_mirna_hits", 0),
        "candidates_analyzed": stats.get("total_candidates", 0),
        "output_dir": str(result_path),
        "evidence": reconciliation.evidence.model_dump(mode="json"),
        "reconciliation_file": str(reconciliation_file),
    }
