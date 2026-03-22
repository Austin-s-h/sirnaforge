"""Typed Python bridge for ZFN Nextflow modules.

These helpers keep domain logic in Python while Nextflow handles process orchestration.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from sirnaforge.models.zfn import (
    DimerMode,
    GenomicAnnotationConfig,
    MatchOrientation,
    Strand,
    ZFNAlgorithm,
    ZFNDesignParameters,
    ZFNDesignResult,
    ZFNHalfSiteConstraints,
    ZFNOffTargetSite,
    ZFNShardingConfig,
    ZFNSpacerConstraints,
)
from sirnaforge.utils.fasta import load_fasta_contig_lengths
from sirnaforge.zfn.design import ZFNDesigner
from sirnaforge.zfn.rank import rank_sites
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher, build_zfn_shard_specs

RegionLiteral = Literal["exon", "promoter", "intron", "intergenic", "unknown"]


@dataclass(slots=True)
class ZFNShardRow:
    """One row in the ZFN shard manifest."""

    shard_id: str
    chrom: str
    core_start_1: int
    core_end_1: int
    scan_start_1: int
    scan_end_1: int
    max_mismatches: int


def _parse_bool(value: str | bool) -> bool:
    """Parse bool-like strings used by workflow parameters."""
    if isinstance(value, bool):
        return value
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _parse_csv_tokens(value: str) -> list[str]:
    """Split comma-delimited values, dropping empty tokens."""
    return [token.strip() for token in value.split(",") if token.strip()]


def make_zfn_shard_manifest(
    *,
    genome_fasta: Path,
    left_half_site: str,
    right_half_site: str,
    spacer_lengths: str,
    max_mismatches: int,
    sharding_enabled: str | bool,
    shard_chunk_mb: float,
    shard_overlap_bp: int,
    shard_chromosomes: str,
    output_tsv: Path,
) -> dict[str, int | bool]:
    """Build shard TSV for ZFN search.

    Sharding is applied whenever enabled, including for single-contig references.
    The direct Python searcher remains authoritative for overlap and chunk planning.
    """
    contig_lengths = load_fasta_contig_lengths(genome_fasta)
    params_obj = ZFNDesignParameters(
        left_half_site=left_half_site,
        right_half_site=right_half_site,
        search_space_fasta=str(genome_fasta),
        search_space_reference=None,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=max_mismatches),
        spacer_constraints=ZFNSpacerConstraints(
            allowed_spacer_lengths=[int(value) for value in _parse_csv_tokens(spacer_lengths)]
        ),
        sharding=ZFNShardingConfig(
            enabled=_parse_bool(sharding_enabled),
            chunk_size_bp=max(1, int(float(shard_chunk_mb) * 1_000_000)),
            overlap_bp=int(shard_overlap_bp),
            chromosomes=_parse_csv_tokens(shard_chromosomes),
        ),
    )
    shard_specs = build_zfn_shard_specs(contig_lengths, params_obj)
    if not shard_specs:
        raise ValueError("No shards generated. Check zfn_shard_chromosomes against FASTA contig names.")

    rows: list[ZFNShardRow] = []
    for shard in shard_specs:
        rows.append(
            ZFNShardRow(
                shard_id=shard.shard_id,
                chrom=shard.chrom,
                core_start_1=shard.core_start0 + 1,
                core_end_1=shard.core_end0,
                scan_start_1=shard.scan_start0 + 1,
                scan_end_1=shard.scan_end0,
                max_mismatches=max_mismatches,
            )
        )

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "shard_id",
                "chrom",
                "core_start_1",
                "core_end_1",
                "scan_start_1",
                "scan_end_1",
                "max_mismatches",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "shard_id": row.shard_id,
                    "chrom": row.chrom,
                    "core_start_1": row.core_start_1,
                    "core_end_1": row.core_end_1,
                    "scan_start_1": row.scan_start_1,
                    "scan_end_1": row.scan_end_1,
                    "max_mismatches": row.max_mismatches,
                }
            )

    return {
        "shards": len(rows),
        "contigs": len({row.chrom for row in rows}),
        "sharding_active": _parse_bool(sharding_enabled),
        "overlap_bp": rows[0].scan_end_1 - rows[0].core_end_1 if rows else 0,
    }


def run_zfn_shard_search(
    *,
    shard_id: str,
    shard_chrom: str,
    scan_start_1: int,
    scan_end_1: int,
    core_start_1: int | None = None,
    core_end_1: int | None = None,
    shard_max_mismatches: int,
    left_half_site: str,
    right_half_site: str,
    genome_fasta: Path,
    algorithm: str,
    dimer_mode: str,
    spacer_lengths: str,
    annotation_file: Path | None,
    output_sites_csv: Path,
    output_summary_json: Path,
) -> dict[str, int | str]:
    """Execute one shard search and persist shard-scoped outputs.

    The search scans ``scan_start_1..scan_end_1`` (1-based, inclusive) for
    sequence context, but only sites whose coordinates fall entirely within the
    *core* window (``core_start_1..core_end_1``) are written to output.  This
    prevents duplicate sites across overlapping shards.  When ``core_start_1``
    or ``core_end_1`` are not provided they default to the scan boundaries.
    """
    effective_core_start = core_start_1 if core_start_1 is not None else scan_start_1
    effective_core_end = core_end_1 if core_end_1 is not None else scan_end_1
    spacer_list = [int(s) for s in _parse_csv_tokens(spacer_lengths)]

    params_obj = ZFNDesignParameters(
        left_half_site=left_half_site,
        right_half_site=right_half_site,
        search_space_fasta=str(genome_fasta),
        search_space_reference=None,
        algorithm=ZFNAlgorithm(algorithm),
        dimer_mode=DimerMode(dimer_mode),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=spacer_list),
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=shard_max_mismatches),
        sharding=ZFNShardingConfig(enabled=False, chromosomes=[shard_chrom]),
    )

    annotation = None
    if annotation_file is not None and annotation_file.name != "NO_ANNOTATION":
        annotation = GenomicAnnotationConfig(annotation_path=str(annotation_file))

    searcher = ExhaustiveZFNOffTargetSearcher()
    designer = ZFNDesigner(searcher=searcher)
    filtered_sites = searcher.search_region(
        params=params_obj,
        chrom=shard_chrom,
        scan_start0=scan_start_1 - 1,
        scan_end0=scan_end_1,
        core_start0=effective_core_start - 1,
        core_end0=effective_core_end,
        annotation=annotation,
        top_n_sites=None,
    )
    candidate = designer.build_candidate(params_obj, filtered_sites)

    filtered_result = ZFNDesignResult(
        parameters=params_obj,
        annotation=annotation,
        candidates=[candidate],
        off_target_sites=filtered_sites,
        processing_time_s=0.0,
        tool_versions=designer.tool_versions(),
    )

    output_sites_csv.parent.mkdir(parents=True, exist_ok=True)
    filtered_result.save_offtargets_csv(str(output_sites_csv))

    summary = filtered_result.get_summary()
    summary["shard_id"] = shard_id
    summary["chrom"] = shard_chrom
    summary["scan_start_1"] = scan_start_1
    summary["scan_end_1"] = scan_end_1
    candidates = [candidate.model_dump(mode="json") for candidate in filtered_result.candidates]

    output_summary_json.parent.mkdir(parents=True, exist_ok=True)
    output_summary_json.write_text(
        json.dumps({"candidates": candidates, "summary": summary}, indent=2, default=str),
        encoding="utf-8",
    )

    return {"shard_id": shard_id, "sites": len(filtered_sites), "chrom": shard_chrom}


def _site_from_row(row: dict[str, str]) -> ZFNOffTargetSite:
    """Deserialize one CSV row into a typed off-target site."""
    nearest_gene_raw = row.get("nearest_gene")
    nearest_gene = nearest_gene_raw if nearest_gene_raw not in {None, "", "nan", "NaN"} else None
    raw_region = (row.get("region") or "unknown").lower()
    if raw_region not in {"exon", "promoter", "intron", "intergenic", "unknown"}:
        raw_region = "unknown"
    typed_region: RegionLiteral = cast(RegionLiteral, raw_region)

    return ZFNOffTargetSite(
        site_id=row["site_id"],
        chrom=row["chrom"],
        start_1based=int(row["start_1based"]),
        end_1based=int(row["end_1based"]),
        strand=Strand(row["strand"]),
        orientation=MatchOrientation(row["orientation"]),
        spacer_len=int(row["spacer_len"]),
        sequence=row["sequence"],
        left_mismatches=int(row["left_mismatches"]),
        right_mismatches=int(row["right_mismatches"]),
        total_mismatches=int(row["total_mismatches"]),
        score=float(row["score"]),
        region=typed_region,
        nearest_gene=nearest_gene,
        left_aligned=row["left_aligned"],
        right_aligned=row["right_aligned"],
    )


def _write_sites_csv(path: Path, sites: list[ZFNOffTargetSite]) -> None:
    """Write ranked sites to canonical CSV schema."""
    fieldnames = [
        "site_id",
        "chrom",
        "start_1based",
        "end_1based",
        "strand",
        "orientation",
        "spacer_len",
        "sequence",
        "left_mismatches",
        "right_mismatches",
        "total_mismatches",
        "score",
        "region",
        "nearest_gene",
        "left_aligned",
        "right_aligned",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for site in sites:
            writer.writerow(
                {
                    "site_id": site.site_id,
                    "chrom": site.chrom,
                    "start_1based": site.start_1based,
                    "end_1based": site.end_1based,
                    "strand": site.strand.value,
                    "orientation": site.orientation.value,
                    "spacer_len": site.spacer_len,
                    "sequence": site.sequence,
                    "left_mismatches": site.left_mismatches,
                    "right_mismatches": site.right_mismatches,
                    "total_mismatches": site.total_mismatches,
                    "score": site.score,
                    "region": site.region,
                    "nearest_gene": site.nearest_gene,
                    "left_aligned": site.left_aligned,
                    "right_aligned": site.right_aligned,
                }
            )


def aggregate_zfn_shard_results(
    *,
    shard_csv_glob: str,
    output_sites_csv: Path,
    output_summary_json: Path,
) -> dict[str, int]:
    """Merge shard CSVs, deduplicate coordinate/orientation collisions, and rank globally."""
    files = sorted(Path().glob(shard_csv_glob))

    all_sites: list[ZFNOffTargetSite] = []
    for file_path in files:
        with file_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                all_sites.append(_site_from_row(row))

    deduped_by_key: dict[tuple[str, int, int, str], ZFNOffTargetSite] = {}
    for site in sorted(all_sites, key=lambda item: item.score, reverse=True):
        key = (site.chrom, site.start_1based, site.end_1based, site.orientation.value)
        if key not in deduped_by_key:
            deduped_by_key[key] = site

    ranked = rank_sites(list(deduped_by_key.values()))

    output_sites_csv.parent.mkdir(parents=True, exist_ok=True)
    _write_sites_csv(output_sites_csv, ranked)

    output_summary_json.parent.mkdir(parents=True, exist_ok=True)
    output_summary_json.write_text(
        json.dumps(
            {
                "candidates": [],
                "summary": {
                    "off_target_sites": len(ranked),
                    "shards": len(files),
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return {"off_target_sites": len(ranked), "shards": len(files)}
