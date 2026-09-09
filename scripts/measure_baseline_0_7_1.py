#!/usr/bin/env python3
"""Recompute every 0.7.1 baseline measurement from a frozen canonical run.

Reads only the frozen artifacts of the run (candidate CSV, per-hit TSVs, aggregate JSON)
plus the screening reference FASTA, and writes ``measurements.json`` next to them. Nothing
here mutates the run, and no siRNAforge source is modified to make a number reachable: gate
verdicts and hit classes that the 0.7.1 writers do not persist are recomputed here with the
repo's own pure functions.

Usage (from the repo root):

    uv run python scripts/measure_baseline_0_7_1.py \
        --run-dir work/baseline_0.7.1/run \
        --reference work/baseline_0.7.1/reference/GRCh38_cdna_primary_1tx_per_gene_plus_TP53.fa \
        --output work/baseline_0.7.1/measurements.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from sirnaforge.core.hit_classification import ClassificationContext, HitClass, classify_hit
from sirnaforge.core.repeat_detection import (
    DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
    RepeatDetector,
    normalize_guide_sequence,
)
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.models.sirna import FilterCriteria, OffTargetFilterCriteria, SiRNACandidate

POSTSCREEN_SIRNA_TERMS = {
    "off_target": "score_off_target",
    "target_accessibility": "score_target_accessibility",
    "asymmetry": "score_asymmetry",
    "gc_content": "score_gc_content",
}
POSTSCREEN_MIRNA_EXTRA = {
    "ago_start": "score_ago_start",
    "pos1_mismatch": "score_pos1_mismatch",
    "supp_13_16": "score_supp_13_16",
}


def sha256(path: Path) -> str:
    """SHA-256 of a file, streamed."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def variance_shares(frame: pd.DataFrame, term_columns: dict[str, str]) -> dict[str, Any]:
    """Per-term cov(contribution, composite)/var(composite), which decomposes exactly."""
    composite = frame["composite_score"].astype(float)
    var = float(composite.var(ddof=1))
    shares: dict[str, float] = {}
    for term, column in term_columns.items():
        contribution = frame[column].astype(float)
        cov = float(contribution.cov(composite))
        shares[term] = 100.0 * cov / var
    return {
        "n": int(len(frame)),
        "composite_variance": var,
        "composite_sd": float(composite.std(ddof=1)),
        "composite_mean": float(composite.mean()),
        "shares_percent": shares,
        "checksum_percent": sum(shares.values()),
    }


def read_hits(path: Path) -> pd.DataFrame:
    """Load a per-hit TSV, or an empty frame when it is absent/empty."""
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


@dataclass
class Reconstruction:
    """Per-candidate counters recomputed from the per-hit table."""

    label: str
    per_query: dict[str, dict[str, int]]
    class_counts: Counter[str]
    symbol_lookup_missing: int
    species_index_missing: int
    resolved_symbols: int
    total_hits: int


def reconstruct(
    hits: pd.DataFrame,
    context: ClassificationContext,
    qname_to_guide: dict[str, str],
    label: str,
) -> Reconstruction:
    """Re-classify every hit and accumulate per-query counters."""
    per_query: dict[str, dict[str, int]] = {}
    class_counts: Counter[str] = Counter()
    symbol_missing = 0
    index_missing = 0
    resolved = 0

    for row in hits.to_dict("records"):
        qname = str(row["qname"])
        guide = qname_to_guide.get(qname, str(row.get("qseq", "")))
        verdict = classify_hit(row, guide, context)
        class_counts[verdict.hit_class.value] += 1
        if verdict.matched_symbol:
            resolved += 1
        if verdict.symbol_lookup_missing:
            symbol_missing += 1
        if verdict.species_index_missing:
            index_missing += 1

        bucket = per_query.setdefault(
            qname,
            {
                "on_target": 0,
                "ortholog": 0,
                "repeat": 0,
                "off_target": 0,
                "off_target_nm_le2": 0,
                "off_target_nm_ge3": 0,
                "human_0mm": 0,
                "human_1mm": 0,
                "human_2mm": 0,
                "seed_0mm": 0,
                "human_total": 0,
            },
        )
        bucket[verdict.hit_class.value] += 1
        if verdict.hit_class is not HitClass.OFF_TARGET:
            continue

        nm = int(row["nm"])
        seed_mm = int(row["seed_mismatches"])
        species = str(row.get("species") or "")
        treated_as_human = species.lower() in {"human", "hsa", "homo_sapiens", ""}
        if nm <= 2:
            bucket["off_target_nm_le2"] += 1
        else:
            bucket["off_target_nm_ge3"] += 1
        if treated_as_human:
            bucket["human_total"] += 1
            if nm in (0, 1, 2):
                bucket[f"human_{nm}mm"] += 1
        if seed_mm == 0:
            bucket["seed_0mm"] += 1

    return Reconstruction(
        label=label,
        per_query=per_query,
        class_counts=class_counts,
        symbol_lookup_missing=symbol_missing,
        species_index_missing=index_missing,
        resolved_symbols=resolved,
        total_hits=int(len(hits)),
    )


def gate_independent_failures(
    candidates: pd.DataFrame,
    recon: Reconstruction,
    mirna_by_query: dict[str, dict[str, int]],
    representative: dict[str, str],
    filters: FilterCriteria,
    offtarget_filters: OffTargetFilterCriteria,
) -> dict[str, dict[str, Any]]:
    """Evaluate every gate independently on every candidate row.

    The shipped writers record only the FIRST failing gate per candidate, and a later stage
    overwrites the verdict, so the per-gate populations cannot be read off the CSV. Each gate
    below is evaluated on its own inputs against its own configured threshold.
    """
    results: dict[str, dict[str, Any]] = {}

    any_failure = pd.Series(False, index=candidates.index)

    def record(name: str, mask: pd.Series, threshold: Any, reads: str) -> None:
        nonlocal any_failure
        results[name] = {
            "independent_failures": int(mask.sum()),
            "fraction_of_rows": float(mask.mean()),
            "threshold": threshold,
            "reads": reads,
            "active": threshold is not None,
        }
        if name != "EXCESS_OFF_TARGETS_nm_le2_scope":
            any_failure = any_failure | mask.fillna(False).astype(bool)

    gc = candidates["gc_content"].astype(float)
    record(
        "GC_OUT_OF_RANGE",
        (gc < filters.gc_min) | (gc > filters.gc_max),
        f"{filters.gc_min}-{filters.gc_max}",
        "gc_content",
    )

    poly = (
        candidates["guide_sequence"]
        .astype(str)
        .apply(lambda seq: any(base * (filters.max_poly_runs + 1) in seq.upper() for base in "ACGU T".replace(" ", "")))
    )
    record("POLY_RUNS", poly, filters.max_poly_runs, "guide_sequence homopolymer run")

    paired = candidates["paired_fraction"].astype(float)
    record("EXCESS_PAIRING", paired > filters.max_paired_fraction, filters.max_paired_fraction, "paired_fraction")

    asym = candidates["asymmetry_score"].astype(float)
    record("LOW_ASYMMETRY", asym < filters.min_asymmetry_score, filters.min_asymmetry_score, "asymmetry_score")

    empirical = candidates["empirical_score"].astype(float)
    record(
        "LOW_EMPIRICAL_SCORE",
        empirical < filters.min_empirical_score,
        filters.min_empirical_score,
        "empirical_score",
    )

    if filters.min_isoform_coverage is None:
        record(
            "LOW_ISOFORM_COVERAGE",
            pd.Series(False, index=candidates.index),
            None,
            "isoform_coverage (gate off by default)",
        )
    else:
        coverage = candidates["isoform_coverage"].astype(float)
        record(
            "LOW_ISOFORM_COVERAGE",
            coverage.notna() & (coverage < filters.min_isoform_coverage),
            filters.min_isoform_coverage,
            "isoform_coverage",
        )

    record(
        "REPEAT_ELEMENT",
        candidates["repeat_flagged"].fillna(False).astype(bool),
        "repeat_transcript_fraction",
        "repeat_flagged",
    )

    def counter(field: str) -> pd.Series:
        return candidates["id"].map(lambda cid: recon.per_query.get(representative.get(cid, cid), {}).get(field, 0))

    def mirna_counter(field: str) -> pd.Series:
        return candidates["id"].map(lambda cid: mirna_by_query.get(representative.get(cid, cid), {}).get(field, 0))

    for name, field, threshold in (
        ("TRANSCRIPTOME_PERFECT_MATCH", "human_0mm", offtarget_filters.max_transcriptome_hits_0mm),
        ("TRANSCRIPTOME_1MM", "human_1mm", offtarget_filters.max_transcriptome_hits_1mm),
        ("TRANSCRIPTOME_2MM", "human_2mm", offtarget_filters.max_transcriptome_hits_2mm),
        ("TRANSCRIPTOME_SEED_PERFECT", "seed_0mm", offtarget_filters.max_transcriptome_seed_perfect),
        ("EXCESS_OFF_TARGETS", "off_target", offtarget_filters.max_off_target_count),
    ):
        values = counter(field)
        mask = pd.Series(False, index=candidates.index) if threshold is None else values > threshold
        record(name, mask, threshold, f"recomputed {field}")

    perfect_seed = mirna_counter("seed_0mm")
    threshold = offtarget_filters.max_mirna_perfect_seed
    record(
        "MIRNA_PERFECT_SEED",
        pd.Series(False, index=candidates.index) if threshold is None else perfect_seed > threshold,
        threshold,
        "recomputed mirna seed_mismatches == 0 (human)",
    )
    record(
        "HIGH_RISK_MIRNA",
        mirna_counter("high_risk") > 0
        if offtarget_filters.fail_on_high_risk_mirna
        else pd.Series(False, index=candidates.index),
        offtarget_filters.fail_on_high_risk_mirna,
        "recomputed mirna seed 0mm and offtarget_score < 5.0 (human)",
    )

    total = counter("human_total") + mirna_counter("human_total")
    threshold = offtarget_filters.max_total_offtarget_hits
    record(
        "TOTAL_OFFTARGETS",
        pd.Series(False, index=candidates.index) if threshold is None else total > threshold,
        threshold,
        "recomputed human transcriptome + human miRNA hits",
    )

    # Same cap, rescoped to genuine nm <= 2 hits (issue #101 item 3).
    cap = offtarget_filters.max_off_target_count
    values = counter("off_target_nm_le2")
    record(
        "EXCESS_OFF_TARGETS_nm_le2_scope",
        pd.Series(False, index=candidates.index) if cap is None else values > cap,
        cap,
        "recomputed off_target hits restricted to nm <= 2",
    )
    results["__union__"] = {
        "rows_failing_at_least_one_gate": int(any_failure.sum()),
        "rows_failing_no_gate": int((~any_failure).sum()),
        "rows": int(len(candidates)),
    }
    return results


def mirna_counters(hits: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Per-query miRNA counters, matching what _integrate_offtarget_results accumulates."""
    out: dict[str, dict[str, int]] = {}
    if hits.empty:
        return out
    for row in hits.to_dict("records"):
        qname = str(row["qname"])
        bucket = out.setdefault(qname, {"total": 0, "human_total": 0, "seed_0mm": 0, "seed_1mm": 0, "high_risk": 0})
        bucket["total"] += 1
        species = str(row.get("species") or "").lower()
        is_human = species in {"human", "hsa", "homo_sapiens", ""}
        if is_human:
            bucket["human_total"] += 1
        seed_mm = int(row["seed_mismatches"])
        if seed_mm == 0:
            bucket["seed_0mm"] += 1
            if float(row["offtarget_score"]) < 5.0:
                bucket["high_risk"] += 1
        elif seed_mm == 1:
            bucket["seed_1mm"] += 1
    return out


def main() -> None:  # noqa: PLR0912
    """Recompute the baseline measurements and write measurements.json."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True, help="Screened reference FASTA")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gc-min", type=float, default=30.0, help="as-run gc_min (workflow CLI default)")
    parser.add_argument("--gc-max", type=float, default=60.0, help="as-run gc_max (workflow CLI default)")
    parser.add_argument("--query-gene-id", default="ENSG00000141510")
    parser.add_argument("--query-gene-symbol", default="TP53")
    args = parser.parse_args()

    run_dir: Path = args.run_dir
    candidates_csv = run_dir / "sirnaforge" / "candidates_all.csv"
    aggregated = run_dir / "off_target" / "results" / "aggregated"
    offtargets_tsv = aggregated / "combined_offtargets.tsv"
    mirna_tsv = aggregated / "combined_mirna_hits.tsv"
    combined_summary = aggregated / "combined_summary.json"
    mirna_summary = aggregated / "combined_mirna_summary.json"
    workflow_summary = run_dir / "logs" / "workflow_summary.json"

    candidates = pd.read_csv(candidates_csv)
    offtargets = read_hits(offtargets_tsv)
    mirna_hits = read_hits(mirna_tsv)

    measurements: dict[str, Any] = {
        "run_dir": str(run_dir),
        "artifacts": {},
        "shape": {},
        "variance_share": {},
        "offtarget_scope": {},
        "gates": {},
        "column_redundancy": {},
    }

    for path in sorted(run_dir.rglob("*")):
        if path.is_file() and path.stat().st_size > 0 and ".nextflow" not in path.parts:
            measurements["artifacts"][str(path.relative_to(run_dir))] = {
                "bytes": path.stat().st_size,
                "sha256": sha256(path) if path.stat().st_size < 200_000_000 else "skipped_large",
            }

    # ---------------- C. shape ----------------
    guides = candidates["guide_sequence"].astype(str).map(normalize_guide_sequence)
    rows_per_guide = guides.value_counts()
    measurements["shape"] = {
        "candidate_rows": int(len(candidates)),
        "distinct_guide_sequences": int(rows_per_guide.size),
        "rows_per_guide_median": float(statistics.median(rows_per_guide.tolist())),
        "rows_per_guide_max": int(rows_per_guide.max()),
        "rows_per_guide_min": int(rows_per_guide.min()),
        "offtarget_alignments": int(len(offtargets)),
        "mirna_hits": int(len(mirna_hits)),
        "distinct_transcripts": int(candidates["transcript_id"].nunique()),
        "passing_rows": int((candidates["passes_filters"].astype(str) == "PASS").sum()),
        "screened_rows": int(candidates["off_target_screened"].fillna(False).astype(bool).sum()),
        "scored_after_screening_rows": int(candidates["scored_after_screening"].fillna(False).astype(bool).sum()),
    }

    # ---------------- A. post-screen variance share ----------------
    scored = candidates[
        candidates["scored_after_screening"].fillna(False).astype(bool) & candidates["composite_score"].notna()
    ].copy()
    vectors = sorted(set(scored["weight_vector"].dropna().astype(str)))
    measurements["variance_share"]["weight_vectors_present"] = vectors
    for vector in vectors:
        subset = scored[scored["weight_vector"].astype(str) == vector]
        columns = dict(POSTSCREEN_SIRNA_TERMS)
        if vector.endswith("mirna_v4"):
            columns.update(POSTSCREEN_MIRNA_EXTRA)
        columns = {term: col for term, col in columns.items() if col in subset.columns}
        measurements["variance_share"][vector] = variance_shares(subset, columns)

    # design-stage comparison, same estimator
    design_scored = candidates[candidates["design_score"].notna()].copy()
    if not design_scored.empty:
        measurements["variance_share"]["design_stage_note"] = (
            "design_score contributions are not persisted per term, so no design-stage decomposition is computed here"
        )

    # ---------------- reconstruction of per-hit classification ----------------
    qname_to_guide = dict(zip(candidates["id"].astype(str), candidates["guide_sequence"].astype(str), strict=True))
    guide_to_ids: dict[str, list[str]] = {}
    for cid, guide in zip(candidates["id"].astype(str), guides, strict=True):
        guide_to_ids.setdefault(guide, []).append(cid)
    representative = {cid: ids[0] for ids in guide_to_ids.values() for cid in ids}

    on_target_transcript_ids = frozenset(
        str(t).split(".")[0].upper() for t in candidates["transcript_id"].dropna().unique()
    )
    query_gene_symbols = frozenset({args.query_gene_symbol.upper()})

    empty_context = ClassificationContext(
        query_gene_ids=frozenset(),
        query_gene_symbols=query_gene_symbols,
        on_target_transcript_ids=on_target_transcript_ids,
        query_species="human",
        index=TranscriptGeneIndex(),
        repeat_flagged_guides=frozenset(),
        requested_species=frozenset({"human"}),
    )
    as_run = reconstruct(offtargets, empty_context, qname_to_guide, "as_run_no_index")

    enriched_index = TranscriptGeneIndex()
    species_index = enriched_index.build("human", args.reference)
    enriched_context = ClassificationContext(
        query_gene_ids=frozenset({args.query_gene_id.upper()}),
        query_gene_symbols=query_gene_symbols,
        on_target_transcript_ids=on_target_transcript_ids,
        query_species="human",
        index=enriched_index,
        repeat_flagged_guides=frozenset(),
        requested_species=frozenset({"human"}),
    )
    enriched = reconstruct(offtargets, enriched_context, qname_to_guide, "with_transcript_index")

    # The run itself could not run repeat detection (it needs the query species' cDNA, which
    # only the transcriptome-materialisation path provides). Scan offline against the same
    # reference the aligner used, so the REPEAT class is measured rather than assumed absent.
    repeat_scan = RepeatDetector(threshold_fraction=DEFAULT_REPEAT_TRANSCRIPT_FRACTION).scan(
        set(guides), args.reference
    )
    repeat_guides = repeat_scan.repeat_sequences
    repeat_context = ClassificationContext(
        query_gene_ids=frozenset({args.query_gene_id.upper()}),
        query_gene_symbols=query_gene_symbols,
        on_target_transcript_ids=on_target_transcript_ids,
        query_species="human",
        index=enriched_index,
        repeat_flagged_guides=repeat_guides,
        requested_species=frozenset({"human"}),
    )
    with_repeats = reconstruct(offtargets, repeat_context, qname_to_guide, "with_index_and_repeats")

    resolvable = 0
    unresolvable = 0
    if not offtargets.empty:
        for rname in offtargets["rname"].astype(str):
            if species_index.symbol_for(rname):
                resolvable += 1
            else:
                unresolvable += 1

    measurements["hit_classes"] = {
        "as_run_no_transcript_index": dict(as_run.class_counts),
        "with_transcript_index": dict(enriched.class_counts),
        "with_transcript_index_and_repeat_scan": dict(with_repeats.class_counts),
        "offline_repeat_scan": {
            "distinct_guides_scanned": len(set(guides)),
            "reference_transcript_count": repeat_scan.reference_transcript_count,
            "threshold_fraction": repeat_scan.threshold_fraction,
            "guides_flagged_as_repeat": len(repeat_guides),
        },
        "as_run_symbol_lookup_missing": as_run.symbol_lookup_missing,
        "as_run_species_index_missing": as_run.species_index_missing,
        "with_index_symbol_lookup_missing": enriched.symbol_lookup_missing,
        "with_index_matched_symbol_hits": enriched.resolved_symbols,
        "symbol_resolution": {
            "hits_with_resolvable_symbol": resolvable,
            "hits_without_resolvable_symbol": unresolvable,
            "miss_rate": (unresolvable / (resolvable + unresolvable)) if (resolvable + unresolvable) else None,
            "reference_transcripts_indexed": species_index.transcript_count,
            "reference_transcripts_without_symbol": species_index.missing_symbol_count,
            "reference_symbol_miss_rate": (
                species_index.missing_symbol_count / species_index.transcript_count
                if species_index.transcript_count
                else None
            ),
        },
    }

    # verify the reconstruction against what the run itself wrote
    recon_off = candidates["id"].map(
        lambda cid: as_run.per_query.get(representative.get(cid, cid), {}).get("off_target", 0)
    )
    measurements["hit_classes"]["reconstruction_check"] = {
        "rows_compared": int(len(candidates)),
        "rows_matching_off_target_count": int((recon_off == candidates["off_target_count"].astype(int)).sum()),
    }

    # ---------------- B. off-target scope deltas ----------------
    nm_series = offtargets["nm"].astype(int) if not offtargets.empty else pd.Series(dtype=int)
    counted = [
        row
        for row in (offtargets.to_dict("records") if not offtargets.empty else [])
        if classify_hit(row, qname_to_guide.get(str(row["qname"]), str(row.get("qseq", ""))), empty_context).hit_class
        is HitClass.OFF_TARGET
    ]
    counted_nm = pd.Series([int(r["nm"]) for r in counted], dtype=int) if counted else pd.Series(dtype=int)
    off_counts = candidates["off_target_count"].astype(int)
    nm_le2 = candidates["id"].map(
        lambda cid: as_run.per_query.get(representative.get(cid, cid), {}).get("off_target_nm_le2", 0)
    )
    cap = OffTargetFilterCriteria().max_off_target_count

    measurements["offtarget_scope"] = {
        "alignments_total": int(len(offtargets)),
        "alignments_nm_distribution": {str(k): int(v) for k, v in sorted(nm_series.value_counts().items())}
        if not offtargets.empty
        else {},
        "counted_as_off_target": len(counted),
        "counted_nm_ge3": int((counted_nm >= 3).sum()) if len(counted_nm) else 0,
        "counted_nm_ge3_fraction": float((counted_nm >= 3).mean()) if len(counted_nm) else None,
        "counted_nm_le2": int((counted_nm <= 2).sum()) if len(counted_nm) else 0,
        "cap": cap,
        "excess_off_targets_current_scope": int((off_counts > cap).sum()) if cap is not None else None,
        "excess_off_targets_nm_le2_scope": int((nm_le2 > cap).sum()) if cap is not None else None,
        "off_target_count_max": int(off_counts.max()),
        "off_target_count_median": float(off_counts.median()),
        "off_target_count_nm_le2_max": int(nm_le2.max()) if len(nm_le2) else 0,
    }

    # ---------------- column redundancy ----------------
    measurements["column_redundancy"] = {
        "transcriptome_hits_total_equals_off_target_count": {
            "rows": int(len(candidates)),
            "equal_rows": int((candidates["transcriptome_hits_total"].astype(int) == off_counts).sum()),
        },
        "mirna_hits_high_risk_equals_mirna_hits_0mm_seed": {
            "rows": int(len(candidates)),
            "equal_rows": int(
                (candidates["mirna_hits_high_risk"].astype(int) == candidates["mirna_hits_0mm_seed"].astype(int)).sum()
            ),
            "nonzero_rows": int((candidates["mirna_hits_0mm_seed"].astype(int) > 0).sum()),
        },
    }

    # ---------------- D. per-gate independent failures ----------------
    mirna_by_query = mirna_counters(mirna_hits)
    gates = gate_independent_failures(
        candidates,
        as_run,
        mirna_by_query,
        representative,
        FilterCriteria(gc_min=args.gc_min, gc_max=args.gc_max),
        OffTargetFilterCriteria(),
    )
    ranked = sorted(
        ((n, v) for n, v in gates.items() if n != "__union__"),
        key=lambda item: -item[1]["independent_failures"],
    )
    measurements["gates"] = {
        "per_gate": gates,
        "ranked": [{"gate": name, **payload} for name, payload in ranked],
        "labelled_failures_in_csv": {
            str(k): int(v) for k, v in candidates["passes_filters"].astype(str).value_counts().items()
        },
    }

    # Same table under the classification a run with a transcript index and repeat detection
    # would have produced, so the share of gate failures that are missing-evidence artifacts
    # is measured rather than argued.
    repeat_flagged_ids = (
        candidates["guide_sequence"].astype(str).map(lambda seq: normalize_guide_sequence(seq) in repeat_guides)
    )
    candidates_repeat_aware = candidates.copy()
    candidates_repeat_aware["repeat_flagged"] = repeat_flagged_ids
    gates_enriched = gate_independent_failures(
        candidates_repeat_aware,
        with_repeats,
        mirna_by_query,
        representative,
        FilterCriteria(gc_min=args.gc_min, gc_max=args.gc_max),
        OffTargetFilterCriteria(),
    )
    measurements["gates"]["per_gate_with_index_and_repeats"] = gates_enriched
    measurements["gates"]["ranked_with_index_and_repeats"] = [
        {"gate": name, **payload}
        for name, payload in sorted(
            ((n, v) for n, v in gates_enriched.items() if n != "__union__"),
            key=lambda item: -item[1]["independent_failures"],
        )
    ]

    # what binds once off_target is rescoped, and again once asymmetry is switched off
    ordering = [name for name, _ in ranked]
    measurements["gates"]["binding_after_rescope"] = [
        name
        for name in ordering
        if name not in {"EXCESS_OFF_TARGETS", "EXCESS_OFF_TARGETS_nm_le2_scope"}
        and gates[name]["independent_failures"] > 0
    ]
    measurements["gates"]["binding_after_rescope_and_asymmetry_off"] = [
        name
        for name in ordering
        if name not in {"EXCESS_OFF_TARGETS", "EXCESS_OFF_TARGETS_nm_le2_scope", "LOW_ASYMMETRY"}
        and gates[name]["independent_failures"] > 0
    ]

    for name, path in (("combined_summary", combined_summary), ("combined_mirna_summary", mirna_summary)):
        if path.exists():
            measurements.setdefault("run_summaries", {})[name] = json.loads(path.read_text())
    if workflow_summary.exists():
        # Only the header sections: this file embeds every hit row and is tens of MB.
        full = json.loads(workflow_summary.read_text())
        measurements.setdefault("run_summaries", {})["workflow_summary"] = {
            key: full[key]
            for key in (
                "workflow_config",
                "transcript_summary",
                "transcript_annotation_summary",
                "orf_summary",
                "design_summary",
                "design_parameters",
                "repeat_summary",
            )
            if key in full
        }

    measurements["filter_status_values"] = [member.value for member in SiRNACandidate.FilterStatus]

    # Contribution-column shape, so a term with a tiny variance share can be told apart from a
    # term that is simply near-constant on this run.
    contribution_stats: dict[str, Any] = {}
    for term, column in {**POSTSCREEN_SIRNA_TERMS, **POSTSCREEN_MIRNA_EXTRA}.items():
        if column not in scored.columns:
            continue
        series = scored[column].astype(float).dropna()
        if series.empty:
            continue
        contribution_stats[term] = {
            "min": float(series.min()),
            "max": float(series.max()),
            "mean": float(series.mean()),
            "sd": float(series.std(ddof=1)),
            "distinct_values": int(series.nunique()),
        }
    measurements["variance_share"]["contribution_stats"] = contribution_stats
    measurements["artifact_totals"] = {
        "files": len(measurements["artifacts"]),
        "total_bytes": sum(item["bytes"] for item in measurements["artifacts"].values()),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(measurements, indent=2, default=str))
    print(
        json.dumps(
            {k: measurements[k] for k in ("shape", "offtarget_scope", "column_redundancy")}, indent=2, default=str
        )
    )
    print(json.dumps(measurements["variance_share"], indent=2, default=str))
    print(json.dumps(measurements["gates"]["ranked"], indent=2, default=str))
    print(json.dumps(measurements["hit_classes"], indent=2, default=str))


if __name__ == "__main__":
    main()
