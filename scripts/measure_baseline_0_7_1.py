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

Multi-species runs name one reference per screened species; ``--ortholog-symbol-alias`` adds an
extra accepted ortholog symbol as its own measured scope (mouse p53 is ``Trp53``, not ``Tp53``,
so symbol equality alone cannot recognise it):

    uv run python scripts/measure_baseline_0_7_1.py \
        --run-dir work/baseline_0.7.1/two_species/run \
        --reference human=work/baseline_0.7.1/reference/GRCh38_cdna_primary_1tx_per_gene_plus_TP53.fa \
        --reference mouse=~/.cache/sirnaforge/transcriptomes/Mus_musculus.GRCm39.cdna.all.fa \
        --ortholog-symbol-alias TRP53 \
        --output work/baseline_0.7.1/two_species/measurements_two_species.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

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
        return cast("pd.DataFrame", pd.DataFrame())
    return pd.read_csv(path, sep="\t")


HUMAN_SPECIES_ALIASES = frozenset({"human", "hsa", "homo_sapiens", ""})

EMPTY_QUERY_BUCKET: dict[str, int] = {
    "on_target": 0,
    "ortholog": 0,
    "repeat": 0,
    "off_target": 0,
    "off_target_nm_le2": 0,
    "off_target_nm_ge3": 0,
    # Species-scoped off-target counters. D3 proposes "human hits at nm <= 2", so both halves of
    # that scope need their own counter as well as the two mixed corners of the 2x2.
    "off_target_human_any_nm": 0,
    "off_target_human_nm_le2": 0,
    "off_target_nonhuman_any_nm": 0,
    "off_target_nonhuman_nm_le2": 0,
    "human_0mm": 0,
    "human_1mm": 0,
    "human_2mm": 0,
    "seed_0mm": 0,
    "human_total": 0,
}


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
    # Per-species views, keyed by the species label written into the hit table. Populated for
    # every run; on a single-species run each carries exactly one key.
    per_species_class_counts: dict[str, Counter[str]] = field(default_factory=dict)
    per_species_alignment_nm: dict[str, Counter[int]] = field(default_factory=dict)
    per_species_offtarget_nm: dict[str, Counter[int]] = field(default_factory=dict)
    per_species_shortfalls: dict[str, Counter[str]] = field(default_factory=dict)
    ortholog_symbols: Counter[str] = field(default_factory=Counter)


def counter_series(
    candidates: pd.DataFrame,
    per_query: Mapping[str, Mapping[str, int]],
    representative: Mapping[str, str],
    field_name: str,
) -> pd.Series:
    """Per-candidate-row value of one recomputed counter, resolved via the guide representative.

    A candidate row carries no hits of its own: the hit table keys on one representative id per
    distinct guide, so every row reads its guide's representative. Named rather than inlined as a
    lambda so the id type is stated once instead of being re-inferred at seven call sites.
    """

    def value(candidate_id: Any) -> int:
        key = str(candidate_id)
        bucket = per_query.get(representative.get(key, key), {})
        return int(bucket.get(field_name, 0))

    return candidates["id"].map(value)


def _accumulate_offtarget(bucket: dict[str, int], *, nm: int, seed_mismatches: int, treated_as_human: bool) -> None:
    """Add one genuine off-target hit to a per-candidate bucket, in every scope it belongs to."""
    bucket["off_target_nm_le2" if nm <= 2 else "off_target_nm_ge3"] += 1
    prefix = "off_target_human" if treated_as_human else "off_target_nonhuman"
    bucket[f"{prefix}_any_nm"] += 1
    if nm <= 2:
        bucket[f"{prefix}_nm_le2"] += 1
    if treated_as_human:
        bucket["human_total"] += 1
        if nm in (0, 1, 2):
            bucket[f"human_{nm}mm"] += 1
    if seed_mismatches == 0:
        bucket["seed_0mm"] += 1


def reconstruct(
    hits: pd.DataFrame,
    context: ClassificationContext,
    qname_to_guide: dict[str, str],
    label: str,
) -> Reconstruction:
    """Re-classify every hit and accumulate per-query and per-species counters."""
    per_query: dict[str, dict[str, int]] = {}
    class_counts: Counter[str] = Counter()
    per_species_class: dict[str, Counter[str]] = {}
    per_species_align_nm: dict[str, Counter[int]] = {}
    per_species_off_nm: dict[str, Counter[int]] = {}
    per_species_shortfall: dict[str, Counter[str]] = {}
    ortholog_symbols: Counter[str] = Counter()
    symbol_missing = 0
    index_missing = 0
    resolved = 0

    for row in hits.to_dict("records"):
        qname = str(row["qname"])
        guide = qname_to_guide.get(qname, str(row.get("qseq", "")))
        verdict = classify_hit(cast("Mapping[str, Any]", row), guide, context)
        class_counts[verdict.hit_class.value] += 1
        if verdict.matched_symbol:
            resolved += 1
        if verdict.symbol_lookup_missing:
            symbol_missing += 1
        if verdict.species_index_missing:
            index_missing += 1

        nm = int(row["nm"])
        seed_mm = int(row["seed_mismatches"])
        species = str(row.get("species") or "")
        species_key = species or "(blank)"
        treated_as_human = species.lower() in HUMAN_SPECIES_ALIASES

        per_species_class.setdefault(species_key, Counter())[verdict.hit_class.value] += 1
        per_species_align_nm.setdefault(species_key, Counter())[nm] += 1
        shortfalls = per_species_shortfall.setdefault(species_key, Counter())
        if verdict.symbol_lookup_missing:
            shortfalls["symbol_lookup_missing"] += 1
        if verdict.species_index_missing:
            shortfalls["species_index_missing"] += 1
        if verdict.hit_class is HitClass.ORTHOLOG and verdict.matched_symbol:
            ortholog_symbols[f"{species_key}:{verdict.matched_symbol}"] += 1

        bucket = per_query.setdefault(qname, dict(EMPTY_QUERY_BUCKET))
        bucket[verdict.hit_class.value] += 1
        if verdict.hit_class is not HitClass.OFF_TARGET:
            continue

        per_species_off_nm.setdefault(species_key, Counter())[nm] += 1
        _accumulate_offtarget(bucket, nm=nm, seed_mismatches=seed_mm, treated_as_human=treated_as_human)

    return Reconstruction(
        label=label,
        per_query=per_query,
        class_counts=class_counts,
        symbol_lookup_missing=symbol_missing,
        species_index_missing=index_missing,
        resolved_symbols=resolved,
        total_hits=int(len(hits)),
        per_species_class_counts=per_species_class,
        per_species_alignment_nm=per_species_align_nm,
        per_species_offtarget_nm=per_species_off_nm,
        per_species_shortfalls=per_species_shortfall,
        ortholog_symbols=ortholog_symbols,
    )


def species_split(recon: Reconstruction) -> dict[str, Any]:
    """Render one Reconstruction's per-species views, with the nm >= 3 share per species.

    The nm >= 3 share is the number the internal all-species run reported as 95.5% and the
    human-only canonical run measured at 6.53%; splitting it by species is what tells those
    two apart.
    """
    out: dict[str, Any] = {}
    for species in sorted(recon.per_species_class_counts):
        classes = recon.per_species_class_counts[species]
        align_nm = recon.per_species_alignment_nm.get(species, Counter())
        off_nm = recon.per_species_offtarget_nm.get(species, Counter())
        align_total = sum(align_nm.values())
        off_total = sum(off_nm.values())
        off_ge3 = sum(count for nm, count in off_nm.items() if nm >= 3)
        align_ge3 = sum(count for nm, count in align_nm.items() if nm >= 3)
        out[species] = {
            "alignments": align_total,
            "alignment_nm_distribution": {str(nm): align_nm[nm] for nm in sorted(align_nm)},
            "alignments_nm_ge3": align_ge3,
            "alignments_nm_ge3_fraction": (align_ge3 / align_total) if align_total else None,
            "hit_classes": dict(classes),
            "counted_off_target": off_total,
            "counted_off_target_nm_distribution": {str(nm): off_nm[nm] for nm in sorted(off_nm)},
            "counted_off_target_nm_ge3": off_ge3,
            "counted_off_target_nm_ge3_fraction": (off_ge3 / off_total) if off_total else None,
            "shortfalls": dict(recon.per_species_shortfalls.get(species, Counter())),
        }
    totals_off = sum(v["counted_off_target"] for v in out.values())
    totals_ge3 = sum(v["counted_off_target_nm_ge3"] for v in out.values())
    return {
        "per_species": out,
        "all_species_counted_off_target": totals_off,
        "all_species_counted_off_target_nm_ge3": totals_ge3,
        "all_species_counted_off_target_nm_ge3_fraction": (totals_ge3 / totals_off) if totals_off else None,
        "ortholog_matched_symbols": dict(recon.ortholog_symbols.most_common(20)),
    }


def excess_offtarget_scopes(
    candidates: pd.DataFrame,
    recon: Reconstruction,
    representative: dict[str, str],
    cap: int | None,
    guides: pd.Series,
) -> dict[str, Any]:
    """The four-way EXCESS_OFF_TARGETS table D3 turns on, at one cap.

    Scopes: (a) all species / any nm — today's default; (b) human only / any nm;
    (c) all species / nm <= 2; (d) human only / nm <= 2 — D3's proposal.

    Reported at row level and at distinct-guide level, because candidate rows are ~14.5x
    redundant per guide on this run and the row count therefore overstates how many
    molecules a scope actually rejects.
    """
    scopes = {
        "a_all_species_any_nm": "off_target",
        "b_human_only_any_nm": "off_target_human_any_nm",
        "c_all_species_nm_le2": "off_target_nm_le2",
        "d_human_only_nm_le2": "off_target_human_nm_le2",
    }
    table: dict[str, Any] = {
        "cap": cap,
        "rows": int(len(candidates)),
        "distinct_guides": int(guides.nunique()),
    }
    for scope, field_name in scopes.items():
        values = counter_series(candidates, recon.per_query, representative, field_name)
        failing = pd.Series(False, index=candidates.index) if cap is None else values > cap
        table[scope] = {
            "counter": field_name,
            "failures": None if cap is None else int(failing.sum()),
            "failures_distinct_guides": None if cap is None else int(guides[failing.to_numpy()].nunique()),
            "max": int(values.max()) if len(values) else 0,
            "median": float(values.median()) if len(values) else None,
            "mean": float(values.mean()) if len(values) else None,
        }
    return table


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

    def counter(field_name: str) -> pd.Series:
        return counter_series(candidates, recon.per_query, representative, field_name)

    def mirna_counter(field_name: str) -> pd.Series:
        return counter_series(candidates, mirna_by_query, representative, field_name)

    for name, gate_field, threshold in (
        ("TRANSCRIPTOME_PERFECT_MATCH", "human_0mm", offtarget_filters.max_transcriptome_hits_0mm),
        ("TRANSCRIPTOME_1MM", "human_1mm", offtarget_filters.max_transcriptome_hits_1mm),
        ("TRANSCRIPTOME_2MM", "human_2mm", offtarget_filters.max_transcriptome_hits_2mm),
        ("TRANSCRIPTOME_SEED_PERFECT", "seed_0mm", offtarget_filters.max_transcriptome_seed_perfect),
        ("EXCESS_OFF_TARGETS", "off_target", offtarget_filters.max_off_target_count),
    ):
        values = counter(gate_field)
        mask = pd.Series(False, index=candidates.index) if threshold is None else values > threshold
        record(name, mask, threshold, f"recomputed {gate_field}")

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


def cross_run_delta(hits: pd.DataFrame, other_run_dir: Path) -> dict[str, Any]:
    """Diff this run's hit table against another frozen run's, on both candidate id and guide.

    Two runs of identical inputs are expected to agree hit-for-hit, but the ``qname`` written
    into the hit table is one arbitrarily chosen candidate id per DISTINCT guide, and that
    choice is not stable across runs. So the id-keyed diff and the guide-keyed diff answer
    different questions, and only the guide-keyed one is a reproducibility statement.
    """
    other_path = other_run_dir / "off_target" / "results" / "aggregated" / "combined_offtargets.tsv"
    other = read_hits(other_path)
    site_key = ["species", "rname", "coord", "strand", "cigar", "nm", "seed_mismatches"]
    out: dict[str, Any] = {
        "other_run_dir": str(other_run_dir),
        "other_alignments": int(len(other)),
        "this_alignments": int(len(hits)),
    }
    if hits.empty or other.empty:
        return out
    # Restrict to species screened by BOTH runs: a species only one run screened is a scope
    # difference, not a disagreement, and mixing the two makes the delta uninterpretable.
    shared_species = sorted(set(hits["species"].astype(str)) & set(other["species"].astype(str)))
    out["species_compared"] = shared_species
    out["species_only_this_run"] = sorted(set(hits["species"].astype(str)) - set(other["species"].astype(str)))
    out["species_only_other_run"] = sorted(set(other["species"].astype(str)) - set(hits["species"].astype(str)))
    hits = hits[hits["species"].astype(str).isin(shared_species)]
    other = other[other["species"].astype(str).isin(shared_species)]
    out["this_alignments_in_shared_species"] = int(len(hits))
    out["other_alignments_in_shared_species"] = int(len(other))
    for label, key in (("guide_keyed", ["qseq", *site_key]), ("candidate_id_keyed", ["qname", *site_key])):
        left = set(hits[key].astype(str).agg("|".join, axis=1))
        right = set(other[key].astype(str).agg("|".join, axis=1))
        out[label] = {
            "join_key": key,
            "shared": len(left & right),
            "only_this_run": len(left - right),
            "only_other_run": len(right - left),
        }
    this_map = dict(zip(hits["qseq"].astype(str), hits["qname"].astype(str), strict=True))
    other_map = dict(zip(other["qseq"].astype(str), other["qname"].astype(str), strict=True))
    shared_guides = set(this_map) & set(other_map)
    changed = {g for g in shared_guides if this_map[g] != other_map[g]}
    out["representative_qname_stability"] = {
        "shared_guides": len(shared_guides),
        "guides_whose_representative_qname_changed": len(changed),
        "fraction_changed": (len(changed) / len(shared_guides)) if shared_guides else None,
    }
    return out


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
        is_human = species in HUMAN_SPECIES_ALIASES
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
    parser.add_argument(
        "--reference",
        action="append",
        required=True,
        metavar="[SPECIES=]FASTA",
        help=(
            "Screened reference FASTA. Repeat once per screened species as 'species=path'; a bare "
            "path is read as 'human=path' so single-species invocations are unchanged."
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gc-min", type=float, default=30.0, help="as-run gc_min (workflow CLI default)")
    parser.add_argument("--gc-max", type=float, default=60.0, help="as-run gc_max (workflow CLI default)")
    parser.add_argument("--query-gene-id", default="ENSG00000141510")
    parser.add_argument("--query-gene-symbol", default="TP53")
    parser.add_argument(
        "--ortholog-symbol-alias",
        action="append",
        default=[],
        metavar="SYMBOL",
        help=(
            "Extra gene symbol accepted as the query gene's ortholog name, e.g. TRP53 for mouse "
            "p53. Measured as an ADDITIONAL scope; the default scopes never use it."
        ),
    )
    parser.add_argument(
        "--compare-run-dir",
        type=Path,
        default=None,
        help=(
            "Another frozen run to diff this one's hit table against, on both candidate `id` and "
            "guide sequence. Phase 4 needs the delta set; the two join keys do not agree."
        ),
    )
    args = parser.parse_args()

    references: dict[str, Path] = {}
    for entry in args.reference:
        text = str(entry)
        species_name, _, path_text = text.partition("=")
        if not path_text:
            species_name, path_text = "human", text
        references[species_name.strip().lower()] = Path(path_text)
    # The query species' reference stays addressable under the historical attribute name, so the
    # single-species code paths below are untouched.
    query_reference = references.get("human") or next(iter(references.values()))

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
    species_index = enriched_index.build("human", query_reference)
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
        set(guides), query_reference
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
    recon_off = counter_series(candidates, as_run.per_query, representative, "off_target")
    measurements["hit_classes"]["reconstruction_check"] = {
        "rows_compared": int(len(candidates)),
        "rows_matching_off_target_count": int((recon_off == candidates["off_target_count"].astype(int)).sum()),
    }

    # ---------------- F. species dimension (plan decision D3) ----------------
    # Every counter split by the species label the aligner wrote, plus the four-way
    # EXCESS_OFF_TARGETS table under three classification scopes: as run (no index for any
    # species), with an index for every screened species, and with the ortholog symbol alias.
    cap_default = OffTargetFilterCriteria().max_off_target_count

    multi_index = TranscriptGeneIndex()
    reference_identity: dict[str, Any] = {}
    for species_name, fasta in sorted(references.items()):
        built = multi_index.build(species_name, fasta)
        reference_identity[species_name] = {
            "fasta": str(fasta),
            "exists": fasta.exists(),
            "transcripts_indexed": built.transcript_count,
            "transcripts_without_symbol": built.missing_symbol_count,
            "query_symbol_transcripts": sorted(built.transcripts_for_symbol(args.query_gene_symbol)),
            "alias_symbol_transcripts": {
                alias.upper(): sorted(built.transcripts_for_symbol(alias)) for alias in args.ortholog_symbol_alias
            },
        }

    multi_context = ClassificationContext(
        query_gene_ids=frozenset({args.query_gene_id.upper()}),
        query_gene_symbols=query_gene_symbols,
        on_target_transcript_ids=on_target_transcript_ids,
        query_species="human",
        index=multi_index,
        repeat_flagged_guides=frozenset(),
        requested_species=frozenset(references),
    )
    multi = reconstruct(offtargets, multi_context, qname_to_guide, "all_species_indexed")

    alias_symbols = frozenset(query_gene_symbols | {alias.upper() for alias in args.ortholog_symbol_alias})
    alias_context = ClassificationContext(
        query_gene_ids=frozenset({args.query_gene_id.upper()}),
        query_gene_symbols=alias_symbols,
        on_target_transcript_ids=on_target_transcript_ids,
        query_species="human",
        index=multi_index,
        repeat_flagged_guides=frozenset(),
        requested_species=frozenset(references),
    )
    alias = reconstruct(offtargets, alias_context, qname_to_guide, "all_species_indexed_with_alias")

    scope_labels = {
        "as_run_no_transcript_index": as_run,
        "all_species_indexed": multi,
        "all_species_indexed_with_alias": alias,
    }
    measurements["species_dimension"] = {
        "references": reference_identity,
        "query_gene_symbols_default": sorted(query_gene_symbols),
        "query_gene_symbols_with_alias": sorted(alias_symbols),
        "species_labels_in_hit_table": sorted(
            {str(s or "(blank)") for s in (offtargets["species"] if not offtargets.empty else [])}
        ),
        "split": {label: species_split(recon) for label, recon in scope_labels.items()},
        "excess_off_targets_four_scopes": {
            label: excess_offtarget_scopes(candidates, recon, representative, cap_default, guides)
            for label, recon in scope_labels.items()
        },
        # Row-level vs hit-level species share. A guide's hits are attributed to every candidate
        # row carrying that guide, so the row-weighted split is not the hit-level split.
        "species_share_of_off_target_count": {
            "hit_level": {
                species: view["counted_off_target"] for species, view in species_split(as_run)["per_species"].items()
            },
            "row_weighted": {
                "note": "sum over candidate rows of that row's per-species counted off-targets",
                **{
                    species: int(counter_series(candidates, as_run.per_query, representative, field_name).sum())
                    for species, field_name in (
                        ("human", "off_target_human_any_nm"),
                        ("non_human", "off_target_nonhuman_any_nm"),
                    )
                },
            },
        },
    }

    # conservation_score, computable for the first time now that a non-query species was screened.
    if "conservation_score" in candidates.columns:
        cons = pd.to_numeric(candidates["conservation_score"], errors="coerce")
        measurements["species_dimension"]["conservation_score"] = {
            "rows": int(len(cons)),
            "null_rows": int(cons.isna().sum()),
            "non_null_rows": int(cons.notna().sum()),
            "distinct_values": {str(k): int(v) for k, v in sorted(cons.dropna().value_counts().items())},
            "min": float(cons.min()) if cons.notna().any() else None,
            "max": float(cons.max()) if cons.notna().any() else None,
            "mean": float(cons.mean()) if cons.notna().any() else None,
        }
    if "ortholog_hits" in candidates.columns:
        orth = pd.to_numeric(candidates["ortholog_hits"], errors="coerce")
        measurements["species_dimension"]["ortholog_hits_as_written_by_the_run"] = {
            "rows": int(len(orth)),
            "nonzero_rows": int((orth > 0).sum()),
            "max": int(orth.max()) if orth.notna().any() else None,
        }
    if "ortholog_species" in candidates.columns:
        measurements["species_dimension"]["ortholog_species_as_written_by_the_run"] = {
            str(k): int(v) for k, v in candidates["ortholog_species"].fillna("(empty)").value_counts().items()
        }

    # ---------------- B. off-target scope deltas ----------------
    nm_series = offtargets["nm"].astype(int) if not offtargets.empty else pd.Series(dtype=int)
    counted = [
        row
        for row in (offtargets.to_dict("records") if not offtargets.empty else [])
        if classify_hit(
            cast("Mapping[str, Any]", row),
            qname_to_guide.get(str(row["qname"]), str(row.get("qseq", ""))),
            empty_context,
        ).hit_class
        is HitClass.OFF_TARGET
    ]
    counted_nm = pd.Series([int(r["nm"]) for r in counted], dtype=int) if counted else pd.Series(dtype=int)
    off_counts = candidates["off_target_count"].astype(int)
    nm_le2 = counter_series(candidates, as_run.per_query, representative, "off_target_nm_le2")
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
        # Issue #100: does any published artifact distinguish a species that was requested and
        # completed from one requested and unavailable? Collect every field that could carry it.
        offtarget_summary = full.get("offtarget_summary") or {}
        filtering_stats = offtarget_summary.get("filtering_stats") or {}
        measurements.setdefault("run_summaries", {})["species_requiredness_evidence"] = {
            "workflow_config_genome_species": (full.get("workflow_config") or {}).get("genome_species"),
            "workflow_config_species_explicitly_requested": (full.get("workflow_config") or {}).get(
                "species_explicitly_requested"
            ),
            "offtarget_status": offtarget_summary.get("status"),
            "offtarget_warnings": offtarget_summary.get("warnings"),
            "filtering_stats_requested_species": filtering_stats.get("requested_species"),
            "filtering_stats_screened_species": filtering_stats.get("screened_species"),
            "filtering_stats_unscreened_species": filtering_stats.get("unscreened_species"),
            "filtering_stats_species_index_misses": filtering_stats.get("species_index_misses"),
            "filtering_stats_hit_classes": filtering_stats.get("hit_classes"),
            "filtering_stats_per_species": filtering_stats.get("per_species"),
        }

    if args.compare_run_dir is not None:
        measurements["cross_run_hit_table_delta"] = cross_run_delta(offtargets, args.compare_run_dir)

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
    print(json.dumps(measurements["species_dimension"], indent=2, default=str))
    print(json.dumps(measurements.get("run_summaries", {}).get("species_requiredness_evidence", {}), indent=2))


if __name__ == "__main__":
    main()
