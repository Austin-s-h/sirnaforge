#!/usr/bin/env python3
"""Run a fixed exploratory siRNAforge feature check on the DSIR 2012 qRT-PCR panel.

The panel contains measured residual target-mRNA expression for guide sequences selected by DSIR.
This analysis checks only the direction of the frozen siRNAforge asymmetry and GC features within
each target, then bootstraps targets as clusters. It does not validate the full default profile or
unintended-repression liability, and it never fits or promotes weights.

Run from the repository root with the project environment:

    .venv/bin/python scripts/validate_dsir2012_independent.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import RNA

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.core.scoring import SCORING_WEIGHT_SET_VERSION
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.scoring_profile import SHIPPED_PROFILE
from sirnaforge.models.sirna import DesignParameters, ScoringWeights, SiRNACandidate

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "tests/data/benchmarks/ml_reviews/independent_panels/dsir_2012.jsonl"
SOURCES = ROOT / "tests/data/benchmarks/ml_reviews/source_extracts"
OLIGOGYM = ROOT / "tests/data/benchmarks/oligogym/records.csv"
NATIVE_OBSERVATIONS = ROOT / "tests/data/benchmarks/oligogym_native_design/observations.csv"
FEATURE_FILES = (
    Path(__file__).resolve(),
    ROOT / "src/sirnaforge/core/design.py",
    ROOT / "src/sirnaforge/core/thermodynamics.py",
    ROOT / "src/sirnaforge/models/scoring_profile.py",
)
BOOTSTRAP_SEED = 102
BOOTSTRAP_DRAWS = 20_000
_RC = str.maketrans("ACGT", "TGCA")


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one input file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical(sequence: str) -> str:
    """Normalize a sequence to DNA alphabet for identity comparison."""
    return "".join(base for base in sequence.upper().replace("U", "T") if base in "ACGT")


def ranks(values: np.ndarray) -> np.ndarray:
    """Return one-based average ranks with stable handling of ties."""
    order = np.argsort(values, kind="mergesort")
    result = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        result[order[start:stop]] = (start + stop - 1) / 2 + 1
        start = stop
    return result


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    """Compute Spearman correlation without adding a SciPy runtime dependency."""
    x = ranks(left)
    y = ranks(right)
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def read_panel() -> list[dict[str, Any]]:
    """Load and validate the normalized S1/S3 table join."""
    rows = [json.loads(line) for line in PANEL.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 129:
        raise ValueError(f"Expected 129 matched table rows, found {len(rows)}")
    identifiers: set[str] = set()
    guides: set[str] = set()
    for row in rows:
        if row["collection_id"] != "dsir_2012":
            raise ValueError(f"Unexpected collection for {row['source_record_id']}")
        if row["endpoint"] != "residual_expression" or row["direction"] != "lower_is_better":
            raise ValueError(f"Unexpected endpoint semantics for {row['source_record_id']}")
        guide = str(row["guide_sequence"]).upper().replace("T", "U")
        passenger = str(row["passenger_sequence"]).upper().replace("T", "U")
        if len(guide) != 21 or len(passenger) != 21 or set(guide + passenger) - set("ACGU"):
            raise ValueError(f"Expected complete 21-nt strands for {row['source_record_id']}")
        guide_core = canonical(guide)[:19]
        passenger_core = canonical(passenger)[:19]
        if guide_core.translate(_RC)[::-1] != passenger_core or row["paired_length"] != 19:
            raise ValueError(f"Expected a complementary 19-bp core for {row['source_record_id']}")
        if row["source_record_id"] in identifiers or canonical(guide) in guides:
            raise ValueError(f"Duplicate source record or guide: {row['source_record_id']}")
        identifiers.add(row["source_record_id"])
        guides.add(canonical(guide))
    return rows


def local_sequence_sets() -> tuple[dict[str, set[str]], dict[str, set[str]], list[Path]]:
    """Collect local 19-nt core and exact 21-nt sequence identities for overlap audit."""
    cores: dict[str, set[str]] = defaultdict(set)
    full: dict[str, set[str]] = defaultdict(set)
    source_files = source_extract_files()
    input_files = [*source_files, OLIGOGYM, NATIVE_OBSERVATIONS]

    def add(source: str, sequence: str) -> None:
        value = canonical(sequence)
        if len(value) < 19:
            return
        core = value[:19]
        cores[source].update((core, core.translate(_RC)[::-1]))
        if len(value) >= 21:
            exact = value[:21]
            full[source].update((exact, exact.translate(_RC)[::-1]))

    for path in source_files:
        source = path.stem
        source_core_count = len(cores[source])
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            for key in ("guide_sequence", "sequence", "paired_guide_sequence", "full_guide_sequence"):
                if row.get(key):
                    add(source, str(row[key]))
        if len(cores[source]) == source_core_count:
            raise ValueError(f"No usable guide sequences found in overlap-audit input: {path}")

    for source, path in (("oligogym", OLIGOGYM), ("oligogym_native", NATIVE_OBSERVATIONS)):
        if not path.exists():
            raise FileNotFoundError(f"Required overlap-audit input is missing: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                for key in ("guide_sequence", "paired_guide_sequence", "full_guide_sequence"):
                    if row.get(key):
                        add(source, row[key])
        if not cores[source]:
            raise ValueError(f"No usable guide sequences found in overlap-audit input: {path}")
    return cores, full, input_files


def source_extract_files() -> list[Path]:
    """Return the tracked source extracts, rejecting a missing audit corpus."""
    files = sorted(SOURCES.glob("*.jsonl"))
    if not files:
        raise ValueError(f"No source-extract JSONL files found in overlap-audit input directory: {SOURCES}")
    return files


def target_cluster_interval(values: np.ndarray) -> tuple[float, float, float]:
    """Return target mean and percentile interval from a fixed cluster bootstrap."""
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(0, len(values), size=(BOOTSTRAP_DRAWS, len(values)))
    means = values[draws].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(values.mean()), float(low), float(high)


def analyze() -> dict[str, Any]:
    """Compute target-stratified sequence-feature associations and audit provenance."""
    rows = read_panel()
    calculator = ThermodynamicCalculator()
    designer = SiRNADesigner(DesignParameters())
    features: dict[str, list[tuple[float, float, float]]] = defaultdict(list)

    for row in rows:
        guide = str(row["guide_sequence"]).upper().replace("T", "U")
        passenger = str(row["passenger_sequence"]).upper().replace("T", "U")
        gc_percent = 100.0 * (guide.count("G") + guide.count("C")) / len(guide)
        # The source duplexes have a 19-bp core plus two 3' overhang bases on
        # each strand. Score asymmetry on the paired core so terminal windows
        # do not include unpaired overhang bases. This is a core-only proxy,
        # not a direct validation of candidate geometry for a shipped 21-nt
        # fully complementary duplex.
        guide_core = guide[:19]
        passenger_core = passenger[:19]
        candidate = SiRNACandidate(
            id=str(row["source_record_id"]),
            transcript_id=str(row["target_gene"]),
            position=1,
            guide_sequence=guide_core,
            passenger_sequence=passenger_core,
            gc_content=100.0 * (guide_core.count("G") + guide_core.count("C")) / len(guide_core),
            length=len(guide_core),
            asymmetry_score=0.0,
        )
        asymmetry = calculator.calculate_asymmetry_score(candidate)[2]
        # The shipped GC term is computed from the full guide strand.
        gc_score = designer._calculate_gc_score(gc_percent)
        residual = float(row["residual_expression_percent"])
        features[str(row["target_gene"])].append((asymmetry, gc_score, residual))

    targets = sorted(features)
    per_target: dict[str, dict[str, float]] = {
        "asymmetry_19nt_core_proxy": {},
        "gc_content_full_guide": {},
    }
    for target in targets:
        values = np.asarray(features[target], dtype=float)
        per_target["asymmetry_19nt_core_proxy"][target] = spearman(values[:, 0], values[:, 2])
        per_target["gc_content_full_guide"][target] = spearman(values[:, 1], values[:, 2])

    cores, exact, overlap_input_files = local_sequence_sets()
    overlap_rows: dict[str, list[str]] = {"core_19": [], "exact_21": []}
    core_union = set().union(*cores.values()) if cores else set()
    exact_union = set().union(*exact.values()) if exact else set()
    for row in rows:
        guide = canonical(str(row["guide_sequence"]))
        if guide[:19] in core_union:
            overlap_rows["core_19"].append(str(row["source_record_id"]))
        if guide[:21] in exact_union:
            overlap_rows["exact_21"].append(str(row["source_record_id"]))

    summaries: dict[str, dict[str, float | int]] = {}
    for term, target_rhos in per_target.items():
        rho_values = np.asarray([target_rhos[target] for target in targets], dtype=float)
        if not np.isfinite(rho_values).all():
            raise ValueError(f"Undefined target-level Spearman value for {term}")
        mean, low, high = target_cluster_interval(rho_values)
        summaries[term] = {
            "equal_target_mean_spearman": mean,
            "target_bootstrap_95pct_low": low,
            "target_bootstrap_95pct_high": high,
            "targets_with_negative_rho": int(np.sum(rho_values < 0)),
            "targets": len(targets),
        }

    weights = ScoringWeights()
    design_weights = weights.design.as_mapping()
    profile_payload = {
        "profile": SHIPPED_PROFILE.model_dump(mode="json"),
        "term_records": [term.model_dump(mode="json") for term in SHIPPED_PROFILE.term_records()],
    }
    profile_json = json.dumps(profile_payload, sort_keys=True, separators=(",", ":"))
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()
    return {
        "analysis_id": "dsir-2012-target-stratified-feature-check-v2",
        "purpose": "exploratory efficacy-term check; no fitting, full-profile claim, or default promotion",
        "source": {
            "doi": "10.1371/journal.pone.0048057",
            "supplements": ["Table S1 (s003)", "Table S3 (s005)"],
            "endpoint": "residual non-cleaved target mRNA percent; lower is greater knockdown",
            "selection": "DSIR-designed and high-score-selected guides; second-round candidates were adaptively selected using first-round results",
            "assay_conditions": "20 nM siRNA in HeLa cells; 72 h endpoint except BCL2L1 at 24 h",
        },
        "inputs": {
            "panel_jsonl_sha256": sha256(PANEL),
            "feature_code_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in FEATURE_FILES},
            "git_revision": revision,
            "profile_id": SHIPPED_PROFILE.profile_id,
            "canonical_profile_terms_sha256": hashlib.sha256(profile_json.encode()).hexdigest(),
            "weight_set_version": SCORING_WEIGHT_SET_VERSION,
            "design_vector_weights": {
                "asymmetry": float(design_weights["asymmetry"]),
                "gc_content": float(design_weights["gc_content"]),
            },
            "viennarna_version": getattr(RNA, "__version__", "unknown"),
            "sequence_overlap_audit_input_sha256": {
                str(path.relative_to(ROOT)): sha256(path) for path in overlap_input_files
            },
        },
        "counts": {
            "measured_table_rows": len(rows),
            "unique_guides": len({canonical(str(row["guide_sequence"])) for row in rows}),
            "target_groups": len(targets),
            "target_group_sizes": {target: len(features[target]) for target in targets},
            "local_sequence_overlap_rows": overlap_rows,
            "local_sequence_overlap_scope": "tracked source extracts and OligoGym sequence files only; not an exhaustive audit of every external training corpus",
        },
        "per_target_spearman_vs_residual_expression": per_target,
        "cluster_summaries": summaries,
        "limitations": [
            "one study and HeLa qRT-PCR assays at 20 nM; 72 h except BCL2L1 at 24 h; target groups, not individual guides, are the bootstrap units",
            "the source used two selection rounds, with the second round informed by first-round findings",
            "asymmetry is evaluated on the verified 19-nt paired core as an applicability-limited proxy; the source has two-base 3-prime overhangs unlike the shipped fully complementary 21-nt duplex geometry",
            "guide selection by DSIR restricts the score range and is not a prospective sample of siRNA designs",
            "the article reports a 128-sequence panel in prose while the S1/S3 ID join contains 129 unique guide rows",
            "the paper analyzed several same-panel features; this post hoc source recovery is exploratory",
            "the data do not provide a broad measured unintended-repression endpoint; three CSNK2A1 guides were checked against CSNK2A2 only",
            "no score weights or thresholds were fit, selected, or promoted",
        ],
    }


def main() -> None:
    """Print the reproducible analysis record as JSON."""
    print(json.dumps(analyze(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
