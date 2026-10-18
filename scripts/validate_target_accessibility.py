#!/usr/bin/env python
r"""Validate the target_accessibility scoring term against measured siRNA knockdown.

The calibration record for issue #95. The tracked regression guard is
`tests/unit/test_target_accessibility.py`, which runs the same three checks against a small
vendored subset; this script runs them against the full third-party benchmark, whose data is
NOT tracked and must be supplied by path.

Three things are checked, in the order they matter:

1. **Scored metric.** Spearman rho of P(seed-end 8-mer of the target site is unpaired) against
   measured inhibition. Expected >= +0.24 (full benchmark, n=2,779, W=150/L=100: +0.267).
2. **Geometry control.** The same statistic taken at the *non-seed* end of the target site is a
   measured near-null (+0.067). If the control ever rises to match the scored metric, the target
   site is being indexed backwards -- the specific bug this guards.
3. **Published-construct sanity.** The widely used anti-TP53 pSUPER site must not fall in the
   bottom decile of the scored term across the transcript's windows.

Geometry, since it is the thing that is easy to get wrong: target site T[1..L] 5'->3' and guide
G[1..L] 5'->3' are antiparallel, so G[i] pairs T[L+1-i]. The guide seed G[2..8] therefore pairs
the *3' end* of the target site.

Data (see tests/unit/data/README.md for the full provenance record):
  --benchmark-csv   columns: efficacy (inhibition, higher = more knockdown), siRNA_seq (guide),
                    B (accession)
  --transcripts     FASTA of the accessions' full-length sequences
  --tp53-fasta      TP53-201 (ENST00000269305.9) cDNA

Measured efficacy is from:

    Huesken D, Lange J, Mickanin C, Weiler J, Asselbergs F, Warner J, Meloon B, Engel S,
    Rosenberg A, Cohen D, Labow M, Reinhardt M, Natt F, Hall J. "Design of a genome-wide siRNA
    library using an artificial neural network." Nat Biotechnol. 2005 Aug;23(8):995-1001.
    doi:10.1038/nbt1118. PMID: 16025102.

obtained via the third-party redistribution https://github.com/apkrfi/unMod-siRNA-Pred
(`raw_data/All_Dataset.csv`), which also contributes 385 rows labelled `Other` whose provenance
could not be established. Those rows are kept here -- this script is a calibration record, not a
shipped artifact -- but are excluded from the tracked test fixture, which cites only Huesken.
The primary paper is not open access, so these values were NOT verified against it.

Transcript sequences are NCBI Nucleotide (`nuccore`) records fetched with `efetch`, not part of the
Huesken paper. TP53-201 is from Ensembl.

Usage:
    uv run python scripts/validate_target_accessibility.py \\
        --benchmark-csv work/sirna_bench.csv \\
        --transcripts work/bench_tx.fa \\
        --tp53-fasta work/tp53.fa
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from sirnaforge.core.scoring import target_accessibility_sub_score
from sirnaforge.core.thermodynamics import SEED_END_WINDOW_NT, TargetAccessibilityProfile
from sirnaforge.models.sirna import (
    DEFAULT_ACCESSIBILITY_LOG_FLOOR,
    DEFAULT_PLFOLD_MAX_BP_SPAN,
    DEFAULT_PLFOLD_WINDOW,
)

# Acceptance thresholds from issue #95.
MIN_SCORED_RHO = 0.24
MAX_CONTROL_RHO = 0.10
MIN_TP53_PERCENTILE = 10.0

# The pSUPER/pSuper-Retro anti-p53 construct: 78 Europe PMC full-text hits, the field's de facto
# positive control for p53 knockdown. Given as the target site (sense), as it occurs in TP53-201.
TP53_PSUPER_SITE = "GACTCCAGTGGTAATCTAC"


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA into {accession without version: sequence}."""
    sequences: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name = line[1:].split()[0].split(".")[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def spearman(a: pd.Series, b: pd.Series) -> float:
    """Spearman rho as Pearson on ranks. Avoids adding a scipy dependency for one number."""
    return float(np.corrcoef(a.rank().to_numpy(), b.rank().to_numpy())[0, 1])


def site_metrics(profile: TargetAccessibilityProfile, start: int, length: int) -> tuple[float | None, float | None]:
    """Return (seed-end 8-mer opening, non-seed 8-mer opening) for a site.

    The seed-end window ends at the site's 3'-most base; the control window ends 8 nt in from the
    site's 5' end, i.e. covers the other end of the same site. Same length, same transcript, same
    fold -- the two differ only in which end of the site they cover.
    """
    scored = profile.site_accessibility(start, length).seed_end_8mer
    control = profile.site_accessibility(start, SEED_END_WINDOW_NT).seed_end_8mer
    return scored, control


def evaluate_benchmark(
    benchmark_csv: Path, transcripts_fasta: Path, window_size: int, max_bp_span: int
) -> pd.DataFrame:
    """Score every locatable benchmark siRNA against its transcript."""
    transcripts = read_fasta(transcripts_fasta)
    table = pd.read_csv(benchmark_csv)

    rows: list[dict[str, float | str]] = []
    for accession, group in table.groupby("B"):
        sequence = transcripts.get(str(accession))
        if not sequence:
            continue
        sequence = sequence.upper()
        guide_length = int(group.siRNA_seq.str.len().max())
        profile = TargetAccessibilityProfile.fold(
            sequence, window_size=window_size, max_bp_span=max_bp_span, u_max=guide_length
        )
        for _, row in group.iterrows():
            site = str(Seq(str(row.siRNA_seq).upper()).reverse_complement())
            start = sequence.find(site)
            if start < 0:
                continue
            scored, control = site_metrics(profile, start, len(site))
            if scored is None or control is None:
                continue
            rows.append(
                {
                    "accession": str(accession),
                    "efficacy": float(row.efficacy),
                    "p_seed_end": scored,
                    "p_non_seed_end": control,
                }
            )
    scored_rows: pd.DataFrame = pd.DataFrame(rows)
    return scored_rows


def tp53_percentile(tp53_fasta: Path, window_size: int, max_bp_span: int) -> tuple[float, int, int]:
    """Percentile of the pSUPER site's scored term among all same-length windows of TP53-201."""
    sequences = read_fasta(tp53_fasta)
    sequence = next(iter(sequences.values())).upper()
    start = sequence.find(TP53_PSUPER_SITE)
    if start < 0:
        raise SystemExit(f"{TP53_PSUPER_SITE} not found in {tp53_fasta}")

    length = len(TP53_PSUPER_SITE)
    profile = TargetAccessibilityProfile.fold(sequence, window_size=window_size, max_bp_span=max_bp_span, u_max=length)
    values = [
        p
        for i in range(len(sequence) - length + 1)
        if (p := profile.site_accessibility(i, length).seed_end_8mer) is not None
    ]
    site_value = profile.site_accessibility(start, length).seed_end_8mer
    assert site_value is not None
    percentile = 100.0 * sum(1 for v in values if v <= site_value) / len(values)
    return percentile, start, len(values)


def main() -> int:
    """Run the three acceptance checks and return a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark-csv", type=Path, required=True, help="Measured-efficacy CSV")
    parser.add_argument("--transcripts", type=Path, required=True, help="FASTA of the benchmark accessions")
    parser.add_argument("--tp53-fasta", type=Path, required=True, help="TP53-201 cDNA FASTA")
    parser.add_argument("--window", type=int, default=DEFAULT_PLFOLD_WINDOW, help="RNAplfold window W")
    parser.add_argument("--max-bp-span", type=int, default=DEFAULT_PLFOLD_MAX_BP_SPAN, help="RNAplfold span L")
    args = parser.parse_args()

    print(f"RNAplfold W={args.window} L={args.max_bp_span}, log floor {DEFAULT_ACCESSIBILITY_LOG_FLOOR}")

    measured = evaluate_benchmark(args.benchmark_csv, args.transcripts, args.window, args.max_bp_span)
    if measured.empty:
        raise SystemExit("no benchmark siRNAs could be located in the supplied transcripts")

    scored_rho = spearman(measured.efficacy, measured.p_seed_end)
    control_rho = spearman(measured.efficacy, measured.p_non_seed_end)
    features = measured.p_seed_end.map(target_accessibility_sub_score)
    percentile, start, windows = tp53_percentile(args.tp53_fasta, args.window, args.max_bp_span)

    print(
        f"located {len(measured)} of {len(pd.read_csv(args.benchmark_csv))} siRNAs "
        f"across {measured.accession.nunique()} transcripts\n"
    )

    checks = [
        (
            f"scored term (seed-end {SEED_END_WINDOW_NT}-mer) rho",
            scored_rho,
            scored_rho >= MIN_SCORED_RHO,
            f">= +{MIN_SCORED_RHO}",
        ),
        (
            f"geometry control (non-seed {SEED_END_WINDOW_NT}-mer) rho",
            control_rho,
            control_rho <= MAX_CONTROL_RHO,
            f"<= +{MAX_CONTROL_RHO}",
        ),
        (
            f"TP53 pSUPER site percentile (offset {start}, {windows} windows)",
            percentile,
            percentile > MIN_TP53_PERCENTILE,
            f"> {MIN_TP53_PERCENTILE} (not bottom decile)",
        ),
    ]
    for label, value, ok, expectation in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {label:56s} {value:+8.3f}  (want {expectation})")

    print(f"\nnormalised feature: min {features.min():.3f} median {features.median():.3f} max {features.max():.3f}")
    return 0 if all(ok for _, _, ok, _ in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
