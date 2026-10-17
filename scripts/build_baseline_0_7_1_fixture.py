#!/usr/bin/env python3
"""Cut a small, deterministic fixture out of the frozen 0.7.1 canonical run.

The canonical run is 131 MB and gitignored, so downstream work packages cannot test against
it directly. This carves out a slice small enough to track: every candidate row and every
off-target/miRNA hit row for a deterministically chosen set of guide sequences, keeping the
column layout byte-identical to the full artifacts.

Guides are selected to exercise the properties the 0.7.1 issues argue about, then sorted so
the output is reproducible:

* guides above ``max_off_target_count`` (the EXCESS_OFF_TARGETS population)
* guides with perfect miRNA seed hits (MIRNA_PERFECT_SEED == HIGH_RISK_MIRNA)
* guides on multiple transcripts (the guide-vs-id join key)
* guides that pass every gate

Usage (from the repo root):

    uv run python scripts/build_baseline_0_7_1_fixture.py \
        --run-dir work/baseline_0.7.1/run \
        --out-dir tests/unit/data/baseline_0_7_1
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

PER_BUCKET = 8


def main() -> None:
    """Write the tracked fixture slice."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    aggregated = args.run_dir / "off_target" / "results" / "aggregated"
    candidates = pd.read_csv(args.run_dir / "sirnaforge" / "candidates_all.csv")
    offtargets = pd.read_csv(aggregated / "combined_offtargets.tsv", sep="\t")
    mirna = pd.read_csv(aggregated / "combined_mirna_hits.tsv", sep="\t")

    rows_per_guide = candidates["guide_sequence"].value_counts()

    def pick(frame: pd.DataFrame) -> list[str]:
        return sorted(frame["guide_sequence"].dropna().unique())[:PER_BUCKET]

    selected: list[str] = []
    selected += pick(candidates[candidates["off_target_count"] > 15])
    selected += pick(candidates[candidates["mirna_hits_0mm_seed"] > 0])
    selected += pick(candidates[candidates["guide_sequence"].map(rows_per_guide) > 5])
    selected += pick(candidates[candidates["passes_filters"].astype(str) == "PASS"])
    guides = sorted(dict.fromkeys(selected))

    slice_candidates = candidates[candidates["guide_sequence"].isin(guides)].sort_values(["guide_sequence", "id"])
    # Hit rows are keyed on the representative candidate id chosen by deduplication, so keep
    # every id present in the slice rather than trying to re-derive which one was submitted.
    ids = set(slice_candidates["id"].astype(str))
    slice_offtargets = offtargets[offtargets["qname"].astype(str).isin(ids)].sort_values(["qname", "rname", "coord"])
    slice_mirna = mirna[mirna["qname"].astype(str).isin(ids)].sort_values(["qname", "mirna_id"])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    slice_candidates.to_csv(args.out_dir / "candidates_sample.csv", index=False)
    slice_offtargets.to_csv(args.out_dir / "offtarget_hits_sample.tsv", sep="\t", index=False)
    slice_mirna.to_csv(args.out_dir / "mirna_hits_sample.tsv", sep="\t", index=False)

    print(
        f"guides={len(guides)} candidate_rows={len(slice_candidates)} "
        f"offtarget_rows={len(slice_offtargets)} mirna_rows={len(slice_mirna)}"
    )


if __name__ == "__main__":
    main()
