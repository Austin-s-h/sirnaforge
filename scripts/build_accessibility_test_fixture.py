#!/usr/bin/env python
r"""Regenerate the tracked target-accessibility regression fixture from the full benchmark.

The fixture (`tests/unit/data/`) is a deterministic subset of a third-party measured-knockdown
dataset, small enough to track and to fold inside `make test-dev`. This script is the definition
of that subset, so the fixture is regenerable rather than mystery data. See
`tests/unit/data/README.md` for provenance.

Selection rule, in order, with no tuning against the resulting correlation:

0. Keep only rows from the citeable source (``--sources``, default ``Huesken``). The redistributed
   dataset also carries 385 rows labelled ``Other`` whose provenance could not be established, so
   they cannot be attributed and are excluded from the tracked artifact. Sources happen to split
   cleanly by accession, so this drops whole transcripts rather than thinning any of them.
1. Keep only rows whose target site (reverse complement of the guide) occurs **exactly once** in
   its full-length transcript, so every row's coordinates are unambiguous.
2. Rank accessions by siRNA density (rows per kb of transcript) descending, ties broken by
   accession, and take them until the cumulative transcript length would exceed
   ``--max-transcript-nt``. Density is the criterion because the fixture's cost is transcript
   length (folding time and tracked bytes) while its statistical power is row count.
3. Order the surviving rows by (accession, target-site start) and keep every k-th row, with
   ``k = max(1, floor(n / target_rows))``, then truncate to ``--target-rows``.

Usage:
    uv run python scripts/build_accessibility_test_fixture.py \\
        --benchmark-csv work/sirna_bench.csv \\
        --transcripts work/bench_tx.fa \\
        --out-dir tests/unit/data
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

FASTA_WRAP = 60


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


def locate_rows(table: pd.DataFrame, transcripts: dict[str, str]) -> pd.DataFrame:
    """Keep rows whose target site occurs exactly once in its transcript, recording the offset."""
    rows: list[dict[str, object]] = []
    for _, row in table.iterrows():
        sequence = transcripts.get(str(row.B))
        if not sequence:
            continue
        sequence = sequence.upper()
        site = str(Seq(str(row.siRNA_seq).upper()).reverse_complement())
        if sequence.count(site) != 1:
            continue
        rows.append(
            {
                "accession": str(row.B),
                "guide_sequence": str(row.siRNA_seq).upper(),
                "site_start": sequence.find(site),
                "efficacy": float(row.efficacy),
            }
        )
    return pd.DataFrame(rows)


def choose_accessions(located: pd.DataFrame, transcripts: dict[str, str], max_nt: int) -> list[str]:
    """Densest accessions first, until the transcript budget is spent."""
    density = (
        located.groupby("accession")
        .size()
        .rename("rows")
        .to_frame()
        .assign(nt=lambda d: [len(transcripts[a]) for a in d.index])
    )
    density["per_kb"] = 1000.0 * density["rows"] / density["nt"]
    ranked = density.sort_values(["per_kb", "rows"], ascending=[False, False], kind="mergesort")

    kept: list[str] = []
    budget = 0
    for accession, record in ranked.iterrows():
        if budget + int(record.nt) > max_nt:
            continue
        kept.append(str(accession))
        budget += int(record.nt)
    return sorted(kept)


def main() -> None:
    """Write the fixture CSV and FASTA."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark-csv", type=Path, required=True)
    parser.add_argument("--transcripts", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--target-rows", type=int, default=180)
    parser.add_argument(
        "--sources",
        default="Huesken",
        help="Comma-separated values of the dataset's `sources` column to keep (default: Huesken, "
        "the only group with an establishable citation)",
    )
    parser.add_argument("--max-transcript-nt", type=int, default=10000)
    args = parser.parse_args()

    transcripts = read_fasta(args.transcripts)
    table = pd.read_csv(args.benchmark_csv)
    wanted = {s.strip() for s in args.sources.split(",") if s.strip()}
    if "sources" not in table.columns:
        raise SystemExit("benchmark CSV has no `sources` column; cannot restrict to a citeable subset")
    table = table[table.sources.isin(wanted)]
    if table.empty:
        raise SystemExit(f"no rows with sources in {sorted(wanted)}")
    located = locate_rows(table, transcripts)
    kept = choose_accessions(located, transcripts, args.max_transcript_nt)

    subset = located[located.accession.isin(kept)].sort_values(
        ["accession", "site_start"], kind="mergesort", ignore_index=True
    )
    stride = max(1, len(subset) // args.target_rows)
    subset = subset.iloc[::stride].head(args.target_rows).reset_index(drop=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "sirna_efficacy_subset.csv"
    subset.to_csv(csv_path, index=False)

    used = sorted(subset.accession.unique())
    fasta_path = args.out_dir / "sirna_efficacy_subset_transcripts.fa"
    with fasta_path.open("w") as handle:
        for accession in used:
            sequence = transcripts[accession].upper()
            handle.write(f">{accession}\n")
            for i in range(0, len(sequence), FASTA_WRAP):
                handle.write(sequence[i : i + FASTA_WRAP] + "\n")

    total_nt = sum(len(transcripts[a]) for a in used)
    print(f"{len(subset)} rows across {len(used)} transcripts ({total_nt} nt), stride {stride}")
    print(f"sources kept: {sorted(wanted)}; accessions: {', '.join(used)}")
    print(f"efficacy range {subset.efficacy.min():.3f}-{subset.efficacy.max():.3f}")
    print(f"wrote {csv_path} ({csv_path.stat().st_size} B) and {fasta_path} ({fasta_path.stat().st_size} B)")


if __name__ == "__main__":
    main()
