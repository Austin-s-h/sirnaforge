# Runnable benchmark inputs

The original synthetic FASTA adapter is retained for design-only smoke tests, but it must not be
used as the target input to the full biological workflow. OligoGym's `fasta` field contains two
oligo strands, not an mRNA transcript; unequal strand lengths encode overhangs or asymmetric
duplex architecture.

For native transcript benchmarking, map the paired region to the downloaded Ensembl human cDNA:

```bash
uv run python scripts/map_oligo_benchmarks_to_ensembl.py \
  --cdna-fasta /path/to/Homo_sapiens.GRCh38.cdna.all.fa \
  --output-dir tests/data/benchmarks/oligogym_native
```

This writes `observations.csv` (one row per assay), `mappings.csv` (one row per native transcript
match), and native target FASTAs for panels compatible with siRNAforge's 19–23 nt designer.
Ichihara uses the 19-nt paired region; Martinelli uses 21 nt. Shmushkovich is a 15-nt/20-nt
asymmetric hsiRNA panel and is mapped for review but not emitted as a design FASTA.

The generated benchmark inputs under `oligogym/` are designed for the existing
`sirnaforge design` command. They are not native transcript benchmarks: OligoGym provides
oligonucleotide duplexes and assay labels, but generally does not provide the full transcript
sequence and site coordinates needed for biologically faithful target-accessibility scoring.

Regenerate them with:

```bash
uv run python scripts/prepare_oligo_benchmarks.py
```

The adapter creates one synthetic transcript per observation:

```text
70 × A + reverse-complement(guide) + 70 × A
```

The known target site begins at position 71 (1-based). The corresponding row in `records.csv`
preserves the measured label, label direction, original dataset row, guide/passenger sequences,
and the synthetic-context warning.

Example runs:

```bash
uv run sirnaforge design \
  tests/data/benchmarks/oligogym/inputs/martinelli_2023_1_len21.fasta \
  --length 21 \
  --output benchmark_results/martinelli_2023_1_len21.csv

uv run sirnaforge design \
  tests/data/benchmarks/oligogym/inputs/shmushkovich_2018_1_len20.fasta \
  --length 20 \
  --output benchmark_results/shmushkovich_2018_1_len20.csv
```

The generated `oligogym/manifest.json` contains commands for every dataset/length combination.
Use the `efficacy_higher_is_better` column for ranking comparisons; Shmushkovich is transformed
from percentage remaining to percentage knockdown for that comparison, while both original labels
remain available.

These files are intended for design-parameter smoke tests and ranking experiments. Do not use
their accessibility values as evidence about real transcript accessibility until native transcript
sequences and coordinates have been joined and version-pinned.
