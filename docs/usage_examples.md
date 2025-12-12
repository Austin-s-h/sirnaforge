# Usage Examples

Scenario-driven snippets that build on the basics from [Getting Started](getting_started.md) and [Workflows](workflows.md). Each example combines Typer commands from `src/sirnaforge/` exactly as they run in production.

## Advanced Workflow Presets

### Publication-Quality Screen

```bash
uv run sirnaforge workflow TP53 \
  --output-dir tp53_publication \
  --top-n 50 \
  --gc-min 35 --gc-max 60 \
  --max-poly-runs 2 \
  --genome-species "human,mouse,rat" \
  --design-mode mirna \
  --mirna-db mirgenedb \
  --mirna-species "human,mouse" \
  --modification-file examples/modification_patterns/standard_2ome.json \
  --verbose
```

- Mirrors the stricter defaults enforced inside `sirnaforge/core/design.py` and `core/off_target.py` (multi-species alignments, miRNA-aware scoring, chemical modifications).
- Use the same flags inside Docker by prefixing with `docker run --rm -v $(pwd):/workspace -w /workspace ghcr.io/austin-s-h/sirnaforge:latest sirnaforge …`.

### Custom FASTA + Transcriptome Override

```bash
uv run sirnaforge workflow TP53 \
  --input-fasta examples/sample_transcripts.fasta \
  --transcriptome-fasta ensembl_mouse_cdna \
  --species "human,mouse" \
  --output-dir tp53_custom_inputs
```

- `--input-fasta` feeds the enumerator exactly the sequences you supply while still naming outputs with the positional gene argument.
- `--transcriptome-fasta` reindexes the transcriptome reference used by off-target analysis (local path, URL, or preset). Presets map to the cache registry defined in `sirnaforge/data/species_registry.py`.
- Pair with `--offtarget-indices human:/refs/hg38 mouse:/refs/mm39` when you maintain your own BWA-MEM2 indices.

## Gene Search & Design Automations

### Curated Transcript Panel

```bash
mkdir -p transcripts
for gene in TP53 BRCA1 EGFR KRAS; do
    uv run sirnaforge search "$gene" \
      --database ensembl \
      --output transcripts/${gene}.fasta \
      --verbose
done

uv run sirnaforge design transcripts/*.fasta \
  --output curated_panel.csv \
  --top-n 25 \
  --gc-min 40 --gc-max 55
```

- `sirnaforge search` talks to the live providers implemented in `sirnaforge/data/gene_search.py`; constraining `--database` avoids redundant API calls.
- `sirnaforge design` accepts globbed FASTA inputs, so you can mix curated transcripts before or after manual review.

## Batch Processing with QC Summaries

```bash
#!/usr/bin/env bash
set -euo pipefail

genes=(TP53 BRCA1 EGFR MYC KRAS)
out_dir="batch_$(date +%Y%m%d)"
mkdir -p "$out_dir"

for gene in "${genes[@]}"; do
  uv run sirnaforge workflow "$gene" \
    --output-dir "${out_dir}/${gene}" \
    --top-n 30 \
    --gc-min 35 --gc-max 60 \
    --species "human,mouse" \
    --verbose

  pass_csv="${out_dir}/${gene}/sirnaforge/${gene}_pass.csv"
  total_csv="${out_dir}/${gene}/sirnaforge/${gene}_all.csv"

  python - <<'PY'
import csv, pathlib, sys
pass_path, total_path = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
if not pass_path.exists() or not total_path.exists():
    raise SystemExit(f"Missing outputs for {pass_path.parent.parent.name}")
with total_path.open() as fh:
    total = sum(1 for _ in csv.DictReader(fh))
with pass_path.open() as fh:
    passing = sum(1 for _ in csv.DictReader(fh))
rate = 0 if total == 0 else passing * 100 / total
print(f"{pass_path.parent.parent.name}: {passing}/{total} ({rate:.1f}%)")
PY "$pass_csv" "$total_csv"
done
```

- Uses the exact CSV schema emitted by `workflow.py` (`sirnaforge/*_pass.csv` and `sirnaforge/*_all.csv`).
- Swap the inline Python block for pandas or DuckDB once you’re ready to integrate with LIMS.

## Design-Only with Chemical Modifications

```bash
uv run sirnaforge design custom_transcripts.fasta \
  --output custom_sirnas.csv \
  --top-n 50 \
  --gc-min 35 --gc-max 60 \
  --max-poly-runs 2 \
  --length 21 \
  --modification-file examples/modification_patterns/fda_approved_onpattro.json
```

- Pattern files live in `examples/modification_patterns/` and are parsed by `sirnaforge/utils/modification_patterns.py`.
- Attach metadata afterwards with:

```bash
sirnaforge sequences annotate custom_sirnas.fasta modifications.json -o custom_sirnas_annotated.fasta
```

## miRNA-Biogenesis Mode

```bash
uv run sirnaforge workflow MIR21 \
  --design-mode mirna \
  --gc-min 32 --gc-max 52 \
  --top-n 60 \
  --mirna-db mirgenedb \
  --mirna-species "human" \
  --output-dir mirna_mode_demo
```

- Activates the alternative scoring weights and seed-cleanliness heuristics defined in `sirnaforge/models/sirna.py` (`DesignMode.MIRNA`).
- Expect additional `quality_flags` for miRNA-specific failures (dirty seed windows, lack of supplementary pairing, etc.).

## Nextflow Off-Target Pipeline Hand-Off

```bash
uv run sirnaforge workflow TP53 --output-dir tp53_nextflow

nextflow run nextflow_pipeline/main.nf \
  --candidates tp53_nextflow/off_target/input_candidates.fasta \
  --outdir tp53_nextflow/off_target_results \
  --genome_species "human,mouse"
```

- `off_target/input_candidates.fasta` contains both high-confidence and dirty-control guides—the exact payload expected by `nextflow_pipeline/modules/local/offtarget_align.nf`.
- Override indices with `--offtarget-indices human:/refs/hg38 mouse:/refs/mm39` for custom BWA-MEM2 builds.

## Docker-Friendly Pattern

Wrap any snippet with:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge <command> <args>
```

- The image caches Ensembl transcriptomes and miRNA databases (see {ref}`Installation → Docker <docker-full-bioinformatics-stack>`).
- Replace `--rm` with `--name sirnaforge_batch` for long-lived runs so you can inspect logs mid-execution.

## See Also

- [Getting Started](getting_started.md) — installation and your first workflow.
- [Workflows](workflows.md) — command-by-command explanations and output diagrams.
- [Scoring Overview](scoring.md) — interpret the metrics in `*_pass.csv`.
- [CLI Reference](cli_reference.md) — authoritative `--help` output captured during the docs build.
