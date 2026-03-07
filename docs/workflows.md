# Workflows

siRNAforge provides flexible workflows from quick gene-to-siRNA analysis to custom multi-step pipelines.

## Complete Workflow

The `workflow` command handles everything: gene search → design → scoring → filtering.

```bash
sirnaforge workflow TP53 --output-dir results/
```

(output-structure)=
### Output Structure

```
results/
├── sirnaforge/
│   ├── TP53_pass.csv          # ✓ Use these for experiments
│   └── TP53_all.csv           # All candidates with scores
├── transcripts/
│   ├── TP53_transcripts.fasta # Retrieved sequences
│   └── TP53_canonical.fasta   # Canonical isoform
├── logs/
│   └── workflow_summary.json  # Analysis statistics
└── off_target/
    └── input_candidates.fasta # For Nextflow pipeline
```

### Common Options

```bash
# Stricter quality filters
sirnaforge workflow TP53 --gc-min 35 --gc-max 55 --top-n 30

# More candidates
sirnaforge workflow BRCA1 --top-n 50

# Verbose output
sirnaforge workflow EGFR --verbose
```

## Step-by-Step Workflow

For more control, run individual commands:

```bash
# 1. Get transcripts
sirnaforge search TP53 -o transcripts.fasta

# 2. Design candidates
sirnaforge design transcripts.fasta -o candidates.csv --top-n 20
```

## Batch Processing

Process multiple genes efficiently:

```bash
for gene in TP53 BRCA1 EGFR KRAS; do
  sirnaforge workflow $gene --output-dir results/$gene
done
```

## ZFN Workflow

Use ZFN mode when you want half-site constrained off-target discovery.

```bash
sirnaforge workflow CCR5_ZFN_RUN \
  --design-mode zfn \
  --zfn-left-half-site GTCATCCTCATC \
  --zfn-right-half-site AAACTGCAAAAG \
  --zfn-search-space ensembl_human_hg38_primary \
  --zfn-spacer-lengths 5,6 \
  --zfn-max-mismatches 2 \
  --zfn-algorithm zfn_v2 \
  --output-dir zfn_output/
```

Expected artifacts:

- `zfn_output/sirnaforge/zfn_candidate_summary.json`
- `zfn_output/sirnaforge/zfn_offtarget_sites.csv`
- `zfn_output/logs/workflow_summary.json`

Advanced runtime controls:

- `SIRNAFORGE_ZFN_USE_NEXTFLOW=1` to force Nextflow-backed ZFN search path.
- `SIRNAFORGE_ZFN_SHARDING_JSON='{"enabled": true, "chunk_size_mb": 3, "overlap_bp": 40000}'` for custom shard/window behavior on large references.

## Variant Targeting (documented runs)

Five offline variant-resolution runs were executed with the bundled FASTA (`examples/sample_transcripts.fasta`) and demo VCF (`examples/variant_demo.vcf`). Minimal outputs live under `docs/documented_workflows/`:

- Summary: `docs/documented_workflows/variant_examples.json` (command, mode, min-af, resolved variants)
- Reports: `docs/documented_workflows/*/logs/resolved_variants.json`

Example commands mirrored in the summary:

```bash
# Avoid mode with default AF threshold
sirnaforge workflow TP53 --input-fasta examples/sample_transcripts.fasta --snp-file examples/variant_demo.vcf --variant-mode avoid --min-af 0.01 --output-dir workflow_output/avoid_baseline

# Target common alleles
sirnaforge workflow TP53 --input-fasta examples/sample_transcripts.fasta --snp-file examples/variant_demo.vcf --variant-mode target --min-af 0.03 --output-dir workflow_output/target_focus
```

## Docker Workflow

```bash
docker run --rm -v $(pwd):/data -w /data \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir results/
```

:::{note}
Inside the container, Nextflow uses the `local` profile automatically—no Docker-in-Docker needed.
:::

(nextflow-pipeline)=
## Nextflow Pipeline

siRNAforge runs off-target analysis via an embedded Nextflow workflow.

```bash
# Full workflow (includes embedded Nextflow off-target analysis)
sirnaforge workflow TP53 --output-dir results/

# Off-target only (runs the embedded Nextflow workflow)
sirnaforge offtarget \
  --input-candidates-fasta results/off_target/input_candidates.fasta \
  --output-dir results/off_target \
  --species "human,mouse"
```

Advanced: discover and run the embedded pipeline directly:

```bash
PIPELINE_NF=$(uv run python -c "from sirnaforge.pipeline.nextflow.runner import NextflowRunner; print(NextflowRunner().get_main_workflow())")
nextflow run "$PIPELINE_NF" --help
```

## Developer Validation Notes

For maintainers and advanced users validating the internal ZFN Nextflow bridge behavior (including `uv run --project ...` usage outside repository root), see:

- [ZFN Nextflow Bridge: Sanity Validation](developer/zfn_nextflow_bridge_validation.md)
