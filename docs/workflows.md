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
