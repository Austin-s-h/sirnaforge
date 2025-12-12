# Workflows

siRNAforge provides flexible workflows from quick gene-to-siRNA analysis to custom multi-step pipelines.

## Complete Workflow

The `workflow` command handles everything: gene search → design → scoring → filtering.

```bash
sirnaforge workflow TP53 --output-dir results/
```

### Output Structure {#output-structure}

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

## Nextflow Pipeline {#nextflow-pipeline}

For large-scale off-target analysis across multiple genomes:

```bash
# Generate candidates first
sirnaforge workflow TP53 --output-dir results/

# Run Nextflow pipeline
nextflow run nextflow_pipeline/main.nf \
  --candidates results/off_target/input_candidates.fasta \
  --outdir nextflow_results/ \
  --genome_species "human,mouse"
```

See [nextflow_pipeline/README.md](https://github.com/austin-s-h/sirnaforge/blob/master/nextflow_pipeline/README.md) for configuration options.
