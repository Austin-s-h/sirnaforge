# siRNA Off-target Analysis Pipeline

A comprehensive Nextflow pipeline for performing off-target analysis of siRNA candidates using the sirnaforge Python package.

## Overview

This pipeline provides:
- **Sequence validation** - Validates siRNA sequences for correct length and nucleotides using sirnaforge utilities
- **Multi-species analysis** - Analyzes off-targets across human, rat, rhesus macaque, and other genomes
- **BWA-MEM2 alignment** - Uses sensitive full-length alignment with sirnaforge BwaAnalyzer
- **Seed-aware scoring** - Prioritizes mismatches in seed region (positions 2-8)
- **Comprehensive reporting** - Generates TSV, JSON, and HTML reports

## Quick Start

### Prerequisites

- Nextflow (≥21.04.0)
- Docker or Conda
- BWA-MEM2 (for alignment)
- Python 3.8+ with sirnaforge package installed

### Basic Usage

```bash
# Run with default parameters
nextflow run nextflow_pipeline/main.nf \
    -profile docker \
    --input nextflow_pipeline/candidates.fasta \
    --outdir results

# Specify custom genome species and parameters
nextflow run nextflow_pipeline/main.nf \
    --input my_sirnas.fasta \
    --genome_species "human,rat" \
    --bwa_k 10 \
    --outdir results
```

### Test Run

Run the pipeline with the bundled test dataset:

```bash
# Local execution
nextflow run nextflow_pipeline/main.nf \
    --input nextflow_pipeline/candidates.fasta \
    --outdir results/test_run \
    --genome_species human,rat \
    -with-trace -with-timeline timeline.html

# Docker execution
nextflow run nextflow_pipeline/main.nf \
    --input nextflow_pipeline/candidates.fasta \
    --outdir results/test_run_docker \
    -profile docker
```

Expected outputs:
- `combined_offtargets.tsv` — combined tabular hits
- `combined_offtargets.json` — structured JSON results
- `final_summary.txt` — analysis summary
- `offtarget_report.html` — interactive HTML report

## Configuration

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input` | None | Input FASTA file with siRNA candidates |
| `outdir` | 'results' | Output directory |
| `genome_species` | 'human,rat,rhesus' | Comma-separated species list |
| `genome_config` | 'genomes.yaml' | Genome configuration file |
| `sirna_length` | 21 | Expected siRNA length |
| `max_hits` | 10000 | Maximum hits per candidate |
| `bwa_k` | 12 | BWA seed length |
| `bwa_T` | 15 | BWA minimum score threshold |
| `seed_start` | 2 | Seed region start (1-based) |
| `seed_end` | 8 | Seed region end (1-based) |
| `download_indexes` | false | Auto-download S3 index prefixes |

### Genome Indices

Configure genome index paths in `genomes.yaml`:

```yaml
human:
  index_prefix: /path/to/human/genome/index

rat:
  index_prefix: /path/to/rat/genome/index

rhesus:
  index_prefix: s3://ngi-igenomes/igenomes/Macaca_mulatta/NCBI/Mmul_10/Sequence/WholeGenomeFasta/
```

Or pass at runtime:
```bash
--genome_indices "human:/abs/path/to/GRCh38,rat:/abs/path/to/Rnor6"
```

### Execution Profiles

```bash
# Local execution (default)
nextflow run main.nf -profile local

# Docker execution
nextflow run main.nf -profile docker

# Test profile (minimal resources)
nextflow run main.nf -profile test
```

## Pipeline Steps

### 1. Sequence Validation (`VALIDATE_INPUT`)

- Validates sequence length matches `sirna_length`
- Checks for valid nucleotides (A,T,C,G) using sirnaforge FastaUtils
- Removes invalid sequences
- Generates validation report

**Outputs:**
- `validation_report.txt` - Summary of validation results
- `validated_sequences.fasta` - Clean sequences for analysis

### 2. Index Resolution (`RESOLVE_GENOME_INDEX`)

- Resolves genome index paths from configuration
- Supports local paths and S3 URLs
- Handles index downloading if enabled
- Validates index existence

**Outputs:**
- Index paths for each species

### 3. Off-target Analysis (`RUN_OFFTARGET_ANALYSIS`)

For each genome species:
- Uses sirnaforge BwaAnalyzer for BWA-MEM2 alignment
- Filters and ranks results by off-target score
- Scores alignments based on seed region mismatches
- Generates comprehensive hit annotations

**Outputs (per species):**
- `{species}_offtargets.tsv` - Tabular results
- `{species}_offtargets.json` - Structured data
- `{species}_summary.txt` - Analysis summary

### 4. Results Aggregation (`AGGREGATE_RESULTS`)

- Combines results from all species using sirnaforge utilities
- Generates cross-species statistics
- Creates HTML report with visualizations

**Final Outputs:**
- `combined_offtargets.tsv` - All hits across species
- `combined_offtargets.json` - Complete structured data
- `final_summary.txt` - Overall summary
- `offtarget_report.html` - Interactive HTML report

## Integration with sirnaforge

The pipeline is tightly integrated with the sirnaforge Python package:

- **FastaUtils**: For sequence parsing and validation
- **BwaAnalyzer**: For BWA-MEM2 alignment and scoring
- **Off-target utilities**: For result aggregation and reporting

This ensures consistency with the broader sirnaforge ecosystem and provides access to advanced scoring algorithms.

### Integration with Python Workflow

The pipeline integrates seamlessly with the sirnaforge Python workflow:

```python
from sirnaforge.workflow import run_sirna_workflow

# Run complete workflow including off-target analysis
results = await run_sirna_workflow(
    gene_query="TP53",
    genome_species=["human", "rat", "rhesus"]
)
```

## Troubleshooting

### Common Issues

1. **Index not found**: Check `genomes.yaml` paths or use `--genome_indices`
2. **BWA-MEM2 not found**: Ensure BWA-MEM2 is installed and in PATH
3. **Empty results**: Check input sequences are valid siRNAs (21 nt, valid nucleotides)
4. **Container issues**: Use `--container` to specify custom Docker image
5. **sirnaforge import errors**: Ensure sirnaforge package is installed in container/environment

### Debug Mode

```bash
nextflow run main.nf --verbose -with-trace -with-timeline timeline.html
```

### Index Setup

The pipeline requires genome indexes for BWA-MEM2. See [INDEX_README.md](INDEX_README.md) for detailed setup instructions.

## Example Output Structure

```
results/
├── validation/
│   ├── validation_report.txt
│   └── validated_sequences.fasta
├── human/
│   ├── human_offtargets.tsv
│   ├── human_offtargets.json
│   └── human_summary.txt
├── rat/
│   ├── rat_offtargets.tsv
│   ├── rat_offtargets.json
│   └── rat_summary.txt
├── combined_offtargets.tsv
├── combined_offtargets.json
├── final_summary.txt
└── offtarget_report.html
```

## Advanced Usage

### Custom Scoring Parameters

```bash
nextflow run main.nf \
    --seed_start 1 \
    --seed_end 7 \
    --bwa_k 8 \
    --bwa_T 20 \
    --max_hits 5000
```

### Multi-species Analysis

```bash
nextflow run main.nf \
    --genome_species "human,rat,mouse,rhesus" \
    --outdir results/multi_species
```

### S3 Index Support

Enable automatic downloading of genome indexes from S3:

```bash
nextflow run main.nf \
    --download_indexes true \
    --genome_config genomes.yaml
```

## Performance Notes

- BWA-MEM2 alignment is the most resource-intensive step
- Memory usage scales with genome size and number of candidates
- Use `--max_hits` to limit memory usage for large genomes
- Docker execution may be slower but provides better reproducibility

## Contact

See the main sirnaforge repository README for support and contribution guidelines.
