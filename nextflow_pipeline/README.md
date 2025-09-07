# siRNA Off-target Analysis Pipeline

A comprehensive Nextflow pipeline for performing off-target analysis of siRNA candidates using multiple alignment tools and genome species.

## Overview

This pipeline provides:
- **Sequence validation** - Validates siRNA sequences for correct length and nucleotides
- **Multi-species analysis** - Analyzes off-targets across human, rat, rhesus macaque, and other genomes
- **Dual alignment strategy** - Uses both Bowtie (seed-based) and BWA-MEM2 (sensitive) alignment
- **Seed-aware scoring** - Prioritizes mismatches in seed region (positions 2-8)
- **Comprehensive reporting** - Generates TSV, JSON, and HTML reports

## Quick Start

### Prerequisites

- Nextflow (≥21.04.0)
- Docker or Conda
- BWA-MEM2
- Bowtie
- Python 3.8+ with BioPython

### Basic Usage

```bash
# Run with default parameters

nextflow run nextflow_pipeline/main.nf -profile docker --input nextflow_pipeline/candidates.fasta --outdir results

# Specify custom genome species and parameters
nextflow run nextflow_pipeline/main.nf \
    --input my_sirnas.fasta \
    --outdir my_results \
    --genome_species "human,mouse,rat" \
    --max_hits 5000 \
    --bwa_k 10
```



### Quickstart: Run the Nextflow pipeline with the bundled test dataset

This quickstart shows how to run the Nextflow workflow included in this repository using the example/test FASTA bundled with the project. The pipeline is written to handle missing genome indexes gracefully, so you can exercise the workflow end-to-end without pre-building indexes.

Prerequisites:
- Install Nextflow (https://www.nextflow.io) and ensure Java 11+ is available on your PATH.
- Docker is optional; use `-profile docker` to run tasks in containers if you have Docker installed.

Run the pipeline (local runner):

```bash
# from the repository root
nextflow run ./sirna_redesigned/nextflow_pipeline/main.nf \
  --input transcripts.fasta \
  --outdir results/test_run \
  --genome_species human,rat \
  -with-trace -with-timeline timeline.html
```

Run the pipeline using Docker (optional):

```bash
nextflow run ./sirna_redesigned/nextflow_pipeline/main.nf \
  --input transcripts.fasta \
  --outdir results/test_run_docker \
  --genome_species human \
  -profile docker
```

What to expect:
- Outputs will be written under the `--outdir` you specify (for example `results/test_run`).
- Important files:
  - `combined_offtargets.tsv` — combined tabular hits
  - `combined_offtargets.json` — combined JSON results
  - `final_summary.txt` — text summary of the run
  - `offtarget_report.html` — simple HTML report you can open in a browser

Quick checks after the run:

```bash
ls -l results/test_run
cat results/test_run/final_summary.txt
xdg-open results/test_run/offtarget_report.html  # or open the HTML file in your browser
```

If you want to provide real genome indexes, set `params.genome_indices` in a Nextflow params file or pass a map/JSON to `--genome_indices` when launching Nextflow; otherwise the pipeline will emit warning messages and create placeholder/empty result files for missing indexes (useful for smoke-testing the workflow).

```

### Integration with Python Workflow

The pipeline is designed to integrate seamlessly with the Python siRNA design workflow:

```python
from sirna_design.workflow import run_sirna_workflow

# Run complete workflow including off-target analysis
results = await run_sirna_workflow(
    gene_query="TP53",
    output_dir="tp53_analysis",
    genome_species=["human", "rat", "rhesus"]
)
```

## Configuration

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input` | None | Input FASTA file with siRNA candidates |
| `outdir` | 'results' | Output directory |
| `genome_species` | 'human,rat,rhesus' | Comma-separated species list |
| `sirna_length` | 21 | Expected siRNA length |
| `max_hits` | 10000 | Maximum hits per candidate |
| `bwa_k` | 12 | BWA seed length |
| `bwa_T` | 15 | BWA minimum score threshold |
| `seed_start` | 2 | Seed region start (1-based) |
| `seed_end` | 8 | Seed region end (1-based) |

### Genome Indices

Configure genome index paths in `nextflow.config`:

```groovy
params {
    genome_indices = [
        'human': '/data/genomes/human/bwa_index/GRCh38',
        'rat': '/data/genomes/rat/bwa_index/Rnor6',
        'rhesus': '/data/genomes/rhesus/bwa_index/Mmul10',
        'mouse': '/data/genomes/mouse/bwa_index/GRCm39'
    ]
}
```

### Execution Profiles

The pipeline includes several execution profiles:

```bash
# Local execution (default)
nextflow run main.nf -profile local

# Cluster execution (SLURM)
nextflow run main.nf -profile cluster

# Cloud execution (AWS Batch)
nextflow run main.nf -profile cloud

# Test profile (minimal resources)
nextflow run main.nf -profile test
```

## Pipeline Steps

### 1. Sequence Validation (`VALIDATE_INPUT`)

- Validates sequence length matches `sirna_length`
- Checks for valid nucleotides (A,T,C,G)
- Removes invalid sequences
- Generates validation report

**Outputs:**
- `validation_report.txt` - Summary of validation results
- `validated_sequences.fasta` - Clean sequences for analysis

### 2. Off-target Analysis (`RUN_OFFTARGET_ANALYSIS`)

For each genome species:

- Runs Bowtie for fast seed-based alignment
- Runs BWA-MEM2 for sensitive full-length alignment
- Combines and filters results
- Scores alignments based on seed region mismatches

**Outputs (per species):**
- `{species}_offtargets.tsv` - Tabular results
- `{species}_offtargets.json` - Structured data
- `{species}_summary.txt` - Analysis summary

### 3. Results Aggregation (`AGGREGATE_RESULTS`)

- Combines results from all species
- Generates cross-species statistics
- Creates HTML report with visualizations

**Final Outputs:**
- `combined_offtargets.tsv` - All hits across species
- `combined_offtargets.json` - Complete structured data
- `final_summary.txt` - Overall summary
- `offtarget_report.html` - Interactive HTML report

## Output Files

### TSV Format

Tab-delimited file with columns:
- `qname` - Query sequence name
- `qseq` - Query sequence
- `rname` - Reference sequence name
- `coord` - Genomic coordinate
- `strand` - Alignment strand
- `cigar` - CIGAR alignment string
- `mapq` - Mapping quality
- `AS` - Alignment score
- `NM` - Number of mismatches
- `mismatch_count` - Total mismatches
- `mismatch_positions` - Positions of mismatches
- `offtarget_score` - Composite off-target score
- `source` - Alignment tool (bowtie/bwa-mem2)
- `species` - Genome species

### JSON Format

Structured format with detailed mismatch information:

```json
[
  {
    "qname": "siRNA_001",
    "qseq": "AAGCUGUGAAGCAACGCGAAG",
    "rname": "chr1",
    "coord": "chr1:1000000",
    "strand": "+",
    "mismatch_positions": [
      {"pos": 3, "is_seed": true},
      {"pos": 15, "is_seed": false}
    ],
    "offtarget_score": 16.0,
    "species": "human"
  }
]
```

### HTML Report

Interactive report featuring:
- Executive summary
- Per-species hit counts
- Parameter details
- Quality metrics
- Download links for all results

## Advanced Usage

### Custom Scoring

Modify the scoring algorithm in `modules/offtarget_wrapper.py`:

```python
def score_alignment(self, alignment):
    """Custom scoring function"""
    mismatch_positions = alignment.get('mismatch_positions', [])
    
    score = 0.0
    for pos in mismatch_positions:
        if self.seed_start <= pos <= self.seed_end:
            # Custom seed penalty
            score += 10.0  
        else:
            # Custom non-seed penalty
            score += 2.0
    
    return score
```

### Adding New Species

1. Add BWA index to `nextflow.config`:
```groovy
params.genome_indices.zebrafish = '/data/genomes/zebrafish/bwa_index/GRCz11'
```

2. Include in species list:
```bash
nextflow run main.nf --genome_species "human,rat,rhesus,zebrafish"
```

### Resource Optimization

Adjust resources in `nextflow.config`:

```groovy
process {
    withName:RUN_OFFTARGET_ANALYSIS {
        cpus = 16
        memory = '32 GB'
        time = '8h'
    }
}
```

## Integration Examples

### With Python Workflow

```python
# Complete siRNA design + off-target workflow
import asyncio
from sirna_design.workflow import run_sirna_workflow

async def analyze_gene():
    results = await run_sirna_workflow(
        gene_query="BRCA1",
        output_dir="brca1_analysis",
        top_n_candidates=50,
        top_n_offtarget=20,
        genome_species=["human", "mouse", "rat"]
    )
    
    print(f"Analysis complete: {results['workflow_config']['processing_time']:.1f}s")
    return results

# Run the analysis
results = asyncio.run(analyze_gene())
```

### Standalone Nextflow

```bash
# Generate siRNA candidates first
uv run sirna design transcripts.fasta --output candidates.csv --top-n 100

# Extract sequences for off-target analysis
awk -F'\t' 'NR>1 {print ">"$1"\\n"$4}' candidates.csv > candidates.fasta

# Run off-target pipeline
nextflow run main.nf \
    --input candidates.fasta \
    --outdir offtarget_results \
    --genome_species "human,rat,rhesus" \
    --max_hits 1000
```

## Troubleshooting

### Common Issues

**1. Index not found**
```
Error: Index path /data/genomes/human/bwa_index does not exist
```
Solution: Update `params.genome_indices` in `nextflow.config`

**2. BWA-MEM2 not found**
```
Error: bwa-mem2: command not found
```
Solution: Install BWA-MEM2 or use Docker profile

**3. Memory issues**
```
Process exceeded available memory
```
Solution: Increase memory in process configuration

### Performance Tips

1. **Use SSD storage** for genome indices
2. **Optimize BWA parameters** for your sequence length
3. **Limit max_hits** for faster processing
4. **Use cluster profile** for large datasets

## Citation

If you use this pipeline, please cite:

```bibtex
@software{sirna_offtarget_pipeline,
  title={siRNA Off-target Analysis Pipeline},
  author={siRNA Design Toolkit Team},
  year={2024},
  url={https://github.com/your-org/sirna-design}
}
```

## Support

For questions and issues:
- Create an issue on GitHub
- Check the documentation
- Contact the development team

## License

This pipeline is released under the MIT License.
