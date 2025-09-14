# Getting Started with siRNAforge

Welcome to siRNAforge! This guide will help you get started with designing high-quality siRNAs.

## Installation

### Prerequisites

- Python 3.9-3.12 (Python 3.13+ not supported due to ViennaRNA compatibility)
- Git
- [uv](https://github.com/astral-sh/uv) (recommended) or pip

### Quick Installation

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge

# One-command setup
make install-dev
```

### Alternative Installation

```bash
# Using pip directly
pip install -e .[dev]

# Or for production only
pip install -e .
```

## First Steps

### 1. Verify Installation

```bash
# Check if siRNAforge is working
uv run sirnaforge --help

# Check version
uv run sirnaforge version
```

### 2. Run Your First Analysis

```bash
# Complete workflow for TP53 gene
uv run sirnaforge workflow TP53 --output-dir my_first_analysis

# This will:
# 1. Search for TP53 transcripts
# 2. Design siRNA candidates
# 3. Score and rank them
# 4. Generate comprehensive results
```

### 3. Explore the Results

After running the workflow, check the `my_first_analysis/` directory:

```
my_first_analysis/
├── transcripts/
│   ├── TP53_transcripts.fasta      # Retrieved transcripts
│   └── TP53_canonical.fasta        # Canonical isoform
├── sirnaforge/
│   ├── TP53_all.csv                # All designed siRNAs (with thermo fields)
│   ├── TP53_pass.csv               # Passing candidates only
│   └── manifest.json               # FAIR metadata for generated files
├── logs/
│   ├── workflow_stream.log         # Console log stream
│   └── workflow_summary.json       # High-level results summary
└── orf_reports/
    └── TP53_orf_validation.txt     # ORF analysis results
```

### 4. Customize Your Analysis

```bash
# High-quality design with strict parameters
uv run sirnaforge workflow BRCA1 \
  --output-dir brca1_strict \
  --gc-min 35 --gc-max 50 \
  --top-n 30 \
  --verbose

# Fast screening with relaxed parameters
uv run sirnaforge workflow EGFR \
  --output-dir egfr_fast \
  --top-n 10 \
  --gc-min 25 --gc-max 65
```

## Basic Workflow Options

### Option 1: Complete Workflow (Recommended)

```bash
# One command does everything
uv run sirnaforge workflow [GENE] --output-dir [OUTPUT]
```

**Advantages:**
- Fully automated
- Consistent parameters
- Complete documentation
- Ready for publication

### Option 2: Step-by-Step Workflow

```bash
# 1. Search for transcripts
uv run sirnaforge search TP53 -o tp53_transcripts.fasta

# 2. Validate input
uv run sirnaforge validate tp53_transcripts.fasta

# 3. Design siRNAs
uv run sirnaforge design tp53_transcripts.fasta -o tp53_results.tsv
```

**Advantages:**
- Full control over each step
- Custom intermediate processing
- Easier debugging

## Understanding the Output

### siRNA Results Tables

Two CSVs are produced in `sirnaforge/`:

- `*_all.csv`: All designed candidates (includes pass/fail) with thermodynamic and structural fields
- `*_pass.csv`: Subset of candidates that pass all filters

Key columns include:

- `id`, `transcript_id`, `position`
- `guide_sequence`, `passenger_sequence`, `gc_content`, `asymmetry_score`
- `structure`, `mfe`, `paired_fraction` (if structure prediction enabled)
- `duplex_stability_dg`, `dg_5p`, `dg_3p`, `delta_dg_end`, `melting_temp_c`
- `transcript_hit_count`, `transcript_hit_fraction`, `off_target_count`
- `composite_score`, `passes_filters`

### Quality Scores

- Composite Score: Overall siRNA quality (higher = better)
- Thermodynamic metrics: ΔG, end-region ΔGs, asymmetry, MFE, approximate Tm
  - **Asymmetry Score**: Critical for guide strand selection (see {doc}`tutorials/custom_scoring` for details)
  - **GC Content**: Optimal 35-60% (ideally 40-55%) for stability/accessibility balance
  - **MFE**: Moderate values (-2 to -8 kcal/mol) preferred for effective RISC processing
  - **End Stability**: Positive delta_dg_end (+2 to +6 kcal/mol) promotes correct strand orientation
- Off-target metrics: Counts/scores added after off-target analysis when available

## Common Parameters

### Quality Control

```bash
--gc-min 30          # Minimum GC content (%)
--gc-max 52          # Maximum GC content (%)
--max-poly-runs 3    # Max consecutive identical nucleotides
--length 21          # siRNA length (nucleotides)
```

### Output Control

```bash
--top-n 20           # Number of candidates to generate
--output-dir results # Output directory
--verbose            # Detailed progress information
```

### Database Options

```bash
--database ensembl   # Use Ensembl (default)
--database refseq    # Use RefSeq
--all               # Search all databases
```

## Troubleshooting

### Common Issues

#### 1. Gene Not Found
```bash
# Error: Gene 'MYGENE' not found
# Solutions:
uv run sirnaforge search MYGENE --all --verbose  # Try all databases
uv run sirnaforge search ENSG00000123456         # Use Ensembl ID directly
```

#### 2. No Valid Transcripts
```bash
# Check what transcripts were found
uv run sirnaforge search MYGENE --no-sequence --verbose

# Try including more transcript types
uv run sirnaforge search MYGENE --types "protein_coding,lncRNA,miRNA"
```

#### 3. Few siRNA Candidates
```bash
# Relax parameters to get more candidates
uv run sirnaforge design input.fasta \
  --gc-min 25 --gc-max 65 \
  --max-poly-runs 4 \
  --top-n 50
```

### Getting Help

```bash
# Command-specific help
uv run sirnaforge workflow --help
uv run sirnaforge design --help

# Show current configuration
uv run sirnaforge config

# Validate your input files
uv run sirnaforge validate your_file.fasta
```

## Next Steps

Once you're comfortable with basic usage:

1. **Explore Advanced Features**: Multi-species off-target analysis
2. **Python API**: Use siRNAforge programmatically
3. **Nextflow Pipeline**: Scale to large datasets
4. **Custom Scoring**: Develop application-specific scoring functions

See the {doc}`tutorials/index` for detailed walkthroughs of advanced features.
