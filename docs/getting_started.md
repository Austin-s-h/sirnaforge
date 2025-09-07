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
├── sirna_design/
│   ├── TP53_sirna_candidates.tsv   # All designed siRNAs
│   ├── TP53_top_candidates.tsv     # Top-ranked siRNAs
│   └── TP53_design_summary.json    # Design parameters & stats
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

### siRNA Results Table

The main results file (`*_sirna_candidates.tsv`) contains:

| Column | Description |
|--------|-------------|
| `sirna_id` | Unique siRNA identifier |
| `guide_sequence` | Guide (antisense) strand sequence |
| `passenger_sequence` | Passenger (sense) strand sequence |
| `position` | Position on transcript |
| `gc_content` | GC content percentage |
| `composite_score` | Overall quality score |
| `thermodynamic_score` | Folding energy score |
| `off_target_score` | Off-target prediction score |

### Quality Scores

- **Composite Score**: Overall siRNA quality (higher = better)
- **Thermodynamic Score**: RNA folding favorability
- **Off-target Score**: Specificity prediction (higher = more specific)

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
