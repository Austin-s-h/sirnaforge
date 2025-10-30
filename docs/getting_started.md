# Getting Started & Quick Reference

Get from installation to your first siRNA analysis in minutes, with essential commands and parameters at your fingertips.

## Installation

**Prerequisites:** Python 3.9-3.12, Git

```bash
# Install uv package manager
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make install-dev
```

**Alternative:** Docker, pip, production setups → [Deployment Guide](deployment.md)

## Your First Analysis

### Verify Installation
```bash
uv run sirnaforge --help
uv run sirnaforge version
```

### Complete Workflow
```bash
# End-to-end analysis for TP53
uv run sirnaforge workflow TP53 --output-dir my_first_analysis
```

**What this does:**
1. Search TP53 transcripts from Ensembl
2. Design siRNA candidates
3. Score using thermodynamic asymmetry
4. Rank and filter best candidates
5. Generate results

### Results Structure
```
my_first_analysis/
├── sirnaforge/
│   ├── TP53_pass.csv        # Best candidates for lab use
│   └── TP53_all.csv         # All candidates with scores
├── transcripts/
│   ├── TP53_transcripts.fasta      # All retrieved transcripts
│   ├── TP53_canonical.fasta        # Canonical isoform only
│   └── temp_for_design.fasta       # Transcripts used for design
├── logs/
│   ├── workflow_stream.log         # Console output
│   └── workflow_summary.json       # High-level summary
└── orf_reports/
    └── TP53_orf_validation.txt     # ORF analysis results
```

### Examine Your Results

```bash
# View top candidates
head -6 my_first_analysis/sirnaforge/TP53_pass.csv

# Check workflow summary
cat my_first_analysis/logs/workflow_summary.json
```

## Understanding Your Results

### Key Output Files

| File | Purpose | Use For |
|------|---------|---------|
| `*_pass.csv` | High-quality candidates that pass all filters | **Laboratory experiments** |
| `*_all.csv` | Complete candidate list with detailed scores | Analysis and custom filtering |
| `workflow_summary.json` | High-level analysis statistics | QC and reporting |

> **Data Models**: Results contain [`SiRNACandidate`](api_reference.rst) objects with full thermodynamic analysis. See [API Reference](api_reference.rst) for field descriptions.

## Quality Indicators

**Key metrics:**
- `asymmetry_score` ≥0.65 (optimal)
- `gc_content` 35-60% (balanced)
- `melting_temp` 55-65°C (effective)

### Interpreting Scores

**Score Ranges:**
- **Composite Score (7-10)**: Higher is better; combines all factors
- **Thermodynamic Score (0-1)**: Measures RNA folding favorability and asymmetry
- **GC Content (35-60%)**: Optimal range 40-55% for stability and accessibility
- **Asymmetry Score (0.65-1.0)**: Higher values indicate better guide strand selection
- **MFE (-2 to -8 kcal/mol)**: Moderate stability preferred for effective processing

## Customizing Your Analysis

### High-Quality Design
```bash
# Stricter quality parameters for research publications
uv run sirnaforge workflow TP53 \
  --output-dir tp53_high_quality \
  --gc-min 35 --gc-max 50 \
  --top-n 30 \
  --verbose
```

### Step-by-Step Workflow
```bash
# 1. Search transcripts
uv run sirnaforge search TP53 -o transcripts.fasta --verbose

# 2. Validate sequences
uv run sirnaforge validate transcripts.fasta

# 3. Design siRNAs
uv run sirnaforge design transcripts.fasta \
  -o results.csv \
  --top-n 25 \
  --verbose
```

### Batch Processing Multiple Genes
```bash
# Process multiple genes
genes=("TP53" "BRCA1" "EGFR")
for gene in "${genes[@]}"; do
    uv run sirnaforge workflow "$gene" \
        --output-dir "analysis_$gene" \
        --gc-min 35 --gc-max 55 \
        --top-n 20
done
```

## Essential Commands & Parameters

### Core Commands

| Command | Purpose | Key Options |
|---------|---------|------------|
| `workflow` | Gene → siRNA pipeline | `--output-dir`, `--top-n` |
| `search` | Find transcripts | `--database`, `--output` |
| `design` | Generate siRNAs | `--top-n`, `--gc-min/max` |
| `validate` | Check FASTA | `--verbose` |
| `version` | Show version | - |

### Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--gc-min` / `--gc-max` | 30 % / 52 % | GC content window |
| `--length` | 21 | Candidate length in nt |
| `--top-n` | 10 (`design`), 20 (`workflow`) | Number of candidates retained |
| `--output-dir` | `sirna_workflow_output` | Output directory for workflows |
| `--output` | `sirna_results.tsv` | Output path for single-step commands |
| `--verbose` | `false` | Emit detailed progress messages |

### Common Usage Patterns

```bash
# High-quality candidates
uv run sirnaforge workflow GENE --gc-min 35 --gc-max 50 --top-n 30

# Fast screening
uv run sirnaforge design input.fasta --skip-structure --top-n 10

# Multi-species analysis
uv run sirnaforge workflow GENE --genome-species "human,mouse,rat"

# Custom parameters
uv run sirnaforge workflow GENE --length 19 --gc-min 40 --max-poly-runs 2
```

### Quality Thresholds

**Recommended for lab experiments:**
- GC Content: 35-60% (optimal: 40-55%)
- Asymmetry Score: ≥0.65
- Melting Temperature: 55-65°C
- End Stability (ΔΔG): +2 to +6 kcal/mol

> **Configuration**: These thresholds are enforced through [`FilterCriteria`](api_reference.rst) in [`DesignParameters`](api_reference.rst). See [API docs](api_reference.rst) for customization options.

## Next Steps

- **� [Usage Examples](usage_examples.md)** - Comprehensive real-world examples
- **⚙️ [CLI Reference](cli_reference.md)** - Complete parameter documentation
- **🧬 [Custom Scoring Guide](tutorials/custom_scoring.md)** - Advanced thermodynamic principles

## Troubleshooting

### Common Issues

**Command not found**: Run `uv run sirnaforge` instead of `sirnaforge`
**No candidates found**: Try relaxing GC content with `--gc-min 25 --gc-max 70`
**Slow performance**: Use `--skip-structure` for faster processing

> **Need help?** Check our [Development Guide](developer/development.md) or open an issue on GitHub.
