# Getting Started

Get from installation to your first siRNA analysis in minutes.

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

## First Analysis

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
├── transcripts/TP53_transcripts.fasta
└── logs/workflow_summary.json
```

## Understanding Your Results

### Key Output Files

| File | Purpose | Use For |
|------|---------|---------|
| `*_pass.csv` | High-quality candidates that pass all filters | **Laboratory experiments** |
| `*_all.csv` | Complete candidate list with detailed scores | Analysis and custom filtering |
| `workflow_summary.json` | High-level analysis statistics | QC and reporting |

## Quality Indicators

**Key metrics:**
- `asymmetry_score` ≥0.65 (optimal)
- `gc_content` 35-60% (balanced)
- `melting_temp` 55-65°C (effective)

## Step-by-Step Alternative

```bash
# 1. Search transcripts
uv run sirnaforge search TP53 -o transcripts.fasta

# 2. Design candidates
uv run sirnaforge design transcripts.fasta -o results.csv

# 3. Validate input
uv run sirnaforge validate transcripts.fasta
```

## Next Steps

- **📖 [Quick Reference](quick_reference.md)** - Essential commands and parameters
- **🔬 [Usage Examples](usage_examples.md)** - Comprehensive real-world examples
- **⚙️ [CLI Reference](cli_reference.md)** - Complete parameter documentation
- **🧬 [Custom Scoring Guide](tutorials/custom_scoring.md)** - Advanced thermodynamic principles

## Troubleshooting

### Common Issues

**Command not found**: Run `uv run sirnaforge` instead of `sirnaforge`
**No candidates found**: Try relaxing GC content with `--gc-min 25 --gc-max 70`
**Slow performance**: Use `--skip-structure` for faster processing

> **Need help?** Check our [Development Guide](development.md) or open an issue on GitHub.
