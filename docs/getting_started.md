# Getting Started with siRNAforge

Welcome to siRNAforge! This guide will get you from installation to your first siRNA analysis in minutes.

## Installation

### Prerequisites
- Python 3.9-3.12
- Git

### Quick Installation

```bash
# Install uv (fast Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge
make install-dev
```

> **Need different installation options?** See [Deployment Guide](deployment.md) for Docker, pip, or production setups.

## Your First Analysis

### 1. Verify Installation
```bash
uv run sirnaforge --help
uv run sirnaforge version
```

### 2. Complete Workflow (Recommended)
```bash
# End-to-end analysis for TP53
uv run sirnaforge workflow TP53 --output-dir my_first_analysis
```

This single command will:
1. **Search** for TP53 transcripts from Ensembl
2. **Design** siRNA candidates
3. **Score** them using thermodynamic asymmetry
4. **Rank** and filter the best candidates
5. **Generate** comprehensive results

### 3. Explore Results
```
my_first_analysis/
├── sirnaforge/
│   ├── TP53_pass.csv        # Best candidates (CSV format)
│   └── TP53_all.csv         # All candidates with scores
├── transcripts/
│   └── TP53_transcripts.fasta
└── logs/
    └── workflow_summary.json # Analysis summary
```

> **Want to customize parameters?** See our [Quick Reference](QUICK_REFERENCE.md) for common options or [Usage Examples](USAGE_EXAMPLES.md) for comprehensive examples.

## Understanding Your Results

### Key Output Files

| File | Purpose | Use For |
|------|---------|---------|
| `*_pass.csv` | High-quality candidates that pass all filters | **Laboratory experiments** |
| `*_all.csv` | Complete candidate list with detailed scores | Analysis and custom filtering |
| `workflow_summary.json` | High-level analysis statistics | QC and reporting |

### siRNA Quality Indicators

Look for these in your results:
- **asymmetry_score**: ≥0.65 (optimal thermodynamic asymmetry)
- **gc_content**: 35-60% (balanced stability)
- **melting_temp**: 55-65°C (effective silencing range)

### Step-by-Step Workflow (Alternative)

For more control over the process:

```bash
# 1. Search for transcripts
uv run sirnaforge search TP53 -o transcripts.fasta

# 2. Design candidates
uv run sirnaforge design transcripts.fasta -o results.csv

# 3. Validate quality
uv run sirnaforge validate transcripts.fasta
```

## Next Steps

- **📖 [Quick Reference](QUICK_REFERENCE.md)** - Essential commands and parameters
- **🔬 [Usage Examples](USAGE_EXAMPLES.md)** - Comprehensive real-world examples
- **⚙️ [CLI Reference](CLI_REFERENCE.md)** - Complete parameter documentation
- **🧬 [Custom Scoring Guide](tutorials/custom_scoring.md)** - Advanced thermodynamic principles

## Troubleshooting

### Common Issues

**Command not found**: Run `uv run sirnaforge` instead of `sirnaforge`
**No candidates found**: Try relaxing GC content with `--gc-min 25 --gc-max 70`
**Slow performance**: Use `--skip-structure` for faster processing

> **Need help?** Check our [Development Guide](development.md) or open an issue on GitHub.
