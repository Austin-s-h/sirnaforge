# 🧬 siRNAforge Quick Reference

**Comprehensive siRNA design toolkit for gene silencing.** — Essential commands at a glance.

## 🚀 Quick Start

```bash
# Complete workflow (recommended)
uv run sirnaforge workflow TP53 --output-dir results

# Step-by-step
uv run sirnaforge search TP53 -o transcripts.fasta
uv run sirnaforge design transcripts.fasta -o results.tsv
```

## 📋 Essential Commands

| Command | Purpose | Example |
|---------|---------|---------|
| `workflow` | Complete gene → siRNA analysis | `sirnaforge workflow TP53` |
| `search` | Find gene transcripts | `sirnaforge search BRCA1` |
| `design` | Generate siRNA candidates | `sirnaforge design input.fasta` |
| `validate` | Check FASTA files | `sirnaforge validate input.fasta` |
| `config` | Show parameters | `sirnaforge config` |
| `version` | Show version | `sirnaforge version` |

## ⚙️ Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--top-n` | 10-20 | Number of candidates |
| `--gc-min` | 30 | Min GC content (%) |
| `--gc-max` | 52 | Max GC content (%) |
| `--length` | 21 | siRNA length (nt) |
| `--verbose` | - | Detailed output |

## 🎯 Common Patterns

```bash
# High-quality design
sirnaforge design input.fasta --gc-min 35 --gc-max 50 --top-n 30

# Fast processing
sirnaforge design input.fasta --skip-structure --skip-off-targets

# Multi-species analysis
sirnaforge workflow TP53 --genome-species "human,rat,mouse"

# Batch processing
for gene in TP53 BRCA1 EGFR; do
    sirnaforge workflow "$gene" --output-dir "analysis_$gene"
done
```

## 🛠️ Development

```bash
# Setup
make install-dev

# Quality checks
make test && make lint

# Documentation
make docs-cli && make docs-examples
```

## 🔗 Quick Links

- [Full CLI Reference](CLI_REFERENCE.md)
- [Usage Examples](USAGE_EXAMPLES.md)
- [Thermodynamic Metrics Guide](THERMODYNAMIC_GUIDE.md)
- [Thermodynamic Asymmetry Scoring](tutorials/custom_scoring.md)
- [Main README](README.md)

---
*siRNAforge — Professional siRNA design for modern researchers*
