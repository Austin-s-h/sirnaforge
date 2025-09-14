# 🧬 siRNAforge Quick Reference

**Essential commands and parameters at a glance.**

## 🚀 Quick Start

```bash
# Complete workflow (most common)
uv run sirnaforge workflow TP53 --output-dir results

# Step-by-step approach
uv run sirnaforge search TP53 -o transcripts.fasta
uv run sirnaforge design transcripts.fasta -o results.csv
```

## 📋 Core Commands

| Command | Purpose | Key Options |
|---------|---------|------------|
| `workflow` | Complete gene → siRNA pipeline | `--output-dir`, `--top-n`, `--gc-min/max` |
| `search` | Find gene transcripts | `--database`, `--output`, `--verbose` |
| `design` | Generate siRNA candidates | `--top-n`, `--gc-min/max`, `--length` |
| `validate` | Check FASTA format | `--verbose` |
| `config` | Show current settings | - |
| `version` | Version information | - |

## ⚙️ Essential Parameters

### Quality Control
| Parameter | Default | Range | Purpose |
|-----------|---------|-------|---------|
| `--gc-min` | 30 | 15-50 | Minimum GC content (%) |
| `--gc-max` | 52 | 50-85 | Maximum GC content (%) |
| `--length` | 21 | 19-23 | siRNA length (nucleotides) |
| `--max-poly-runs` | 3 | 2-5 | Max consecutive identical bases |

### Output Control
| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--top-n` | 10-20 | Number of candidates to generate |
| `--output-dir` | - | Output directory for workflow |
| `--output` | stdout | Output file for individual commands |
| `--verbose` | false | Detailed progress information |

## 🎯 Common Usage Patterns

```bash
# High-quality candidates
sirnaforge workflow GENE --gc-min 35 --gc-max 50 --top-n 30

# Fast screening
sirnaforge design input.fasta --skip-structure --top-n 10

# Multi-species analysis
sirnaforge workflow GENE --genome-species "human,mouse,rat"

# Custom parameters
sirnaforge workflow GENE --length 19 --gc-min 40 --max-poly-runs 2
```

## 🔬 Quality Thresholds

**Recommended for lab experiments:**
- GC Content: 35-60% (optimal: 40-55%)
- Asymmetry Score: ≥0.65
- Melting Temperature: 55-65°C
- End Stability (ΔΔG): +2 to +6 kcal/mol

## � More Information

> **Need examples?** See [Usage Examples](USAGE_EXAMPLES.md)
> **Full parameters?** See [CLI Reference](CLI_REFERENCE.md)
> **Development setup?** See [Development Guide](development.md)
> **Understanding results?** See [Thermodynamic Guide](THERMODYNAMIC_GUIDE.md)
