# Quick Reference

Essential commands and parameters.

## Quick Start

```bash
# Complete workflow
uv run sirnaforge workflow TP53 --output-dir results

# Step-by-step
uv run sirnaforge search TP53 -o transcripts.fasta
uv run sirnaforge design transcripts.fasta -o results.csv
```

## Commands

| Command | Purpose | Key Options |
|---------|---------|------------|
| `workflow` | Gene → siRNA pipeline | `--output-dir`, `--top-n` |
| `search` | Find transcripts | `--database`, `--output` |
| `design` | Generate siRNAs | `--top-n`, `--gc-min/max` |
| `validate` | Check FASTA | `--verbose` |
| `version` | Show version | - |

## Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--gc-min` / `--gc-max` | 30 % / 52 % | GC content window |
| `--length` (`design`) / `--length` (`workflow`) | 21 | Candidate length in nt |
| `--top-n` | 10 (`design`), 20 (`workflow`) | Number of candidates retained |
| `--output-dir` | varies (`sirna_workflow_output`) | Output directory for workflows |
| `--output` | `sirna_results.tsv` | Output path for single-step commands |
| `--verbose` | `false` | Emit detailed progress messages |

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

## ℹ️ More Information

> **Need examples?** See [Usage Examples](usage_examples.md)
> **Full parameters?** See [CLI Reference](cli_reference.md)
> **Development setup?** See [Development Guide](development.md)
> **Understanding results?** See [Thermodynamic Guide](thermodynamic_guide.md)
