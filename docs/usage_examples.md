# Usage Examples

Real-world examples demonstrating siRNAforge's core features and advanced capabilities.

> **New to siRNAforge?** Start with [Getting Started](getting_started.md) for installation and quick start.

## Core Workflows

### Minimal Example

`````{tab-set}

````{tab-item} uv
```bash
# Simplest workflow - uses all defaults
uv run sirnaforge workflow TP53
```
````

````{tab-item} Docker
```bash
# Simplest workflow - uses all defaults
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53
```
````

`````

**Defaults:** GC 30-60%, length 21nt, top 20 candidates, human genome, miRNA disabled

### Comprehensive Example

`````{tab-set}

````{tab-item} uv
```bash
# Publication-quality analysis with all major features
uv run sirnaforge workflow TP53 \
  --output-dir tp53_analysis \
  --top-n 50 \
  --gc-min 35 --gc-max 60 \
  --max-poly-runs 2 \
  --length 21 \
  --genome-species "human,mouse,rat" \
  --design-mode mirna \
  --mirna-db mirgenedb \
  --mirna-species "human,mouse" \
  --modification-file examples/modification_patterns/standard_2ome.json \
  --verbose
```
````

````{tab-item} Docker
```bash
# Publication-quality analysis with all major features
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 \
    --output-dir tp53_analysis \
    --top-n 50 \
    --gc-min 35 --gc-max 60 \
    --max-poly-runs 2 \
    --length 21 \
    --genome-species "human,mouse,rat" \
    --design-mode mirna \
    --mirna-db mirgenedb \
    --mirna-species "human,mouse" \
    --modification-file examples/modification_patterns/standard_2ome.json \
    --verbose
```
````

`````

**Features:** Multi-species off-target, miRNA seed analysis, chemical modifications, verbose logging

### Custom FASTA + Transcriptome Override

`````{tab-set}

````{tab-item} uv
```bash
# Design directly from a provided FASTA while aligning against mouse transcriptome
uv run sirnaforge workflow TP53 \
  --input-fasta examples/sample_transcripts.fasta \
  --transcriptome-fasta ensembl_mouse_cdna \
  --output-dir tp53_custom_inputs
```
````

````{tab-item} Docker
```bash
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 \
    --input-fasta examples/sample_transcripts.fasta \
    --transcriptome-fasta ensembl_mouse_cdna \
    --output-dir tp53_custom_inputs
```
````

`````

**Tips:**
- Combine `--input-fasta` with a gene query to keep familiar file names while sourcing transcripts from your own assemblies.
- `--transcriptome-fasta` accepts preset identifiers (`ensembl_*`), HTTP(S) URLs, or local files; siRNAforge caches indexes automatically so repeated runs are fast.
- Leave the transcriptome flag unset to fall back to `ensembl_human_cdna`.

## Gene Search

### Minimal Example

`````{tab-set}

````{tab-item} uv
```bash
# Search gene symbol across all databases
uv run sirnaforge search BRCA1

# Save to file
uv run sirnaforge search TP53 --output tp53.fasta
```
````

````{tab-item} Docker
```bash
# Search gene symbol across all databases
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search BRCA1

# Save to file
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search TP53 --output tp53.fasta
```
````

`````

### Comprehensive Example

`````{tab-set}

````{tab-item} uv
```bash
# Multi-database search with filters
uv run sirnaforge search HOTAIR \
  --database ensembl \
  --types "lncRNA,antisense" \
  --all \
  --verbose

# Batch search multiple genes
for gene in TP53 BRCA1 EGFR KRAS; do
    uv run sirnaforge search $gene \
      --output transcripts/${gene}.fasta \
      --database ensembl
done
```
````

````{tab-item} Docker
```bash
# Multi-database search with filters
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search HOTAIR \
    --database ensembl \
    --types "lncRNA,antisense" \
    --all \
    --verbose

# Batch search multiple genes
for gene in TP53 BRCA1 EGFR KRAS; do
    docker run --rm -v $(pwd):/workspace -w /workspace \
      ghcr.io/austin-s-h/sirnaforge:latest \
      sirnaforge search $gene \
        --output transcripts/${gene}.fasta \
        --database ensembl
done
```
````

`````

**Options:** `--database` (ensembl/refseq/gencode), `--output`, `--types`, `--all`, `--verbose`

## Design Only (No Off-Target)

### Minimal Example

`````{tab-set}

````{tab-item} uv
```bash
# Design from custom FASTA
uv run sirnaforge design examples/sample_transcripts.fasta
```
````

````{tab-item} Docker
```bash
# Design from custom FASTA
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge design examples/sample_transcripts.fasta
```
````

`````

### Comprehensive Example

`````{tab-set}

````{tab-item} uv
```bash
# Custom parameters with modifications
uv run sirnaforge design custom_transcripts.fasta \
  --output custom_sirnas.csv \
  --top-n 50 \
  --gc-min 35 --gc-max 60 \
  --max-poly-runs 2 \
  --length 21 \
  --modification-file examples/modification_patterns/standard_2ome.json \
  --verbose
```
````

````{tab-item} Docker
```bash
# Custom parameters with modifications
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge design custom_transcripts.fasta \
    --output custom_sirnas.csv \
    --top-n 50 \
    --gc-min 35 --gc-max 60 \
    --max-poly-runs 2 \
    --length 21 \
    --modification-file examples/modification_patterns/standard_2ome.json \
    --verbose
```
````

`````

## Batch Processing

### Minimal Example

`````{tab-set}

````{tab-item} uv
```bash
# Process multiple genes with defaults
genes=("TP53" "BRCA1" "EGFR")
for gene in "${genes[@]}"; do
    uv run sirnaforge workflow "$gene" \
        --output-dir "batch/$gene"
done
```
````

````{tab-item} Docker
```bash
# Process multiple genes with defaults
genes=("TP53" "BRCA1" "EGFR")
for gene in "${genes[@]}"; do
    docker run --rm -v $(pwd):/workspace -w /workspace \
      ghcr.io/austin-s-h/sirnaforge:latest \
      sirnaforge workflow "$gene" \
        --output-dir "batch/$gene"
done
```
````

`````

### Comprehensive Example

`````{tab-set}

````{tab-item} uv
```bash
#!/bin/bash
# Automated batch with quality control

genes=("TP53" "BRCA1" "EGFR" "MYC" "KRAS")
output="batch_$(date +%Y%m%d)"
qc="${output}/qc_summary.txt"

echo "Batch Analysis - $(date)" > "$qc"

for gene in "${genes[@]}"; do
    echo "Processing $gene..."

    uv run sirnaforge workflow "$gene" \
        --output-dir "${output}/${gene}" \
        --top-n 30 \
        --gc-min 35 --gc-max 60 \
        --genome-species "human,mouse" \
        --verbose

    # Quality check
    pass="${output}/${gene}/sirna_design/${gene}_pass.csv"
    all="${output}/${gene}/sirna_design/${gene}_all.csv"

    if [[ -f "$pass" && -f "$all" ]]; then
        total=$(tail -n +2 "$all" | wc -l)
        passing=$(tail -n +2 "$pass" | wc -l)
        rate=$(( passing * 100 / total ))
        echo "$gene: $passing/$total ($rate%)" >> "$qc"
    else
        echo "$gene: FAILED" >> "$qc"
    fi
done

cat "$qc"
```
````

````{tab-item} Docker
```bash
#!/bin/bash
# Automated batch with quality control

genes=("TP53" "BRCA1" "EGFR" "MYC" "KRAS")
output="batch_$(date +%Y%m%d)"
qc="${output}/qc_summary.txt"

echo "Batch Analysis - $(date)" > "$qc"

for gene in "${genes[@]}"; do
    echo "Processing $gene..."

    docker run --rm -v $(pwd):/workspace -w /workspace \
      ghcr.io/austin-s-h/sirnaforge:latest \
      sirnaforge workflow "$gene" \
        --output-dir "${output}/${gene}" \
        --top-n 30 \
        --gc-min 35 --gc-max 60 \
        --genome-species "human,mouse" \
        --verbose

    # Quality check
    pass="${output}/${gene}/sirna_design/${gene}_pass.csv"
    all="${output}/${gene}/sirna_design/${gene}_all.csv"

    if [[ -f "$pass" && -f "$all" ]]; then
        total=$(tail -n +2 "$all" | wc -l)
        passing=$(tail -n +2 "$pass" | wc -l)
        rate=$(( passing * 100 / total ))
        echo "$gene: $passing/$total ($rate%)" >> "$qc"
    else
        echo "$gene: FAILED" >> "$qc"
    fi
done

cat "$qc"
```
````

`````

## Chemical Modifications

### Minimal Example

`````{tab-set}

````{tab-item} uv
```bash
# Apply built-in pattern
uv run sirnaforge design examples/sample_transcripts.fasta \
  --modification-file examples/modification_patterns/standard_2ome.json
```
````

````{tab-item} Docker
```bash
# Apply built-in pattern
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge design examples/sample_transcripts.fasta \
    --modification-file examples/modification_patterns/standard_2ome.json
```
````

`````

**Built-in patterns:** `standard_2ome.json`, `fda_approved_onpattro.json`, `maximal_stability.json`

### Comprehensive Example

`````{tab-set}

````{tab-item} uv
```bash
# Create custom pattern
cat > custom_mods.json << 'EOFJSON'
{
  "name": "Custom Therapeutic",
  "description": "2'-OMe terminal + PS backbone",
  "modifications": [
    {"position": 1, "strand": "guide", "type": "2'-O-methyl", "base": "*"},
    {"position": 19, "strand": "guide", "type": "2'-O-methyl", "base": "*"},
    {"position": 1, "strand": "passenger", "type": "phosphorothioate", "base": "*"}
  ]
}
EOFJSON

# Apply with thermodynamic recalculation
uv run sirnaforge design examples/sample_transcripts.fasta \
  --modification-file custom_mods.json \
  --output modified.csv \
  --verbose
```
````

````{tab-item} Docker
```bash
# Create custom pattern
cat > custom_mods.json << 'EOFJSON'
{
  "name": "Custom Therapeutic",
  "description": "2'-OMe terminal + PS backbone",
  "modifications": [
    {"position": 1, "strand": "guide", "type": "2'-O-methyl", "base": "*"},
    {"position": 19, "strand": "guide", "type": "2'-O-methyl", "base": "*"},
    {"position": 1, "strand": "passenger", "type": "phosphorothioate", "base": "*"}
  ]
}
EOFJSON

# Apply with thermodynamic recalculation
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge design examples/sample_transcripts.fasta \
    --modification-file custom_mods.json \
    --output modified.csv \
    --verbose
```
````

`````

See [Modification Integration Guide](modification_integration_guide.md) for full specification.

## Python API

### Minimal Example

```python
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters

# Simple programmatic design
params = DesignParameters(sirna_length=21, top_candidates=10)
designer = SiRNADesigner(params)

results = designer.design_from_fasta("examples/sample_transcripts.fasta")
for candidate in results:
    print(f"{candidate.id}: GC={candidate.gc_content:.1f}%")
```

### Comprehensive Example

```python
#!/usr/bin/env python3
"""Automated batch analysis with QC"""

import subprocess
import pandas as pd
from pathlib import Path

def run_batch(genes, output_dir, **params):
    """Run workflow on multiple genes with QC"""
    results = {}

    for gene in genes:
        gene_dir = output_dir / gene
        cmd = [
            "uv", "run", "sirnaforge", "workflow", gene,
            "--output-dir", str(gene_dir),
            "--top-n", str(params.get("top_n", 30)),
            "--gc-min", str(params.get("gc_min", 35)),
            "--gc-max", str(params.get("gc_max", 60)),
            "--verbose"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            pass_file = gene_dir / "sirna_design" / f"{gene}_pass.csv"
            if pass_file.exists():
                df = pd.read_csv(pass_file)
                results[gene] = {
                    "count": len(df),
                    "mean_gc": df["gc_content"].mean(),
                    "mean_score": df["composite_score"].mean()
                }
        else:
            results[gene] = {"error": result.stderr}

    return results

# Run analysis
genes = ["TP53", "BRCA1", "EGFR"]
results = run_batch(genes, Path("batch"), top_n=30, gc_min=35, gc_max=60)

# Summary
summary = pd.DataFrame.from_dict(results, orient="index")
summary.to_csv("batch_summary.csv")
print(summary)
```

See [Python API Tutorial](tutorials/python_api.md) for detailed examples.

## Specialized Scenarios

### Difficult Targets

```bash
# Relaxed parameters for low-GC or high-structure sequences
uv run sirnaforge design difficult.fasta \
  --gc-min 20 --gc-max 75 \
  --max-poly-runs 5 \
  --skip-structure \
  --top-n 40
```

### High-Throughput Screening

```bash
# Fast mode - skip expensive computations
uv run sirnaforge design large_set.fasta \
  --skip-structure \
  --top-n 10 \
  --output screening.csv
```

### Parameter Comparison

```bash
# Compare stringency levels
GENE="TP53"

uv run sirnaforge workflow $GENE \
  --output-dir "${GENE}_conservative" \
  --gc-min 40 --gc-max 55 --max-poly-runs 2 --top-n 20

uv run sirnaforge workflow $GENE \
  --output-dir "${GENE}_moderate" \
  --gc-min 35 --gc-max 60 --max-poly-runs 3 --top-n 30

uv run sirnaforge workflow $GENE \
  --output-dir "${GENE}_permissive" \
  --gc-min 25 --gc-max 70 --max-poly-runs 4 --top-n 50
```

## Troubleshooting

### Input Validation

```bash
# Check FASTA format
uv run sirnaforge validate custom.fasta

# Debug gene search
uv run sirnaforge search RARE_GENE --all --verbose

# Enable verbose logging
uv run sirnaforge workflow GENE --verbose 2>&1 | tee debug.log
```

### Common Issues

**Low passing candidates:**
```bash
# Relax filters
uv run sirnaforge design input.fasta \
  --gc-min 25 --gc-max 70 \
  --max-poly-runs 4 \
  --top-n 50
```

**Performance optimization:**
```bash
# Skip structure prediction
uv run sirnaforge design large.fasta --skip-structure --top-n 15
```

**Multi-species requires Docker:**
```bash
# Build container (one-time)
make docker

# Run with BWA-MEM2 enabled
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --genome-species "human,mouse"
```

## Best Practices

### Recommendations

1. **Start Simple:** Use `workflow` command with gene symbols
2. **Validate First:** Run `validate` on custom FASTA files
3. **Check Quality:** Review pass rates before selecting candidates
4. **Iterate:** Adjust parameters based on initial results
5. **Document:** Save commands/scripts for reproducibility

### Parameter Guidelines

| Use Case | GC Range | Poly-runs | Top-N | Notes |
|----------|----------|-----------|-------|-------|
| Research | 35-60% | 2 | 30-50 | Conservative, high quality |
| Screening | 40-60% | 2-3 | 20-30 | Balanced |
| Initial | 30-70% | 3-4 | 30-50 | Permissive |
| Difficult | 20-75% | 4-5 | 40-60 | Relaxed for low-GC |

### Quality Thresholds

Typical pass rates: 60-80%. If < 30%:
- Check input quality (`validate`)
- Relax GC range (±5-10%)
- Increase `--max-poly-runs`
- Review sequence complexity

---

## See Also

- [Getting Started](getting_started.md) - Installation and quick start
- [CLI Reference](cli_reference.md) - Complete parameter docs
- [Gene Search Guide](gene_search.md) - Database access
- [Python API Tutorial](tutorials/python_api.md) - Programmatic usage
- [Modification Guide](modification_integration_guide.md) - Chemical mods
