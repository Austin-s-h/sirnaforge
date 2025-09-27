# 🧬 siRNAforge Usage Examples

Comprehensive real-world examples for advanced siRNAforge usage patterns.

> **New to siRNAforge?** Start with [Getting Started](getting_started.md) for basic installation and first analysis.

## Production Workflows

### High-Quality Research Design
```bash
# Comprehensive analysis with strict quality control
uv run sirnaforge workflow TP53 \
  --output-dir tp53_publication \
  --top-n 50 \
  --gc-min 35 --gc-max 55 \
  --max-poly-runs 2 \
  --verbose
```

### Clinical Screening Pipeline
```bash
# Multi-gene analysis for therapeutic targets
genes=("TP53" "BRCA1" "EGFR" "KRAS" "PIK3CA")
for gene in "${genes[@]}"; do
    uv run sirnaforge workflow "$gene" \
        --output-dir "clinical_screen/$gene" \
        --top-n 30 \
        --gc-min 40 --gc-max 60 \
        --genome-species "human,mouse" \
        --verbose
done
```

## Advanced Gene Search Scenarios

### Multi-Database Search Strategy
```bash
# Comprehensive search across all databases
uv run sirnaforge search RARE_GENE --all --verbose

# Species-specific searches
uv run sirnaforge search BRCA1 --database ensembl --species homo_sapiens
uv run sirnaforge search BRCA1 --database ensembl --species mus_musculus

# Transcript type filtering for lncRNA analysis
uv run sirnaforge search HOTAIR \
  --types "lncRNA,antisense" \
  --database ensembl \
  --verbose
```

### Custom Transcript Handling
```bash
# Working with your own transcript sequences
# 1. Validate format first
uv run sirnaforge validate custom_transcripts.fasta

# 2. Design with custom sequences
uv run sirnaforge design custom_transcripts.fasta \
  --top-n 30 \
  --gc-min 35 --gc-max 55 \
  --output custom_sirnas.csv
```

## Specialized Design Workflows

### Publication-Quality Analysis
```bash
# Research-grade parameters with comprehensive validation
uv run sirnaforge workflow GENE_OF_INTEREST \
  --output-dir publication_analysis \
  --top-n 50 \
  --gc-min 35 --gc-max 55 \
  --max-poly-runs 2 \
  --length 21 \
  --genome-species "human,mouse,rat" \
  --verbose

# Generate summary report
echo "Analysis completed for publication dataset" > analysis_notes.txt
cat publication_analysis/logs/workflow_summary.json >> analysis_notes.txt
```

### High-Throughput Screening
```bash
# Fast screening with relaxed parameters
uv run sirnaforge design large_gene_set.fasta \
  --skip-structure \
  --top-n 10 \
  --gc-min 25 --gc-max 70 \
  --max-poly-runs 4 \
  --output screening_results.csv
```

### Species-Specific Off-Target Analysis
```bash
# Multi-species safety screening
uv run sirnaforge workflow THERAPEUTIC_TARGET \
  --output-dir safety_analysis \
  --genome-species "human,macaque,mouse" \
  --top-n 25 \
  --verbose
```

### Difficult Targets (Low GC, High Structure)
```bash
# Relaxed parameters for challenging sequences
uv run sirnaforge design difficult_targets.fasta \
  --gc-min 20 --gc-max 75 \
  --max-poly-runs 5 \
  --skip-structure \
  --top-n 40 \
  --verbose
```
## Batch Processing & Automation

### Batch Gene Analysis
```bash
#!/bin/bash
# Process multiple genes with consistent parameters

genes=("TP53" "BRCA1" "EGFR" "MYC" "KRAS" "PIK3CA")
output_base="batch_analysis_$(date +%Y%m%d)"

for gene in "${genes[@]}"; do
    echo "🧬 Processing $gene..."
    uv run sirnaforge workflow "$gene" \
        --output-dir "${output_base}/${gene}" \
        --top-n 25 \
        --gc-min 35 --gc-max 55 \
        --verbose
    echo "✅ $gene completed"
done

# Generate batch summary
echo "Batch analysis completed: $(date)" > "${output_base}/batch_summary.txt"
```

### Parameter Set Comparisons
```bash
# Compare different stringency levels
GENE="TP53"

# High stringency
uv run sirnaforge workflow $GENE \
  --output-dir "${GENE}_high_stringency" \
  --gc-min 40 --gc-max 50 \
  --max-poly-runs 2 \
  --top-n 20

# Medium stringency
uv run sirnaforge workflow $GENE \
  --output-dir "${GENE}_medium_stringency" \
  --gc-min 35 --gc-max 60 \
  --max-poly-runs 3 \
  --top-n 30

# Low stringency (more candidates)
uv run sirnaforge workflow $GENE \
  --output-dir "${GENE}_low_stringency" \
  --gc-min 25 --gc-max 70 \
  --max-poly-runs 4 \
  --top-n 50
```

### Automated Quality Assessment
```bash
#!/bin/bash
# Automated pipeline with quality checks

GENE=$1
OUTPUT_DIR="analysis_${GENE}_$(date +%Y%m%d)"

echo "🔍 Starting analysis for $GENE"

# Run workflow
uv run sirnaforge workflow $GENE \
  --output-dir $OUTPUT_DIR \
  --top-n 30 \
  --verbose

# Check results quality
PASS_COUNT=$(tail -n +2 "$OUTPUT_DIR/sirnaforge/${GENE}_pass.csv" | wc -l)
ALL_COUNT=$(tail -n +2 "$OUTPUT_DIR/sirnaforge/${GENE}_all.csv" | wc -l)

echo "📊 Quality Report:"
echo "  - Total candidates: $ALL_COUNT"
echo "  - Passing filters: $PASS_COUNT"
echo "  - Pass rate: $(( PASS_COUNT * 100 / ALL_COUNT ))%"

if [ $PASS_COUNT -lt 5 ]; then
    echo "⚠️  Warning: Low number of passing candidates"
    echo "💡 Consider relaxing parameters or checking input quality"
fi
```

## Integration & Automation

### Makefile Integration
```makefile
# Project Makefile with siRNAforge targets
GENE ?= TP53
OUTPUT_DIR = analysis_$(shell date +%Y%m%d)

.PHONY: sirna-design sirna-batch sirna-qc

sirna-design:
	@echo "🧬 Designing siRNAs for $(GENE)"
	uv run sirnaforge workflow $(GENE) \
		--output-dir $(OUTPUT_DIR)/$(GENE) \
		--top-n 25 \
		--verbose

sirna-batch:
	@for gene in TP53 BRCA1 EGFR MYC KRAS; do \
		echo "Processing $$gene..."; \
		$(MAKE) sirna-design GENE=$$gene; \
	done

sirna-qc:
	@echo "📊 Generating QC report"
	@find $(OUTPUT_DIR) -name "*_pass.csv" -exec wc -l {} + > qc_summary.txt
```

### Python API Integration
```python
#!/usr/bin/env python3
"""
Advanced siRNAforge automation with Python
"""
import subprocess
import json
import pandas as pd
from pathlib import Path

def run_sirnaforge_workflow(gene, output_dir, **kwargs):
    """Run siRNAforge workflow with custom parameters"""
    cmd = [
        "uv", "run", "sirnaforge", "workflow", gene,
        "--output-dir", str(output_dir)
    ]

    # Add optional parameters
    for key, value in kwargs.items():
        cmd.extend([f"--{key.replace('_', '-')}", str(value)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0

def analyze_results(output_dir, gene):
    """Analyze siRNAforge results"""
    pass_file = output_dir / "sirnaforge" / f"{gene}_pass.csv"
    all_file = output_dir / "sirnaforge" / f"{gene}_all.csv"

    if pass_file.exists() and all_file.exists():
        pass_df = pd.read_csv(pass_file)
        all_df = pd.read_csv(all_file)

        return {
            "gene": gene,
            "total_candidates": len(all_df),
            "passing_candidates": len(pass_df),
            "pass_rate": len(pass_df) / len(all_df) * 100,
            "mean_asymmetry": pass_df['asymmetry_score'].mean(),
            "mean_gc": pass_df['gc_content'].mean()
        }
    return None

# Example usage
genes = ["TP53", "BRCA1", "EGFR"]
results = []

for gene in genes:
    output_dir = Path(f"analysis_{gene}")

    # Run analysis
    success = run_sirnaforge_workflow(
        gene, output_dir,
        top_n=30,
        gc_min=35,
        gc_max=55,
        verbose=True
    )

    if success:
        analysis = analyze_results(output_dir, gene)
        if analysis:
            results.append(analysis)

# Generate summary report
summary_df = pd.DataFrame(results)
summary_df.to_csv("batch_analysis_summary.csv", index=False)
print("📊 Analysis completed. Summary saved to batch_analysis_summary.csv")
```

## Troubleshooting & Optimization

### Performance Optimization
```bash
# For large datasets - skip expensive computations
uv run sirnaforge design large_input.fasta \
  --skip-structure \
  --skip-off-targets \
  --top-n 10 \
  --output fast_results.csv

# Memory-efficient processing
uv run sirnaforge workflow GENE \
  --output-dir results \
  --top-n 15 \
  --max-poly-runs 3
```

### Debugging Failed Analyses
```bash
# Verbose mode for troubleshooting
uv run sirnaforge workflow PROBLEMATIC_GENE \
  --output-dir debug_output \
  --verbose 2>&1 | tee debug.log

# Check intermediate files
uv run sirnaforge search PROBLEMATIC_GENE \
  --output debug_transcripts.fasta \
  --verbose

# Validate input quality
uv run sirnaforge validate debug_transcripts.fasta
```

### Quality Control Scripts
```bash
#!/bin/bash
# qc_analysis.sh - Quality control for siRNAforge results

analysis_dir=$1
if [ -z "$analysis_dir" ]; then
    echo "Usage: $0 <analysis_directory>"
    exit 1
fi

echo "🔍 Quality Control Report for $analysis_dir"
echo "=============================================="

# Count total and passing candidates
for csv_file in "$analysis_dir"/sirnaforge/*_all.csv; do
    if [ -f "$csv_file" ]; then
        gene=$(basename "$csv_file" _all.csv)
        total=$(tail -n +2 "$csv_file" | wc -l)
        pass_file="$analysis_dir/sirnaforge/${gene}_pass.csv"

        if [ -f "$pass_file" ]; then
            passing=$(tail -n +2 "$pass_file" | wc -l)
            pass_rate=$(( passing * 100 / total ))
            echo "📊 $gene: $passing/$total candidates pass filters ($pass_rate%)"
        fi
    fi
done

# Check for common issues
echo ""
echo "🔧 Diagnostic Checks:"
if [ -f "$analysis_dir/logs/workflow_summary.json" ]; then
    echo "✅ Workflow completed successfully"
else
    echo "❌ No workflow summary found - check for errors"
fi

# Size checks
echo "📁 Output sizes:"
du -sh "$analysis_dir"/* 2>/dev/null | head -5
```

## Best Practices Summary

### Parameter Selection Guidelines
```bash
# Research/Publication quality
--gc-min 35 --gc-max 55 --max-poly-runs 2 --top-n 30

# High-throughput screening
--gc-min 25 --gc-max 70 --max-poly-runs 4 --top-n 15

# Difficult targets
--gc-min 20 --gc-max 75 --skip-structure --top-n 50

# Quick testing
--top-n 5 --skip-structure --skip-off-targets
```

### Workflow Recommendations
1. **Start simple**: Use `workflow` command for most use cases
2. **Validate inputs**: Always run `validate` on custom FASTA files
3. **Check results**: Review pass rates and quality metrics
4. **Iterate parameters**: Adjust stringency based on results
5. **Document settings**: Save configurations for reproducibility

> **More Resources:**
> - [CLI Reference](cli_reference.md) - Complete parameter documentation
> - [Custom Scoring Guide](tutorials/custom_scoring.md) - Advanced thermodynamic principles
> - [Development Guide](development.md) - Contributing and extending siRNAforge
