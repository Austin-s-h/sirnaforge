# 🧬 siRNAforge Usage Examples

Comprehensive examples demonstrating real-world usage patterns for siRNAforge — Comprehensive siRNA design toolkit for gene silencing.

## Quick Start Examples

### Complete Workflow (Recommended)
```bash
# End-to-end analysis for a single gene
uv run sirnaforge workflow TP53 --output-dir tp53_analysis

# With custom parameters
uv run sirnaforge workflow BRCA1 \
  --output-dir brca1_results \
  --top-n 30 \
  --gc-min 35 --gc-max 55 \
  --verbose
```

### Step-by-Step Workflow
```bash
# 1. Search for gene transcripts
uv run sirnaforge search TP53 -o tp53_transcripts.fasta

# 2. Validate the input
uv run sirnaforge validate tp53_transcripts.fasta

# 3. Design siRNA candidates
uv run sirnaforge design tp53_transcripts.fasta \
  -o tp53_sirnas.tsv \
  -n 20
```

## Gene Search Examples

### Basic Searches
```bash
# Simple gene search
uv run sirnaforge search TP53

# Search with custom output
uv run sirnaforge search BRCA1 -o brca1_sequences.fasta

# Verbose output for debugging
uv run sirnaforge search EGFR --verbose
```

### Advanced Searches
```bash
# Search multiple databases
uv run sirnaforge search MYC --all --verbose

# Only canonical isoforms
uv run sirnaforge search KRAS --canonical-only

# Specific transcript types
uv run sirnaforge search PIK3CA \
  --types "protein_coding,lncRNA,miRNA" \
  --verbose

# Metadata only (no sequences)
uv run sirnaforge search AKT1 --no-sequence
```

### Database-Specific Searches
```bash
# Ensembl (default)
uv run sirnaforge search TP53 --database ensembl

# RefSeq
uv run sirnaforge search BRCA1 --database refseq

# GENCODE
uv run sirnaforge search EGFR --database gencode
```

## siRNA Design Examples

### Basic Design
```bash
# Design from FASTA file
uv run sirnaforge design transcripts.fasta

# Custom number of candidates
uv run sirnaforge design transcripts.fasta -n 50

# Custom output file
uv run sirnaforge design transcripts.fasta -o my_sirnas.tsv
```

### High-Quality Design
```bash
# Strict parameters for high-quality siRNAs
uv run sirnaforge design transcripts.fasta \
  --gc-min 35 --gc-max 50 \
  --max-poly-runs 2 \
  --top-n 30 \
  --verbose
```

### Fast Design (Large Datasets)
```bash
# Skip computationally expensive steps
uv run sirnaforge design large_transcripts.fasta \
  --skip-structure \
  --skip-off-targets \
  --top-n 10
```

### Design with Off-Target Analysis
```bash
# Include genome index for off-target screening
uv run sirnaforge design transcripts.fasta \
  --genome-index /path/to/human_genome.idx \
  --top-n 20 \
  --verbose
```

### Avoiding SNPs
```bash
# Exclude regions with common SNPs
uv run sirnaforge design transcripts.fasta \
  --snp-file common_variants.vcf \
  --verbose
```

## Advanced Workflow Examples

### Multi-Gene Analysis
```bash
#!/bin/bash
# Batch process multiple genes

genes=("TP53" "BRCA1" "EGFR" "MYC" "KRAS")

for gene in "${genes[@]}"; do
    echo "🧬 Processing $gene..."
    uv run sirnaforge workflow "$gene" \
      --output-dir "analysis_${gene}" \
      --top-n 25 \
      --verbose
    echo "✅ $gene complete"
done
```

### High-Throughput Analysis
```bash
# Process many genes from a list
while IFS= read -r gene; do
    echo "Processing $gene..."
    uv run sirnaforge workflow "$gene" \
      --output-dir "batch_analysis/${gene}" \
      --top-n 15 \
      --gc-min 30 --gc-max 60
done < gene_list.txt
```

### Custom Parameter Sets
```bash
# Define parameter sets for different use cases

# Strict parameters
STRICT="--gc-min 35 --gc-max 50 --max-poly-runs 2"

# Relaxed parameters
RELAXED="--gc-min 25 --gc-max 65 --max-poly-runs 4"

# Fast parameters
FAST="--skip-structure --skip-off-targets --top-n 5"

# Use parameter sets
uv run sirnaforge design input.fasta $STRICT -o strict_results.tsv
uv run sirnaforge design input.fasta $RELAXED -o relaxed_results.tsv
uv run sirnaforge design input.fasta $FAST -o fast_results.tsv
```

## Multi-Species Off-Target Analysis

### Basic Multi-Species
```bash
# Analyze against human, rat, and rhesus genomes
uv run sirnaforge workflow TP53 \
  --genome-species "human,rat,rhesus" \
    --top-n 15
```

### Extended Species Panel
```bash
# Include additional model organisms
uv run sirnaforge workflow BRCA1 \
  --genome-species "human,mouse,rat,rhesus,dog" \
  --top-n 20 \
  --verbose
```

### Species-Specific Analysis
```bash
# Focus on specific species
uv run sirnaforge workflow EGFR \
  --genome-species "human,mouse" \
  --top-n 25
```

## File Validation Examples

### Basic Validation
```bash
# Check a single file
uv run sirnaforge validate input.fasta

# Validate before processing
uv run sirnaforge validate transcripts.fasta && \
uv run sirnaforge design transcripts.fasta
```

### Batch Validation
```bash
# Validate multiple files
for file in *.fasta; do
    echo "Validating $file..."
    if uv run sirnaforge validate "$file"; then
        echo "✅ $file is valid"
    else
        echo "❌ $file has issues"
    fi
done
```

## Configuration and Version Examples

### Check Configuration
```bash
# View current settings
uv run sirnaforge config

# Save configuration for reference
uv run sirnaforge config > my_sirnaforge_config.txt

# Check version
uv run sirnaforge version
```

## Integration Examples

### Makefile Integration
```makefile
# Add to your Makefile
GENE ?= TP53
OUTPUT_DIR = analysis_$(shell date +%Y%m%d)

sirna-analysis:
	uv run sirnaforge workflow $(GENE) \
		--output-dir $(OUTPUT_DIR) \
		--top-n 20 \
		--verbose

sirna-batch:
	@for gene in TP53 BRCA1 EGFR; do \
		echo "Processing $$gene..."; \
		uv run sirnaforge workflow "$$gene" \
			--output-dir "$(OUTPUT_DIR)/$$gene"; \
	done
```

### Shell Script Integration
```bash
#!/bin/bash
# comprehensive_analysis.sh

set -euo pipefail

GENE=${1:-"TP53"}
OUTPUT_BASE="sirna_analysis_$(date +%Y%m%d)"
OUTPUT_DIR="$OUTPUT_BASE/$GENE"

echo "🧬 siRNAforge Analysis Pipeline"
echo "Gene: $GENE"
echo "Output: $OUTPUT_DIR"
echo "=========================="

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run analysis
echo "1. Running complete workflow..."
uv run sirnaforge workflow "$GENE" \
    --output-dir "$OUTPUT_DIR" \
    --top-n 30 \
    --verbose

# Generate summary
echo "2. Generating summary..."
echo "Gene: $GENE" > "$OUTPUT_DIR/summary.txt"
echo "Date: $(date)" >> "$OUTPUT_DIR/summary.txt"
echo "siRNAforge version: $(uv run sirnaforge version)" >> "$OUTPUT_DIR/summary.txt"

echo "✅ Analysis complete!"
echo "Results available in: $OUTPUT_DIR"
```

### Python Integration
```python
#!/usr/bin/env python3
"""
Python wrapper for siRNAforge workflows
"""

import subprocess
import sys
from pathlib import Path

def run_sirnaforge_workflow(gene: str, output_dir: str = None, **kwargs):
    """Run siRNAforge workflow from Python"""

    cmd = ["uv", "run", "sirnaforge", "workflow", gene]

    if output_dir:
        cmd.extend(["--output-dir", output_dir])

    # Add additional parameters
    for key, value in kwargs.items():
        cmd.extend([f"--{key.replace('_', '-')}", str(value)])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✅ Workflow completed for {gene}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"❌ Workflow failed for {gene}: {e.stderr}")
        return None

# Usage example
if __name__ == "__main__":
    genes = ["TP53", "BRCA1", "EGFR"]

    for gene in genes:
        output = run_sirnaforge_workflow(
            gene=gene,
            output_dir=f"python_analysis_{gene}",
            top_n=25,
            gc_min=35,
            gc_max=55,
            verbose=True
        )

        if output:
            print(f"Analysis for {gene} completed successfully")
```

## Performance Optimization Examples

### Large Dataset Processing
```bash
# Split large FASTA files for parallel processing
split -l 2000 large_transcripts.fasta batch_

# Process batches in parallel
for batch in batch_*; do
    (
        echo "Processing $batch..."
        uv run sirnaforge design "$batch" \
          --skip-structure \
          --top-n 5 \
          -o "${batch}_results.tsv"
    ) &
done
wait

# Combine results
cat batch_*_results.tsv > combined_results.tsv
```

### Memory-Efficient Processing
```bash
# For very large datasets, process one at a time
find . -name "*.fasta" -type f | while read -r file; do
    echo "Processing $file..."
    uv run sirnaforge design "$file" \
      --skip-structure \
      --skip-off-targets \
      --top-n 3 \
      -o "${file%.fasta}_results.tsv"
done
```

## Error Handling Examples

### Robust Processing
```bash
#!/bin/bash
# Robust processing with error handling

process_gene() {
    local gene=$1
    local max_retries=3
    local attempt=1

    while [ $attempt -le $max_retries ]; do
        echo "Attempt $attempt for $gene..."

        if uv run sirnaforge workflow "$gene" \
            --output-dir "analysis_$gene" \
            --verbose 2>"error_$gene.log"; then
            echo "✅ $gene completed successfully"
            return 0
        else
            echo "❌ $gene failed (attempt $attempt)"
            ((attempt++))
            sleep 5
        fi
    done

    echo "❌ $gene failed after $max_retries attempts"
    return 1
}

# Process genes with error handling
genes=("TP53" "BRCA1" "EGFR")
for gene in "${genes[@]}"; do
    process_gene "$gene"
done
```

### Validation Pipeline
```bash
#!/bin/bash
# Complete validation pipeline

validate_and_process() {
    local input_file=$1
    local output_prefix=${input_file%.fasta}

    echo "🔍 Validating $input_file..."
    if ! uv run sirnaforge validate "$input_file"; then
        echo "❌ Validation failed for $input_file"
        return 1
    fi

    echo "🎯 Processing $input_file..."
    if uv run sirnaforge design "$input_file" \
        -o "${output_prefix}_sirnas.tsv" \
        --verbose; then
        echo "✅ Successfully processed $input_file"
        return 0
    else
        echo "❌ Processing failed for $input_file"
        return 1
    fi
}

# Process all FASTA files
for fasta_file in *.fasta; do
    validate_and_process "$fasta_file"
done
```

## Quality Control Examples

### Comprehensive QC Pipeline
```bash
#!/bin/bash
# Quality control pipeline

run_qc_pipeline() {
    local gene=$1
    local output_dir="qc_${gene}_$(date +%Y%m%d)"

    mkdir -p "$output_dir"

    echo "🧬 QC Pipeline for $gene"
    echo "======================="

    # 1. Search and validate
    echo "1. Searching for transcripts..."
    uv run sirnaforge search "$gene" \
        -o "$output_dir/${gene}_transcripts.fasta" \
        --verbose

    echo "2. Validating sequences..."
    uv run sirnaforge validate "$output_dir/${gene}_transcripts.fasta"

    # 2. Design with multiple parameter sets
    echo "3. Designing with strict parameters..."
    uv run sirnaforge design "$output_dir/${gene}_transcripts.fasta" \
        --gc-min 35 --gc-max 50 \
        --max-poly-runs 2 \
        --top-n 50 \
        -o "$output_dir/${gene}_strict.tsv"

    echo "4. Designing with relaxed parameters..."
    uv run sirnaforge design "$output_dir/${gene}_transcripts.fasta" \
        --gc-min 25 --gc-max 65 \
        --max-poly-runs 4 \
        --top-n 50 \
        -o "$output_dir/${gene}_relaxed.tsv"

    # 3. Generate summary
    echo "5. Generating QC summary..."
    {
        echo "QC Summary for $gene"
        echo "===================="
        echo "Date: $(date)"
        echo "siRNAforge version: $(uv run sirnaforge version)"
        echo ""
        echo "Files generated:"
        ls -la "$output_dir"
        echo ""
        echo "Strict parameters results:"
        wc -l "$output_dir/${gene}_strict.tsv"
        echo ""
        echo "Relaxed parameters results:"
        wc -l "$output_dir/${gene}_relaxed.tsv"
    } > "$output_dir/qc_summary.txt"

    echo "✅ QC complete for $gene. Results in $output_dir"
}

# Run QC for multiple genes
for gene in TP53 BRCA1 EGFR; do
    run_qc_pipeline "$gene"
done
```

---

## Tips and Best Practices

### 🎯 **Design Quality**
- Use `--gc-min 35 --gc-max 50` for high-quality siRNAs
- Set `--max-poly-runs 2` to avoid off-targets
- Include `--verbose` for detailed progress information

### 🚀 **Performance**
- Use `--skip-structure --skip-off-targets` for large datasets
- Process files in batches for very large analyses
- Consider `--top-n 5` for initial screening

### 🔍 **Troubleshooting**
- Always run `validate` before `design`
- Check `config` if results seem unexpected
- Use `--verbose` to diagnose issues
- Try different databases if search fails

### 📁 **File Management**
- Use descriptive output directory names with dates
- Keep parameter sets consistent across analyses
- Archive results with metadata files

---

*For more examples and use cases, see the [examples directory](examples/USAGE_EXAMPLES.md) and the [CLI Reference](CLI_REFERENCE.md).*
