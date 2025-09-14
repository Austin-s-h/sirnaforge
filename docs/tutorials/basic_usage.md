# Basic Usage Tutorial

This tutorial walks through your first siRNA design project using siRNAforge.

## Learning Objectives

By the end of this tutorial, you will:
- Understand the basic siRNAforge workflow
- Successfully design siRNAs for a gene of interest
- Interpret siRNA scoring and ranking
- Customize design parameters for your needs

## Prerequisites

- siRNAforge installed (see {doc}`../getting_started`)
- Basic familiarity with command-line interfaces
- Understanding of siRNA biology (helpful but not required)

## Installation and Setup

### Step 1: Verify Your Installation

First, let's make sure siRNAforge is working correctly:

```bash
# Check that siRNAforge is available
uv run sirnaforge --help

# Verify version and dependencies
uv run sirnaforge version
```

**Expected Output:**
```
╭──── Version Info ─────╮
│ 🧬 siRNAforge Toolkit │
│ Version: 0.1.0        │
│ Author: Austin S. H.  │
╰───────────────────────╯
```

If you see this output, you're ready to proceed!

## Basic Workflow

### Step 2: Your First siRNA Design

Let's design siRNAs for the well-known tumor suppressor gene TP53:

```bash
# Create a directory for our tutorial
mkdir sirnaforge_tutorial
cd sirnaforge_tutorial

# Run complete workflow for TP53
uv run sirnaforge workflow TP53 --output-dir tp53_analysis --verbose
```

**Expected Structure:**
```
tp53_analysis/
├── transcripts/
│   ├── TP53_transcripts.fasta      # All retrieved transcripts
│   ├── TP53_canonical.fasta        # Canonical isoform only
│   └── temp_for_design.fasta       # Transcripts used for design
├── sirnaforge/
│   ├── TP53_all.csv                # All designed candidates
│   ├── TP53_pass.csv               # Passing candidates only
│   └── manifest.json               # FAIR metadata
├── logs/
│   ├── workflow_stream.log         # Console output
│   └── workflow_summary.json       # High-level summary
└── orf_reports/
  └── TP53_orf_validation.txt     # ORF analysis results
```

### Step 3: Explore the Results

After the workflow completes, examine the output structure:

```bash
# Inspect top rows of passing candidates CSV (header + 5 rows)
head -6 tp53_analysis/sirnaforge/TP53_pass.csv
```

### Step 4: Examine Top Candidates

Look at the top-ranked siRNA candidates:

```bash
# View top 10 candidates
head -11 tp53_analysis/sirna_design/TP53_top_candidates.tsv
```

**Key Columns to Note:**
- `sirna_id`: Unique identifier

```bash
# Sort by composite score to see the best candidates
sort -k7 -nr tp53_analysis/sirna_design/TP53_top_candidates.tsv | head -5
```

### Step 5: Interpret the Scores

**Interpreting Scores:**
- **Composite Score (7-10)**: Higher is better; combines all factors
- **Thermodynamic Score (0-1)**: Measures RNA folding favorability and asymmetry (see {doc}`custom_scoring` for detailed explanation of thermodynamic asymmetry)
- **GC Content (35-60%)**: Optimal range 40-55% for balance of stability and accessibility
- **Asymmetry Score (0.65-1.0)**: Higher values indicate better guide strand selection into RISC
- **MFE (-2 to -8 kcal/mol)**: Moderate stability preferred for effective RISC processing

## Advanced Techniques

### Step 6: Customize Parameters

Now let's try a design with stricter quality parameters:

```bash
# High-quality design with strict GC content
uv run sirnaforge workflow TP53 \
  --output-dir tp53_strict \
  --gc-min 35 \
  --gc-max 50 \
  --top-n 30 \
  --verbose
```

**Compare Results:**
```bash
# Compare number of passing candidates
wcl_pass_a=$(($(wc -l < tp53_analysis/sirnaforge/TP53_pass.csv)-1))
wcl_pass_b=$(($(wc -l < tp53_strict/sirnaforge/TP53_pass.csv)-1))
echo "PASS count (analysis): $wcl_pass_a"
echo "PASS count (strict):   $wcl_pass_b"

# Compare average GC content in PASS sets (CSV; gc_content column)
awk -F',' 'NR>1 {sum+=$6; n++} END {if(n>0) printf("Average GC: %.2f%%\n", sum/n); else print "No rows"}' \
  tp53_analysis/sirnaforge/TP53_pass.csv

awk -F',' 'NR>1 {sum+=$6; n++} END {if(n>0) printf("Average GC: %.2f%%\n", sum/n); else print "No rows"}' \
  tp53_strict/sirnaforge/TP53_pass.csv
```

### Step 7: Step-by-Step Workflow

Let's break down the workflow into individual steps for better understanding:

```bash
# Create a new directory for step-by-step analysis
mkdir tp53_stepwise
cd tp53_stepwise

# Step 1: Search for transcripts
uv run sirnaforge search TP53 -o tp53_transcripts.fasta --verbose

# Step 2: Validate the sequences
uv run sirnaforge validate tp53_transcripts.fasta

# Step 3: Design siRNAs
uv run sirnaforge design tp53_transcripts.fasta \
  -o tp53_sirnas.tsv \
  --top-n 25 \
  --verbose
```

**Advantages of Step-by-Step:**
- Full control over each step
- Can modify intermediate files
- Easier to debug issues
- Can apply custom processing

### Step 8: Different Genes and Parameters

Try designing for different genes with various parameters:

```bash
# Cancer-related genes with different stringencies
uv run sirnaforge workflow BRCA1 --output-dir brca1_analysis --gc-min 40 --gc-max 55
uv run sirnaforge workflow EGFR --output-dir egfr_analysis --top-n 50 --gc-min 25
uv run sirnaforge workflow MYC --output-dir myc_analysis --length 21
```

### Step 9: Batch Processing

Process multiple genes systematically:

```bash
#!/bin/bash
# Create batch_analysis.sh

genes=("TP53" "BRCA1" "BRCA2" "EGFR" "MYC")

for gene in "${genes[@]}"; do
  echo "🧬 Processing $gene..."
  uv run sirnaforge workflow "$gene" \
    --output-dir "batch_${gene}" \
    --gc-min 35 \
    --gc-max 50 \
    --top-n 20 \
    --verbose
  echo "✅ $gene completed"
done

# Make it executable and run
chmod +x batch_analysis.sh
./batch_analysis.sh
```

## Practical Exercise

**Challenge**: Design siRNAs for your gene of interest

1. Choose a gene relevant to your research
2. Run the complete workflow with default parameters
3. Examine the results and identify the top 3 candidates
4. Rerun with parameters optimized for your application
5. Compare and select final candidates

**Questions to Consider:**
- How many candidates were generated?
- What's the GC content range of top candidates?
- Are there any obvious off-target concerns?
- How do composite scores compare across candidates?

## Common Issues and Solutions

### Issue: Gene Not Found
```bash
# Try different search strategies
uv run sirnaforge search MYGENE --all --verbose     # All databases
uv run sirnaforge search ENSG00000123456            # Direct Ensembl ID
```

### Issue: Few Candidates Generated
```bash
# Relax parameters
uv run sirnaforge workflow MYGENE \
  --gc-min 25 --gc-max 65 \
  --top-n 50
```

### Issue: Understanding Output
```bash
# Check configuration used
uv run sirnaforge config

# Validate input files
uv run sirnaforge validate input.fasta

# Get help for specific commands
uv run sirnaforge design --help
```

## What's Next?

Now that you understand basic usage:

1. **{doc}`advanced_workflows`** - Multi-species analysis and optimization
2. **{doc}`python_api`** - Automate workflows programmatically
3. **{doc}`pipeline_integration`** - Scale to large datasets
4. **{doc}`custom_scoring`** - Develop specialized scoring functions

## Key Takeaways

✅ **Complete Workflow**: `sirnaforge workflow` handles everything automatically
✅ **Quality Control**: Use GC content and composite scores to select candidates
✅ **Customization**: Adjust parameters based on your specific requirements
✅ **Validation**: Always validate inputs and examine outputs carefully
✅ **Batch Processing**: Automate analysis of multiple genes efficiently
