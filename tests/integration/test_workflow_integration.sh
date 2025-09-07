#!/bin/bash
"""
Test script for the complete siRNA workflow integration.

This script demonstrates:
1. Gene search and transcript retrieval
2. ORF validation 
3. siRNA design
4. Off-target analysis with Nextflow
5. Comprehensive reporting

Usage:
    ./test_workflow_integration.sh TP53
    ./test_workflow_integration.sh BRCA1 --genome-species human,mouse
"""

set -e  # Exit on any error

# Default parameters
GENE_QUERY="${1:-TP53}"
OUTPUT_DIR="test_workflow_$(date +%Y%m%d_%H%M%S)"
GENOME_SPECIES="${2:-human,rat}"

echo "🧬 Testing Complete siRNA Workflow Integration"
echo "=============================================="
echo "Gene Query: $GENE_QUERY"
echo "Output Directory: $OUTPUT_DIR"
echo "Genome Species: $GENOME_SPECIES"
echo ""

# Check if uv is available
if ! command -v uv &> /dev/null; then
    echo "❌ Error: uv is not installed. Please install uv first."
    exit 1
fi

# Check if nextflow is available
if ! command -v nextflow &> /dev/null; then
    echo "⚠️  Warning: nextflow is not installed. Off-target analysis will be skipped."
    SKIP_NEXTFLOW=true
else
    SKIP_NEXTFLOW=false
fi

echo "📋 Step 1: Running complete workflow..."
echo "---------------------------------------"

# Run the complete workflow
uv run sirna workflow \
    "$GENE_QUERY" \
    --output-dir "$OUTPUT_DIR" \
    --genome-species "$GENOME_SPECIES" \
    --top-n 20 \
    --offtarget-n 10 \
    --verbose

echo ""
echo "📊 Step 2: Analyzing results..."
echo "-------------------------------"

# Check if results were generated
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "❌ Error: Output directory was not created"
    exit 1
fi

# Count output files
TRANSCRIPT_FILES=$(find "$OUTPUT_DIR/transcripts" -name "*.fasta" 2>/dev/null | wc -l)
ORF_REPORTS=$(find "$OUTPUT_DIR/orf_reports" -name "*.txt" 2>/dev/null | wc -l)
SIRNA_RESULTS=$(find "$OUTPUT_DIR/sirna_design" -name "*.csv" 2>/dev/null | wc -l)
OFFTARGET_FILES=$(find "$OUTPUT_DIR/off_target" -name "*.tsv" 2>/dev/null | wc -l)

echo "Results Summary:"
echo "  📄 Transcript files: $TRANSCRIPT_FILES"
echo "  🔍 ORF reports: $ORF_REPORTS"
echo "  🎯 siRNA results: $SIRNA_RESULTS"  
echo "  🚫 Off-target files: $OFFTARGET_FILES"

# Show file sizes
echo ""
echo "File Details:"
if [ -f "$OUTPUT_DIR/workflow_summary.json" ]; then
    echo "  📊 Workflow summary: $(wc -c < "$OUTPUT_DIR/workflow_summary.json") bytes"
fi

if [ -d "$OUTPUT_DIR/transcripts" ]; then
    for file in "$OUTPUT_DIR/transcripts"/*.fasta; do
        if [ -f "$file" ]; then
            sequences=$(grep -c "^>" "$file" 2>/dev/null || echo "0")
            echo "  📄 $(basename "$file"): $sequences sequences"
        fi
    done
fi

if [ -d "$OUTPUT_DIR/sirna_design" ]; then
    for file in "$OUTPUT_DIR/sirna_design"/*.csv; do
        if [ -f "$file" ]; then
            candidates=$(tail -n +2 "$file" 2>/dev/null | wc -l || echo "0")
            echo "  🎯 $(basename "$file"): $candidates candidates"
        fi
    done
fi

echo ""
echo "📈 Step 3: Validating workflow integration..."
echo "--------------------------------------------"

# Test individual components
echo "Testing individual CLI commands:"

echo "  🔍 Gene search..."
SEARCH_OUTPUT="${OUTPUT_DIR}/test_search.fasta"
uv run sirna search "$GENE_QUERY" -o "$SEARCH_OUTPUT" --canonical-only

if [ -f "$SEARCH_OUTPUT" ]; then
    SEARCH_SEQS=$(grep -c "^>" "$SEARCH_OUTPUT" 2>/dev/null || echo "0")
    echo "    ✅ Found $SEARCH_SEQS canonical sequences"
else
    echo "    ❌ Search failed"
fi

echo "  🎯 siRNA design..."
if [ -f "$SEARCH_OUTPUT" ] && [ "$SEARCH_SEQS" -gt 0 ]; then
    DESIGN_OUTPUT="${OUTPUT_DIR}/test_design.csv"
    uv run sirna design "$SEARCH_OUTPUT" -o "$DESIGN_OUTPUT" --top-n 10
    
    if [ -f "$DESIGN_OUTPUT" ]; then
        DESIGN_CANDIDATES=$(tail -n +2 "$DESIGN_OUTPUT" 2>/dev/null | wc -l || echo "0")
        echo "    ✅ Generated $DESIGN_CANDIDATES candidates"
    else
        echo "    ❌ Design failed"
    fi
else
    echo "    ⏭️  Skipped (no input sequences)"
fi

echo "  📋 FASTA validation..."
if [ -f "$SEARCH_OUTPUT" ]; then
    uv run sirna validate "$SEARCH_OUTPUT" > /dev/null 2>&1
    echo "    ✅ Validation completed"
else
    echo "    ⏭️  Skipped (no FASTA file)"
fi

# Test Nextflow integration if available
if [ "$SKIP_NEXTFLOW" = false ] && [ -f "$OUTPUT_DIR/sirna_design/${GENE_QUERY}_top_candidates.fasta" ]; then
    echo ""
    echo "🚀 Step 4: Testing Nextflow integration..."
    echo "-----------------------------------------"
    
    CANDIDATES_FASTA="$OUTPUT_DIR/sirna_design/${GENE_QUERY}_top_candidates.fasta"
    NEXTFLOW_OUTPUT="$OUTPUT_DIR/nextflow_test"
    
    # Check if nextflow pipeline exists
    PIPELINE_PATH="$(dirname "$0")/main.nf"
    if [ ! -f "$PIPELINE_PATH" ]; then
        PIPELINE_PATH="../nextflow_pipeline/main.nf"
    fi
    
    if [ -f "$PIPELINE_PATH" ]; then
        echo "  🧪 Running Nextflow pipeline..."
        
        # Run with test profile for quick execution
        nextflow run "$PIPELINE_PATH" \
            --input "$CANDIDATES_FASTA" \
            --outdir "$NEXTFLOW_OUTPUT" \
            --genome_species "human" \
            --max_hits 100 \
            -profile test \
            -resume || echo "    ⚠️  Nextflow execution failed (indices may not be available)"
        
        if [ -d "$NEXTFLOW_OUTPUT" ]; then
            NF_FILES=$(find "$NEXTFLOW_OUTPUT" -name "*.tsv" -o -name "*.json" 2>/dev/null | wc -l)
            echo "    ✅ Nextflow generated $NF_FILES result files"
        fi
    else
        echo "    ❌ Nextflow pipeline not found at $PIPELINE_PATH"
    fi
else
    echo ""
    echo "⏭️  Step 4: Nextflow integration skipped"
    if [ "$SKIP_NEXTFLOW" = true ]; then
        echo "    Reason: Nextflow not installed"
    else
        echo "    Reason: No candidate sequences available"
    fi
fi

echo ""
echo "🎉 Integration Test Complete!"
echo "============================"
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Recommended next steps:"
echo "1. Review the ORF validation report in orf_reports/"
echo "2. Examine top siRNA candidates in sirna_design/"
echo "3. Check off-target analysis results in off_target/"
echo "4. Open workflow_summary.json for complete metrics"

if [ -f "$OUTPUT_DIR/off_target/offtarget_report.html" ]; then
    echo "5. View the interactive HTML report: $OUTPUT_DIR/off_target/offtarget_report.html"
fi

echo ""
echo "To run the workflow again with different parameters:"
echo "uv run sirna workflow $GENE_QUERY --output-dir my_analysis --genome-species human,mouse,rat"
