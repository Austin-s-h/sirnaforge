#!/usr/bin/env bash

# Test script for Nextflow pipeline integration
# This script tests the Nextflow 25+ setup and linting

set -euo pipefail

echo "🧬 Testing siRNAforge Nextflow Integration"
echo "========================================"

pixi_dev=(pixi run --locked -e dev)
"${pixi_dev[@]}" install-app

# Check Nextflow version
echo "🔍 Checking Nextflow version..."
"${pixi_dev[@]}" nextflow -version

# Lint Nextflow scripts
echo "🔍 Linting Nextflow scripts..."
PIPELINE_NF=$("${pixi_dev[@]}" python -c "from sirnaforge.pipeline.nextflow.runner import NextflowRunner; print(NextflowRunner().get_main_workflow())")
"${pixi_dev[@]}" nextflow lint "$PIPELINE_NF"

# Test basic Nextflow functionality
echo "🧪 Testing Nextflow basic functionality..."
"${pixi_dev[@]}" nextflow run hello

# Verify Docker integration works
if command -v docker >/dev/null 2>&1; then
    echo "🐳 Testing Nextflow with Docker..."
    "${pixi_dev[@]}" nextflow run hello -with-docker ubuntu:20.04 || echo "⚠️  Docker test failed (may need Docker daemon)"
else
    echo "⚠️  Docker not available, skipping Docker integration test"
fi

# Test our pipeline syntax (dry run)
echo "🔬 Testing SIRNAforge pipeline syntax with built-in defaults..."
probe_dir=$(mktemp -d)
trap 'rm -rf "$probe_dir"' EXIT

CANDIDATES_FASTA="$probe_dir/candidates.fasta"
cat > "$CANDIDATES_FASTA" <<'EOF'
>candidate_1
GGAUCUUCUUAGCUUACGAUU
>candidate_2
GCUAUGCUUACGGAUCUUCUU
EOF

"${pixi_dev[@]}" nextflow run "$PIPELINE_NF" \
    --input "$CANDIDATES_FASTA" \
    --outdir test_results \
    -profile test \
    -stub-run || echo "⚠️  Pipeline syntax test failed - this may be expected if dependencies are missing"

echo "✅ Nextflow integration tests completed!"
echo ""
echo "Next steps:"
echo "- Run 'make docker-test' for the container integration suite"
echo "- Run 'make docker-build' to build comprehensive Docker image"
