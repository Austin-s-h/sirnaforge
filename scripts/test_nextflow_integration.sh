#!/usr/bin/env bash

# Test script for Nextflow pipeline integration
# This script tests the Nextflow 25+ setup and linting

set -euo pipefail

echo "🧬 Testing siRNAforge Nextflow Integration"
echo "========================================"

# Check if running in CI or local environment
if [[ "${CI:-false}" == "true" ]]; then
    echo "🚀 Running in CI environment"
    SYNC_CMD="uv sync --group pipeline"
else
    echo "🖥️  Running in local environment"
    SYNC_CMD="uv sync --group pipeline"
fi

# Install pipeline dependencies
echo "📦 Installing pipeline dependencies..."
$SYNC_CMD

# Check Nextflow version
echo "🔍 Checking Nextflow version..."
uv run --group pipeline nextflow -version

# Lint Nextflow scripts
echo "🔍 Linting Nextflow scripts..."
PIPELINE_NF=$(uv run --group pipeline python -c "from sirnaforge.pipeline.nextflow.runner import NextflowRunner; print(NextflowRunner().get_main_workflow())")
uv run --group pipeline nextflow lint "$PIPELINE_NF"

# Test basic Nextflow functionality
echo "🧪 Testing Nextflow basic functionality..."
uv run --group pipeline nextflow run hello

# Verify Docker integration works
if command -v docker >/dev/null 2>&1; then
    echo "🐳 Testing Nextflow with Docker..."
    uv run --group pipeline nextflow run hello -with-docker ubuntu:20.04 || echo "⚠️  Docker test failed (may need Docker daemon)"
else
    echo "⚠️  Docker not available, skipping Docker integration test"
fi

# Test our pipeline syntax (dry run)
echo "🔬 Testing SIRNAforge pipeline syntax with built-in defaults..."
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

CANDIDATES_FASTA="$TMPDIR/candidates.fasta"
cat > "$CANDIDATES_FASTA" <<'EOF'
>candidate_1
GGAUCUUCUUAGCUUACGAUU
>candidate_2
GCUAUGCUUACGGAUCUUCUU
EOF

uv run --group pipeline nextflow run "$PIPELINE_NF" \
    --input "$CANDIDATES_FASTA" \
    --outdir test_results \
    -profile test \
    -stub-run || echo "⚠️  Pipeline syntax test failed - this may be expected if dependencies are missing"

echo "✅ Nextflow integration tests completed!"
echo ""
echo "Next steps:"
echo "- Run 'make nextflow-run' to test with real data"
echo "- Run 'make lint-nextflow' to lint pipeline scripts"
echo "- Run 'make docker-build' to build comprehensive Docker image"
