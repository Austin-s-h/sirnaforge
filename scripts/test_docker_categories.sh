#!/usr/bin/env bash
set -euo pipefail

echo "🧪 Testing Docker Test Categories"
echo "================================"
echo ""

# Test 1: Smoke tests (should be very few, very fast)
echo "1. 🔥 Smoke Tests (Must always pass):"
echo "   These are the bare minimum tests that validate Docker image functionality"
make docker-test-smoke || {
    echo "❌ SMOKE TESTS FAILED - Docker image is not functional!"
    echo "   This should never happen and blocks all releases."
    exit 1
}
echo "   ✅ Smoke tests passed"
echo ""

# Test 2: Integration tests (can fail in pre-release)
echo "2. 🔧 Integration Tests (Can fail in pre-release):"
echo "   These test complex workflows that might have issues in development versions"
set +e  # Don't exit on failure
make docker-test-integration
integration_exit_code=$?
set -e

if [ $integration_exit_code -eq 0 ]; then
    echo "   ✅ Integration tests passed"
else
    echo "   ⚠️  Integration tests failed (exit code: $integration_exit_code)"
    echo "   This is acceptable for pre-releases but should be fixed before full release"
fi
echo ""

# Test 3: Show the difference
echo "3. 📊 Test Category Summary:"
echo "   - Smoke tests: Focus on basic Docker/CLI functionality"
echo "   - Integration tests: Focus on complex workflows and external dependencies"
echo ""

echo "🎯 Test categorization complete!"
echo "For CI/CD:"
echo "  - Pre-release: Smoke tests must pass, integration can fail"
echo "  - Full release: Both smoke and integration tests must pass"
