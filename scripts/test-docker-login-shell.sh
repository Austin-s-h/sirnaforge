#!/bin/bash
# Test script to verify Docker login shell PATH fix
# This script verifies that both login and non-login shells can access sirnaforge and nextflow

set -e

IMAGE="${1:-sirnaforge:latest}"

echo "Testing Docker image: $IMAGE"
echo "========================================"
echo ""

# Test 1: Non-login shell (should work before and after fix)
echo "Test 1: Non-login shell (bash -c)"
echo "-----------------------------------"
if docker run --rm "$IMAGE" /bin/bash -c 'command -v sirnaforge && command -v nextflow' > /dev/null 2>&1; then
    echo "✓ PASS: Non-login shell finds sirnaforge and nextflow"
else
    echo "✗ FAIL: Non-login shell cannot find sirnaforge or nextflow"
    exit 1
fi
echo ""

# Test 2: Login shell (the fix target)
echo "Test 2: Login shell (bash -lc)"
echo "-------------------------------"
if docker run --rm "$IMAGE" /bin/bash -lc 'command -v sirnaforge && command -v nextflow' > /dev/null 2>&1; then
    echo "✓ PASS: Login shell finds sirnaforge and nextflow"
else
    echo "✗ FAIL: Login shell cannot find sirnaforge or nextflow"
    docker run --rm "$IMAGE" /bin/bash -lc 'echo "PATH=$PATH"'
    exit 1
fi
echo ""

# Test 3: Verify PATH contains conda directories in login shell
echo "Test 3: Verify PATH in login shell"
echo "-----------------------------------"
PATH_OUTPUT=$(docker run --rm "$IMAGE" /bin/bash -lc 'echo "$PATH"')
if echo "$PATH_OUTPUT" | grep -q "/opt/conda/bin"; then
    echo "✓ PASS: Login shell PATH contains /opt/conda/bin"
    echo "   PATH: $PATH_OUTPUT"
else
    echo "✗ FAIL: Login shell PATH missing /opt/conda/bin"
    echo "   PATH: $PATH_OUTPUT"
    exit 1
fi
echo ""

echo "========================================"
echo "All tests passed! ✓"
echo "Both login and non-login shells work correctly."
