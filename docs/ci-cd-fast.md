# Fast CI/CD with Toy Data

This document describes the improved CI/CD GitHub Actions workflow using minimal toy data for quick validation.

## Overview

The new fast CI/CD workflow addresses the previous issue where Docker processes requested far too many resources in tests. The solution provides:

- **Ultra-fast execution**: < 15 minutes total runtime
- **Minimal resources**: 256MB memory, 0.5 CPU cores for Docker tests
- **Toy data**: Minimal FASTA files (< 500 bytes) for quick testing
- **Smoke tests**: Essential functionality validation without full integration

## Workflow Structure

### 1. Quick Validation (< 5 minutes)
- Basic Python setup and dependency installation
- Critical lint checks (syntax errors only)
- Unit tests with toy data using `ci` marker

### 2. Docker Smoke Tests (< 10 minutes)
- Minimal Docker image build
- Ultra-lightweight container tests (256MB memory)
- Basic CLI functionality validation
- Toy workflow execution

### 3. Quick Integration (< 5 minutes)
- Package build validation
- Quick test coverage analysis
- Basic functionality verification

## Test Profiles

### Smoke Profile
```python
TestProfile(
    name="smoke",
    description="Ultra-minimal smoke tests for CI/CD - toy data only",
    docker_cpus=0.5,
    docker_memory="256m", 
    docker_memory_swap="512m",
    timeout=60
)
```

### Usage

#### GitHub Actions
The workflow is triggered on:
- Push to `master` or `dev` branches
- Pull requests to `master`
- Manual dispatch

#### Local Testing
```bash
# Run smoke tests locally
make docker-test-smoke

# Run with minimal resources
docker run --rm \
    --cpus=0.5 \
    --memory=256m \
    --memory-swap=512m \
    -v $(pwd):/workspace -w /workspace sirnaforge:latest \
    pytest -m "smoke" --maxfail=1
```

## Toy Data Files

### Location
- `tests/unit/data/toy_transcripts.fasta` - Minimal transcript sequences
- `tests/unit/data/toy_candidates.fasta` - Minimal siRNA candidates

### Characteristics
- **Size**: < 500 bytes total
- **Sequences**: 20-70 nucleotides (optimal for siRNA design)
- **Bases**: Valid RNA bases only (A, U, G, C)
- **Purpose**: Fast validation, not biological accuracy

## Comparison with Full CI/CD

| Aspect | Fast CI/CD | Full CI/CD |
|--------|------------|------------|
| Runtime | < 15 min | 30-60 min |
| Memory | 256MB-4GB | 4GB-12GB |
| CPU | 0.5-2 cores | 2-8 cores |
| Tests | Smoke only | Full suite |
| Data | Toy data | Real datasets |
| Purpose | Quick feedback | Comprehensive validation |

## When to Use

### Fast CI/CD
- Pull request validation
- Development branch pushes
- Quick feedback loops
- Resource-constrained environments

### Full CI/CD 
- Release preparation
- Master branch validation
- Comprehensive testing
- Performance validation

## Integration with Existing Workflows

The fast CI/CD complements the existing comprehensive CI/CD workflow:

```yaml
# Fast workflow (.github/workflows/ci-fast.yml)
name: 🚀 Fast CI/CD with Toy Data
on:
  push: [dev, feature/*]
  pull_request: [master]

# Comprehensive workflow (.github/workflows/ci.yml)  
name: 🧬 siRNA Design Toolkit CI/CD
on:
  push: [master]
  release: [published]
```

## Troubleshooting

### Common Issues

1. **Memory allocation failures**
   - Reduce Docker memory limits further if needed
   - Use `--memory=128m` for ultra-constrained environments

2. **Timeout issues**
   - Increase timeout values in test profile
   - Optimize toy data to be even smaller

3. **CLI command failures**
   - Check that `sirnaforge version` works
   - Verify Docker container has all dependencies

### Performance Tuning

```bash
# Monitor resource usage
docker stats

# Profile test execution time
pytest --durations=10 -m "smoke"

# Validate toy data size
find tests/unit/data -name "toy_*" -exec wc -c {} \;
```