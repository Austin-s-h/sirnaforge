# 🧪 siRNAforge Testing Guide

Complete guide for running tests correctly with the Make-based workflow system.

## Quick Reference Card

### ⚡ Development Testing (Python-only)
```bash
make test-local-python       # Fastest validation (3-4s, 30 tests)
make test-unit              # Unit tests (30-35s, 31 tests)
make lint                   # Code quality check (3-5s)
make check                  # Combined lint + fast tests
```

### 🐳 Container Testing (Full environment)
```bash
make docker-test-smoke      # CI/CD minimal (256MB RAM)
make docker-test-fast       # Development (2GB RAM)
make docker-test-full       # Comprehensive (8GB RAM)
```

## Testing Strategy Overview

siRNAforge uses a **tiered testing approach** optimized for different development phases and resource constraints.

### Test Categories

| Category | Purpose | Time | Tests | Resources |
|----------|---------|------|-------|-----------|
| `test-local-python` | Development iteration | 12-15s | 30 | Minimal |
| `test-unit` | Unit validation | 30-35s | 31 | Low |
| `test-fast` | Pre-commit checks | 25-30s | 53+ | Low |
| `test-ci` | CI/CD pipeline | 40-60s | All | Medium |
| `docker-test-*` | Container validation | Varies | All | High |

## Local Development Testing

### 1. Initial Setup (Required - Run Once)
```bash
# Install all development dependencies
make install-dev
# Expected: 60-120 seconds, creates uv.lock, installs pre-commit hooks
# ✅ Success indicator: "Development environment ready!"
```

### 2. Fast Iteration Cycle
```bash
# Fastest validation (recommended for active development)
make test-local-python
# Expected: 12-15 seconds, 30 tests
# ✅ Success: All tests pass, no Docker required

# Alternative: Unit tests only
make test-unit
# Expected: 30-35 seconds, 31 tests
```

### 3. Code Quality Checks
```bash
# Quick linting (fast)
make lint
# Expected: 3-5 seconds
# Tools: ruff check, ruff format --check, mypy

# Auto-fix linting issues
make format
make lint-fix
# Expected: 5-10 seconds, auto-fixes code style issues

# Combined quality + fast tests
make check
# Expected: 35-40 seconds, runs lint-fix + test-fast
```

### 4. Pre-Commit Validation
```bash
# Run all fast tests (excludes slow integration tests)
make test-fast
# Expected: 25-30 seconds, 53+ tests
# Excludes: Docker, external service tests

# Full local test suite (may include some failures)
make test
# Expected: 60+ seconds, includes all test categories
# Note: Some Docker integration tests may fail without Docker setup
```

## Docker Testing (Comprehensive Validation)

### Prerequisites
- Docker installed and running
- 4GB+ RAM available to Docker
- Built Docker image: `make docker`

### 1. Build Docker Image
```bash
make docker
# Expected: ~19 minutes first time, creates sirnaforge:0.1.2
# ✅ Success: "Docker image built: sirnaforge:0.1.2"
# Image size: ~2GB (includes all bioinformatics tools)
```

### 2. Tiered Docker Testing

#### Ultra-Fast Smoke Tests (CI/CD)
```bash
make docker-test-smoke
# Resources: 256MB RAM, 0.5 CPU
# Purpose: Minimal validation for CI pipelines
# Expected: <30 seconds
```

#### Fast Development Tests
```bash
make docker-test-fast
# Resources: 2GB RAM, 1 CPU
# Purpose: Quick validation without resource strain
# Expected: 1-2 minutes
# Runs: Fast tests only, minimal resource usage
```

#### Standard Development Tests
```bash
make docker-test
# Resources: 4GB RAM, 2 CPU
# Purpose: Standard development validation
# Expected: 2-5 minutes
# Includes: Most test categories with resource limits
```

#### Full Integration Tests
```bash
make docker-test-full
# Resources: 8GB RAM, 4 CPU
# Purpose: Comprehensive validation for releases
# Expected: 5-10 minutes
# Includes: All tests with high resource allocation
```

### 3. Manual Docker Verification

#### Basic Functionality
```bash
# Version check (✅ Verified working)
docker run --rm sirnaforge:0.1.2 sirnaforge version
# Expected output: Version info box with 0.1.2

# Help system
docker run --rm sirnaforge:0.1.2 sirnaforge --help
docker run --rm sirnaforge:0.1.2 sirnaforge design --help
```

#### Workflow Testing
```bash
# Test with sample data (✅ Verified working)
docker run --rm -v $(pwd)/examples:/data sirnaforge:0.1.2 \
  sirnaforge design /data/sample_transcripts.fasta \
  -o /tmp/results.tsv --top-n 5 --skip-structure --skip-off-targets

# Expected:
# - Configuration display
# - Processing 3 sequences → 2 candidates
# - Processing time: ~0.02s
# - Results summary table
# - "Results saved" message
```

#### Interactive Development
```bash
# Development shell
make docker-dev
# Opens interactive bash session in container
# All tools available: sirnaforge, nextflow, bwa-mem2, samtools, ViennaRNA
```

## CI/CD Integration

### GitHub Actions / GitLab CI
```bash
# Generate artifacts for CI systems
make test-ci
# Creates: pytest-report.xml, coverage.xml
# Runs: CI-optimized test suite with proper markers
```

### Resource-Constrained Environments
```bash
# For CI systems with limited resources
make docker-test-smoke
# Minimal resource usage, essential validation only

# For development CI with moderate resources
make docker-test-fast
# Balanced testing with reasonable resource requirements
```

## Troubleshooting

### Common Issues

#### 1. Setup Issues
```bash
# Issue: "uv sync" fails or takes too long
# Solution: Ensure network connectivity, increase timeout
# Command: uv sync --timeout 300

# Issue: Pre-commit hooks not installed
# Solution: Re-run setup
make install-dev
```

#### 2. Test Failures
```bash
# Issue: Docker tests fail with "No such option"
# Solution: Check CLI syntax with --help
docker run --rm sirnaforge:0.1.2 sirnaforge design --help

# Issue: Resource exhaustion in Docker tests
# Solution: Use appropriate test level for your system
make docker-test-fast  # Instead of docker-test-full

# Issue: Some integration tests fail
# Solution: This is normal without full Docker/Nextflow setup
# Use: make test-fast  # Excludes problematic integration tests
```

#### 3. Performance Issues
```bash
# Issue: Tests run slowly
# Solutions:
make test-local-python  # Fastest option (12-15s)
pytest -v -n auto       # Parallel execution if supported
pytest -k "not slow"    # Skip slow tests manually

# Issue: Docker build is slow
# Solutions:
# - Ensure 8GB+ RAM allocated to Docker
# - Use BuildKit: export DOCKER_BUILDKIT=1
# - Subsequent builds are faster (layer caching)
```

## Best Practices

### Development Workflow
1. **Setup once**: `make install-dev`
2. **Fast iteration**: `make test-local-python` after changes
3. **Quality check**: `make lint` before commits
4. **Pre-commit**: `make check` before pushing
5. **Validation**: `make docker-test-fast` before releases

### Resource Management
- **Local development**: Use `test-local-python` and `test-unit`
- **CI/CD**: Use `test-ci` with artifacts or `docker-test-smoke`
- **Release validation**: Use `docker-test-full` with adequate resources
- **Quick validation**: Use `make check` for lint + fast tests

### Timeouts and Expectations
- **Never cancel** `uv sync --dev` (can take 60-120s first time)
- **Docker builds** take ~19 minutes first time, much faster subsequently
- **Unit tests** should complete in 30-35 seconds
- **Fast tests** should complete in 12-15 seconds
- **CI tests** may take 60+ seconds but generate proper artifacts

## Quick Health Checks

```bash
# Local installation verification
uv run sirnaforge version
uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test.tsv

# Docker environment verification
docker run --rm sirnaforge:0.1.2 sirnaforge version
```

**📋 For Docker operations and deployment:** See the Docker documentation in the `docker/` directory

This guide focuses on testing workflows across development phases and CI/CD environments.
