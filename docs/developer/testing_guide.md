# Testing Guide

Tiered testing approach for different development phases and resources.

## Quick Commands

### Development (Python-only)
```bash
make test-dev           # Fast unit tests (~15s) - development iteration
make test-unit          # Unit tests only (~30s)
make lint               # Code quality checks (~5s)
make check              # format + test-dev (~40s)
```

### Docker (Full environment)
```bash
make docker-test        # Tests INSIDE container
make docker-build       # Build Docker image
make docker-shell       # Interactive debugging
```

### Test Categories by Tier

| Category | Purpose | Time | Tests | Resources |
|----------|---------|------|-------|-----------|
| `test-dev` | Development iteration | ~15s | 30 | Minimal |
| `test-ci` | CI/CD validation | ~40s | 50+ | Low |
| `test-release` | Full validation | ~60s | All | Medium |
| `test-unit` | Unit tests only | ~30s | 31 | Low |
| `test-integration` | Integration tests | ~30s | 20+ | Low |
| `docker-test` | Container tests | ~60s | All | High |

## Local Development Testing

### 1. Initial Setup (Required - Run Once)
```bash
# Install all development dependencies
make dev
# Expected: 60-120 seconds, installs deps + pre-commit hooks
# ✅ Success indicator: "Ready for development!"
```

### 2. Fast Iteration Cycle
```bash
# Fastest validation (recommended for active development)
make test-dev
# Expected: ~15 seconds, 30 tests
# ✅ Success: All tests pass, no Docker required

# Alternative: Unit tests only
make test-unit
# Expected: ~30 seconds, 31 tests
```

### 3. Code Quality Checks
```bash
# Quick linting (fast)
make lint
# Expected: ~5 seconds
# Tools: ruff check, ruff format --check, mypy

# Auto-fix linting issues
make format
# Expected: ~5-10 seconds, auto-fixes code style issues

# Combined quality + fast tests
make check
# Expected: ~40 seconds, runs format + lint + test-dev
```

### 4. Pre-Commit Validation
```bash
# Run CI-tier tests (quick smoke tests for CI/CD)
make test-ci
# Expected: ~40 seconds
# Includes smoke tests with coverage reports

# Full release validation
make test-release
# Expected: ~60 seconds, includes all tests with coverage
# Note: Some tests may require Docker or network access

# Full local test suite (all tests, may have skips/failures)
make test
# Expected: 60+ seconds, includes all test categories
# Note: Some Docker integration tests may skip without Docker setup
```

## Docker Testing (Comprehensive Validation)

### Prerequisites
- Docker installed and running
- 4GB+ RAM available to Docker
- Image built with: `make docker-build`

### 1. Build Docker Image
```bash
make docker-build
# Expected: ~15-20 minutes first time, creates sirnaforge:latest
# ✅ Success: "Docker image: sirnaforge:latest"
# Image size: ~2.5GB (includes all bioinformatics tools)
```

### 2. Run Tests in Container
```bash
# Run tests INSIDE Docker container (validates image setup)
make docker-test
# Expected: ~60 seconds
# Tests all container-based functionality
# ✅ Success: All tests pass, verifying Docker environment

# Enter interactive shell for debugging
make docker-shell
# Expected: Interactive bash prompt inside container
# Useful for: Debugging, manual testing, tool validation
```

### 3. Manual Docker Verification

#### Basic Functionality
```bash
# Version check
docker run --rm sirnaforge:latest sirnaforge version
# Expected output: Version information

# Help system
docker run --rm sirnaforge:latest sirnaforge --help
docker run --rm sirnaforge:latest sirnaforge design --help
```

#### Workflow Testing
```bash
# Test with sample data
docker run --rm -v $(pwd)/examples:/data sirnaforge:latest \
  sirnaforge design /data/sample_transcripts.fasta \
  -o /tmp/results.csv --top-n 5

# Expected: Results file created with siRNA candidates
```

## Best Practices

### Development Workflow
1. **Setup once**: `make dev`
2. **Fast iteration**: `make test-dev` after changes
3. **Quality check**: `make lint` before commits
4. **Pre-commit**: `make check` before pushing
5. **Validation**: `make test-release` before releases

### Resource Management
- **Local development**: Use `test-dev` and `test-unit`
- **CI/CD**: Use `test-ci` with artifacts
- **Release validation**: Use `test-release` with full coverage
- **Quick validation**: Use `make check` for lint + fast tests

### Timeouts and Expectations
- **Never cancel** `uv sync --dev` (can take 60-120s first time)
- **Docker builds** take ~15-20 minutes first time, much faster subsequently
- **Unit tests** should complete in ~30 seconds
- **Fast tests** should complete in ~15 seconds
- **CI tests** may take 40+ seconds but generate proper artifacts

## Quick Health Checks

```bash
# Local installation verification
uv run sirnaforge version
uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test.csv

# Docker environment verification
docker run --rm sirnaforge:latest sirnaforge version
```

**📋 For Docker operations and deployment:** See the Docker documentation in the `docker/` directory

This guide focuses on testing workflows across development phases and CI/CD environments.
