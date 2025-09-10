# GitHub Copilot Instructions for siRNAforge

**ALWAYS follow these instructions first.** Only fallback to additional search and context gathering if the information in these instructions is incomplete or found to be in error.

## Repository Overview

siRNAforge is a comprehensive siRNA design toolkit built with modern Python (3.9-3.12). It combines core siRNA design algorithms with off-target analysis using bioinformatics tools (BWA-MEM2, Nextflow, ViennaRNA). The repository uses uv for dependency management, has extensive test coverage with multiple test categories, and provides both CLI and Python API interfaces.

## Essential Setup Commands

### Installation (VALIDATED)
```bash
# Install uv package manager (if not available)
pip install uv

# Install development environment - NEVER CANCEL: takes 60-120 seconds
uv sync --dev
# Timeout: Set 180+ seconds for uv sync commands

# Verify installation
uv run sirnaforge --help
```

### Build Process (CRITICAL LIMITATIONS)
```bash
# Clean before building (REQUIRED due to cached temp files)
make clean

# Package build - FAILS due to Nextflow temp file symlinks
# DO NOT USE: make build  
# Known issue: Workflow tests create symlinks that break tar packaging
# Error: "symlink destination for [...] is outside of the target directory"

# Use alternative verification instead of building
uv run sirnaforge version  # Verifies package installation works
```

## Testing (EXHAUSTIVELY VALIDATED)

### Unit Tests (FASTEST - Use for Development)
```bash
# Unit tests only - NEVER CANCEL: takes 30-35 seconds
make test-unit
# Timeout: Set 60+ seconds
# Expected: 31 tests pass, ~30 seconds
```

### Local Python Development Tests
```bash
# Local Python tests - NEVER CANCEL: takes 12-15 seconds  
make test-local-python
# Timeout: Set 30+ seconds
# Expected: 30 tests pass, ~12 seconds
```

### Fast Tests (Excludes Slow/External Dependencies)
```bash
# Fast tests excluding slow ones - NEVER CANCEL: takes 25-30 seconds
make test-fast
# Timeout: Set 60+ seconds
# Expected: 53+ tests pass, 3+ skipped, possible 1 failure (pipeline integration)
```

### Full Test Suite
```bash
# All tests - NEVER CANCEL: takes 60+ seconds, some will fail without Docker/Nextflow
make test
# Timeout: Set 120+ seconds
# Expected: Most pass, some skip/fail due to missing external tools
```

## Code Quality (FAST)

### Linting and Formatting (VALIDATED)
```bash
# Format code - takes 3-5 seconds
make format

# Run all quality checks - NEVER CANCEL: takes 3-5 seconds
make lint
# Timeout: Set 30+ seconds
# Expected: "All checks passed!", "36 files already formatted", "Success: no issues found"

# Combined quality + fast tests - takes 30-40 seconds total
make check
```

## Documentation (VALIDATED)

### CLI Documentation Generation
```bash
# Generate CLI reference - NEVER CANCEL: takes 6-8 seconds
make docs-cli
# Timeout: Set 30+ seconds
# Creates: docs/CLI_REFERENCE.md
```

### Full Documentation Build
```bash
# Build Sphinx documentation - NEVER CANCEL: takes 8-12 seconds
make docs
# Timeout: Set 60+ seconds
# Expected: Some SSL warnings for external inventories (normal)
# Creates: docs/_build/html/
```

### Serve Documentation Locally
```bash
# Serve docs on http://localhost:8000
make docs-serve
```

## Application Testing (MANUALLY VALIDATED)

### CLI Functionality Verification
```bash
# Test basic CLI commands - takes 1-2 seconds each
uv run sirnaforge version
uv run sirnaforge config
uv run sirnaforge --help

# Test basic design functionality - NEVER CANCEL: takes 1-2 seconds
mkdir -p /tmp/test_output
uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test_output/results.csv --top-n 10
# Expected: Completes in ~1 second, generates CSV with siRNA candidates
```

### Manual Validation Scenarios (ALWAYS RUN AFTER CHANGES)
```bash
# Scenario 1: Basic siRNA design workflow
uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test_design.csv
# Expected: Success message, CSV file created with candidates, completes in ~1 second

# Scenario 2: Validate input files  
uv run sirnaforge validate examples/sample_transcripts.fasta
# Expected: Validation passes, shows sequence statistics, completes in <1 second

# Scenario 3: Complete workflow with local FASTA
mkdir -p /tmp/workflow_test
uv run sirnaforge workflow --input-fasta examples/sample_transcripts.fasta test_gene --output-dir /tmp/workflow_test
# Expected: Full workflow completes in 12-15 seconds, creates complete output structure

# Scenario 4: Test CLI help system
uv run sirnaforge workflow --help
# Expected: Detailed help output with all options

# Scenario 5: Network timeout test (will fail gracefully)
timeout 10 uv run sirnaforge search TP53 --output /tmp/test_search.fasta
# Expected: Times out after 10 seconds, this is normal behavior
```

## External Dependencies (IMPORTANT LIMITATIONS)

### What Works in Base Environment
- ✅ Python siRNA design algorithms
- ✅ Basic CLI functionality  
- ✅ Unit and integration tests (Python-only)
- ✅ Documentation generation
- ✅ Code quality tools

### What Requires Docker/Conda Environment
- ❌ Nextflow pipeline execution (`nextflow` command not found)
- ❌ BWA-MEM2 alignment (`bwa-mem2` command not found)
- ❌ SAMtools processing (`samtools` command not found)
- ❌ ViennaRNA folding (`RNAfold` command not found)
- ❌ Full off-target analysis pipeline

### Docker Usage (Alternative for Full Pipeline)
```bash
# Build Docker image - NEVER CANCEL: expected 10-20 minutes
make docker
# Timeout: Set 30+ minutes

# Run in Docker with full environment
docker run -v $(pwd):/workspace -w /workspace sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir results
```

## Critical Timing Information

| Command | Expected Time | Timeout Setting | Notes |
|---------|---------------|-----------------|-------|
| `uv sync --dev` | 60-120 seconds | 180+ seconds | NEVER CANCEL - downloads many packages |
| `make test-unit` | 30-35 seconds | 60+ seconds | 31 tests, fastest validation |
| `make test-local-python` | 12-15 seconds | 30+ seconds | Pure Python, no external deps |
| `make test-fast` | 25-30 seconds | 60+ seconds | May have 1 pipeline failure |
| `make lint` | 3-5 seconds | 30+ seconds | Fast quality check |
| `make docs` | 8-12 seconds | 60+ seconds | SSL warnings normal |
| `make docs-cli` | 6-8 seconds | 30+ seconds | CLI reference generation |
| `sirnaforge design` | 1-2 seconds | 30+ seconds | Basic functionality test |
| `sirnaforge validate` | <1 second | 30+ seconds | Input validation test |
| `sirnaforge workflow` (local) | 12-15 seconds | 60+ seconds | Complete workflow with local FASTA |

## Development Workflow (VALIDATED)

### Making Changes
```bash
# 1. Create feature branch
git checkout -b feature/my-feature

# 2. Make changes and validate immediately
make lint                    # 3-5 seconds
make test-local-python      # 12-15 seconds - fastest validation

# 3. If working on core functionality, run unit tests
make test-unit              # 30-35 seconds

# 4. Before committing, run quality checks
make check                  # 35-40 seconds total (lint + fast tests)
```

### CI Validation Commands (Match CI Pipeline)
```bash
# These match the CI pipeline and should pass before PR:
make lint                   # Code quality
make test-fast              # Most tests (excludes Docker/Nextflow)
make docs                   # Documentation build
```

## Common Issues and Solutions

### Common Issues and Solutions

### Build Failures
- **Issue**: `uv build` fails with "symlink destination [...] is outside of the target directory"
- **Cause**: Workflow tests create symlinks that break tar packaging 
- **Solution**: No current workaround for package building
- **Alternative**: Use `uv run sirnaforge version` to verify package functionality

### Test Failures
- **Issue**: Pipeline integration tests fail without Docker
- **Expected**: 1 failure in `test_offtarget_analysis_docker_integration` is normal
- **Solution**: Run `make test-local-python` for development validation

### Missing External Tools
- **Issue**: Commands like `nextflow`, `bwa-mem2` not found
- **Expected**: These require Docker/conda environment from `docker/environment-nextflow.yml`
- **Solution**: Use Docker image or focus on Python-only development

### SSL Certificate Warnings in Documentation
- **Issue**: SSL warnings when building docs
- **Expected**: External inventory access issues are normal in sandboxed environments
- **Solution**: Warnings can be ignored, docs build successfully

### Network Timeouts
- **Issue**: `sirnaforge search` commands hang or timeout
- **Expected**: External API access may be blocked in sandboxed environments
- **Solution**: Use `--input-fasta` with local files for testing workflows

## File Structure Reference

### Key Directories
- `src/sirnaforge/` - Main package source code
- `tests/unit/` - Fast Python-only tests (30 tests, ~12 seconds)
- `tests/integration/` - Cross-component tests (require Docker)
- `tests/pipeline/` - Nextflow pipeline tests (require Docker + Nextflow)
- `nextflow_pipeline/` - Nextflow DSL2 workflow definitions
- `docs/` - Sphinx documentation source
- `examples/` - Working examples and test data

### Important Files
- `pyproject.toml` - Python packaging and tool configuration
- `Makefile` - Development workflow automation (VALIDATED commands)
- `uv.lock` - Reproducible dependency resolution
- `docker/environment-nextflow.yml` - Conda environment with bioinformatics tools

## Repository-Specific Notes

- Uses modern `src/` layout with `uv` package manager
- Multiple test categories with different resource requirements
- External tools integration requires containerized environment
- Build process has known issues with Nextflow temp file caching
- Documentation generation works despite SSL warnings for external inventories

Always run `make test-local-python` after making changes for fastest validation (12 seconds).
Always run `make lint` before committing (3 seconds).
For comprehensive validation before PR: `make check` (35-40 seconds total).