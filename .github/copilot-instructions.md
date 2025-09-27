# GitHub Copilot Instructions for siRNAforge

**Modern Python bioinformatics toolkit using uv package manager and make-based workflows.**

## Quick Reference

**Package Manager:** `uv` (fast Python package management)
**Python:** 3.9-3.12 with `src/` source layout
**Build Tool:** `make` commands for all workflows
**Testing:** Multiple test categories with pytest
**Packaging:** Modern `pyproject.toml` with hatchling

## Essential Commands

### Setup (REQUIRED - 60-120 seconds)
```bash
uv sync --dev  # Install all dependencies - NEVER CANCEL, timeout 180s
```

### Development Cycle (FAST)
```bash
make lint                    # Code quality (3-5s)
make test-local-python      # Fastest tests (12-15s)
make test-unit              # Unit tests (30-35s)
make check                  # Lint + fast tests (35-40s)
```

### Package Commands
```bash
uv run sirnaforge --help    # Verify installation
uv run sirnaforge version   # Check version (alt to make build)
make clean                  # Clean artifacts before building
# Note: make build fails due to symlink issues
```

## File Structure

```
src/sirnaforge/          # Main package (modern src/ layout)
├── cli.py              # CLI entry point
├── core/               # Core algorithms
├── models/             # Data models
└── workflow.py         # Main workflows

tests/
├── unit/               # Fast Python tests (30 tests, 12s)
├── integration/        # Cross-component tests
└── pipeline/           # Nextflow pipeline tests

pyproject.toml          # Modern Python packaging config
Makefile               # All development workflows
uv.lock                # Dependency lock file
```

## Make Commands (Organized by Speed/Purpose)

### Fast Commands (< 5 seconds)
```bash
make lint                   # Code quality (ruff, mypy)
make format                 # Auto-format code
make version                # Show version
make clean                  # Clean artifacts
```

### Development Commands (10-40 seconds)
```bash
make test-local-python      # Fastest tests (12-15s)
make test-unit              # Unit tests (30-35s)
make check                  # Lint + fast tests (35-40s)
make docs                   # Build documentation (8-12s)
```

### Comprehensive Commands (60+ seconds)
```bash
make test                   # All tests (60s+, some may fail)
make docker                 # Build Docker image (10-20 min)
```

## uv Package Manager

**Why uv?** 10-100x faster than pip, unified dependency resolution, modern lockfile format.

### Key uv Commands
```bash
uv sync --dev               # Install dev dependencies (first run: 60-120s)
uv sync --no-dev           # Production dependencies only
uv add package_name        # Add new dependency
uv run command             # Run command in virtual environment
uv build                   # Build package (currently fails - use version check)
```

### Dependency Management
- `pyproject.toml` - Project config and dependencies
- `uv.lock` - Exact dependency versions (committed)
- `[dependency-groups]` - Separate dev dependencies in pyproject.toml

## Python Packaging (Modern Style)

### Key Files
- `pyproject.toml` - Single config file (replaces setup.py, setup.cfg, requirements.txt)
- `src/sirnaforge/` - Source layout (prevents import issues)
- `uv.lock` - Reproducible builds

### Package Structure
```
src/sirnaforge/
├── __init__.py          # Version and package info
├── cli.py              # Typer-based CLI entry point
├── core/               # Core siRNA algorithms
├── models/             # Pydantic data models
├── workflow.py         # Main workflow orchestration
└── utils/              # Helper functions
```

### Build System
- **Backend:** `hatchling` (modern, simple)
- **Build fails:** Nextflow symlinks break tar packaging
- **Alternative:** `uv run sirnaforge version` to verify functionality

## Testing Strategy (4 Categories)

### 1. Unit Tests (`make test-unit`)
- **Purpose:** Test individual functions/classes
- **Speed:** 30-35 seconds (31 tests)
- **Use for:** Core algorithm development

### 2. Local Python (`make test-local-python`)
- **Purpose:** Python-only integration tests
- **Speed:** 12-15 seconds (30 tests) - FASTEST
- **Use for:** Development iteration

### 3. Fast Tests (`make test-fast`)
- **Purpose:** All tests except slow external ones
- **Speed:** 25-30 seconds (53+ tests)
- **Use for:** Pre-commit validation

### 4. Full Tests (`make test`)
- **Purpose:** Complete test suite including Docker/Nextflow
- **Speed:** 60+ seconds, some failures expected
- **Use for:** Release validation

## External Dependencies

### ✅ Works in Base Environment
- Python siRNA algorithms
- CLI functionality
- Unit/integration tests
- Documentation builds
- Code quality tools

### ❌ Requires Docker
- Nextflow pipeline execution
- BWA-MEM2 alignment
- SAMtools processing
- ViennaRNA folding
- Full off-target analysis

### Docker Alternative
```bash
make docker                 # Build image (10-20 min)
make docker-run GENE=TP53   # Run workflow in container
```

## Best Practices for AI Coding Agents

### Development Workflow
1. **Setup:** `uv sync --dev` (one-time, 60-120s)
2. **Iterate:** `make test-local-python` (12s validation)
3. **Quality:** `make lint` (3s) before any commit
4. **Validate:** `make check` (35s) before push

### Timeout Guidelines
- `uv sync --dev`: 180+ seconds (NEVER cancel)
- `make test-*`: 60+ seconds
- `make lint`: 30+ seconds
- CLI commands: 30+ seconds

### Manual Verification After Changes
```bash
uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test.csv
# Expected: Success in 1-2 seconds, CSV output created
```

### Common Gotchas
- **Build failure:** Use `uv run sirnaforge version` instead of `make build`
- **Test failures:** 1 Docker integration failure is normal
- **SSL warnings:** Normal in docs builds, can ignore
- **Network timeouts:** Use `--input-fasta` with local files

### Quick Health Check
```bash
make lint && make test-local-python
# Should complete in ~15-20 seconds total
```
