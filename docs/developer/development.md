# Development Guide

This guide covers development setup, contribution guidelines, and best practices for siRNAforge development.

> **New to siRNAforge?** Start with [Getting Started](../getting_started.md) for basic installation and usage.

## Development Environment Setup

### Prerequisites

- Python 3.9-3.12
- Git
- `uv` package manager
- `make` (optional but recommended)

### Quick Start

```bash
# Clone repository
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge

# One-command setup (installs deps + pre-commit hooks)
make dev
```

This runs:
1. `uv sync` - Installs all dependencies with dev extras
2. `uv run pre-commit install` - Sets up git hooks for code quality

### Manual Setup

If you prefer not to use `make`:

```bash
# Install with development dependencies
uv sync

# Install pre-commit hooks
uv run pre-commit install

# Verify installation
uv run sirnaforge --help
```

## Make Commands Reference

siRNAforge uses a comprehensive `Makefile` for common development tasks. Run `make help` to see all available commands.

### Testing Commands

#### By Tier (Recommended)
```bash
make test-dev        # Fast unit tests (~15s) - development iteration
make test-ci         # Smoke tests for CI/CD with coverage
make test-release    # Comprehensive validation (all tests + coverage)
make test            # All tests (may have skips/failures)
```

These tiers match the pytest marker structure:
- `dev` tier: Fast unit tests for rapid iteration
- `ci` tier: Smoke tests for CI/CD validation
- `release` tier: Full integration and heavy tests

#### By Requirements
```bash
make test-requires-docker    # Tests requiring Docker daemon
make test-requires-network   # Tests requiring network access
make test-requires-nextflow  # Tests requiring Nextflow
```

### Code Quality Commands

```bash
make lint       # Check code quality (ruff + mypy)
make format     # Auto-format code (ruff format + fixes)
make check      # Quick validation: format + test-dev
make pre-commit # Run pre-commit hooks on all files
make security   # Run security scans (bandit + safety)
```

### Docker Commands

```bash
make docker-build    # Build Docker image
make docker-test     # Run tests INSIDE container
make docker-shell    # Interactive shell in container
make docker-run      # Run workflow in Docker (GENE=TP53)
make docker-ensure   # Ensure image exists (build if needed)
```

**Example:**
```bash
# Run workflow for specific gene in Docker
make docker-run GENE=BRCA1
```

### Documentation Commands

```bash
make docs        # Build HTML documentation
make docs-serve  # Serve docs at localhost:8000
```

### Utility Commands

```bash
make clean      # Clean build and cache artifacts
make version    # Show current version
make build      # Build package distribution
make example    # Run basic example workflow
make help       # Show all available commands
```

## Development Workflow

## Branching & Release Strategy

siRNAforge follows a simple two-branch model that keeps `master` stable while enabling fast iteration on `dev`.

- `dev`: active integration branch (unstable, where PRs land)
- `master`: stable production branch (release-ready only)

This structure helps prevent accidental “production releases” from everyday feature merges while still allowing frequent Docker image publishing for preview/testing.

```{raw} html
<div class="mermaid">
flowchart LR
    subgraph Work[Daily Development]
        F1[feature/*] --> PR1[PR to dev]
        PR1 --> D[dev]
        D --> CI1[CI checks]
    end

    subgraph Release[Stable Release]
        D --> PR2[Release PR: dev → master]
        PR2 --> M[master]
        M --> CI2[CI + Release workflow]
        CI2 --> DOCKER_STABLE[Publish Docker (version + latest)]
    end

    subgraph Pre[Pre-release / Preview]
        D --> CI3[CI + Dev publish]
        CI3 --> DOCKER_DEV[Publish Docker (dev tag)]
    end
</div>
```

### How to Work Day-to-Day

1. Branch from `dev`: `feature/<name>` or `fix/<name>`
2. Open PR into `dev` (CI must pass)
3. Cut a stable release by opening a PR from `dev` into `master`

### Recommended GitHub Branch Protection

For `master` (stable):
- Require pull requests (no direct pushes)
- Require status checks (CI must pass)
- Require at least one review
- Block force pushes

For `dev` (integration):
- Require pull requests
- Require status checks

### Versioning Guidance (PEP 440)

- Stable releases: `X.Y.Z` (published to PyPI)
- Pre-releases: `X.Y.ZrcN` / `X.Y.ZaN` / `X.Y.ZbN` (published to TestPyPI by default)

### Docker Image Tags (Recommended)

- `ghcr.io/<owner>/<repo>:dev`: latest image from the `dev` branch (unstable)
- `ghcr.io/<owner>/<repo>:sha-<shortsha>`: immutable image for a specific commit (useful for debugging)
- `ghcr.io/<owner>/<repo>:<version>` and `:latest`: stable release images published from `master`

### Typical Development Cycle

```bash
# 1. Create feature branch from dev
git checkout dev
git pull
git checkout -b feature/my-feature

# 2. Make changes and iterate quickly
make test-dev          # Fast feedback (~15s)

# 3. Format and check quality
make format            # Auto-format code
make lint              # Check for issues

# 4. Run comprehensive tests before commit
make check             # format + test-dev

# 5. Commit changes (pre-commit hooks run automatically)
git commit -m "feat: add my feature"

# 6. Run full validation before PR
make test-release      # All tests with coverage
```

### Working with Tests

#### Running Specific Tests

```bash
# Run specific test file
uv run pytest tests/unit/test_design.py -v

# Run specific test function
uv run pytest tests/unit/test_design.py::test_basic_design -v

# Run tests matching pattern
uv run pytest -k "test_gc_content" -v

# Run with live output
uv run pytest tests/unit/ -v -s
```

#### Test Markers

Tests are organized by markers for flexible execution:

```bash
# By tier (auto-assigned)
uv run pytest -m "dev"        # Fast unit tests
uv run pytest -m "ci"         # Smoke tests
uv run pytest -m "release"    # Integration tests

# By type
uv run pytest -m "unit"          # Unit tests
uv run pytest -m "integration"   # Integration tests
uv run pytest -m "smoke"         # Quick gate tests

# By requirements
uv run pytest -m "requires_docker"    # Needs Docker
uv run pytest -m "requires_network"   # Needs internet
uv run pytest -m "requires_nextflow"  # Needs Nextflow

# By execution environment
uv run pytest -m "runs_in_container"  # Runs inside Docker

# Combine markers
uv run pytest -m "unit and not requires_network"
```

### Docker Development

```bash
# Build image
make docker-build

# Test the image
make docker-test

# Debug in container
make docker-shell
# Then inside: pytest tests/container/ -v

# Run workflow in container
make docker-run GENE=TP53

# Check if image exists (builds if missing)
make docker-ensure
```

## Code Style and Standards

### Python Code Style

We use several tools to maintain code quality:

#### Ruff Formatter
```bash
# Format codebase (same engine as `make format`)
uv run ruff format src tests

# Check formatting without modifying files
uv run ruff format --check src tests
```

#### Ruff (Linting and Import Sorting)
```bash
# Lint code and apply fixes
uv run ruff check --fix src tests

# Lint without changing files (used by `make lint`)
uv run ruff check src tests
```

#### MyPy (Type Checking)
```bash
# Type check the codebase
uv run mypy src

# Type check specific files
uv run mypy src/sirnaforge/core/design.py
```

### Code Organization Principles

#### 1. Type Annotations
All functions should include type annotations:

```python
from typing import List, Optional
from pathlib import Path

def process_sequences(
    input_file: Path,
    output_file: Optional[Path] = None,
    max_candidates: int = 10
) -> List[SiRNACandidate]:
    """Process sequences with full type annotations."""
    pass
```

#### 2. Pydantic Models
Use Pydantic for data validation and serialization:

```python
from pydantic import BaseModel, Field, validator

class DesignParameters(BaseModel):
    """Design parameters with validation."""

    sirna_length: int = Field(21, ge=19, le=23)
    gc_min: float = Field(30.0, ge=0, le=100)
    gc_max: float = Field(60.0, ge=0, le=100)

    @validator('gc_max')
    def gc_max_greater_than_min(cls, v, values):
        if 'gc_min' in values and v <= values['gc_min']:
            raise ValueError('gc_max must be greater than gc_min')
        return v
```

#### 3. Error Handling
Use specific exceptions with clear messages:

```python
class SiRNAForgeException(Exception):
    """Base exception for siRNAforge errors."""

class DesignException(SiRNAForgeException):
    """siRNA design specific errors."""

def design_candidates(sequences: List[str]) -> List[SiRNACandidate]:
    if not sequences:
        raise DesignException("No sequences provided for design")

    try:
        # Design logic
        return candidates
    except Exception as e:
        raise DesignException(f"Design failed: {e}") from e
```

#### 4. Documentation
Use comprehensive docstrings:

```python
def calculate_thermodynamics(sequence: str, temperature: float = 37.0) -> ThermodynamicResult:
    """Calculate thermodynamic properties of RNA sequence.

    Args:
        sequence: RNA sequence (A, U, G, C only)
        temperature: Temperature in Celsius for calculations

    Returns:
        ThermodynamicResult containing folding energy, structure, and stability metrics

    Raises:
        ValidationException: If sequence contains invalid characters
        CalculationException: If thermodynamic calculation fails

    Example:
        >>> result = calculate_thermodynamics("AUGCAUGCAUGC")
        >>> print(f"ΔG = {result.free_energy:.2f} kcal/mol")
    """
    pass
```

## Testing Best Practices

### 1. Test Organization

```
tests/
├── unit/                    # Fast, isolated tests
│   ├── test_design.py      # Core algorithm tests
│   ├── test_models.py      # Data model tests
│   └── test_utils.py       # Utility function tests
├── integration/             # Component interaction tests
│   ├── test_workflow_integration.py
│   └── test_api_integration.py
└── pipeline/               # Pipeline and external tool tests
    ├── test_nextflow.py
    └── test_docker.py
```

### 2. Test Fixtures

Use pytest fixtures for common test data:

```python
# conftest.py
import pytest
from pathlib import Path

@pytest.fixture
def sample_sequences():
    """Sample RNA sequences for testing."""
    return [
        "AUGCGCGAUCUCGAUGCAUGU",
        "GCCAUGCGAUCGAUGCGUAUC",
        "AUCGAUGCGCGAUGCUGUGAU"
    ]

@pytest.fixture
def temp_fasta_file(tmp_path, sample_sequences):
    """Temporary FASTA file with sample sequences."""
    fasta_file = tmp_path / "test.fasta"
    with open(fasta_file, 'w') as f:
        for i, seq in enumerate(sample_sequences):
            f.write(f">seq_{i}\n{seq}\n")
    return fasta_file
```

### 3. Parametrized Tests

Test multiple scenarios efficiently:

```python
@pytest.mark.parametrize("gc_min,gc_max,expected_count", [
    (30, 70, 10),  # Relaxed parameters
    (40, 60, 5),   # Strict parameters
    (45, 55, 2),   # Very strict parameters
])
def test_gc_filtering(sample_sequences, gc_min, gc_max, expected_count):
    """Test GC content filtering with various parameters."""
    designer = SiRNADesigner(DesignParameters(gc_min=gc_min, gc_max=gc_max))
    candidates = designer.design_from_sequences(sample_sequences)
    assert len(candidates) == expected_count
```

### 4. Mock External Dependencies

```python
from unittest.mock import patch, MagicMock

@patch('sirnaforge.data.gene_search.httpx.AsyncClient')
async def test_gene_search_api_error(mock_client):
    """Test handling of API errors during gene search."""
    # Setup mock to simulate API failure
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = httpx.HTTPStatusError("API Error", request=None, response=None)
    mock_client.return_value.__aenter__.return_value.get.return_value = mock_response

    # Test that error is handled gracefully
    searcher = GeneSearcher()
    result = await searcher.search_gene("TESTGENE")
    assert not result.success
    assert "API Error" in result.error
```

## Contribution Guidelines

### 1. Pull Request Process

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature/amazing-feature`
3. **Make** your changes with appropriate tests
4. **Run** quality checks: `make check`
5. **Commit** with clear messages: `git commit -m 'Add amazing feature'`
6. **Push** to your fork: `git push origin feature/amazing-feature`
7. **Create** a Pull Request with description of changes

### 2. Commit Message Format

Use conventional commit format:

```
type(scope): description

feat(design): add thermodynamic scoring algorithm
fix(cli): handle missing input file gracefully
docs(api): add examples for gene search functions
test(unit): increase coverage for design module
refactor(models): simplify sirna candidate validation
```

### 3. Pull Request Guidelines

**Required Information:**
- Clear description of the change
- Rationale for the change
- Any breaking changes
- Testing approach
- Documentation updates

**Quality Checklist:**
- [ ] Tests pass (`make test`)
- [ ] Code is formatted (`make format`)
- [ ] No linting errors (`make lint`)
- [ ] Documentation updated
- [ ] Changelog updated (for significant changes)

### 4. Code Review Process

**Review Criteria:**
- Code correctness and efficiency
- Test coverage and quality
- Documentation completeness
- Adherence to project standards
- Backward compatibility

## Advanced Development

### 1. Adding New Algorithms

```python
# src/sirnaforge/core/custom_algorithm.py
from abc import ABC, abstractmethod
from .design import BaseDesigner

class CustomAlgorithm(BaseDesigner):
    """Custom siRNA design algorithm."""

    def __init__(self, parameters: DesignParameters):
        super().__init__(parameters)

    def design_candidates(self, sequence: str) -> List[SiRNACandidate]:
        """Implement your custom design logic."""
        candidates = []
        # Your algorithm here
        return candidates

# Register the algorithm
from sirnaforge.core.design import DesignRegistry
DesignRegistry.register("custom", CustomAlgorithm)
```

### 2. Adding New Data Sources

```python
# src/sirnaforge/data/custom_provider.py
from .base import BaseDataProvider

class CustomDataProvider(BaseDataProvider):
    """Custom gene data provider."""

    async def search_gene(self, query: str) -> SearchResult:
        """Implement custom gene search."""
        # Your search logic here
        return SearchResult(success=True, transcripts=transcripts)

# Register the provider
from sirnaforge.data.gene_search import DataProviderRegistry
DataProviderRegistry.register("custom", CustomDataProvider)
```

### 3. Performance Profiling

```bash
# Profile specific functions
uv run python -m cProfile -s cumulative scripts/profile_design.py

# Memory profiling
uv run python -m memory_profiler scripts/profile_memory.py

# Line-by-line profiling
uv run kernprof -l -v scripts/profile_lines.py
```

## Debugging

### 1. Logging

Enable detailed logging during development:

```bash
# Set logging level
export SIRNAFORGE_LOG_LEVEL=DEBUG

# Run with verbose output
uv run sirnaforge workflow TP53 --verbose --output-dir debug_output
```

### 2. Interactive Debugging

```python
# Add breakpoints in code
import pdb; pdb.set_trace()

# Or use ipdb for better interface
import ipdb; ipdb.set_trace()

# Run tests with debugger
uv run pytest tests/unit/test_design.py::test_function --pdb
```

### 3. Testing with Real Data

```bash
# Use examples directory for testing
uv run sirnaforge design examples/sample_transcripts.fasta --verbose

# Test with different parameters
uv run sirnaforge workflow TP53 --output-dir test_output --gc-min 35 --verbose
```

## Release Process

### 1. Version Management

```bash
# Update version in src/sirnaforge/__init__.py
__version__ = "0.2.0"

# Tag the release
git tag -a v0.2.0 -m "Release version 0.2.0"
git push origin v0.2.0
```

### 2. Building and Distribution

```bash
# Build package
make build

# Check build artifacts
ls -la dist/

# Test installation from built package
pip install dist/sirnaforge-0.2.0-py3-none-any.whl
```

### 3. Documentation Deployment

```bash
# Build documentation
make docs

# Publish (copy `docs/_build/html` to hosting of choice)
# Deployment is handled via CI workflows (see .github/workflows)
```

## Getting Help

- **Discord/Slack**: Development discussions
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and ideas
- **Email**: Direct contact with maintainers

Remember: Quality over speed, tests are essential, and clear communication helps everyone!
