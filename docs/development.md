# Development Guide

This guide covers development setup, contribution guidelines, and best practices for siRNAforge development.

## Development Environment Setup

### Prerequisites

- Python 3.9-3.12
- Git
- [uv](https://github.com/astral-sh/uv) (recommended)
- Make (optional, for convenience commands)

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge

# One-command development setup
make install-dev

# Verify installation
uv run pytest
uv run sirnaforge --help
```

### Manual Setup

```bash
# Create virtual environment and install dependencies
uv sync --group dev

# Install pre-commit hooks
uv run pre-commit install

# Verify installation
uv run pytest -x
```

## Development Workflow

### 1. Code Changes

```bash
# Create feature branch
git checkout -b feature/my-new-feature

# Make your changes
# ... edit files ...

# Run tests continuously during development
uv run pytest tests/unit/ --watch

# Run specific test files
uv run pytest tests/unit/test_design.py -v
```

### 2. Quality Checks

```bash
# Format code
make format

# Run all quality checks
make lint

# Run tests with coverage
make test-cov

# Full quality check
make check  # Equivalent to: lint + test
```

### 3. Testing

#### Unit Tests
```bash
# Run all unit tests
uv run pytest tests/unit/

# Run with coverage report
uv run pytest tests/unit/ --cov=sirnaforge --cov-report=html

# Run specific test class/function
uv run pytest tests/unit/test_design.py::TestSiRNADesigner::test_basic_design
```

#### Integration Tests
```bash
# Run integration tests (may be slower)
uv run pytest tests/integration/ -v

# Skip slow tests during development
uv run pytest -m "not slow"
```

#### Pipeline Tests
```bash
# Test Nextflow pipeline components
uv run pytest tests/pipeline/

# Run pipeline integration test
bash tests/integration/test_workflow_integration.sh
```

### 4. Documentation

```bash
# Generate CLI documentation
make docs-cli

# Generate usage examples
make docs-examples

# Build complete documentation
make docs-full

# Build Sphinx documentation (when available)
make docs

# Serve documentation locally
make docs-serve
```

## Code Style and Standards

### Python Code Style

We use several tools to maintain code quality:

#### Black (Code Formatting)
```bash
# Format all code
uv run black src tests

# Check formatting without changes
uv run black --check src tests
```

#### Ruff (Linting and Import Sorting)
```bash
# Lint code and fix automatically
uv run ruff check --fix src tests

# Check without fixing
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
    gc_max: float = Field(52.0, ge=0, le=100)

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
make docs-full

# Deploy to GitHub Pages (if configured)
make docs-deploy
```

## Getting Help

- **Discord/Slack**: Development discussions
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and ideas
- **Email**: Direct contact with maintainers

Remember: Quality over speed, tests are essential, and clear communication helps everyone!
