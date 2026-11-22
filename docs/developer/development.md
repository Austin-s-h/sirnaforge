# Development Guide

This guide covers development setup, contribution guidelines, and best practices for siRNAforge development.

> **New to siRNAforge?** Start with [Getting Started](../getting_started.md) for basic installation and usage.

## Development Environment Setup

### Locked Pixi setup

Install Pixi 0.76.1, Git and Make on Linux x86_64 or x86_64
WSL2. Run the commands from the source checkout.

```bash
make dev
make lint
make test-dev
make docs
```

Point the editor at `.pixi/envs/dev/bin/python`. Use
`pixi run --locked -e dev` for individual commands, or
`pixi shell --locked -e dev` to activate a terminal. Use `make install` for
the runtime environment and `pixi run --locked -e runtime sirnaforge --help`
to check it. The setup tasks install the local package in editable mode with
dependency resolution and build isolation disabled. In a linked Git worktree,
use `make pre-commit` to run the Pixi hooks. Hook installation from the main
checkout applies to the repository's worktrees.

The development and runtime environments share the Python 3.12 solve. The
docs environment uses a separate Python 3.12 solve for the Sphinx toolchain.
The `py310` and `py311` environments cover wheel compatibility. Python 3.10
is a compatibility lane since it has reached upstream end of life. Python 3.10 uses fm-index 3.0.2 and PyArrow 25.0.1
from the existing wheel-consumer lock. The full-tool pilot retains fm-index
4.0.0, PyArrow 26.0.0 and Nextflow 26.10.3 from the validated prototype.

Conda owns Python, Java, Nextflow and its JAR, both BWA variants, samtools,
ViennaRNA and the native Python stack. PyPI supplies fm-index and PyArrow.
The Pixi manifest is separate from PEP 517 project metadata so Pixi does not
automatically import the wheel's PyPI Nextflow and ViennaRNA requirements
into a prefix that already has their Conda providers. The CI ownership
check validates the installed wheel requirements and the dev/docs constraints.
The Conda Nextflow and ViennaRNA records supply version evidence because
their packages lack Python distribution metadata. Use the ownership check
to verify those providers. The Python package audit covers Python
distribution metadata. Pixi pins Nextflow updates through the lock and
disables the launcher's separate release-version check.

Disposable lock solves succeeded for macOS x86_64, macOS ARM64 and Linux
ARM64, including the pinned BWA packages. Native installs and tests on
those platforms are still pending. Full tooling is supported on Linux x86_64 and x86_64 WSL2. For the full tools on other hosts, use an x86_64 Docker engine or
`--platform linux/amd64` on an emulator that exposes AVX2 for BWA-MEM3.
The emulated fallback still needs a platform acceptance receipt.
Standard wheel users can install with `pip install sirnaforge`. Workflow
execution also needs Java and external tools on PATH. RNAfold and the
aligners are provided by the full environment.

#### Update and audit

Use the committed lock for setup and CI. In a dedicated branch or disposable
worktree, edit `pixi.toml` and run `pixi lock`, or run
`pixi update PACKAGE` for a targeted update. Review the complete lock diff,
including transitive changes. Run `make test-release`, `make security`, and
`pixi run --locked -e dev python scripts/check_pixi_environment.py`. Record the impact on
index and cache compatibility. Restore the previous manifest and lock with
Git and reinstall using `pixi install --locked` to verify rollback.

The release workflows check the wheel and sdist with locked Twine, install
the wheel in a separate virtual environment, and pass the tested bytes to
`make docker-build`. The image verifies the wheel hash and package identity.
Run `make test-release` locally for the host and offline container gates.

Reverting the consolidation commits restores the prior release definitions.
The pull request description records the rollback and acceptance criteria.
[Environment consolidation PRD](https://github.com/Austin-s-h/sirnaforge/pull/128)

#### Pilot evidence

The clean WSL2 x86_64 install with Pixi 0.76.1 used a lock with SHA-256
`fbc8defeaeb5910ce2fb0513b5708935c873b28be4e0f0d9b7750792139f586f`.
The host release selection passed 2,217 tests on the first run. A repeat
run passed 2,216 tests and measured 2.276 seconds of CPU time on the
2 MB repeat-scan guard, which has a 2.0-second limit. Python 3.10 and 3.11
each passed 41 smoke tests. The wheel passed an isolated installation and
the sdist rebuilt the identical wheel.

A disposable full-lock update passed 2,217 tests in the host selection.
The offline update checks passed under UID 12345, including the Conda
process and staged reference workflow. Restoring the original lock
reinstalled the previous closure and passed the ownership checks.
The updated closure kept the pinned native tools and changed Python from
3.12.15 to 3.12.16, pandas from 3.1.0 to 3.1.1, and other transitive
dependencies. The root lock retains the evaluated pilot baseline.

The pilot audit flagged urllib3 2.6.3. The bundled AWS CLI capped urllib3
at that version. The release environment drops the bundled AWS CLI and
requires urllib3 2.8 or later. Custom workflows that call the AWS CLI
can supply it externally. The lock also updates Python and the CA bundles.

The initial matched-direct-version image comparison measured 1,744,268,058
bytes for Pixi and 2,254,129,853 bytes for the legacy image. The same wheel
was used in both. Across three alternating fresh container launches, CLI
help medians were 2.621 and 1.919 seconds, respectively. Other validation
jobs shared the host during this trial. The PRD calls for maintainer
acceptance of a median startup regression above 20%. The maintainer accepted
the measured startup cost for the smaller image and single lock.

On the same WSL2 host, one cold build took 98.67 seconds with Pixi and
228.32 seconds with the legacy recipe. Warm builds took 0.64 and 1.20
seconds. The dependency and wheel install stages received 538,153,653
bytes for Pixi and 673,276,745 bytes for the legacy recipe. Their installed
prefixes occupied 1,711,272 and 2,135,796 KiB. These runs used matched direct
versions and different transitive closures. Network byte counts include
transport overhead and exclude the operating-system and base-image steps.

## Make Commands Reference

siRNAforge uses a comprehensive `Makefile` for common development tasks. Run `make help` to see all available commands.

### Testing Commands

#### By Tier (Recommended)
```bash
make test-dev        # Fast unit tests (1,788 tests, ~36s) - development iteration
make test-ci         # Smoke tests for CI/CD with coverage (40 tests, ~13s)
make test-release    # Comprehensive validation, host + container (~6.5min)
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

### Variant mode validation & reporting (new)

- Variant avoid/target/both modes now always write `logs/resolved_variants.json`, even when no variants pass filters, so you can audit filters that excluded everything.
- Variant report rows are validated with Pandera to ensure required fields (`chr`, `pos`, `ref`, `alt`, `vcf_style`, AF) are present and typed. Candidate CSVs now include variant annotations (`variant_mode`, `allele_specific`, `targeted_alleles`, `overlapped_variants`) and a dedicated `logs/candidate_variants.json` captures the candidate↔variant links for downstream consumers.
- Local VCF parsing is covered by `tests/unit/test_variant_resolver.py` using a small fixture (`tests/unit/data/test_variants.vcf`) to keep regressions visible.

### Code Quality Commands

```bash
make lint       # Check code quality (ruff + mypy)
make format     # Auto-format code (ruff format + fixes)
make check      # Quick validation: format + test-dev
make pre-commit # Run pre-commit hooks on all files
make security   # Run security scans (bandit + pip-audit; advisory, never fails the build)
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
        M --> TAG[Reviewed version tag]
        TAG --> CI2[CI + release gate]
        CI2 --> DOCKER_STABLE[Publish PyPI / Docker / docs / GitHub release]
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
3. Open a release PR from `dev` into `master`; its release-candidate gate must pass before merging.

### Hygiene Rules (Recommended)

These conventions keep the repo clean and make it obvious what is safe to deploy.

- **Do not use version branches** (e.g. `v0.3.3`). Use Git tags (`v0.3.3`) for versions; use branches for ongoing work.
- **Branch naming**:
    - `feature/<short-description>` for new capabilities
    - `fix/<short-description>` for bug fixes
    - `chore/<short-description>` for maintenance (deps, CI, docs)
- **PR targets**:
    - Most PRs merge into `dev`
    - Only “release PRs” merge `dev` → `master`
- **Keep branches fresh**: before starting work, update `dev` (`git checkout dev && git pull`). Rebase or merge `dev` into your feature branch if it drifts.
- **Delete merged branches**: once a PR is merged, delete the source branch (feature/fix/chore) to avoid clutter.

### Version & Tag Hygiene

- Treat `master` as “release-ready”: changes merged to `master` should be intentional and reviewed.
- Select the next version by SemVer: compatible fixes increment the patch, new backward-compatible
  behavior increments the minor, and incompatible changes require an explicit compatibility review.
  In the pre-1.0 series, do not mistake an unreleased Git tag for a published version.
- Set the version in `pyproject.toml` and move the release notes out of
  `[Unreleased]` into the corresponding changelog heading. On release day replace `Pending release`
  with the actual `YYYY-MM-DD` date in a reviewed PR. A pending heading cannot be tagged.
- Merge the release PR after Python CI, docs, wheel and the host/container release-candidate gate pass.
  Create an annotated `vX.Y.Z` tag on that reviewed `master` commit and push it. The tag workflow
  checks package/changelog/wheel identity and reruns the full gate before publishing anything.
- If you need pre-release identifiers, prefer PEP 440 formats like `X.Y.Zrc1` (avoid hyphenated variants like `X.Y.Z-rc1`).
- The tag workflow publishes a draft GitHub release, the tested wheel to PyPI, and a versioned
  GHCR image. It checks the published wheel hash and image identity, deploys release docs, then
  publishes the GitHub release. `:latest` advances only after a stable release is verified.
- **Recovery:** `v0.7.1` exists but its workflow failed; do not move or delete that tag. The public
  GitHub release and PyPI 0.7.1 do not exist. Use a new version. If a future tag workflow fails
  after publication begins, preserve its draft and artifacts, inspect each registry and the failed
  job, and resolve recovery with a new reviewed workflow change or new version. Do not rerun an
  already published PyPI version or force-push a release tag: immutable artifacts cannot be rolled back.

### Recommended GitHub Branch Protection

For `master` (stable):
- Require pull requests (no direct pushes)
- Require the CI summary and **Host, container, wheel and docs release gate** checks on release PRs
  (GitHub branch protection must be configured by a maintainer; a workflow alone cannot require itself)
- Require at least one review
- Block force pushes

For `dev` (integration):
- Require pull requests
- Require status checks

### Versioning Guidance (PEP 440)

- Stable releases: `X.Y.Z` (published to PyPI)
- Pre-releases: `X.Y.ZrcN` (published to PyPI and GHCR without moving `:latest`)

### Docker Image Tags (Recommended)

- `ghcr.io/<owner>/<repo>:dev`: latest image from the `dev` branch (unstable)
- `ghcr.io/<owner>/<repo>:sha-<shortsha>`: immutable image for a specific commit (useful for debugging)
- `ghcr.io/<owner>/<repo>:<version>` and `:latest`: stable release images published from `master`

**Consumption guidance**:
- For active development/testing: use `:dev` (fast-moving)
- For reproducible experiments: use `:sha-<shortsha>`
- For production: use `:<version>` (preferred) or `:latest` (convenience)

### Typical Development Cycle

```bash
# 1. Create feature branch from dev
git checkout dev
git pull
git checkout -b feature/my-feature

# 2. Make changes and iterate quickly
make test-dev          # Fast feedback (~36s)

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
pixi run --locked -e dev pytest tests/unit/test_design.py -v

# Run specific test function
pixi run --locked -e dev pytest tests/unit/test_design.py::test_basic_design -v

# Run tests matching pattern
pixi run --locked -e dev pytest -k "test_gc_content" -v

# Run with live output
pixi run --locked -e dev pytest tests/unit/ -v -s
```

#### Test Markers

Tests are organized by markers for flexible execution:

```bash
# By tier (auto-assigned)
pixi run --locked -e dev pytest -m "dev"        # Fast unit tests
pixi run --locked -e dev pytest -m "ci"         # Smoke tests
pixi run --locked -e dev pytest -m "release"    # Integration tests

# By type
pixi run --locked -e dev pytest -m "unit"          # Unit tests
pixi run --locked -e dev pytest -m "integration"   # Integration tests
pixi run --locked -e dev pytest -m "smoke"         # Quick gate tests

# By requirements
pixi run --locked -e dev pytest -m "requires_network"   # Needs internet
pixi run --locked -e dev pytest -m "requires_nextflow"  # Needs Nextflow

# By execution environment
pixi run --locked -e dev pytest -m "runs_in_container"  # Runs inside Docker

# Combine markers
pixi run --locked -e dev pytest -m "unit and not requires_network"
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
pixi run --locked -e dev ruff format src tests

# Check formatting without modifying files
pixi run --locked -e dev ruff format --check src tests
```

#### Ruff (Linting and Import Sorting)
```bash
# Lint code and apply fixes
pixi run --locked -e dev ruff check --fix src tests

# Lint without changing files (used by `make lint`)
pixi run --locked -e dev ruff check src tests
```

#### MyPy (Type Checking)
```bash
# Type check the codebase
pixi run --locked -e dev mypy src

# Type check specific files
pixi run --locked -e dev mypy src/sirnaforge/core/design.py
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
pixi run --locked -e dev python -m cProfile -s cumulative scripts/profile_design.py

# Memory profiling
pixi run --locked -e dev python -m memory_profiler scripts/profile_memory.py

# Line-by-line profiling
pixi run --locked -e dev kernprof -l -v scripts/profile_lines.py
```

## Debugging

### 1. Logging

Enable detailed logging during development:

```bash
# Set logging level
export SIRNAFORGE_LOG_LEVEL=DEBUG

# Run with verbose output
pixi run --locked -e dev sirnaforge workflow TP53 --verbose --output-dir debug_output
```

### 2. Interactive Debugging

```python
# Add breakpoints in code
import pdb; pdb.set_trace()

# Or use ipdb for better interface
import ipdb; ipdb.set_trace()

# Run tests with debugger
pixi run --locked -e dev pytest tests/unit/test_design.py::test_function --pdb
```

### 3. Testing with Real Data

```bash
# Use examples directory for testing
pixi run --locked -e dev sirnaforge design examples/sample_transcripts.fasta --verbose

# Test with different parameters
pixi run --locked -e dev sirnaforge workflow TP53 --output-dir test_output --gc-min 35 --verbose
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
