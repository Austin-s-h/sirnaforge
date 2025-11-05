# Docker Test Categories

siRNAforge uses a tier-based testing approach integrated with Docker for comprehensive validation.

## Test Tiers (Auto-assigned by conftest.py)

### Development Tier (`dev` marker)
**Purpose:** Fast validation for active development
**Command:** `make test-dev` or `make docker-test` (in container)
**Expected:** ~15 seconds, ~30 tests

- Unit tests only
- No external dependencies
- No Docker required (unless in container)
- Fast feedback during development

### CI Tier (`ci` marker)
**Purpose:** Quick validation for CI/CD pipelines
**Command:** `make test-ci`
**Expected:** ~40 seconds, with coverage reports

- Smoke tests (quick sanity checks)
- Basic Docker functionality
- CLI availability and help
- Must always pass in CI

### Release Tier (`release` marker)
**Purpose:** Comprehensive validation before release
**Command:** `make test-release`
**Expected:** ~60 seconds, all tests + coverage

- Full integration tests
- Off-target analysis workflows
- Nextflow pipeline execution
- Bioinformatics tool integration
- Must pass before release

## Docker Testing

### Build Docker Image
```bash
make docker-build
# Expected: ~15-20 minutes first time
# Creates sirnaforge:latest image with all tools
```

### Run Tests in Container
```bash
# Run Docker-specific tests inside container
make docker-test
# Expected: ~60 seconds
# Validates: Container environment, all tools, integrations

# Interactive debugging
make docker-shell
# Interactive bash prompt inside container
```

## 🎯 Pytest Markers

Tests are organized using pytest markers for flexible execution:

```python
# Tier markers (auto-assigned)
@pytest.mark.dev          # Development tier - fast unit tests
@pytest.mark.ci           # CI tier - smoke tests
@pytest.mark.release      # Release tier - full integration

# Type markers
@pytest.mark.unit              # Unit tests
@pytest.mark.integration       # Integration tests  
@pytest.mark.smoke             # Quick sanity checks

# Requirement markers
@pytest.mark.requires_docker    # Needs Docker daemon
@pytest.mark.requires_network   # Needs internet
@pytest.mark.requires_nextflow  # Needs Nextflow
@pytest.mark.runs_in_container  # Runs inside Docker
```

### Running Tests by Marker
```bash
# By tier
pytest -m "dev"      # Fast development tests
pytest -m "ci"       # CI smoke tests
pytest -m "release"  # Full release validation

# By type
pytest -m "unit"          # Unit tests only
pytest -m "integration"   # Integration tests only

# By requirements
pytest -m "requires_docker"    # Docker tests only
pytest -m "not requires_network"  # Skip network tests

# Combinations
pytest -m "unit and not requires_network"
```
