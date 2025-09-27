# Docker Test Categories

Test categorization for CI/CD pre-release handling.

## Test Types

### Smoke Tests (Must Always Pass)
**Purpose:** Basic Docker functionality
**Command:** `make docker-test-smoke`
**CI:** Must pass for ALL releases

- CLI availability and help
- Package imports
- Basic version commands
- Test data validation

### Integration Tests (Can Fail in Pre-release)
**Purpose:** Complex workflows
**Command:** `make docker-test-integration`
**CI:** Can fail in pre-release, must pass for full release

- Off-target analysis workflows
- Nextflow pipeline execution
- Bioinformatics tool integration
- Production configuration

## 🚦 CI/CD Workflow Logic

### Pre-release Mode (`prerelease: true`)
```
✅ Smoke tests pass → Docker image published
⚠️  Integration tests can fail → Release still created with warning
```

### Full Release Mode (`prerelease: false`)
```
✅ Smoke tests pass → Continue to integration
✅ Integration tests pass → Docker image published
❌ Either test category fails → Release blocked
```

## 📋 Make Commands

```bash
# Test individual categories
make docker-test-smoke        # Ultra-fast, must always pass
make docker-test-integration  # Complex workflows, can fail in pre-release

# Combined testing
make docker-test-fast         # Both categories combined
make docker-test-categories   # Test the categorization setup

# Resource-aware testing
make docker-test-smoke        # 256MB RAM, <30s
make docker-test-integration  # 2GB RAM, 1-2min
make docker-test-full         # 8GB RAM, 5-10min
```

## 🧪 Implementation Details

### Pytest Markers
```python
@pytest.mark.smoke           # Basic Docker functionality - MUST pass
@pytest.mark.docker_integration  # Complex workflows - can fail in pre-release
@pytest.mark.docker          # All Docker-based tests
```

### Marker Combinations
```bash
# Smoke tests only
pytest -m 'docker and smoke'

# Integration tests only
pytest -m 'docker and (docker_integration or (not smoke))'

# All Docker tests
pytest -m 'docker'
```

### GitHub Actions Jobs
- `test-docker-smoke`: Always required, blocks release if failed
- `test-docker-integration`: Uses `continue-on-error` in pre-release mode

## 🔍 Validation Script

Use the validation script to test your setup:
```bash
./scripts/test_docker_categories.sh
```

This script runs both test categories and shows how they behave differently in pre-release vs full release scenarios.

---

**Related Documentation:**
- [Testing Guide](testing_guide.md) - Complete testing documentation
- Release workflow in `.github/workflows/release.yml`
- Docker documentation in `docker/` directory
