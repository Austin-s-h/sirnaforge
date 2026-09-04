# Testing Guide

Tiered testing approach for different development phases and resources.

## Quick Commands

### Development (Python-only)

```bash
make test-dev           # Fast marker-based tests (~15s) for iteration
make test               # Full pytest run on host (may include skips)
make lint               # Ruff + mypy checks (~5s)
make format             # Auto-format & autofix style issues
make check              # format + test-dev (mutating quick gate)
```

> **Note:** `make check` runs `make format` first, so it will modify files to enforce style before executing tests.

### Docker (Full environment)

```bash
make docker-build       # Build Docker image
make docker-build-test  # and test
make docker-shell       # Interactive debugging
```

### Test Categories by Tier

| Target                   | Purpose                                            | Time | Scope                           | Resources |
| ------------------------ | -------------------------------------------------- | ---- | ------------------------------- | --------- |
| `test-dev`               | Fast development iteration (pytest `-m dev`)       | ~15s | Fastest unit-style set          | Minimal   |
| `test-ci`                | CI/CD smoke with coverage                          | ~40s | `ci` markers + coverage XML     | Low       |
| `test-release`           | Host + container validation with combined coverage | ~60s | dev+ci+release markers          | Medium    |
| `test`                   | Full pytest run (allows skips/failures)            | 60s+ | Entire suite on host            | Medium    |
| `test-requires-docker`   | Host tests needing Docker daemon                   | ~45s | `requires_docker` marker subset | Medium    |
| `test-requires-network`  | Network-access-required subset                     | ~30s | `requires_network` marker       | Low       |
| `test-requires-nextflow` | Nextflow-specific subset                           | ~45s | `requires_nextflow` marker      | Medium    |
| `docker-test`            | Container-only tests (`runs_in_container`)         | ~60s | tests/container suite           | High      |

### Environment Requirement Markers

Three markers describe an environment a test cannot create for itself, and `tests/conftest.py`
skips rather than fails when it is absent. The shared principle: each probe asks whether the
dependency _works_, not whether it is nominally present — a resolvable hostname and an
executable on `PATH` both lie.

- **`runs_in_container`** — skipped unless the interpreter is genuinely inside a container
  (Docker's `/.dockerenv` or Podman's `/run/.containerenv`). These tests need the image's
  bioinformatics tooling (`bwa-mem2`, `RNAfold`, Nextflow), so run them with `make docker-test`.
  On the host they are expected to skip; `make test-ci` and `make test-release-host` filter
  them out entirely.
- **`requires_nextflow`** — skipped unless `nextflow -version` actually exits 0. The launcher on
  `PATH` is only a shim that fetches the framework JAR into `NXF_HOME` on first use, so
  `shutil.which("nextflow")` returns a path even on an image where no real invocation can
  succeed. This probe mirrors `NextflowRunner.validate_installation`, so a test skips exactly
  when the pipeline itself would report `nextflow_unavailable`.
- **`requires_network`** — skipped unless a _verified_ TLS handshake with `rest.ensembl.org`
  succeeds. The probe handshakes through Python's own trust store rather than merely opening
  a TCP socket, because interception proxies accept the connection and then fail verification
  on every request.

### Running network tests behind a TLS-intercepting proxy

Corporate proxies (Zscaler, Netskope and similar) re-sign HTTPS with a private root that lives
in the OS keychain, which Python and `aiohttp` do not read — so `requires_network` tests skip
even though `curl` works. Point Python at a bundle that includes that root:

```bash
# macOS: append the keychain's proxy root to the system bundle
cat /etc/ssl/cert.pem > /tmp/sf_ca_bundle.pem
security find-certificate -a -c Zscaler -p /Library/Keychains/System.keychain >> /tmp/sf_ca_bundle.pem

SSL_CERT_FILE=/tmp/sf_ca_bundle.pem make test-requires-network
```

Export `SSL_CERT_FILE` (and `REQUESTS_CA_BUNDLE` for `requests`-based tooling) in your shell
profile to make this permanent. Substitute your provider's certificate name for `Zscaler`.

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

- **Local development**: Use `test-dev` for iteration and `test` before commits
- **CI/CD**: Use `test-ci` with artifacts
- **Release validation**: Use `test-release` with full coverage
- **Quick validation**: Use `make check` for lint + fast tests

## ZFN Manual Validation

> **⚠️ EXPERIMENTAL — the two runbooks below validate runtime, not correctness.** The ZFN arm ships
> experimental in 0.6.0 with known unfixed defects tracked in
> [#82](https://github.com/Austin-s-h/sirnaforge/issues/82). Every run recorded in those notes used
> the published CCR5 half-site pair verbatim, which under the default `require_opposite_strands=True`
> matches **no site at all** — not even its own on-target locus — because both published sequences
> occur on the hg38 plus strand. So the recorded site counts, recovery figures and backend
> conclusions are timing and plumbing evidence only. **Do not sign off a ZFN change on them, do not
> cite them as validation, and do not use ZFN output for any decision without independent
> validation.** Re-derive them with `--zfn-right-half-site CTTTTGCAGTTT` (the reverse complement of
> the published `AAACTGCAAAAG`) once #82's orientation defect is resolved.

The heavy ZFN benchmarking and real-reference checks are documented as technical validation notes rather than notebook-management policy.

- use [zfn_backend_tuning.md](zfn_backend_tuning.md) for the backend selection rationale — runtime ordering only, per the caveat above
- use [zfn_hg38_primary_test_commands.md](zfn_hg38_primary_test_commands.md) for full hg38 primary reruns — that page carries the reverse-complement workaround needed to make its commands match
- keep chr3 and hg38 durable behavior in `tests/integration/` so reference resolution, search execution, and annotation are exercised together
- when adding ZFN half-site fixtures, build them from the **published** genomic text rather than from `reverse_complement(published)`. The existing fixtures do the latter, which is why #82's orientation defect survived a file named `test_zfn_realworld_ccr5_data.py`.

## miRNA default-backend rollout validation

Use this sequence when touching the internal miRNA backend seam or the Nextflow batch path.

The operational default for miRNA seed analysis is `pyahocorasick`.
Keep `exhaustive_python` as the correctness oracle for parity checks and treat the BWA path as the semantic comparison baseline in environments where BWA exists.
Backend selection remains intentionally internal for this slice; do not add a public CLI or workflow knob unless the product surface is being widened on purpose.

```bash
make lint

docker run --rm \
  -v $(pwd):/workspace \
  -w /workspace \
  -v ~/.cache/sirnaforge:/home/sirnauser/.cache/sirnaforge \
  -e CI \
  -e GITHUB_ACTIONS \
  -e PYTEST_ADDOPTS= \
  -e SIRNAFORGE_CACHE_DIR=/home/sirnauser/.cache/sirnaforge \
  -e NXF_HOME=/home/sirnauser/.cache/sirnaforge/nextflow/home \
  sirnaforge:latest \
  bash -lc 'export PYTHONPATH=/workspace/.pip:/workspace/src && /opt/conda/bin/python -m pytest -n 0 -v \
    tests/container/test_toy_databases_integration.py::test_toy_mirna_seed_backend_matches_bwa_semantic_hits \
    tests/container/test_workflow_modes.py::test_nextflow_mirna_batch_path_uses_default_backend \
    --override-ini="addopts=-ra -q --strict-markers --strict-config --color=yes"'
```

Treat the default backend as rollout-ready only when all of the following hold:

- `tests/unit/test_mirna_seed_backends.py` continues to protect exhaustive-oracle parity and schema compatibility
- `test_toy_mirna_seed_backend_matches_bwa_semantic_hits` passes in the container environment
- `test_nextflow_mirna_batch_path_uses_default_backend` passes and emits aggregated miRNA artifacts through the embedded Nextflow path
- no new public CLI or workflow backend-selection surface is introduced

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
