# Testing Guide

Tiered testing approach for different development phases and resources.

## Quick Commands

### Development (Python-only)

```bash
make test-dev           # Fast marker-based tests (1,788 tests, ~36s) for iteration
make test               # Full pytest run on host (~2min, may include skips)
make lint               # Ruff + mypy checks (~5s warm, ~30s cold mypy cache)
make format             # Auto-format & autofix style issues
make check              # format + test-dev (mutating quick gate)
```

> **Timings below are measured on a 10-core arm64 laptop with warm caches** (uv, mypy, `docs/_build`,
> and a current `sirnaforge:latest`). They were previously inherited from a 30-test suite and were
> wrong by up to 9x; if you change a tier's population, re-measure rather than editing the number.

> **Note:** `make check` runs `make format` first, so it will modify files to enforce style before executing tests.

### Docker (Full environment)

```bash
make docker-build       # Build Docker image
make docker-build-test  # and test
make docker-shell       # Interactive debugging
```

### Test Categories by Tier

| Target                   | Purpose                                            | Time     | Scope                                      | Resources |
| ------------------------ | -------------------------------------------------- | -------- | ------------------------------------------ | --------- |
| `test-dev`               | Fast development iteration (pytest `-m dev`)       | ~36s     | 1,788 tests — unit-style, hermetic         | Minimal   |
| `test-ci`                | CI/CD smoke with coverage                          | ~13s     | 40 `ci`-marked tests + coverage XML        | Low       |
| `test-release`           | Host + container validation with combined coverage | ~6.5min  | 1,822 host + 46 container, merged coverage | Medium    |
| `test`                   | Full pytest run (allows skips/failures)            | ~2min    | Entire 1,868-test suite on host            | Medium    |
| `test-requires-network`  | Network-access-required subset                     | 15 tests | `requires_network` marker                  | Low       |
| `test-requires-nextflow` | Nextflow-specific subset                           | 2 tests  | `requires_nextflow` marker                 | Medium    |
| `docker-test`            | Container-only tests (`runs_in_container`)         | ~4.7min  | tests/container suite                      | High      |

`test-release` also pays a rebuild whenever the image is stale, because `docker-ensure` rebuilds when
`sirnaforge:latest` was built from a different version _or_ a different source fingerprint. Two very
different costs hide behind that:

- **Source-only change** (anything under `src/`, which is the usual case): **~3 minutes**. The conda and
  bioinformatics-tool layers are cache hits; only the `COPY src` layers and the export re-run. Measured
  at 3m11s.
- **Cold cache** (no Docker layer cache, changed `Dockerfile`, changed `uv.lock`): **15-20 minutes**,
  because the tool stack is rebuilt.

So a full `make test-release` after an ordinary source edit is ~9.5 minutes, not ~25. The two
`requires_*` subsets are sized in tests rather than seconds because both skip entirely without the
capability, and neither is a meaningful timing.

There is no `test-requires-docker` target. `-m requires_docker` matches nothing, and pytest exits 0
on a fully deselected run, so the target was green over an empty set — the worst possible thing to
read during release triage. The marker stays declared for the day a host-side Docker test exists.

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

### What each `test-release` stage is actually for

`test-release` is `docs`, then the host suite, then the container suite, then a coverage rollup.
The two test stages answer different questions and are not interchangeable:

- **Host stage** — every `dev`/`ci`/`release` test that is not `runs_in_container`, serial, with
  coverage. 1,822 selected, ~83s. This is where line coverage comes from: 83% overall, and 90-100% on
  everything 0.7.1 added. It reaches that number with the screen stubbed and Nextflow monkeypatched.
- **Container stage** — 46 tests that drive the installed CLI with `subprocess.run` inside the
  image, ~4.7min. It answers "does the artifact work", which coverage cannot measure: coverage does not
  follow a subprocess, so these tests contribute **0% of the new modules** to the rollup while
  genuinely executing them. That is expected. Do not plumb `COVERAGE_PROCESS_START` into the
  container to raise the number — it would tell you nothing the host stage has not already said.
  Assert the artifacts instead (`tests/container/release_artifacts.py`), which is how the tier
  distinguishes "screened clean" from "never screened".

Expect **8 of the 46 container tests to skip** without a verified TLS route, and know which 8: the TP53
full workflow, the five variant/SNP modes, the gene-symbol search and the custom-transcriptome
off-target run. They are the most realistic tests in the repository, so a green container stage on a
proxied network is a statement about the packaged artifact and the toy-reference screen, not about the
gene name → Ensembl → transcripts → screen path.

The container stage is the only place the report (#103), the `.nf` evidence emission (#100) and the
installed `benchmark`/`report` commands run for real, so it is the only place a packaging or
integration defect in them can be caught. `workflow.py` turns a failed render and a failed manifest
write into `logger.warning`, so the CLI exits 0 either way: the assertions in
`release_artifacts.py` are what make those failures visible.

> **Do not trust the coverage number if anything else touched this working tree.** There are two
> independent ways to get a wrong percentage, and the second is the one that actually bites.
>
> 1. **A second writer.** `[tool.coverage.run] data_file` is the fixed path `.coverage`, and
>    `make test-ci`, a bare `uv run pytest --cov` and the container stage's `--cov-append` all write it.
>    A concurrent run in the same checkout rewrites the host stage's data between step 1 and step 3.
> 2. **A source edit between measurement and report.** `coverage report` stores bare line numbers and
>    re-parses the source file when it renders, so an edit landing after step 1 makes the recorded lines
>    land on the wrong statements. Measured on this repository: the rollup printed **83%** with
>    `workflow.py` at **87%**, and re-reading the _same, untouched_ `.coverage` twelve minutes later
>    printed **78%** with `workflow.py` at **56%** — one concurrent commit to that file, no test re-run,
>    no change to the database. Only the edited file drifted; every other module read identically.
>
> Test counts are unaffected either way, because each pytest run reports its own. If the total looks
> wrong, re-run `make test-release-host` on a quiescent tree and report on that, and set `COVERAGE_FILE`
> to a scratch path for anything else that needs coverage meanwhile.

### A `dev` test must never reach the network

An unmarked test that calls out is worse than a slow test: it passes or fails on the environment.
The two ways it has happened here, both fixed, both cheap to reintroduce:

- **Orthology.** `_process_nextflow_results` resolves cross-species orthologues through Ensembl
  Compara whenever a screen declares more than one species and a hit row comes from a non-query
  species. Behind a TLS-intercepting proxy each call cost ~25s — three attempts per route with 2s + 4s
  backoff — and 20 tests paid it, which was 90% of the `dev` tier. Set
  `WorkflowConfig(ortholog_mapping_file=...)` to state the orthologues in a file instead
  (`tests/unit/data/ortholog_mapping_synthetic.json` is the fixture); that is also the supported
  air-gapped path for a real run.
- **Nextflow.** Driving `step5_offtarget_analysis` runs a real `nextflow run -profile docker`, whose
  failure step 5 swallows — so the test still passes, 23s later. Monkeypatch
  `SiRNAWorkflow._run_nextflow_offtarget_analysis` unless the test is marked `requires_nextflow`.

`make test-dev` should stay under ~45s with no test above ~3s; `-m dev --durations=10` is the check.

The old wording said 30s and 2s, and both clauses had quietly gone false — the tier grew from 30
tests to 1,788 while the numbers stayed. A tripwire that is already tripped is one nobody reads, so
these are re-baselined against measurement: 1,788 tests in ~36s on a 10-core laptop, slowest
item 2.27s. Both ceilings were re-checked after the release run, at load average 13 with another
session working in the same checkout: 36.22s, slowest 2.27s — so the headroom is real and not an idle-machine
artifact. Expect up to 2-3x on a badly contended machine (the same tests have been measured at 54s and
6.1s). If a change pushes either past its ceiling, the question to ask is whether a test acquired a
dependency — those two failure modes are documented above and both have happened.

Do not read the host stage's `--durations` as this ceiling. Under coverage and `-n 0` the same slowest
test measures 5.1s rather than 2.27s; that is instrumentation, not a regression.

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
# Expected: ~36 seconds, 1,788 tests
# ✅ Success: All tests pass, no Docker and no network required
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
# Expected: ~45 seconds, runs format + lint + test-dev
```

### 4. Pre-Commit Validation

```bash
# Run CI-tier tests (quick smoke tests for CI/CD)
make test-ci
# Expected: ~13 seconds, 40 tests
# Includes smoke tests with coverage reports
# Note: this is the historical PR gate and it selects no test added in 0.7.1;
#       `make test-dev` now runs alongside it in CI for exactly that reason.

# Full release validation
make test-release
# Expected: ~6.5 minutes -- docs (~20s), host suite serial (~83s), container suite (~4.7min),
#           coverage rollup (~6s). Add ~3 minutes for a source-only image rebuild, or 15-20
#           on a cold Docker cache. Measured end to end at 9m38s including a source rebuild.
# Note: exactly 15 tests skip without a verified TLS route -- 7 host, 8 container (Ensembl,
#       miRBase); that is honest,
#       and the whole gene-name -> Ensembl -> off-target path is what goes untested with them.

# Full local test suite (all tests, may have skips/failures)
make test
# Expected: ~2 minutes, 1,868 tests, includes all test categories
# Note: the 46 container tests skip on the host by design (`runs_in_container`)
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
# Expected: ~4.7 minutes (46 tests, 8 of them skipping without Ensembl; rebuilds a stale image first)
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
- **Docker builds** take ~15-20 minutes first time, much faster subsequently. `docker-ensure`
  triggers one whenever `src/`, `pyproject.toml`, `uv.lock` or the Dockerfile has changed since the
  image was built, so `make test-release` after a source edit pays that cost once.
- **`make test-dev`** should complete in ~36s (1,788 tests) — see the ceiling above
- **`make test-ci`** completes in ~13s (40 tests) and writes coverage.xml
- **`make test-release`** completes in ~6.5 minutes when the image is current, ~9.5 with a source rebuild

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
