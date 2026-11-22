# Makefile for siRNAforge
# Development and release commands use the root Pixi lock.

.PHONY: help install test lint format build clean docker docs all

# Variables
DOCKER_IMAGE = sirnaforge
VERSION := $(shell sed -n '/^\[project\]/,/^\[/s/^version = "\([^"]*\)".*/\1/p' pyproject.toml)
# Digest of everything the image bakes in. VERSION cannot answer "is this image current?" -- it stays
# 0.7.1 across every commit of an unreleased version, so a version-only check passes an image built
# before the code it is meant to validate. Reads the working tree, not git, so uncommitted edits count.
fingerprint = $(shell find $(1) -type f \
	! -name '*.pyc' ! -path '*/__pycache__/*' -print0 2>/dev/null \
	| sort -z | xargs -0 shasum -a 256 2>/dev/null | shasum -a 256 | cut -c1-16)
WHEEL_INPUTS = src pyproject.toml README.md LICENSE
WHEEL_SRC_FINGERPRINT = $(call fingerprint,$(WHEEL_INPUTS))
SRC_FINGERPRINT = $(call fingerprint,$(WHEEL_INPUTS) pixi.toml pixi.lock scripts/check_pixi_environment.py docker .dockerignore)
NEXTFLOW_IMAGE ?= $(DOCKER_IMAGE):$(VERSION)
SIRNAFORGE_NEXTFLOW_IMAGE ?= $(NEXTFLOW_IMAGE)
export SIRNAFORGE_NEXTFLOW_IMAGE

# Host user mapping (prevents root-owned outputs on bind mounts)
HOST_UID = $(shell id -u)
HOST_GID = $(shell id -g)
# Docker daemons with userns-remap cannot represent directory-service IDs
# unless the container opts into the host user namespace.
DOCKER_HOST_NAMESPACE ?= --userns=host
DOCKER_HOST_USER = $(DOCKER_HOST_NAMESPACE) --user $(HOST_UID):$(HOST_GID)

# SiRNAforge data cache mounts (transcriptomes and miRNA databases)
SIRNAFORGE_CACHE_DIR = $(shell echo "$$HOME/.cache/sirnaforge")
SIRNAFORGE_CACHE_MOUNT = -v $(SIRNAFORGE_CACHE_DIR):/home/sirnauser/.cache/sirnaforge

DOCKER_MOUNT_FLAGS = -v $$(pwd):/workspace -w /workspace $(SIRNAFORGE_CACHE_MOUNT)
# Propagate CI-related env vars into the container so tests can reliably
# skip known-flaky network flows in CI (e.g., Ensembl blocks runner IPs).
DOCKER_TEST_ENV = -e CI -e GITHUB_ACTIONS -e SIRNAFORGE_IN_CONTAINER=1 -e PYTEST_ADDOPTS='--basetemp=/workspace/.pytest_tmp' -e SIRNAFORGE_CACHE_DIR=/home/sirnauser/.cache/sirnaforge -e NXF_HOME=/home/sirnauser/.cache/sirnaforge/nextflow/home -e SIRNAFORGE_NEXTFLOW_IMAGE
DOCKER_RUN = docker run --rm $(DOCKER_HOST_USER) $(DOCKER_MOUNT_FLAGS) $(DOCKER_TEST_ENV) $(DOCKER_IMAGE):latest

# GitHub Actions checkouts (and many local workspaces) are often owned by a
# different UID than the container's default non-root user. Use the host
# user mapping when running tests in the container so files created by the
# container are owned by the host user (avoids root-owned debug artifacts).
DOCKER_TEST_USER = $(DOCKER_HOST_USER)
# Optional Docker resource constraints for reproducible small-runner validation.
DOCKER_TEST_MARKERS ?= runs_in_container and not requires_network
DOCKER_TEST_PYTEST_ARGS ?=
DOCKER_TEST_RESOURCE_FLAGS ?=

PIXI_ENV ?= dev
RUN = pixi run --locked -e $(PIXI_ENV)
RUN_DOCS = pixi run --locked -e docs
RUN_RUNTIME = pixi run --locked -e runtime

# Pytest command shortcuts
PYTEST = $(RUN) pytest
PYTEST_V = $(PYTEST) -v
PYTEST_Q = $(PYTEST) -q

#==============================================================================
# HELP
#==============================================================================

# Default target
help: ## Show available commands
	@echo "siRNAforge - Modern Python siRNA Design Toolkit"
	@echo ""
	@echo "Setup & Installation"
	@echo "  make install          Install production dependencies"
	@echo "  make dev              Quick dev setup (install + pre-commit)"
	@echo ""
	@echo "Testing - By Tier (matches marker structure)"
	@echo "  make test-dev         Fast unit tests for dev iteration (1,788 tests, ~36s)"
	@echo "                        Whole suite is 1,868 tests: 1,788 dev + 46 container + release/ci."
	@echo "  make test-ci          Smoke tests for CI/CD (40 tests, ~13s, writes coverage.xml)"
	@echo "  make test-release     Complete release validation, host + container, combined coverage (~6.5min)"
	@echo "  make test-release-host      Host-only release suite (1,822 tests serial, ~83s, coverage base)"
	@echo "  make test-release-container Container release suite (~4.7min, expects host coverage)"
	@echo "  make test             All tests (may have skips/failures)"
	@echo "                        Rebuild adds ~3min for a source change, 15-20min on a cold cache."
	@echo ""
	@echo "Docker Testing"
	@echo "  make docker-build-test Clean + build + test Docker image (all-in-one)"
	@echo "  make docker-test      Run container validation tests INSIDE Docker"
	@echo "  make docker-build     Build Docker image"
	@echo "  make docker-shell     Interactive shell in Docker"
	@echo "  make docker-nextflow-help Show embedded Nextflow pipeline help"
	@echo "  make cache-info       Show data cache locations and status"
	@echo ""
	@echo "Code Quality"
	@echo "  make lint             Check code quality (ruff + mypy)"
	@echo "  make format           Auto-format code"
	@echo "  make check            lint + format + test-dev"
	@echo ""
	@echo "Other"
	@echo "  make docs             Build documentation"
	@echo "  make clean            Clean build artifacts"
	@echo "  make build            Build package"
	@echo "  make version          Show version"

#==============================================================================
# INSTALLATION
#==============================================================================

install: ## Install locked runtime dependencies and the local package
	$(RUN_RUNTIME) install-app

dev: ## Install locked dev and docs environments
	$(RUN) install-app
	$(RUN_DOCS) install-app
	@if [ "$$(git rev-parse --git-dir)" = "$$(git rev-parse --git-common-dir)" ]; then \
		$(RUN) pre-commit install; \
	else \
		echo "Linked worktree: run make pre-commit to run the hooks."; \
	fi
	@echo "Ready for development!"

#==============================================================================
# TESTING - BY TIER (Matches marker structure)
#==============================================================================

test-dev: ## Development tier - fast unit tests (1,788 tests, ~36s)
	$(PYTEST_V) -m "dev"

test-ci: ## CI tier - smoke tests for CI/CD (host-only, skip Docker/Nextflow suites)
	$(PYTEST_V) -m "ci and not runs_in_container and not requires_docker and not requires_tools and not requires_network" -n 0 --junitxml=pytest-report.xml \
		--cov=sirnaforge --cov-report=xml:coverage.xml --cov-report=term-missing

test-release: ## Run the release gate in order
	$(MAKE) docs
	$(MAKE) test-release-host
	$(MAKE) test-release-container
	$(MAKE) test-release-report

# `-n 0` is deliberate: it overrides the repo-wide `-n 2` in pyproject.toml so this stage runs serial.
# The reason is one wall-clock assertion, not coverage (pytest-cov combines xdist workers fine, and the
# container stage already proves cross-process merging works here). Keeping it serial also keeps the
# release log readable. Retained after `test_performance_guard_2mb_reference` moved to CPU time,
# because serial is what makes stage timings comparable between runs.
test-release-host: ## Host-only release suite, serial (1,822 tests, ~83s, produces base coverage database)
	@echo "Step 1/3: Running host-based tests with coverage..."
	@rm -f .coverage coverage*.xml pytest-*.xml 2>/dev/null || true
	$(PYTEST_V) -m "(dev or ci or release) and not runs_in_container and not requires_network" \
		-n 0 \
		--cov=sirnaforge --cov-report= \
		--junitxml=pytest-host-report.xml
	@echo ""

test-release-container: docker-smoke ## Run offline container tests and combine their coverage with host data
	@test -f .coverage || { echo "Run make test-release-host first or provide host coverage."; exit 1; }
	COVERAGE_FILE=.coverage.pixi-container $(MAKE) docker-test \
		DOCKER_TEST_RESOURCE_FLAGS="$(DOCKER_TEST_RESOURCE_FLAGS) --network none" \
		DOCKER_TEST_PYTEST_ARGS="--cov=sirnaforge --cov-report= --junitxml=pytest-container-report.xml"
	$(RUN) coverage combine --append .coverage.pixi-container

test-release-report: ## Generate coverage rollups used for release verification artifacts
	@echo "Step 3/3: Generating combined coverage reports..."
	@$(RUN) coverage report -m
	@$(RUN) coverage xml -o coverage.xml
	@$(RUN) coverage html -d htmlcov
	@echo ""
	@echo "Release validation complete!"
	@echo "   Test Reports: pytest-host-report.xml, pytest-container-report.xml"
	@echo "   Coverage: coverage.xml (HTML: htmlcov/index.html)"
	@echo ""
	@echo "   Summary:"
	@$(RUN) coverage report --format=total 2>/dev/null | awk '{printf "      Total Coverage: %.1f%%\n", $$1}' || echo "      (see above)"

test: ## Run all tests (shows what passes/skips/fails)
	$(PYTEST_V)

#==============================================================================
# TESTING - SPECIAL CATEGORIES
#==============================================================================

# There is no `test-requires-docker` target. `-m requires_docker` selects 0 of 1,864 tests and pytest
# exits 0 on a fully deselected run, so the target was a green gate over an empty set - advertised in
# docs/developer/testing_guide.md as a ~45s subset, and read as evidence during release triage. The
# marker itself stays declared in pyproject.toml and filtered out of `test-ci`, so the day a
# host-side Docker test is written it lands in a tier that already excludes it from CI.

test-requires-network: ## Explicit local-only live service checks (requires working remote APIs)
	$(PYTEST_V) --run-network -m "requires_network"

test-requires-nextflow: ## Tests requiring Nextflow (2 tests, both container-tier)
	$(PYTEST_V) -m "requires_nextflow"

#==============================================================================
# DOCKER
#==============================================================================

print-src-fingerprint: ## Print the source fingerprint an image must carry to be considered current
	@echo "$(SRC_FINGERPRINT)"

DOCKER_TEST_IMAGE ?= $(DOCKER_IMAGE):$(VERSION)-test

.PHONY: docker-smoke
docker-smoke: docker-ensure ## Smoke-test the runtime image offline as an arbitrary UID
	docker run --rm $(DOCKER_TEST_RESOURCE_FLAGS) --network none --user 12345:12345 $(DOCKER_IMAGE):latest \
		bash -ec 'nextflow -version; bwa-mem2 version; bwa-mem3 version; samtools --version; RNAfold --version; \
		python -c "import Bio, RNA, fm_index, numpy, pandas, pyarrow, pysam, sirnaforge"; \
		for tool in pytest jupyter sphinx-build; do if command -v "$$tool"; then exit 1; fi; done; sirnaforge version'

.PHONY: check-wheel ensure-wheel print-wheel-sha docker-build-locked
check-wheel:
	@set -eu; set -- dist/*.whl; test "$$#" -eq 1 && test -f "$$1"; \
	test "$$(head -n 1 dist/.pixi-build)" = "$(WHEEL_SRC_FINGERPRINT)"; \
	tail -n +2 dist/.pixi-build | sha256sum --check --status

ensure-wheel:
	@if ! $(MAKE) --no-print-directory -s check-wheel >/dev/null 2>&1; then $(MAKE) build; fi

print-wheel-sha: check-wheel ## Print the verified wheel digest
	@cut -d ' ' -f 1 < dist/.pixi-build | tail -n 1

DOCKER_BUILD_ENV ?= runtime
DOCKER_BUILD_IMAGE ?= $(DOCKER_IMAGE):$(VERSION)
docker-build-locked: check-wheel
	docker build --platform linux/amd64 -f docker/Dockerfile \
		--build-arg VERSION=$(VERSION) --build-arg SRC_FINGERPRINT=$(SRC_FINGERPRINT) \
		--build-arg PIXI_ENV=$(DOCKER_BUILD_ENV) --build-arg WHEEL_SHA256="$$(tail -n 1 dist/.pixi-build | cut -d ' ' -f 1)" \
		-t "$(DOCKER_BUILD_IMAGE)" .

docker-build: ensure-wheel ## Build the runtime image from the verified wheel and locked runtime
	$(MAKE) docker-build-locked
	docker tag "$(DOCKER_IMAGE):$(VERSION)" "$(DOCKER_IMAGE):latest"
	@echo "Docker image: $(DOCKER_IMAGE):$(VERSION) (sources $(SRC_FINGERPRINT))"

# Rebuild when :latest is absent, was built from a different version, or was built from different
# sources. Container tests exercise the image's installed package, so an image that predates the code
# under test reports green about code it does not contain.
#
# The version check alone was not enough and had already failed in practice: VERSION stays 0.7.1 across
# every commit of an unreleased version, so a three-day-old image passed the guard while ~180 commits
# landed, and the container tier certified stale code. SRC_FINGERPRINT is the content answer.
docker-ensure: ensure-wheel ## Ensure the image matches the project version AND its sources (build if missing or stale)
	@env=$$(docker image inspect $(DOCKER_IMAGE):latest \
		--format '{{range .Config.Env}}{{println .}}{{end}}' 2>/dev/null); \
	built=$$(printf '%s\n' "$$env" | sed -n 's/^BUILD_VERSION=//p'); \
	fp=$$(printf '%s\n' "$$env" | sed -n 's/^BUILD_SRC_FINGERPRINT=//p'); \
	if [ -z "$$built" ]; then \
		echo "Image $(DOCKER_IMAGE):latest missing or unversioned - building $(VERSION)..."; \
		$(MAKE) docker-build; \
	elif [ "$$built" != "$(VERSION)" ]; then \
		echo "Image $(DOCKER_IMAGE):latest was built from $$built but project is $(VERSION) - rebuilding..."; \
		$(MAKE) docker-build; \
	elif [ -z "$$fp" ] || [ "$$fp" = "unknown" ]; then \
		echo "Image $(DOCKER_IMAGE):latest carries no source fingerprint - rebuilding to establish one..."; \
		$(MAKE) docker-build; \
	elif [ "$$fp" != "$(SRC_FINGERPRINT)" ]; then \
		echo "Image $(DOCKER_IMAGE):latest was built from sources $$fp but the tree is $(SRC_FINGERPRINT) - rebuilding..."; \
		$(MAKE) docker-build; \
	elif [ "$$(printf '%s\n' "$$env" | sed -n 's/^BUILD_WHEEL_SHA256=//p')" != "$$(tail -n 1 dist/.pixi-build | cut -d ' ' -f 1)" ]; then \
		echo "Image $(DOCKER_IMAGE):latest has a different wheel, rebuilding..."; \
		$(MAKE) docker-build; \
	else \
		echo "Image $(DOCKER_IMAGE):latest matches $(VERSION) and sources $(SRC_FINGERPRINT)"; \
	fi

cache-ensure: ## Ensure the host cache directory exists
	@mkdir -p "$(SIRNAFORGE_CACHE_DIR)"

docker-test: cache-ensure docker-ensure ## Run container tests with the same wheel and locked dev tools
	$(MAKE) docker-build-locked DOCKER_BUILD_ENV=dev DOCKER_BUILD_IMAGE=$(DOCKER_TEST_IMAGE)
	docker run --rm $(DOCKER_TEST_RESOURCE_FLAGS) $(DOCKER_TEST_USER) \
		-v "$$(pwd):/workspace" -w /workspace -e CI -e GITHUB_ACTIONS -e COVERAGE_FILE \
		-e SIRNAFORGE_IN_CONTAINER=1 -e PYTHONDONTWRITEBYTECODE=1 \
		-e SIRNAFORGE_CACHE_DIR=/tmp/sirnaforge-cache \
		-e SIRNAFORGE_NEXTFLOW_IMAGE=$(NEXTFLOW_IMAGE) $(DOCKER_TEST_IMAGE) \
		python -m pytest tests/container/ -o pythonpath= -n 0 \
		-m '$(DOCKER_TEST_MARKERS)' $(DOCKER_TEST_PYTEST_ARGS)

docker-build-test: ## Clean debug folder, build Docker image, and run tests
	@echo "Cleaning debug folders..."
	@rm -rf tp53_workflow_debug workflow_output workflow_test_debug_*
	@echo "Building Docker image..."
	@$(MAKE) docker-build
	@echo "Running container tests..."
	@$(MAKE) docker-test

docker-shell: docker-ensure ## Interactive shell in Docker
	docker run -it $(DOCKER_HOST_USER) $(DOCKER_MOUNT_FLAGS) $(DOCKER_IMAGE):latest bash

docker-run: GENE ?= TP53  ## Run workflow in Docker (usage: make docker-run GENE=TP53)
docker-run: cache-ensure docker-ensure
	$(DOCKER_RUN) sirnaforge workflow $(GENE) --output-dir docker_results

docker-nextflow-help: cache-ensure docker-ensure ## Show embedded Nextflow pipeline help inside the container
	docker run --rm $(DOCKER_HOST_USER) $(DOCKER_MOUNT_FLAGS) $(DOCKER_TEST_ENV) $(DOCKER_IMAGE):latest bash -c \
		"PIPELINE_NF=\$$(python -c 'from sirnaforge.pipeline.nextflow.runner import NextflowRunner; print(NextflowRunner().get_main_workflow())') && \
		echo \"Embedded pipeline: \$$PIPELINE_NF\" && \
		nextflow run \"\$$PIPELINE_NF\" --help"

# Aliases
docker: docker-build

#==============================================================================
# CODE QUALITY
#==============================================================================

# scripts/ is tracked, load-bearing tooling (baseline measurement, fixture builders) and was
# outside every gate until 0.7.1, so it is now fully linted and type-checked with no exclusions.
# tests/ stays out of mypy deliberately: `mypy tests` reports 125 pre-existing errors across 14
# files, which is its own work package, not this one.
RUFF_PATHS := src tests scripts
MYPY_PATHS := src scripts

lint: ## Check code quality
	$(RUN) ruff check $(RUFF_PATHS)
	$(RUN) ruff format --check $(RUFF_PATHS)
	$(RUN) mypy $(MYPY_PATHS)
	@echo "Code quality checks passed!"

format: ## Auto-format code
	$(RUN) ruff format $(RUFF_PATHS)
	$(RUN) ruff check --fix $(RUFF_PATHS)
	@echo "Code formatted!"

check: format test-dev ## Quick check: format + fast tests

#==============================================================================
# BUILD & RELEASE
#==============================================================================

build: ## Build package
	$(RUN) python -m hatchling build
	@set -eu; set -- dist/*.whl; test "$$#" -eq 1 && test -f "$$1"; \
	{ printf '%s\n' "$(WHEEL_SRC_FINGERPRINT)"; sha256sum "$$1"; } > dist/.pixi-build
	@echo "Package built in dist/"

clean: ## Clean build and cache artifacts
	rm -rf dist/ build/ src/*.egg-info/ .pytest_cache/ .pytest_tmp/
	rm -rf .coverage .coverage.* htmlcov/ .mypy_cache/ .ruff_cache/
	rm -rf coverage*.xml pytest-*.xml
	rm -rf docs/_build/ work/ .nextflow* nextflow_results/
	rm -rf tp53_workflow_debug workflow_output workflow_test_debug_* docker_results
	find . \( -name .pixi -o -name .venv \) -prune -o -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . \( -name .pixi -o -name .venv \) -prune -o -type f -name "*.pyc" -exec rm -f {} +
	@echo "Cleaned!"

.PHONY: print-version
print-version: ## Print the package version from project metadata
	@echo "$(VERSION)"

version: ## Show version
	@echo "siRNAforge version: $(VERSION)"

#==============================================================================
# DOCUMENTATION
#==============================================================================

docs: ## Build documentation (warnings fail the build)
	$(RUN_DOCS) sphinx-build -W --keep-going -b html docs docs/_build/html
	@echo "Docs: docs/_build/html/index.html"

docs-serve: docs ## Serve docs locally
	@echo "http://localhost:8000 (Ctrl+C to stop)"
	@$(RUN_DOCS) python -m http.server 8000 --directory docs/_build/html

#==============================================================================
# UTILITIES
#==============================================================================

example: ## Run basic example
	@mkdir -p examples/output
	$(RUN) sirnaforge design examples/sample_transcripts.fasta -o examples/output/results.tsv

pre-commit: ## Run pre-commit hooks
	$(RUN) pre-commit run --all-files

nextflow-check: ## Check Nextflow installation
	@$(RUN_RUNTIME) nextflow -version || echo "Nextflow not available"

# Advisory, not a gate: both legs end in `|| true` on purpose, because bandit's 8 pre-existing HIGH
# findings are all B324 weak-MD5 on cache-key derivation (not security decisions) and pip-audit needs
# a network route CI has and a laptop behind an interception proxy may not.
#
# The dependency leg was `safety check`, pinned <3.0.0, which needs pkg_resources and so has emitted
# `Unhandled exception ... No module named 'pkg_resources'` into safety-report.json on Python 3.12 --
# not JSON, no result -- while CI uploaded the file as `security-reports-<sha>`. A green target with
# an empty scan reads as "scanned, clean". pip-audit replaces it, and the summary below states
# outright whether a dependency scan actually happened.
security: ## Run security checks (advisory: reports findings, does not fail the build)
	@echo "Running security scans..."
	@$(RUN) bandit -r src/ -f json -o bandit-report.json || true
	@$(RUN) bandit -r src/ -q || true
	@rm -f pip-audit-report.json
	@$(RUN) pip-audit --progress-spinner=off --format=json --output=pip-audit-report.json || true
	@if [ -s pip-audit-report.json ]; then \
		echo "Dependency scan: pip-audit wrote pip-audit-report.json"; \
	else \
		echo "Dependency scan: NO RESULT - pip-audit produced nothing (offline, or proxy TLS). Treat as unscanned."; \
	fi
	@echo "Security scan complete (advisory reports: bandit-report.json, pip-audit-report.json)"

cache-info: ## Show data cache locations and status
	@echo "SiRNAforge Data Cache Information"
	@echo "=================================="
	@echo ""
	@echo "Cache Directories:"
	@echo "  Transcriptomes: $(SIRNAFORGE_CACHE_DIR)/transcriptomes"
	@echo "  miRNA DBs:      $(SIRNAFORGE_CACHE_DIR)/mirna"
	@echo ""
	@if [ -d "$(SIRNAFORGE_CACHE_DIR)" ]; then \
		echo "Status: ✅ Cache directory exists"; \
		echo ""; \
		echo "Disk Usage:"; \
		du -sh "$(SIRNAFORGE_CACHE_DIR)" 2>/dev/null || echo "  (unable to calculate)"; \
		if [ -d "$(SIRNAFORGE_CACHE_DIR)/transcriptomes" ]; then \
			echo "  Transcriptomes: $$(du -sh $(SIRNAFORGE_CACHE_DIR)/transcriptomes 2>/dev/null | cut -f1)"; \
			echo "  Files: $$(find $(SIRNAFORGE_CACHE_DIR)/transcriptomes -type f 2>/dev/null | wc -l)"; \
		fi; \
		if [ -d "$(SIRNAFORGE_CACHE_DIR)/mirna" ]; then \
			echo "  miRNA DBs: $$(du -sh $(SIRNAFORGE_CACHE_DIR)/mirna 2>/dev/null | cut -f1)"; \
			echo "  Files: $$(find $(SIRNAFORGE_CACHE_DIR)/mirna -type f 2>/dev/null | wc -l)"; \
		fi; \
		echo ""; \
		echo "Docker Mount: ✅ Will be mounted to container at /home/sirnauser/.cache/sirnaforge"; \
	else \
		echo "Status: ⚠️  Cache directory does not exist yet"; \
		echo "        It will be created automatically when downloading data."; \
		echo ""; \
		echo "Docker Mount: ⏭️  No mount (cache empty)"; \
	fi
	@echo ""
	@echo "To populate cache, run workflows or use CLI commands like:"
	@echo "  $(RUN_RUNTIME) sirnaforge workflow <gene> --species human"
