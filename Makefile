# Makefile for siRNAforge
# Uses uv for Python package management

.PHONY: help install test lint format build clean docker docs all

# Variables
DOCKER_IMAGE = sirnaforge
VERSION = $(shell uv run python -c "from sirnaforge import __version__; print(__version__)" 2>/dev/null || echo "0.1.0")

# Docker configuration
UV_CACHE_MOUNT = $(shell \
	if [ -n "$$CI" ] || [ -n "$$GITHUB_ACTIONS" ]; then echo ""; \
	elif [ -d "$$(uv cache dir 2>/dev/null)" ] && [ -w "$$(uv cache dir 2>/dev/null)" ]; then \
		echo "-v $$(uv cache dir):/home/sirnauser/.cache/uv"; \
	else echo ""; fi)

# SiRNAforge data cache mounts (transcriptomes and miRNA databases)
SIRNAFORGE_CACHE_DIR = $(shell echo "$$HOME/.cache/sirnaforge")
SIRNAFORGE_CACHE_MOUNT = $(shell \
	if [ -d "$(SIRNAFORGE_CACHE_DIR)" ]; then \
		echo "-v $(SIRNAFORGE_CACHE_DIR):/home/sirnauser/.cache/sirnaforge"; \
	else echo ""; fi)

DOCKER_MOUNT_FLAGS = -v $$(pwd):/workspace -w /workspace $(UV_CACHE_MOUNT) $(SIRNAFORGE_CACHE_MOUNT)
# Propagate CI-related env vars into the container so tests can reliably
# skip known-flaky network flows in CI (e.g., Ensembl blocks runner IPs).
DOCKER_TEST_ENV = -e UV_LINK_MODE=copy -e CI -e GITHUB_ACTIONS -e PYTEST_ADDOPTS='--basetemp=/workspace/.pytest_tmp'
DOCKER_RUN = docker run --rm $(DOCKER_MOUNT_FLAGS) $(DOCKER_TEST_ENV) $(DOCKER_IMAGE):latest

# GitHub Actions checkouts (and many local workspaces) are often owned by a
# different UID than the container's default non-root user, which can prevent
# writing coverage/junit outputs into /workspace. Run test containers as root.
DOCKER_TEST_USER = --user 0:0

# Pytest command shortcuts
PYTEST = uv run pytest
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
	@echo "  make test-dev         Fast unit tests for dev iteration (~15s)"
	@echo "  make test-ci          Smoke tests for CI/CD"
	@echo "  make test-release     Complete release validation (host + container tests with combined coverage)"
	@echo "  make test-release-host      Host-only release suite (generates coverage base)"
	@echo "  make test-release-container Container release suite (expects host coverage)"
	@echo "  make test             All tests (may have skips/failures)"
	@echo ""
	@echo "Docker Testing"
	@echo "  make docker-build-test Clean + build + test Docker image (all-in-one)"
	@echo "  make docker-test      Run container validation tests INSIDE Docker"
	@echo "  make docker-build     Build Docker image"
	@echo "  make docker-shell     Interactive shell in Docker"
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

install: ## Install production dependencies
	uv sync --no-dev

dev: ## Quick dev setup (install + pre-commit)
	uv sync
	uv run pre-commit install
	@echo "Ready for development!"

#==============================================================================
# TESTING - BY TIER (Matches marker structure)
#==============================================================================

test-dev: ## Development tier - fast unit tests (~15s)
	$(PYTEST_V) -m "dev"

test-ci: ## CI tier - smoke tests for CI/CD (host-only, skip Docker/Nextflow suites)
	$(PYTEST_V) -m "ci and not runs_in_container and not requires_docker and not requires_tools" -n 0 --junitxml=pytest-report.xml \
		--cov=sirnaforge --cov-report=xml:coverage.xml --cov-report=term-missing

test-release: docs test-release-host test-release-container test-release-report ## Release tier - comprehensive validation (host + container tests with combined coverage)

test-release-host: ## Host-only release suite (produces base coverage database)
	@echo "Step 1/3: Running host-based tests with coverage..."
	@rm -f .coverage coverage*.xml pytest-*.xml 2>/dev/null || true
	$(PYTEST_V) -m "(dev or ci or release) and not runs_in_container" \
		-n 0 \
		--cov=sirnaforge --cov-report= \
		--junitxml=pytest-host-report.xml
	@echo ""

test-release-container: docker-ensure ## Container release suite (expects .coverage from host stage)
	@if [ ! -f ".coverage" ]; then \
		echo "Missing .coverage from host tests. Run 'make test-release-host' first or provide the artifact before running container tests."; \
		exit 1; \
	fi
	@echo "Step 2/3: Running container tests (appending coverage)..."
	@mkdir -p .pytest_tmp && chmod 777 .pytest_tmp 2>/dev/null || true
	docker run --rm $(DOCKER_TEST_USER) $(DOCKER_MOUNT_FLAGS) -e CI -e GITHUB_ACTIONS -e PYTEST_ADDOPTS='' $(DOCKER_IMAGE):latest bash -c \
		"pip install --quiet pytest pytest-cov && \
		/opt/conda/bin/pytest tests/container/ -v -m 'runs_in_container' \
		--cov=sirnaforge --cov-append --cov-report= \
		--junitxml=/workspace/pytest-container-report.xml \
		--override-ini='addopts=-ra -q --strict-markers --strict-config --color=yes'"
	@echo ""

test-release-report: ## Generate coverage rollups used for release verification artifacts
	@echo "Step 3/3: Generating combined coverage reports..."
	@uv run coverage report -m
	@uv run coverage xml -o coverage.xml
	@uv run coverage html -d htmlcov
	@echo ""
	@echo "Release validation complete!"
	@echo "   Test Reports: pytest-host-report.xml, pytest-container-report.xml"
	@echo "   Coverage: coverage.xml (HTML: htmlcov/index.html)"
	@echo ""
	@echo "   Summary:"
	@uv run coverage report --format=total 2>/dev/null | awk '{printf "      Total Coverage: %.1f%%\n", $$1}' || echo "      (see above)"

test: ## Run all tests (shows what passes/skips/fails)
	$(PYTEST_V) || true

#==============================================================================
# TESTING - SPECIAL CATEGORIES
#==============================================================================

test-requires-docker: ## Tests requiring Docker daemon (run on host)
	$(PYTEST_V) -m "requires_docker"

test-requires-network: ## Tests requiring network access
	$(PYTEST_V) -m "requires_network"

test-requires-nextflow: ## Tests requiring Nextflow
	$(PYTEST_V) -m "requires_nextflow"

#==============================================================================
# DOCKER
#==============================================================================

docker-build: ## Build Docker image
	docker build -f docker/Dockerfile -t $(DOCKER_IMAGE):$(VERSION) -t $(DOCKER_IMAGE):latest .
	@echo "Docker image: $(DOCKER_IMAGE):$(VERSION)"

docker-ensure: ## Ensure Docker image exists (build if missing)
	@docker image inspect $(DOCKER_IMAGE):latest >/dev/null 2>&1 || $(MAKE) docker-build

docker-test: docker-ensure ## Run tests INSIDE Docker container (validates image)
	@mkdir -p .pytest_tmp && chmod 777 .pytest_tmp 2>/dev/null || true
	docker run --rm $(DOCKER_TEST_USER) $(DOCKER_MOUNT_FLAGS) -e CI -e GITHUB_ACTIONS -e PYTEST_ADDOPTS='' $(DOCKER_IMAGE):latest bash -c "pip install --quiet pytest && /opt/conda/bin/pytest tests/container/ -v -m 'runs_in_container' --override-ini='addopts=-ra -q --strict-markers --strict-config --color=yes'"

docker-build-test: ## Clean debug folder, build Docker image, and run tests
	@echo "Cleaning debug folders..."
	@rm -rf tp53_workflow_debug workflow_output workflow_test_debug_*
	@echo "Building Docker image..."
	@$(MAKE) docker-build
	@echo "Running container tests..."
	@$(MAKE) docker-test

docker-shell: docker-ensure ## Interactive shell in Docker
	docker run -it $(DOCKER_MOUNT_FLAGS) $(DOCKER_IMAGE):latest bash

docker-run: GENE ?= TP53  ## Run workflow in Docker (usage: make docker-run GENE=TP53)
docker-run: docker-ensure
	$(DOCKER_RUN) sirnaforge workflow $(GENE) --output-dir docker_results

# Aliases
docker: docker-build

#==============================================================================
# CODE QUALITY
#==============================================================================

lint: ## Check code quality
	uv run ruff check src tests
	uv run ruff format --check src tests
	uv run mypy src
	@echo "Code quality checks passed!"

format: ## Auto-format code
	uv run ruff format src tests
	uv run ruff check --fix src tests
	@echo "Code formatted!"

check: format test-dev ## Quick check: format + fast tests

#==============================================================================
# BUILD & RELEASE
#==============================================================================

build: ## Build package
	uv build
	@echo "Package built in dist/"

clean: ## Clean build and cache artifacts
	rm -rf dist/ build/ src/*.egg-info/ .pytest_cache/ .pytest_tmp/
	rm -rf .coverage .coverage.* htmlcov/ .mypy_cache/ .ruff_cache/
	rm -rf coverage*.xml pytest-*.xml
	rm -rf docs/_build/ work/ .nextflow* nextflow_results/
	rm -rf tp53_workflow_debug workflow_output workflow_test_debug_* docker_results
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	@echo "Cleaned!"

version: ## Show version
	@echo "siRNAforge version: $(VERSION)"

#==============================================================================
# DOCUMENTATION
#==============================================================================

docs: ## Build documentation (warnings fail the build)
	uv run sphinx-build -W --keep-going -b html docs docs/_build/html
	@echo "Docs: docs/_build/html/index.html"

docs-serve: docs ## Serve docs locally
	@echo "http://localhost:8000 (Ctrl+C to stop)"
	@cd docs/_build/html && uv run python -m http.server 8000

#==============================================================================
# UTILITIES
#==============================================================================

example: ## Run basic example
	@mkdir -p examples/output
	uv run sirnaforge design examples/sample_transcripts.fasta -o examples/output/results.tsv

pre-commit: ## Run pre-commit hooks
	uv run pre-commit run --all-files

nextflow-check: ## Check Nextflow installation
	@uv run nextflow -version || echo "Nextflow not available"

security: ## Run security checks
	@echo "Running security scans..."
	@uv run bandit -r src/ -f json -o bandit-report.json || true
	@uv run bandit -r src/ -q || true
	@(uv run safety check --output json 2>&1 | grep -v "UserWarning" > safety-report.json) || echo '{"vulnerabilities": [], "scan_failed": true}' > safety-report.json
	@echo "Security scan complete (reports: bandit-report.json, safety-report.json)"

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
	@echo "  uv run sirnaforge workflow <gene> --species human"
