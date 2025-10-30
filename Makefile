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

DOCKER_MOUNT_FLAGS = -v $$(pwd):/workspace -w /workspace $(UV_CACHE_MOUNT)
DOCKER_TEST_ENV = -e UV_LINK_MODE=copy -e PYTEST_ADDOPTS='--basetemp=/workspace/.pytest_tmp'
DOCKER_RUN = docker run --rm $(DOCKER_MOUNT_FLAGS) $(DOCKER_TEST_ENV) $(DOCKER_IMAGE):latest

# Pytest command shortcuts
PYTEST = uv run pytest
PYTEST_V = $(PYTEST) -v
PYTEST_Q = $(PYTEST) -q

#==============================================================================
# HELP
#==============================================================================

# Default target
help: ## Show available commands
	@echo "🧬 siRNAforge - Modern Python siRNA Design Toolkit"
	@echo ""
	@echo "📦 Setup & Installation"
	@echo "  make install          Install production dependencies"
	@echo "  make install-dev      Install with dev dependencies (default)"
	@echo "  make dev              Quick dev setup (install-dev + pre-commit)"
	@echo ""
	@echo "🧪 Testing - By Tier (matches marker structure)"
	@echo "  make test-dev         Fast unit tests for dev iteration (~15s)"
	@echo "  make test-ci          Smoke tests for CI/CD"
	@echo "  make test-release     Full integration + heavy tests"
	@echo "  make test             All tests (may have skips/failures)"
	@echo "  make all              Same as 'make test'"
	@echo ""
	@echo "🧪 Testing - By Type"
	@echo "  make test-unit        Unit tests only"
	@echo "  make test-integration Integration tests only"
	@echo ""
	@echo "🐳 Docker Testing"
	@echo "  make docker-test      Run container validation tests INSIDE Docker"
	@echo "  make docker-build     Build Docker image"
	@echo "  make docker-shell     Interactive shell in Docker"
	@echo ""
	@echo "🔧 Code Quality"
	@echo "  make lint             Check code quality (ruff + mypy)"
	@echo "  make format           Auto-format code"
	@echo "  make check            lint + format + test-dev"
	@echo ""
	@echo "🚀 Other"
	@echo "  make docs             Build documentation"
	@echo "  make clean            Clean build artifacts"
	@echo "  make build            Build package"
	@echo "  make version          Show version"

#==============================================================================
# INSTALLATION
#==============================================================================

install: ## Install production dependencies
	uv sync --no-dev

install-dev: ## Install with dev dependencies
	uv sync
	@echo "✅ Development environment ready!"

dev: install-dev ## Quick dev setup with pre-commit hooks
	uv run pre-commit install
	@echo "🚀 Ready for development!"

#==============================================================================
# TESTING - BY TIER (Matches marker structure)
#==============================================================================

test-dev: ## Development tier - fast unit tests (~15s)
	$(PYTEST_V) -m "dev"

test-ci: ## CI tier - smoke tests for CI/CD
	$(PYTEST_V) -m "ci" --junitxml=pytest-report.xml \
		--cov=sirnaforge --cov-report=xml:coverage.xml --cov-report=term-missing

test-release: ## Release tier - comprehensive validation (all tests with coverage)
	$(PYTEST_V) -m "dev or ci or release" --junitxml=pytest-report.xml \
		--cov=sirnaforge --cov-report=xml:coverage.xml --cov-report=html --cov-report=term-missing

test: ## Run all tests (shows what passes/skips/fails)
	$(PYTEST_V) || true

all: test  ## Alias for 'make test'

#==============================================================================
# TESTING - BY TYPE
#==============================================================================

test-unit: ## Unit tests only
	$(PYTEST_V) -m "unit"

test-integration: ## Integration tests only
	$(PYTEST_V) -m "integration"

test-cov: ## Tests with coverage report
	$(PYTEST) --cov=sirnaforge --cov-report=html --cov-report=term-missing

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
	@echo "✅ Docker image: $(DOCKER_IMAGE):$(VERSION)"

docker-ensure: ## Ensure Docker image exists (build if missing)
	@docker image inspect $(DOCKER_IMAGE):latest >/dev/null 2>&1 || $(MAKE) docker-build

docker-test: docker-ensure ## Run tests INSIDE Docker container (validates image)
	@mkdir -p .pytest_tmp && chmod 777 .pytest_tmp 2>/dev/null || true
	docker run --rm $(DOCKER_MOUNT_FLAGS) -e PYTEST_ADDOPTS='' $(DOCKER_IMAGE):latest bash -c "pip install --quiet pytest && /opt/conda/bin/pytest tests/container/ -v -m 'runs_in_container' --override-ini='addopts=-ra -q --strict-markers --strict-config --color=yes'"

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
	@echo "✅ Code quality checks passed!"

format: ## Auto-format code
	uv run ruff format src tests
	uv run ruff check --fix src tests
	@echo "✅ Code formatted!"

check: format test-dev ## Quick check: format + fast tests

#==============================================================================
# BUILD & RELEASE
#==============================================================================

build: ## Build package
	uv build
	@echo "✅ Package built in dist/"

clean: ## Clean build and cache artifacts
	rm -rf dist/ build/ src/*.egg-info/ .pytest_cache/ .pytest_tmp/
	rm -rf .coverage htmlcov/ .mypy_cache/ .ruff_cache/
	rm -rf docs/_build/ work/ .nextflow* nextflow_results/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	@echo "✅ Cleaned!"

version: ## Show version
	@echo "siRNAforge version: $(VERSION)"

#==============================================================================
# DOCUMENTATION
#==============================================================================

docs: ## Build documentation
	uv run sphinx-build -b html docs docs/_build/html
	@echo "✅ Docs: docs/_build/html/index.html"

docs-serve: docs ## Serve docs locally
	@echo "🌐 http://localhost:8000 (Ctrl+C to stop)"
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
	@uv run nextflow -version || echo "❌ Nextflow not available"

security: ## Run security checks
	@echo "🔐 Running security scans..."
	@uv run bandit -r src/ -q || true
	@echo "✅ Security scan complete"
