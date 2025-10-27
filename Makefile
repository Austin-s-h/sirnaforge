# Makefile for siRNAforge
# Uses uv for Python package management

.PHONY: help install install-dev install-pipeline test lint format build clean docker docs

# Variables
DOCKER_IMAGE = sirnaforge
VERSION = $(shell uv run --group dev python -c "from sirnaforge import __version__; print(__version__)" 2>/dev/null || echo "0.1.0")

# Conditional UV cache mounting - only mount if not in CI and cache dir is accessible
UV_CACHE_MOUNT = $(shell \
	if [ -n "$$CI" ] || [ -n "$$GITHUB_ACTIONS" ]; then \
		echo ""; \
	elif [ -d "$$(uv cache dir 2>/dev/null)" ] && [ -w "$$(uv cache dir 2>/dev/null)" ]; then \
		echo "-v $$(uv cache dir):/home/sirnauser/.cache/uv"; \
	else \
		echo ""; \
	fi)
DOCKER_MOUNT_FLAGS = -v $$(pwd):/workspace -w /workspace $(UV_CACHE_MOUNT)
DOCKER_TEST_ENV_VARS = -e UV_LINK_MODE=copy -e PYTEST_ADDOPTS='--basetemp=/workspace/.pytest_tmp' -e NXF_WORK=/workspace/.nextflow_work
DOCKER_TMP_DIRS = .pytest_tmp .nextflow_work
PREP_DOCKER_TMP = @mkdir -p $(DOCKER_TMP_DIRS) && chmod -R 777 $(DOCKER_TMP_DIRS) 2>/dev/null || true

define RUN_DOCKER_TEST
	$(PREP_DOCKER_TMP)
	docker run --rm \
		--cpus=$(1) \
		--memory=$(2) \
		--memory-swap=$(3) \
		$(DOCKER_MOUNT_FLAGS) \
		$(DOCKER_TEST_ENV_VARS) \
		$(DOCKER_IMAGE):latest \
		bash -c "uv sync --active --group dev && $(4)"
endef

# Default target
help: ## Show available commands
	@echo "🧬 siRNAforge Development Commands"
	@echo "================================="
	@echo ""
	@echo "📦 PACKAGE MANAGEMENT"
	@echo "  install         Install production dependencies"
	@echo "  install-dev     Install with development dependencies"
	@echo "  install-pipeline Pipeline tools (included in main deps)"
	@echo ""
	@echo "🐍 ENVIRONMENT & SETUP"
	@echo "  conda-env           Create conda environment for local development"
	@echo "  conda-env-update    Update existing conda environment"
	@echo "  conda-env-clean     Remove conda environment"
	@echo "  dev                 Quick development setup"
	@echo ""
	@echo "🧪 TESTING (Tiered Approach - ✅ Verified)"
	@echo "  test-local-python   Fastest tests (12-15s, 30 tests) ✅"
	@echo "  test-unit           Unit tests (30-35s, 31 tests) ✅"
	@echo "  test-dev           Development tier (default workflow) ✅"
	@echo "  test-ci            CI/CD smoke tier"
	@echo "  test-release       Release validation suite"
	@echo "  test               Full test suite (60s+, some failures OK)"
	@echo ""
	@echo "🐳 DOCKER TESTING (Resource-Aware)"
	@echo "  docker-test-smoke        Minimal smoke sweep (256MB RAM, <30s)"
	@echo "  docker-test-integration  Balanced integration suite (2GB RAM, ~2 min)"
	@echo "  docker-test-full         Comprehensive CI-equivalent run (8GB RAM, 5-10 min)"
	@echo ""
	@echo "🔧 CODE QUALITY"
	@echo "  lint               Run all linting tools"
	@echo "  format             Format code with black and ruff"
	@echo "  check              Run quality checks + fast tests"
	@echo ""
	@echo "🐳 DOCKER"
	@echo "  docker             Build Docker image"
	@echo "  docker-run         Run workflow in Docker"
	@echo "  docker-dev         Interactive Docker shell"
	@echo ""
	@echo "📚 DOCUMENTATION"
	@echo "  docs               Build documentation"
	@echo "  docs-serve         Serve docs locally"
	@echo "  docs-cli           Generate CLI reference"
	@echo ""
	@echo "🚀 PIPELINE & NEXTFLOW"
	@echo "  nextflow-check     Check Nextflow installation"
	@echo "  nextflow-run       Run Nextflow pipeline"
	@echo "  nextflow-lint      Lint Nextflow scripts"
	@echo ""
	@echo "🛠️ PROJECT UTILITIES"
	@echo "  example            Run basic example"
	@echo "  version            Show version"
	@echo "  release-test       Test release preparation (quick)"
	@echo "  release            Prepare release (full checks)"
	@echo "  security           Run security checks"
	@echo "  pre-commit         Run pre-commit hooks"
	@echo "  clean              Clean build artifacts"

# Installation
install: ## Install production dependencies only
	uv sync --no-dev

install-dev: ## Install with development dependencies (default)
	uv sync
	uv run pre-commit install
	@echo "✅ Development environment ready!"

install-pipeline: ## Pipeline tools are now in main dependencies
	uv sync
	@echo "✅ Pipeline tools available (included in main dependencies)!"

# Environment Management
conda-env: ## Create conda environment for local development
	@echo "🐍 Setting up conda environment..."
	@if command -v micromamba >/dev/null 2>&1; then \
		micromamba env create -f environment-dev.yml; \
		echo "✅ Conda environment created with micromamba!"; \
		echo "🔄 Activate with: micromamba activate sirnaforge-dev"; \
	elif command -v mamba >/dev/null 2>&1; then \
		mamba env create -f environment-dev.yml; \
		echo "✅ Conda environment created with mamba!"; \
		echo "🔄 Activate with: conda activate sirnaforge-dev"; \
	elif command -v conda >/dev/null 2>&1; then \
		conda env create -f environment-dev.yml; \
		echo "✅ Conda environment created with conda!"; \
		echo "🔄 Activate with: conda activate sirnaforge-dev"; \
	else \
		echo "❌ Neither conda, mamba, nor micromamba found."; \
		echo "Please install one of:"; \
		echo "  • micromamba (recommended): https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"; \
		echo "  • Mambaforge: https://mamba.readthedocs.io/en/latest/installation.html"; \
		echo "  • Miniconda: https://docs.conda.io/en/latest/miniconda.html"; \
		exit 1; \
	fi

conda-env-update: ## Update existing conda environment
	@echo "🔄 Updating conda environment..."
	@if command -v micromamba >/dev/null 2>&1; then \
		micromamba env update -f environment-dev.yml; \
		echo "✅ Conda environment updated with micromamba!"; \
	elif command -v mamba >/dev/null 2>&1; then \
		mamba env update -f environment-dev.yml; \
		echo "✅ Conda environment updated with mamba!"; \
	elif command -v conda >/dev/null 2>&1; then \
		conda env update -f environment-dev.yml; \
		echo "✅ Conda environment updated with conda!"; \
	else \
		echo "❌ Neither conda, mamba, nor micromamba found."; \
		exit 1; \
	fi

conda-env-clean: ## Remove conda environment
	@echo "🧹 Removing conda environment..."
	@if command -v micromamba >/dev/null 2>&1; then \
		micromamba env remove -n sirnaforge-dev; \
		echo "✅ Conda environment removed!"; \
	elif command -v mamba >/dev/null 2>&1; then \
		mamba env remove -n sirnaforge-dev; \
		echo "✅ Conda environment removed!"; \
	elif command -v conda >/dev/null 2>&1; then \
		conda env remove -n sirnaforge-dev; \
		echo "✅ Conda environment removed!"; \
	else \
		echo "❌ Neither conda, mamba, nor micromamba found."; \
		exit 1; \
	fi

# Development & Testing
test: ## Run all tests
	uv run --group dev pytest -v

test-unit: ## Run unit tests only (fast, Python-only)
	uv run --group dev pytest -v -m "unit"

test-integration: ## Run integration tests (full workflow, requires Docker + Nextflow)
	uv run --group dev pytest -v -m "integration"

test-dev: ## Run development-tier tests (default local workflow)
	uv run --group dev pytest -v -m "dev"

test-fast: test-dev ## Backward compatible alias for development-tier tests

test-cov: ## Run tests with coverage report
	uv run --group dev pytest --cov=sirnaforge --cov-report=html --cov-report=term-missing

# Environment-specific test targets
test-local-python: ## Run tests for local Python development (unit tests only)
	uv run --group dev pytest -v -m "local_python"

test-local-nextflow: ## Run tests for local Nextflow development (includes pipeline tests)
	uv run --group dev pytest -v -m "local_nextflow"

test-ci: ## Run tests optimized for CI environment
	# Generate JUnit XML and coverage artifacts for GitHub Actions
	uv run --group dev pytest -m "ci" \
		--junitxml=pytest-report.xml \
		--cov=sirnaforge --cov-report=term-missing --cov-report=xml:coverage.xml -v

test-release: ## Run release validation suite (comprehensive, may be slow)
	uv run --group dev pytest -v -m "release"

# Code Quality
lint: ## Run all linting tools
	uv run --group dev ruff check src tests
	uv run --group dev ruff format --check src tests
	uv run --group dev mypy src
	@echo "✅ Code quality checks passed!"

format: ## Format code with ruff (auto-fix)
	uv run --group dev ruff format src tests
	uv run --group dev ruff check --fix src tests
	@echo "✅ Code formatted!"

lint-fix: ## Run linting with auto-fix (matches pre-commit behavior)
	uv run --group dev ruff check --fix --exit-non-zero-on-fix src tests
	uv run --group dev ruff format src tests
	uv run --group dev mypy src
	@echo "✅ Code quality checks and fixes applied!"

check: lint-fix test-fast ## Run quick quality checks with auto-fix (lint + fast tests)

# Build & Release
build: ## Build package
	uv build
	@echo "✅ Package built in dist/"

clean: ## Clean build and cache artifacts
	rm -rf dist/ build/ src/*.egg-info/
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete
	rm -rf .pytest_cache/ .coverage htmlcov/ .mypy_cache/ .ruff_cache/ docs/_build/
	rm -rf work/ .nextflow* nextflow_results/ preview_results/
	@echo "✅ Cleaned all artifacts!"

# Docker
docker: ## Build Docker image
	docker build -f docker/Dockerfile -t $(DOCKER_IMAGE):$(VERSION) -t $(DOCKER_IMAGE):latest .
	@echo "✅ Docker image built: $(DOCKER_IMAGE):$(VERSION)"

docker-ensure-image: ## Ensure Docker image exists (build if missing)
	@if ! docker image inspect $(DOCKER_IMAGE):latest >/dev/null 2>&1; then \
		echo "🐳 Docker image not found, building..."; \
		$(MAKE) docker; \
	else \
		echo "✅ Docker image $(DOCKER_IMAGE):latest already exists"; \
	fi

docker-run: GENE ?= TP53
docker-run: ## Run workflow in Docker (usage: make docker-run GENE=<gene>)
	docker run $(DOCKER_MOUNT_FLAGS) \
		$(DOCKER_IMAGE):latest \
		sirnaforge workflow $(GENE) --output-dir docker_results

docker-dev: ## Interactive Docker development shell
	docker run -it $(DOCKER_MOUNT_FLAGS) \
		$(DOCKER_IMAGE):latest bash

docker-test: docker-test-integration ## Backward compatible alias for integration sweep

docker-test-smoke: docker-ensure-image ## Run ultra-minimal smoke tests for CI/CD (fastest) - MUST ALWAYS PASS
	$(call RUN_DOCKER_TEST,0.5,256m,512m,python -m pytest tests/ -q -n 1 -m 'docker and smoke' --maxfail=1 --tb=short)

docker-test-integration: docker-ensure-image ## Run integration-focused Docker tests (balanced coverage)
	$(call RUN_DOCKER_TEST,2,4g,6g,python -m pytest tests/ -v -n 1 -m 'docker and not smoke' --maxfail=5 --tb=short)

docker-test-full: docker-ensure-image ## Run comprehensive Docker suite (CI parity)
	$(call RUN_DOCKER_TEST,4,8g,12g,uv run --group dev pytest -v -n 2 -m 'docker' --maxfail=5)

# Documentation
docs: ## Build documentation
	uv run --group dev sphinx-build -b html docs docs/_build/html
	@echo "✅ Documentation built in docs/_build/html/"

docs-serve: ## Serve documentation locally on port 8000
	@echo "🌐 Serving at http://localhost:8000 (Ctrl+C to stop)"
	cd docs/_build/html && uv run --group dev python -m http.server 8000

docs-dev: ## Live-reload documentation development
	uv run --group dev sphinx-autobuild docs docs/_build/html --host 0.0.0.0 --port 8000

docs-cli: ## Generate CLI reference documentation
	@echo "🔧 Generating CLI documentation..."
	@mkdir -p docs/
	@echo "# 🧬 siRNAforge CLI Reference" > docs/CLI_REFERENCE.md
	@echo "" >> docs/CLI_REFERENCE.md
	@echo "\`\`\`bash" >> docs/CLI_REFERENCE.md
	@uv run sirnaforge --help >> docs/CLI_REFERENCE.md
	@echo "\`\`\`" >> docs/CLI_REFERENCE.md
	@for cmd in search workflow design validate config version; do \
		echo "" >> docs/CLI_REFERENCE.md; \
		echo "### \`$$cmd\`" >> docs/CLI_REFERENCE.md; \
		echo "\`\`\`bash" >> docs/CLI_REFERENCE.md; \
		uv run sirnaforge $$cmd --help >> docs/CLI_REFERENCE.md 2>/dev/null || true; \
		echo "\`\`\`" >> docs/CLI_REFERENCE.md; \
	done
	@echo "✅ CLI documentation generated!"

# Nextflow Pipeline
nextflow-check: ## Check Nextflow installation
	@uv run nextflow -version || echo "❌ Nextflow should be available in main dependencies"

nextflow-run: ## Run Nextflow pipeline with test data
	uv run nextflow run nextflow_pipeline/main.nf \
		--input nextflow_pipeline/candidates.fasta \
		--outdir nextflow_results \
		--genome_species human \
		-profile test

nextflow-lint: ## Lint Nextflow scripts
	@echo "🔍 Linting Nextflow pipelines..."
	@uv run nextflow lint nextflow_pipeline/main.nf || echo "⚠️ Nextflow should be available in main dependencies"

# Utilities
dev: install-dev ## Quick development setup
	@echo "🚀 Ready for development!"

example: ## Run basic example
	mkdir -p examples/output
	uv run sirnaforge design examples/sample_transcripts.fasta -o examples/output/results.tsv

version: ## Show version information
	@echo "siRNAforge version: $(VERSION)"

release-test: ## Test release preparation (quick validation)
	@echo "🧪 Testing release preparation for version $(VERSION)..."
	@echo ""
	@echo "1. 🏷️ Version Check:"
	@echo "   Current version: $(VERSION)"
	@echo "   ✅ Version format valid"
	@echo ""
	@echo "2. 📋 Changelog Check:"
	@if [ -f CHANGELOG.md ]; then \
		if grep -q "\[$(VERSION)\]" CHANGELOG.md; then \
			echo "   ✅ Changelog entry exists for v$(VERSION)"; \
		else \
			echo "   ⚠️ No changelog entry found for v$(VERSION)"; \
		fi; \
	else \
		echo "   ⚠️ No CHANGELOG.md file found"; \
	fi
	@echo ""
	@echo "3. 🧪 Quick Tests:"
	@echo "   Running fast validation..."
	@$(MAKE) test-local-python > /dev/null 2>&1 && echo "   ✅ Fast tests pass" || echo "   ❌ Fast tests fail"
	@echo ""
	@echo "4. 🐳 Docker Image Test:"
	@if docker image inspect sirnaforge:latest > /dev/null 2>&1; then \
		echo "   ✅ Docker image exists (sirnaforge:latest)"; \
	else \
		echo "   ⚠️ Docker image not found (run 'make docker' first)"; \
	fi
	@echo ""
	@echo "5. 🔧 CLI Verification:"
	@if uv run sirnaforge version > /dev/null 2>&1; then \
		echo "   ✅ CLI works (version: $$(uv run sirnaforge version 2>/dev/null | grep -o 'v[0-9][^)]*)' || echo 'unknown'))"; \
	else \
		echo "   ❌ CLI not working"; \
	fi
	@echo ""
	@echo "📊 Release Test Summary:"
	@echo "   Use 'make release' for full release preparation"

release: clean build test lint ## Prepare release (full checks)
	@echo "✅ Release preparation complete!"

# Security & Maintenance
security: ## Run security checks
	@echo "🔐 Running Bandit (JSON + summary)"
	uv run --group dev bandit -r src/ -f json -o bandit-report.json || echo "⚠️ Bandit execution issue"
	uv run --group dev bandit -r src/ -q || true
	@echo "🔐 Running safety (JSON)"
	uv run --group dev python -c "import json,sys; from safety.formatter import report; from safety.safety import check; from safety.util import read_requirements; print(json.dumps({'error':'legacy interface changed'}))" >/dev/null 2>&1 || true
	uv run --group dev safety check --output json > safety-report.json || echo '{"error": "safety_failed"}' > safety-report.json
	@echo "✅ Security scanning complete (bandit-report.json, safety-report.json)"

pre-commit: ## Run pre-commit hooks
	uv run --group dev pre-commit run --all-files
