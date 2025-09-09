# Makefile for siRNAforge
# Uses uv for Python package management

.PHONY: help install install-dev install-pipeline test lint format build clean docker docs

# Variables
DOCKER_IMAGE = sirnaforge
VERSION = $(shell uv run python -c "from sirnaforge import __version__; print(__version__)" 2>/dev/null || echo "0.1.0")

# Default target
help: ## Show available commands
	@echo "🧬 siRNAforge Development Commands"
	@echo "================================="
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

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

# Development & Testing
test: ## Run all tests
	uv run pytest -v

test-unit: ## Run unit tests only
	uv run pytest -v -m "unit"

test-integration: ## Run integration tests
	uv run pytest -v -m "integration"

test-fast: ## Run tests excluding slow ones
	uv run pytest -v -m "not slow"

test-cov: ## Run tests with coverage report
	uv run pytest --cov=sirnaforge --cov-report=html --cov-report=term-missing

# Code Quality
lint: ## Run all linting tools
	uv run ruff check src tests
	uv run black --check src tests
	uv run mypy src
	@echo "✅ Code quality checks passed!"

format: ## Format code with black and ruff
	uv run black src tests
	uv run ruff check --fix src tests
	@echo "✅ Code formatted!"

check: lint test-fast ## Run quick quality checks (lint + fast tests)

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

docker-run: GENE ?= TP53
docker-run: ## Run workflow in Docker (usage: make docker-run GENE=<gene>)
	docker run -v $$(pwd):/workspace -w /workspace $(DOCKER_IMAGE):latest \
		sirnaforge workflow $(GENE) --output-dir docker_results

docker-dev: ## Interactive Docker development shell
	docker run -it -v $$(pwd):/workspace -w /workspace $(DOCKER_IMAGE):latest bash

docker-test: ## Run tests in Docker
	docker run --rm -v $$(pwd):/workspace -w /workspace $(DOCKER_IMAGE):latest \
		uv run pytest -v

# Documentation
docs: ## Build documentation
	uv sync --group dev
	uv run sphinx-build -b html docs docs/_build/html
	@echo "✅ Documentation built in docs/_build/html/"

docs-serve: ## Serve documentation locally on port 8000
	uv sync --group dev
	@echo "🌐 Serving at http://localhost:8000 (Ctrl+C to stop)"
	cd docs/_build/html && uv run python -m http.server 8000

docs-dev: ## Live-reload documentation development
	uv sync --group dev
	uv run sphinx-autobuild docs docs/_build/html --host 0.0.0.0 --port 8000

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

release: clean build test lint ## Prepare release (full checks)
	@echo "✅ Release preparation complete!"

# Security & Maintenance
security: ## Run security checks
	uv run bandit -r src/ || echo "⚠️ Install dev dependencies: make install-dev"
	uv run safety check || echo "⚠️ Install dev dependencies: make install-dev"

pre-commit: ## Run pre-commit hooks
	uv run pre-commit run --all-files
