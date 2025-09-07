# Makefile for siRNAforge
# Uses uv for Python package management

.PHONY: help install install-dev test lint format build clean docker-build docker-run docker-dev

# Default target
help: ## Show this help message
	@echo "🧬 siRNAforge - Development Commands"
	@echo "===================================="
	@echo ""
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
	@echo ""
	@echo "🧪 Test Commands:"
	@echo "  test                 Run all tests"
	@echo "  test-unit            Run only unit tests"
	@echo "  test-integration     Run integration tests"
	@echo "  test-pipeline        Run pipeline tests (requires tools)"
	@echo "  test-docker          Run Docker integration tests"
	@echo "  test-fast            Run tests excluding slow ones"
	@echo "  test-all             Run all tests including slow ones"
	@echo "  test-cov             Run tests with coverage report"
	@echo ""
	@echo "📚 Documentation Commands:"
	@echo "  docs                 Build Sphinx documentation"
	@echo "  docs-serve           Serve documentation locally on :8000"
	@echo "  docs-dev             Build and serve docs with auto-reload"
	@echo "  docs-clean           Clean documentation build artifacts"
	@echo "  docs-check           Quick validation (for CI)"
	@echo "  docs-outdated        Check if generated docs need updating"
	@echo "  docs-cli             Generate CLI reference documentation"
	@echo "  docs-examples        Generate usage examples"
	@echo "  docs-full           Generate all documentation"
	@echo "  docs-full-ci        Generate all documentation (CI mode)"
	@echo "  docs-all            Complete docs rebuild (clean + generate + build)"
	@echo ""
	@echo "🐳 Docker Commands:"
	@echo "  docker-build         Build comprehensive Docker image"
	@echo "  docker-workflow      Run workflow in Docker (GENE=<name>)"
	@echo "  docker-dev           Interactive Docker development"
	@echo "  docker-compose-run   Run with docker-compose (GENE=<name>)"
	@echo ""
	@echo "🔄 Nextflow Commands:"
	@echo "  nextflow-check       Check Nextflow installation"
	@echo "  nextflow-run         Run pipeline with test data"
	@echo "  nextflow-run-docker  Run pipeline with Docker profile"
	@echo "  nextflow-preview     Preview pipeline execution (syntax check)"
	@echo "  lint-nextflow        Lint Nextflow scripts"
	@echo ""
	@echo "🎯 Quick Start:"
	@echo "  make install-dev     Setup development environment"
	@echo "  make docker-build    Build Docker image (recommended)"
	@echo "  make test            Run test suite"
	@echo "  make example-basic   Run basic example"
	@echo ""

# Environment setup
install: ## Install package for production
	uv sync --no-dev

install-dev: ## Install package with development dependencies (default groups)
	uv sync
	uv run pre-commit install
	@echo "✅ Development environment ready!"

install-all: ## Install all dependency groups
	uv sync --all-groups
	uv run pre-commit install
	@echo "✅ Complete environment with all groups ready!"

install-analysis: ## Install with analysis tools (Jupyter, plotting)
	uv sync --group analysis
	@echo "✅ Analysis environment ready!"

install-pipeline: ## Install with pipeline tools (Nextflow, Docker)
	uv sync --group pipeline
	@echo "✅ Pipeline environment ready!"

# Development
test: ## Run tests
	uv run pytest -v

test-unit: ## Run only unit tests
	uv run pytest -v -m "unit"

test-integration: ## Run integration tests
	uv run pytest -v -m "integration"

test-pipeline: ## Run pipeline tests (requires pipeline tools)
	uv run pytest -v -m "pipeline"

test-docker: ## Run Docker integration tests
	uv run pytest -v -m "docker"

test-fast: ## Run tests excluding slow ones
	uv run pytest -v -m "not slow"

test-all: ## Run all tests including slow ones
	uv run pytest -v

test-cov: ## Run tests with coverage
	uv run pytest --cov=sirnaforge --cov-report=html --cov-report=term-missing

lint: ## Run linting tools
	uv run ruff check src tests
	uv run black --check src tests
	uv run mypy src
	@echo "🔍 Running Nextflow lint..."
	@if command -v nextflow >/dev/null 2>&1; then \
		find nextflow_pipeline -name "*.nf" -exec nextflow lint {} \; 2>/dev/null || true; \
	else \
		uv run --group pipeline nextflow lint nextflow_pipeline/main.nf || echo "⚠️  Install Nextflow for pipeline linting: uv sync --group pipeline"; \
	fi

lint-nextflow: ## Run Nextflow-specific linting
	@echo "🔍 Linting Nextflow pipelines..."
	uv run --group pipeline nextflow lint nextflow_pipeline/main.nf
	@find nextflow_pipeline -name "*.nf" -not -path "*/modules/*" | while read -r file; do \
		echo "Linting $$file..."; \
		uv run --group pipeline nextflow lint "$$file" || true; \
	done

format: ## Format code
	uv run black src tests
	uv run ruff check --fix src tests
	@echo "✅ Code formatted!"

# Build
build: ## Build package
	uv build
	@echo "✅ Package built in dist/"

clean: ## Clean build artifacts
	rm -rf dist/ build/ src/*.egg-info/
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete
	rm -rf .pytest_cache/ .coverage htmlcov/ .mypy_cache/ .ruff_cache/
	@echo "✅ Cleaned build artifacts"

# Docker targets
docker-build: ## Build comprehensive Docker image with all dependencies
	docker build -f docker/Dockerfile -t sirnaforge:latest .
	@echo "✅ Docker image built: sirnaforge:latest"
	# Extract version from pyproject.toml (Python 3.11+ tomllib) and pass into docker build
	VERSION=$$(python -c "import tomllib,sys;data=tomllib.loads(open('pyproject.toml','rb').read());print(data.get('project',{}).get('version','0.0.0+local'))")
	docker build -f docker/Dockerfile --build-arg VERSION=$$VERSION -t sirnaforge:$$VERSION -t sirnaforge:latest .
	@echo "✅ Docker image built: sirnaforge:$$VERSION (also tagged as latest)"

docker-run: ## Run workflow in Docker container
	# Try to determine the versioned image tag, fall back to latest
	VERSION=$$(python -c "import tomllib,sys;data=tomllib.loads(open('pyproject.toml','rb').read());print(data.get('project',{}).get('version',''))")
	IMG=sirnaforge:latest
	if [ -n "$$VERSION" ]; then IMG=sirnaforge:$$VERSION; fi
	docker run -v $$(pwd):/workspace -w /workspace $$IMG sirnaforge --help
	
	IMG=sirnaforge:latest
	if [ -n "$$VERSION" ]; then IMG=sirnaforge:$$VERSION; fi
	docker run -v $$(pwd):/workspace -w /workspace $$IMG sirnaforge --help

docker-workflow: ## Run example workflow in Docker (usage: make docker-workflow GENE=TP53)
	@if [ -z "$(GENE)" ]; then echo "Usage: make docker-workflow GENE=<gene_name>"; exit 1; fi
	docker run -v $$(pwd):/workspace -w /workspace sirnaforge:latest \
		sirnaforge workflow $(GENE) --output-dir docker_results

docker-dev: ## Run development Docker image with bash
	docker run -it -v $$(pwd):/workspace -w /workspace sirnaforge:latest bash

docker-compose-up: ## Start services with docker-compose
	docker-compose up --build

docker-compose-run: ## Run workflow with docker-compose (usage: make docker-compose-run GENE=TP53)
	@if [ -z "$(GENE)" ]; then echo "Usage: make docker-compose-run GENE=<gene_name>"; exit 1; fi
	docker-compose run --rm sirnaforge sirnaforge workflow $(GENE) --output-dir compose_results

docker-test: ## Run tests in Docker container
	docker run -v $$(pwd):/workspace -w /workspace sirnaforge:latest \
		uv run pytest tests/

# Nextflow Pipeline Development
nextflow-check: ## Check Nextflow installation and version
	@echo "🔍 Checking Nextflow installation..."
	@uv run --group pipeline nextflow -version || echo "❌ Nextflow not found. Run: make install-pipeline"

nextflow-run: ## Run Nextflow pipeline with test data
	uv run --group pipeline nextflow run nextflow_pipeline/main.nf \
		--input nextflow_pipeline/candidates.fasta \
		--outdir nextflow_results \
		--genome_species human \
		-profile test

nextflow-run-docker: ## Run Nextflow pipeline with Docker profile  
	uv run --group pipeline nextflow run nextflow_pipeline/main.nf \
		--input nextflow_pipeline/candidates.fasta \
		--outdir nextflow_results \
		--genome_species human \
		-profile docker

nextflow-preview: ## Preview pipeline execution (syntax check)
	uv run --group pipeline nextflow run nextflow_pipeline/main.nf \
		--input nextflow_pipeline/candidates.fasta \
		--outdir preview_results \
		--genome_species human \
		-profile test \
		-preview

nextflow-clean: ## Clean Nextflow work directory
	rm -rf work/ .nextflow* nextflow_results/ preview_results/
	@echo "✅ Nextflow work directory cleaned"

nextflow-report: ## Generate Nextflow execution report
	uv run --group pipeline nextflow run nextflow_pipeline/main.nf \
		--input nextflow_pipeline/candidates.fasta \
		--outdir nextflow_results \
		-with-report nextflow_report.html \
		-with-timeline timeline.html \
		-with-dag flowchart.html

nextflow-test: ## Test Nextflow integration
	@bash scripts/test_nextflow_integration.sh

# Documentation
docs: ## Build Sphinx documentation
	uv sync --group docs
	uv run sphinx-build -b html docs docs/_build/html
	@echo "✅ Documentation built in docs/_build/html/"

docs-ci: ## Build documentation for CI/CD (with stricter error handling)
	uv sync --group docs
	uv run sphinx-build -W --keep-going -b html docs docs/_build/html
	@echo "✅ Documentation built for CI/CD in docs/_build/html/"

docs-serve: ## Serve documentation locally
	uv sync --group docs
	@echo "🌐 Serving documentation at http://localhost:8000"
	@echo "Press Ctrl+C to stop the server"
	cd docs/_build/html && uv run python -m http.server 8000

docs-dev: ## Build docs and serve with auto-reload (requires sphinx-autobuild)
	@echo "🔄 Starting documentation development server with auto-reload..."
	uv sync --group docs
	uv run sphinx-autobuild docs docs/_build/html --host 0.0.0.0 --port 8000

docs-clean: ## Clean documentation build artifacts
	rm -rf docs/_build/
	@echo "✅ Documentation build artifacts cleaned"

docs-check: ## Quick documentation validation (for CI lint phase)
	@echo "🔍 Validating documentation..."
	uv sync --group docs
	# Generate CLI docs to ensure CLI is working
	@mkdir -p docs/
	@uv run sirnaforge --help > /dev/null || (echo "❌ CLI not working" && exit 1)
	# Quick syntax check of Sphinx docs (allowing warnings for quick validation)
	uv run sphinx-build -q -E -b html docs docs/_build/html
	@echo "✅ Documentation validation passed!"

docs-outdated: ## Check if generated docs are outdated
	@echo "🔍 Checking if generated docs are outdated..."
	@if [ ! -f docs/CLI_REFERENCE.md ]; then \
		echo "❌ CLI_REFERENCE.md missing - run 'make docs-cli'"; \
		exit 1; \
	fi
	@if [ ! -f docs/USAGE_EXAMPLES.md ]; then \
		echo "❌ USAGE_EXAMPLES.md missing - run 'make docs-examples'"; \
		exit 1; \
	fi
	@echo "✅ Generated docs are present!"

docs-cli: ## Generate CLI reference documentation
	@echo "🔧 Generating CLI documentation..."
	@mkdir -p docs/
	@echo "# 🧬 siRNAforge CLI Reference" > docs/CLI_REFERENCE.md
	@echo "" >> docs/CLI_REFERENCE.md
	@echo "Complete command-line interface documentation for siRNAforge — Comprehensive siRNA design toolkit for gene silencing." >> docs/CLI_REFERENCE.md
	@echo "" >> docs/CLI_REFERENCE.md
	@echo "## Main Commands" >> docs/CLI_REFERENCE.md
	@echo "" >> docs/CLI_REFERENCE.md
	@echo "\`\`\`bash" >> docs/CLI_REFERENCE.md
	@uv run sirnaforge --help >> docs/CLI_REFERENCE.md
	@echo "\`\`\`" >> docs/CLI_REFERENCE.md
	@echo "" >> docs/CLI_REFERENCE.md
	@echo "## Command Details" >> docs/CLI_REFERENCE.md
	@for cmd in search workflow design validate config version; do \
		echo "" >> docs/CLI_REFERENCE.md; \
		echo "### \`$$cmd\`" >> docs/CLI_REFERENCE.md; \
		echo "" >> docs/CLI_REFERENCE.md; \
		echo "\`\`\`bash" >> docs/CLI_REFERENCE.md; \
		uv run sirnaforge $$cmd --help >> docs/CLI_REFERENCE.md; \
		echo "\`\`\`" >> docs/CLI_REFERENCE.md; \
	done
	@echo "✅ CLI documentation generated in docs/CLI_REFERENCE.md"

docs-examples: ## Generate usage examples
	@echo "🧬 Generating usage examples..."
	@mkdir -p docs/examples/
	@echo "# siRNAforge Usage Examples" > docs/examples/USAGE_EXAMPLES.md
	@echo "" >> docs/examples/USAGE_EXAMPLES.md
	@echo "## Basic Examples" >> docs/examples/USAGE_EXAMPLES.md
	@echo "" >> docs/examples/USAGE_EXAMPLES.md
	@echo "\`\`\`bash" >> docs/examples/USAGE_EXAMPLES.md
	@echo "# Complete workflow" >> docs/examples/USAGE_EXAMPLES.md
	@echo "uv run sirnaforge workflow TP53 --output-dir results" >> docs/examples/USAGE_EXAMPLES.md
	@echo "" >> docs/examples/USAGE_EXAMPLES.md
	@echo "# Design from FASTA" >> docs/examples/USAGE_EXAMPLES.md
	@echo "uv run sirnaforge design transcripts.fasta -o results.tsv" >> docs/examples/USAGE_EXAMPLES.md
	@echo "" >> docs/examples/USAGE_EXAMPLES.md
	@echo "# Search for transcripts" >> docs/examples/USAGE_EXAMPLES.md
	@echo "uv run sirnaforge search BRCA1 --all --verbose" >> docs/examples/USAGE_EXAMPLES.md
	@echo "\`\`\`" >> docs/examples/USAGE_EXAMPLES.md
	@echo "✅ Usage examples generated in docs/examples/"

docs-full: docs-cli docs-examples docs-api ## Generate complete documentation set
	@echo "📚 Generating complete documentation..."
	@mkdir -p docs/
	@echo "# siRNAforge Documentation" > docs/README.md
	@echo "" >> docs/README.md
	@echo "Complete documentation for siRNAforge — Comprehensive siRNA design toolkit for gene silencing." >> docs/README.md
	@echo "" >> docs/README.md
	@echo "## Available Documentation" >> docs/README.md
	@echo "" >> docs/README.md
	@echo "- [Getting Started](getting_started.md) - Installation and first steps" >> docs/README.md
	@echo "- [CLI Reference](CLI_REFERENCE.md) - Complete command documentation" >> docs/README.md
	@echo "- [Usage Examples](USAGE_EXAMPLES.md) - Real-world usage patterns" >> docs/README.md
	@echo "- [Gene Search Guide](gene_search.md) - Gene search functionality" >> docs/README.md
	@echo "- [Architecture](architecture.md) - System design and components" >> docs/README.md
	@echo "- [API Reference](api_reference.rst) - Python API documentation" >> docs/README.md
	@echo "- [Tutorials](tutorials/index.md) - Step-by-step guides" >> docs/README.md
	@echo "- [Development Guide](development.md) - Contributing and development" >> docs/README.md
	@echo "" >> docs/README.md
	@echo "Generated on: $$(date)" >> docs/README.md
	@echo "siRNAforge version: $$(uv run sirnaforge version 2>/dev/null || echo 'Development')" >> docs/README.md
	uv sync --group docs
	uv run sphinx-build -b html docs docs/_build/html
	@echo "✅ Complete documentation generated in docs/ and docs/_build/html/"

docs-full-ci: docs-cli docs-examples docs-api ## Generate complete documentation for CI/CD with strict error handling
	@echo "📚 Generating complete documentation for CI/CD..."
	@mkdir -p docs/
	@echo "# siRNAforge Documentation" > docs/README.md
	@echo "" >> docs/README.md
	@echo "Complete documentation for siRNAforge — Comprehensive siRNA design toolkit for gene silencing." >> docs/README.md
	@echo "" >> docs/README.md
	@echo "## Available Documentation" >> docs/README.md
	@echo "" >> docs/README.md
	@echo "- [Getting Started](getting_started.md) - Installation and first steps" >> docs/README.md
	@echo "- [CLI Reference](CLI_REFERENCE.md) - Complete command documentation" >> docs/README.md
	@echo "- [Usage Examples](USAGE_EXAMPLES.md) - Real-world usage patterns" >> docs/README.md
	@echo "- [Gene Search Guide](gene_search.md) - Gene search functionality" >> docs/README.md
	@echo "- [Architecture](architecture.md) - System design and components" >> docs/README.md
	@echo "- [API Reference](api_reference.rst) - Python API documentation" >> docs/README.md
	@echo "- [Tutorials](tutorials/index.md) - Step-by-step guides" >> docs/README.md
	@echo "- [Development Guide](development.md) - Contributing and development" >> docs/README.md
	@echo "" >> docs/README.md
	@echo "Generated on: $$(date)" >> docs/README.md
	@echo "siRNAforge version: $$(uv run sirnaforge version 2>/dev/null || echo 'Development')" >> docs/README.md
	uv sync --group docs
	uv run sphinx-build -W --keep-going -b html docs docs/_build/html
	@echo "✅ Complete documentation generated for CI/CD in docs/ and docs/_build/html/"

docs-api: ## Generate API documentation from docstrings
	@echo "🔧 Generating API documentation..."
	uv sync --group docs
	@echo "✅ API documentation updated in docs/api_reference.rst"

docs-all: docs-clean docs-full ## Complete docs rebuild (clean + generate + build)
	@echo "🚀 Complete documentation rebuild finished!"

# Utilities
check: lint test ## Run all checks (lint + test)

release-check: clean build test lint ## Full release preparation check
	@echo "✅ Release checks passed!"

# Release management
prepare-release: ## Prepare for release (run tests, build, check)
	@echo "🚀 Preparing release..."
	make clean
	make lint
	make test-fast
	make build
	make docs-check
	@echo "✅ Release preparation complete!"

release-check: clean build test lint ## Full release preparation check
	@echo "✅ Release checks passed!"

# Generate test data for CI/Docker
generate-test-data: ## Generate minimal test data for CI/Docker testing
	@echo "📦 Generating test data..."
	./scripts/generate_test_data.sh

version: ## Show version information
	uv run python -c "from sirnaforge import __version__; print(f'Version: {__version__}')"

# Quick development workflow
dev: install-dev ## Quick setup for development

# Pipeline testing
test-nextflow: ## Test Nextflow pipeline (requires Nextflow)
	nextflow run pipeline/main.nf -profile test

# Examples
example-basic: ## Run basic example
	mkdir -p examples/output
	uv run sirnaforge design examples/sample_transcripts.fasta -o examples/output/results.tsv

example-advanced: ## Run advanced example with all features
	mkdir -p examples/output
	uv run sirnaforge design examples/sample_transcripts.fasta \
		-o examples/output/advanced_results.tsv \
		--gc-min 35 --gc-max 55 \
		--top-n 20 \
		--verbose

# Maintenance
update-deps: ## Update dependencies
	uv pip compile pyproject.toml --upgrade
	@echo "✅ Dependencies updated!"

security-check: ## Run security checks
	uv run bandit -r src/
	uv run safety check

pre-commit: ## Run pre-commit hooks
	uv run pre-commit run --all-files

pre-commit-docs: ## Update docs and run pre-commit (ensures docs are current)
	make docs-cli docs-examples
	uv run pre-commit run --all-files

# CI/CD helpers
ci-setup: ## Setup for CI/CD environment
	uv pip install -e .[dev]

ci-test: ## Run tests in CI/CD format
	uv run pytest --junitxml=junit.xml --cov=sirnaforge --cov-report=xml

ci-test-docker: ## Run Docker tests in CI
	docker run --rm -v $$(pwd):/workspace -w /workspace sirnaforge:latest \
		uv run pytest -v -m "docker" --junitxml=junit-docker.xml

# Local registry
registry-run: ## Run local Docker registry
	docker run -d -p 5000:5000 --name registry registry:2

registry-push: ## Push to local registry
	docker tag sirnaforge:latest localhost:5000/sirnaforge:latest
	docker push localhost:5000/sirnaforge:latest
