# 🧬 siRNAforge Test Organization

This document outlines the three main groups of tests in siRNAforge, designed to support different development and deployment scenarios.

## Test Groups Overview

### 1. Unit Tests (`unit`, `local_python`, `ci`)
**Purpose**: Fast, Python-only validation of individual components
**Environment**: Local Python development, CI/CD pipelines
**Requirements**: Python dependencies only (no Docker/Nextflow)
**Location**: `tests/unit/`

**What they test**:
- Model validation and data structures
- CLI command functionality
- Individual algorithm components
- Fast-running validation checks
- Gene search and database interactions
- ORF analysis and sequence processing

**Run commands**:
```bash
# Fastest Make shortcut focused on dev markers
make test-dev

# Full unit marker set via pytest
uv run pytest -m "unit"

# CI-optimized
make test-ci
# or
uv run pytest -m "ci"
```

### 2. Integration Tests (`integration`, `docker`)
**Purpose**: End-to-end workflow validation
**Environment**: Docker containers with bioinformatics tools
**Requirements**: Docker + Nextflow + bioinformatics tools (samtools, bwa-mem2, RNAfold)
**Location**: `tests/docker/`, `tests/pipeline/`, `tests/integration/`

**What they test**:
- Complete workflow from gene input to off-target analysis
- Docker container functionality
- Nextflow pipeline integration
- Bioinformatics tool integration
- Multi-species off-target analysis
- Full siRNA design pipeline

**Run commands**:
```bash
# Container validation suite (runs tests/container)
make docker-test

# Marker-based integration tests on host
uv run pytest -m "integration"
```

### 3. Local Nextflow Tests (`local_nextflow`)
**Purpose**: Development testing with Nextflow workflows
**Environment**: Local machine with Nextflow installed
**Requirements**: Nextflow + bioinformatics tools
**Location**: `tests/pipeline/`, integration scripts

**What they test**:
- Nextflow pipeline execution
- Off-target analysis workflows
- Pipeline configuration and validation
- Local development workflow testing

**Run commands**:
```bash
# Marker-based Nextflow subset (requires local Nextflow install)
uv run pytest -m "local_nextflow"

# Host target that isolates Nextflow requirements
make test-requires-nextflow

# Full integration script
./tests/integration/test_workflow_integration.sh TP53
```

## Test Markers Reference

| Marker | Purpose | Environment |
|--------|---------|-------------|
| `unit` | Individual component tests | Python-only |
| `integration` | Full workflow tests | Docker + Nextflow |
| `docker` | Docker container tests | Container environment |
| `pipeline` | Pipeline-specific tests | Nextflow environment |
| `local_python` | Local Python development | Development machine |
| `local_nextflow` | Local Nextflow development | Development machine |
| `ci` | CI/CD optimized tests | Automated environments |
| `lightweight` | Fast, minimal resource tests | Any environment |
| `slow` | Resource-intensive tests | High-resource environments |

## Development Workflow

### For Python Development
```bash
# Quick validation
make test-dev

# Full local testing
make test
```

### For Pipeline Development
```bash
# Test Nextflow-heavy markers on host
make test-requires-nextflow

# Full Docker integration
make docker-test
```

### For CI/CD
```bash
# Fast CI tests
make test-ci

# Release validation pipeline
make test-release
```

## Resource Requirements

### Unit Tests
- **CPU**: Minimal (1 core)
- **Memory**: < 512MB
- **Time**: < 30 seconds
- **Dependencies**: Python packages only

### Integration Tests
- **CPU**: 2-4 cores
- **Memory**: 2-8GB
- **Time**: 2-15 minutes
- **Dependencies**: Docker, Nextflow, bioinformatics tools

### Local Nextflow Tests
- **CPU**: 2+ cores
- **Memory**: 4+ GB
- **Time**: 5-30 minutes
- **Dependencies**: Nextflow, bioinformatics tools

## File Organization

```
tests/
├── unit/                    # Unit tests (Python-only)
│   ├── test_models.py      # Data model validation
│   ├── test_gene_search.py # Gene search functionality
│   ├── test_fixed_vienna.py # Thermodynamics
│   └── ...
├── docker/                 # Docker integration tests
│   ├── test_container_integration.py
│   └── test_profiles.py
├── pipeline/               # Pipeline integration tests
│   └── test_embedded_pipeline.py
└── integration/            # Full workflow scripts
    └── test_workflow_integration.sh
```</content>
<parameter name="filePath">/home/hovland/sirnaforge/sirnaforge/tests/README.md
