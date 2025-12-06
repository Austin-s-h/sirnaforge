# 🧬 siRNAforge Test Organization

This document outlines the three main groups of tests in siRNAforge, designed to support different development and deployment scenarios.

## Test Output & Cleanup Behavior

### Container Test Persistence (for Failure Inspection)

All container tests use a **dual-strategy** for output management:

**When running in Docker** (`make docker-test`, `make docker-build-test`):
- Output saved to `/workspace/workflow_test_debug_*` directories
- These are mounted from host via `-v $(pwd):/workspace`
- **Persists after test completion** for manual inspection
- Clean with: `make docker-build-test` (removes all `workflow_test_debug_*` dirs)

**When running on host** (`pytest tests/container/`):
- Output saved to pytest's `tmp_path` (`.pytest_tmp/`)
- **Auto-cleanup**: Removed after successful tests
- **Auto-retention**: pytest keeps last 3 failed test directories
- Temporary directories include full output for debugging failures

**Failure inspection workflow:**
```bash
# After a container test failure:
ls -la workflow_test_debug_*/       # Check persistent output
cat workflow_test_debug_*/logs/*.log  # View logs

# Clean before next run:
make docker-build-test  # Cleans all debug dirs + rebuilds + tests
```

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
├── container/              # Docker container integration tests
│   ├── test_container_integration.py  # Basic container validation
│   └── test_workflow_modes.py         # Comprehensive workflow coverage
├── docker/                 # Docker integration tests
│   ├── test_container_integration.py
│   └── test_profiles.py
├── pipeline/               # Pipeline integration tests
│   └── test_embedded_pipeline.py
└── integration/            # Full workflow scripts
    └── test_workflow_integration.sh
```

## Container Test Coverage

The `tests/container/` directory contains comprehensive end-to-end tests that validate different workflow operation modes inside Docker:

### `test_container_integration.py` - Basic Container Validation
- CLI help and version commands
- Python environment and package availability
- Bioinformatics tool availability (BWA-MEM2, SAMtools, ViennaRNA)
- Minimal workflow smoke test
- Full TP53 workflow with real data

### `test_workflow_modes.py` - Comprehensive Workflow Coverage
Tests all major workflow configurations non-redundantly with existing unit tests and test_container_integration.py:

1. **Full Gene Search Workflow** (`test_full_workflow_with_gene_search`)
   - Gene symbol (ACTB) → transcript retrieval → ORF validation → design → off-target
   - Validates Ensembl API integration
   - Checks complete output structure
   - **Non-redundant**: Uses gene search (vs TP53 test which uses --input-fasta)

2. **Custom Transcriptome Off-Target** (`test_custom_transcriptome_offtarget`)
   - Uses `--input-fasta` + `--transcriptome-fasta`
   - Validates transcriptome indexing
   - Verifies off-target analysis runs

3. **Genome Index Override** (`test_genome_index_override`)
   - Tests `--offtarget-indices` flag
   - Validates custom genome reference paths
   - Checks species derivation from overrides

4. **miRNA Design Mode** (`test_mirna_design_mode`)
   - Uses `--design-mode mirna`
   - Validates miRNA-specific scoring columns
   - Checks position-1 nucleotide preference

5. **Chemical Modifications** (`test_modification_pattern_application`)
   - Tests `--modifications` and `--overhang` flags
   - Validates modification annotations in FASTA headers
   - Checks manifest includes modification info

6. **Multi-Species Off-Target** (`test_multi_species_offtarget`)
   - Tests `--species` with multiple values
   - Validates cross-species miRNA seed checks
   - Verifies off-target across genomes

7. **GC Content Filtering** (`test_gc_range_filtering`)
   - Tests `--gc-min` and `--gc-max` filters
   - Validates candidate filtering logic
   - Checks QC reason documentation

8. **Toy Database Workflow** (`test_minimal_toy_workflow`)
    - Fast sanity check using toy databases
    - Validates basic pipeline flow
    - Quick validation for container builds

**Note**: Design-only mode (--input-fasta without transcriptome) is already comprehensively tested in `test_container_integration.py::test_docker_full_tp53_workflow`, so it's intentionally not duplicated here. Output directory structure validation is also covered by the TP53 test's `_verify_workflow_outputs()` helper.

**Run container tests**:
```bash
# Build and test container
make docker-build-test

# Run container tests directly
make docker-test

# Run specific workflow mode test
uv run pytest tests/container/test_workflow_modes.py::test_design_only_mode_no_offtarget -v
```
