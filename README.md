# 🧬 siRNAforge — Comprehensive siRNA Design Tool

<div align="center">
  <img src="docs/branding/sirnaforge_logo_3.svg" alt="siRNAforge Logo" width="200"/>

  **Multi-species gene to siRNA design, off-target prediction, and ranking**
  [![🚀 Release siRNAforge](https://github.com/austin-s-h/sirnaforge/actions/workflows/release.yml/badge.svg?branch=master)](https://github.com/austin-s-h/sirnaforge/actions/workflows/release.yml)
  [![Python 3.9–3.12](https://img.shields.io/badge/python-3.9--3.12-blue.svg)](https://www.python.org/downloads/)
  [![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
  [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
  [![Docker](https://img.shields.io/badge/docker-available-blue?logo=docker)](https://github.com/users/Austin-s-h/packages/container/package/sirnaforge)
  [![Nextflow](https://img.shields.io/badge/nextflow-pipeline-brightgreen?logo=nextflow)](https://github.com/austin-s-h/sirnaforge/tree/main/nextflow_pipeline)
  [![Tests](https://img.shields.io/badge/tests-passing-brightgreen?logo=pytest)](https://github.com/austin-s-h/sirnaforge/actions)
  [![Coverage](https://img.shields.io/badge/coverage-%3E90%25-brightgreen?logo=codecov)](https://github.com/austin-s-h/sirnaforge)
  [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/austin-s-h/sirnaforge/blob/main/LICENSE)
</div>


siRNAforge is a modern, comprehensive toolkit for designing high-quality siRNAs with integrated off-target analysis. Built with Python 3.9-3.12, it combines cutting-edge bioinformatics algorithms with robust software engineering practices to provide a complete gene silencing solution for researchers and biotechnology applications.

## ✨ Key Features

- 🎯 **Algorithm-driven design** - Comprehensive siRNA design with multi-component thermodynamic scoring
- 🔍 **Multi-species off-target analysis** - BWA-MEM2 alignment (transcriptome + miRNA seed modes) across human, rat, rhesus genomes
- 📊 **Advanced scoring system** - Composite scoring with seed-region specificity and secondary structure prediction
- 🧪 **ViennaRNA integration** - Secondary structure prediction for enhanced design accuracy
- 🧬 **Chemical modifications metadata** - Track 2'-O-methyl, 2'-fluoro, PS linkages, overhangs, and provenance
- 🔬 **Nextflow pipeline integration** - Scalable, containerized workflow execution with automatic parallelization
- 🐍 **Modern Python architecture** - Type-safe code with Pydantic models, async/await support, and rich CLI
- ⚡ **Lightning-fast dependency management** - Built with `uv` for sub-second installs and virtual environment management
- 🐳 **Fully containerized** - Docker images with all bioinformatics dependencies pre-installed
- 🧬 **Multi-database support** - Ensembl, RefSeq, GENCODE integration for comprehensive transcript retrieval

> **Note:** Supports Python 3.9-3.12. Python 3.13+ not yet supported due to ViennaRNA dependency constraints.

## 🚀 Quick Start

### Installation Options

**🐳 Docker (Recommended - Complete Environment):**
```bash
# Pull the pre-built image with all dependencies
docker pull ghcr.io/austin-s-h/sirnaforge:latest

# Quick workflow example
docker run -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir results --genome-species human

# With custom parameters
docker run -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow BRCA1 --gc-min 40 --gc-max 60 --sirna-length 21 --top-n 50
```

**🐍 Conda Environment (Alternative - Local Development):**
```bash
# Install micromamba (recommended - fastest), Mambaforge, or Miniconda
# micromamba (fastest option):
curl -LsSf https://micro.mamba.pm/install.sh | bash

# Or Mambaforge:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# Create siRNAforge development environment
make conda-env

# Activate the environment
micromamba activate sirnaforge-dev  # or conda activate sirnaforge-dev

# Install Python dependencies
make install-dev

# Run tests to verify installation
make test-local-python
```

**🖥️ Local Development Installation:**
```bash
# Install uv (lightning-fast Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup with development dependencies
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make install-dev

# Run sanity checks to verify installation
make test-local-python
```

### Essential Dependencies for Off-target Analysis

The Docker image includes all bioinformatics dependencies via conda environment (`docker/environment-nextflow.yml`):

- ✅ **Nextflow** (≥25.04.0) - Workflow orchestration and parallelization
- ✅ **BWA-MEM2** (≥2.2.1) - High-performance genome alignment (transcriptome + miRNA seed analysis)
- ✅ **SAMtools** (≥1.19.2) - SAM/BAM file processing and indexing
- ✅ **ViennaRNA** (≥2.7.0) - RNA secondary structure prediction
- ✅ **AWS CLI** (≥2.0) - Automated genome reference downloads
- ✅ **Java 17** - Nextflow runtime environment

**For local development without Docker:**
```bash
# Option 1: Use conda environment (includes all tools)
make conda-env
micromamba activate sirnaforge-dev  # or conda activate sirnaforge-dev

# Option 2: Install bioinformatics tools via micromamba
curl -LsSf https://micro.mamba.pm/install.sh | bash
micromamba env create -f docker/environment-nextflow.yml
micromamba activate sirnaforge-env
```

### Usage Examples

**🎯 Complete Workflow (Gene Query to Results):**
```bash
# Basic workflow with default parameters
uv run sirnaforge workflow TP53 --output-dir results

# Advanced workflow with custom parameters
uv run sirnaforge workflow BRCA1 \
  --genome-species "human,rat,rhesus" \
  --gc-min 40 --gc-max 60 \
  --sirna-length 21 \
  --top-n 50 \
  --output-dir brca1_analysis

# Workflow from a pre-existing FASTA file (local path or remote URL)
uv run sirnaforge workflow --input-fasta transcripts.fasta \
  --output-dir custom_analysis \
  --offtarget-n 25 \
  custom_gene_name

# Remote FASTA example
uv run sirnaforge workflow --input-fasta https://example.org/transcripts.fasta \
  --output-dir remote_input_run \
  remote_dataset
```

**🔍 Individual Component Usage:**
```bash
# Search for gene transcripts across databases
uv run sirnaforge search TP53 --output transcripts.fasta --database ensembl

# Design siRNAs from transcript sequences
uv run sirnaforge design transcripts.fasta --output results.csv --top-n 20

# Validate input files before processing
uv run sirnaforge validate candidates.fasta

# Display configuration and system information
uv run sirnaforge config

# Show detailed help for any command
uv run sirnaforge --help
uv run sirnaforge workflow --help
```

### Python API

**🔧 Programmatic Access for Custom Workflows:**
```python
import asyncio
from pathlib import Path
from sirnaforge.workflow import run_sirna_workflow
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters, FilterCriteria
from sirnaforge.data.gene_search import search_gene_sync

# Complete async workflow with custom parameters
async def design_sirnas_custom():
    results = await run_sirna_workflow(
        gene_query="TP53",
        output_dir="results",
        database="ensembl",
        top_n_candidates=50,
        top_n_offtarget=15,
        genome_species=["human", "rat", "rhesus"],
        gc_min=40.0,
        gc_max=60.0,
        sirna_length=21,
    )
    return results

# Run the workflow
results = asyncio.run(design_sirnas_custom())
print(f"✅ Designed {len(results.get('top_candidates', []))} siRNA candidates")

# Individual component usage for custom pipelines
def custom_design_pipeline():
    # 1. Search for gene transcripts
    transcripts = search_gene_sync(
        gene_query="BRCA1",
        database="ensembl",
        output_file="transcripts.fasta"
    )

    # 2. Configure design parameters
    design_params = DesignParameters(
        sirna_length=21,
        filters=FilterCriteria(
            gc_min=40,
            gc_max=60,
            avoid_patterns=["AAAA", "TTTT", "GGGG", "CCCC"]
        )
    )

    # 3. Initialize designer and generate candidates
    designer = SiRNADesigner(design_params)
    design_results = designer.design_from_file("transcripts.fasta")

    # 4. Process results
    for candidate in design_results.top_candidates[:10]:
        print(f"Candidate {candidate.id}:")
        print(f"  Guide: {candidate.guide_sequence}")
        print(f"  Score: {candidate.composite_score:.2f}")
        print(f"  GC%: {candidate.gc_content:.1f}")
        print(f"  Transcripts: {len(candidate.transcript_ids)}")
        print()

    return design_results

# Example: Batch processing multiple genes
async def batch_design_genes(genes: list[str]):
    results = {}
    for gene in genes:
        print(f"Processing {gene}...")
        gene_results = await run_sirna_workflow(
            gene_query=gene,
            output_dir=f"results_{gene.lower()}",
            top_n_candidates=20
        )
        results[gene] = gene_results
    return results

# Process multiple cancer-related genes
cancer_genes = ["TP53", "BRCA1", "BRCA2", "EGFR", "MYC"]
batch_results = asyncio.run(batch_design_genes(cancer_genes))
```

## 🏗️ Architecture & Workflow

### Complete Pipeline Overview

```
Gene Query → Transcript Search → ORF Validation → siRNA Design → Off-target Analysis → Ranked Results
     ↓              ↓                ↓               ↓               ↓                    ↓
Multi-database   Canonical       Coding Frame   Thermodynamic   Multi-species BWA    Scored & Filtered
Gene Search      Isoform         Validation     + Structure     Alignment (seed &    siRNA Candidates
(Ensembl/        Selection                      Scoring         transcriptome)       with Off-target
RefSeq/GENCODE)                                                                    Predictions
```

### Core Components

**🔍 Gene Search & Data Layer** (`sirnaforge.data.*`)
- **Multi-database integration**: Ensembl, RefSeq, GENCODE APIs with automatic fallback
- **Canonical transcript selection**: Prioritizes protein-coding, longest transcripts
- **Robust error handling**: Network timeouts, API rate limiting, malformed responses
- **Async/await support**: Non-blocking I/O for improved performance

**🧬 ORF Analysis** (`sirnaforge.data.orf_analysis`)
- **Reading frame validation**: Ensures proper coding sequence targeting
- **Quality control reporting**: Detailed validation logs and metrics
- **Multi-transcript support**: Handles gene isoforms and splice variants

**🎯 siRNA Design Engine** (`sirnaforge.core.design`)
- **Algorithm-based candidate generation**: Systematic 19-23 nucleotide window scanning
- **Multi-component scoring system**:
  - **Thermodynamic properties**: GC content (30-60%), melting temperature optimization
  - **Secondary structure prediction**: ViennaRNA integration for accessibility scoring
  - **Position-specific penalties**: 5' and 3' end optimization
  - **Off-target risk assessment**: Simplified seed-region analysis
- **Composite scoring**: Weighted combination of all scoring components
- **Transcript consolidation**: Deduplicates guide sequences across multiple transcript isoforms

- **🔍 Off-target Analysis** (`sirnaforge.core.off_target`)
  - **Adaptive BWA-MEM2 modes**: Sensitive genome-wide alignment plus ultra-short miRNA seed analysis using tuned parameters
- **Multi-species support**: Human, rat, rhesus macaque genome analysis
- **Advanced scoring**: Position-weighted mismatch penalties with seed-region emphasis
- **Scalable processing**: Batch candidate analysis with parallel execution

**🔬 Nextflow Pipeline Integration** (`nextflow_pipeline/`)
- **Containerized execution**: Docker/Singularity support with pre-built environments
- **Automatic resource management**: Dynamic CPU/memory allocation based on workload
- **Cloud-ready**: AWS S3 genome reference integration with automatic downloading
- **Fault tolerance**: Resume capability and error recovery mechanisms
- **Parallel processing**: Multi-genome, multi-candidate simultaneous analysis

**⚡ Modern Python Architecture**
- **Type safety**: Full mypy compliance with Pydantic models for data validation
- **Async/await**: Non-blocking I/O throughout the pipeline for improved throughput
- **Rich CLI**: Beautiful terminal interface with progress bars, tables, and error formatting
- **Comprehensive testing**: Unit, integration, and pipeline tests with pytest
- **Developer experience**: Pre-commit hooks, automated formatting (black), linting (ruff)

### Repository Structure

```
sirnaforge/
├── 📦 src/sirnaforge/              # Main package (modern src-layout)
│   ├── 🎯 core/                   # Core algorithms and analysis engines
│   │   ├── design.py              # siRNA design, scoring, and candidate generation
│   │   ├── off_target.py          # BWA-MEM2 off-target analysis (transcriptome + miRNA seed)
│   │   └── thermodynamics.py     # ViennaRNA integration & structure prediction
│   ├── 📊 models/                 # Type-safe Pydantic data models
│   │   ├── sirna.py              # siRNA candidates, parameters, results
│   │   └── transcript.py         # Transcript and gene representations
│   ├── 💾 data/                   # Data access and integration layer
│   │   ├── gene_search.py        # Multi-database API integration
│   │   ├── orf_analysis.py       # Reading frame and coding validation
│   │   └── base.py               # Common utilities (FASTA parsing, etc.)
│   ├── 🔧 pipeline/               # Nextflow workflow integration
│   │   ├── nextflow/             # Nextflow execution and config management
│   │   └── resources.py          # Resource and test data management
│   ├── 🛠️ utils/                  # Cross-cutting utilities
│   │   └── logging_utils.py      # Structured logging configuration
│   ├── 📟 cli.py                  # Rich CLI interface with Typer
│   └── workflow.py               # High-level workflow orchestration
├── 🧪 tests/                      # Comprehensive test suite
│   ├── unit/                     # Component-specific unit tests
│   ├── integration/              # Cross-component integration tests
│   ├── pipeline/                 # Nextflow pipeline validation tests
│   └── docker/                   # Container integration tests
├── 🌊 nextflow_pipeline/          # Nextflow DSL2 workflow
│   ├── main.nf                   # Main workflow orchestration
│   ├── nextflow.config           # Execution and resource configuration
│   ├── modules/local/            # Custom process definitions
│   └── subworkflows/local/       # Reusable workflow components
├── 🐳 docker/                     # Container definitions and environments
│   ├── Dockerfile                # Multi-stage production image
│   └── environment-nextflow.yml  # Conda environment specification
├── 📚 docs/                       # Documentation and examples
│   ├── api_reference.rst         # API documentation
│   ├── tutorials/                # Step-by-step guides
│   └── examples/                 # Working code examples
└── 🔧 Configuration files
    ├── pyproject.toml            # Python packaging and tool configuration
    ├── Makefile                  # Development workflow automation
    └── uv.lock                   # Reproducible dependency resolution
## 📊 Output Formats & Results

siRNAforge generates comprehensive, structured outputs for downstream analysis and experimental validation:

### Workflow Output Structure

```
output_directory/
├── 📁 transcripts/                # Retrieved transcript sequences
│   ├── {gene}_transcripts.fasta   # All retrieved transcript isoforms
│   └── temp_for_design.fasta      # Filtered sequences for design
├── 📁 orf_reports/               # Open reading frame validation
│   └── {gene}_orf_validation.txt  # Coding sequence quality report
├── 📁 sirnaforge/                # Core siRNA design results
│   ├── {gene}_sirna_results.csv   # Complete candidate table
│   ├── {gene}_top_candidates.fasta # Top-ranked sequences for validation
│   └── {gene}_candidate_summary.txt # Human-readable summary
├── 📁 off_target/                # Off-target analysis results
│   ├── basic_analysis.json        # Simplified off-target metrics
│   ├── input_candidates.fasta     # Candidates sent for analysis
│   └── results/                   # Detailed Nextflow pipeline outputs
│       ├── aggregated/            # Combined multi-species results
│       └── individual_results/    # Per-candidate detailed analysis
├── 📄 workflow_manifest.json      # Complete workflow configuration
└── 📄 workflow_summary.json       # High-level results summary
```

### Key Output Files

**🎯 `{gene}_sirna_results.csv`** - Complete candidate table with all scoring metrics:
```csv
id,guide_sequence,antisense_sequence,transcript_ids,position,gc_content,melting_temp,thermodynamic_score,secondary_structure_score,off_target_score,composite_score
TP53_001,GUAACAUUUGAGCCUUCUGA,UCAGAAGGCUCAAAUGUUAC,"ENST00000269305;ENST00000455263",245,47.6,52.3,0.85,0.92,0.78,4.22
TP53_002,CAUCAACUGAUUGUGCUGC,GCAGCACAAUCAGUUGAUG,"ENST00000269305",512,52.6,54.1,0.91,0.88,0.82,4.45
...
```

**🧬 `{gene}_top_candidates.fasta`** - Ready-to-order sequences for experimental validation:
```fasta
>TP53_001 score=4.22 gc=47.6% transcripts=2
GUAACAUUUGAGCCUUCUGA
>TP53_002 score=4.45 gc=52.6% transcripts=1
CAUCAACUGAUUGUGCUGC
```

**📋 `{gene}_candidate_summary.txt`** - Human-readable summary report:
```
siRNAforge Design Summary for TP53
Generated: 2025-09-08 14:30:22
=================================

Input Statistics:
- Transcripts processed: 3
- Total sequence length: 2,847 bp
- Coding sequences: 1,182 bp

Design Results:
- Candidates generated: 1,156
- Passed filters: 234
- Top candidates selected: 50

Top 5 Candidates:
1. TP53_001: GUAACAUUUGAGCCUUCUGA (Score: 4.22, GC: 47.6%)
2. TP53_002: CAUCAACUGAUUGUGCUGC (Score: 4.45, GC: 52.6%)
...
```

**🔍 Off-target Analysis Outputs:**
```json
{
  "analysis_summary": {
    "candidates_analyzed": 10,
    "total_off_targets": 15,
    "high_confidence_hits": 3
  },
  "by_species": {
    "human": {"transcriptome_hits": 8, "mirna_hits": 2},
    "rat": {"transcriptome_hits": 3, "mirna_hits": 1},
    "rhesus": {"transcriptome_hits": 1, "mirna_hits": 0}
  },
  "candidates": [
    {
      "candidate_id": "TP53_001",
      "guide_sequence": "GUAACAUUUGAGCCUUCUGA",
      "off_target_score": 0.78,
      "species_analysis": {
        "human": {"hits": 5, "seed_matches": 2},
        "rat": {"hits": 2, "seed_matches": 0}
      }
    }
  ]
}
```

### Integration with Analysis Tools

**🔬 For Laboratory Validation:**
- FASTA files can be directly submitted to oligonucleotide synthesis providers
- CSV files import into Excel/R/Python for further analysis
- Candidate rankings support experimental prioritization

**🖥️ For Computational Analysis:**
- JSON outputs enable programmatic result processing
- Structured CSV format supports statistical analysis and machine learning
- Off-target data facilitates safety assessment and regulatory compliance

**📊 For Visualization and Reporting:**
- Summary reports provide publication-ready candidate lists
- Score distributions support quality control assessment
- Multi-species comparisons enable cross-species research applications

## 🔬 Nextflow Pipeline Integration

The integrated Nextflow pipeline provides scalable, containerized off-target analysis:

### Pipeline Features

- **Multi-Species Analysis** - Human, rat, rhesus macaque genomes
- **Parallel Processing** - Each siRNA candidate processed independently
- **Auto Index Management** - Downloads and builds BWA indices on demand
- **Cloud Ready** - AWS Batch, Kubernetes, SLURM support
- **Comprehensive Results** - TSV, JSON, and HTML outputs

### Usage Examples

```bash
# Standalone pipeline execution
nextflow run nextflow_pipeline/main.nf \
  --input candidates.fasta \
  --genome_species "human,rat,rhesus" \
  --outdir results

# With custom genome indices
nextflow run nextflow_pipeline/main.nf \
  --input candidates.fasta \
  --genome_indices "human:/path/to/human/index" \
  --profile docker

# Using S3-hosted indices
nextflow run nextflow_pipeline/main.nf \
  --input candidates.fasta \
  --download_indexes true \
  --profile aws
```

### Pipeline Output Structure

```
results/
├── aggregated/                    # Final combined results
│   ├── combined_mirna_analysis.tsv
│   ├── combined_transcriptome_analysis.tsv
│   ├── combined_summary.json
│   └── analysis_report.html
└── individual_results/            # Per-candidate results
    ├── candidate_0001/
    ├── candidate_0002/
    └── ...
```

## 🛠️ Development & Quality Assurance

### Modern Development Environment with uv

siRNAforge leverages `uv` for lightning-fast dependency management and development workflows:

```bash
# Complete development setup (recommended)
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make install-dev  # Installs all dev dependencies

# Core development commands
make test-local-python  # Fastest Python-only tests (markers=local_python)
make test-fast          # Quick pytest suite excluding slow markers
make lint               # Ruff (lint + format --check) and mypy
make check              # lint-fix + test-fast for pre-commit parity
make docs               # Build Sphinx documentation
make docker             # Build the production Docker image

# Selective dependency installation
uv sync --group analysis    # Jupyter, plotting, pandas extras
uv sync --group pipeline    # Nextflow, Docker integration
uv sync --group docs        # Sphinx documentation tools
uv sync --group lint        # Pre-commit, mypy, ruff, black

# Production deployment (minimal dependencies)
uv sync --no-dev
```

### Conda Environment Management

For local development with bioinformatics tools, siRNAforge provides conda environment management:

```bash
# Create complete development environment
make conda-env

# Update existing environment with new dependencies
make conda-env-update

# Remove environment (cleanup)
make conda-env-clean

# Activate environment for development
conda activate sirnaforge-dev

# Deactivate when done
conda deactivate
```

The conda environment includes all bioinformatics tools (BWA-MEM2, SAMtools, ViennaRNA, etc.) plus Python development dependencies, providing a complete local development setup without Docker.

### Quality Assurance & Testing

**🧪 Comprehensive Test Suite:**
```bash
# Run all tests with coverage reporting
make test
# Output: >95% code coverage across all modules

# Fast development testing (unit tests only)
make test-fast

# Integration tests (includes external APIs)
uv run pytest tests/integration/ -v

# Pipeline tests (requires Docker/Nextflow)
uv run pytest tests/pipeline/ -v

# Specific test categories
uv run pytest tests/unit/test_design.py::test_scoring_algorithm -v
```

**🔍 Code Quality Tools:**
```bash
# Type checking with mypy (strict mode)
uv run mypy src/
# Result: Success: no issues found in 20 source files

# Code formatting with black
uv run black src tests
make format

# Linting with ruff (fast Python linter)
uv run ruff check src tests
make lint

# All quality checks together
make lint  # Includes ruff, black, mypy, nextflow lint
```

### Available Dependency Groups

| Group | Purpose | Key Tools |
|-------|---------|-----------|
| `dev` | Core development (auto-installed) | pytest, black, ruff |
| `test` | Testing frameworks | pytest-cov, pytest-xdist |
| `lint` | Code quality | mypy, ruff, black |
| `analysis` | Data science workflows | jupyter, matplotlib, pandas |
| `pipeline` | Nextflow integration | workflow tools, containers |
| `docs` | Documentation generation | sphinx, sphinx-rtd-theme |

### Code Quality Standards

- **Type Safety**: Full mypy coverage with Pydantic models
- **Formatting**: Black + Ruff for consistent style
- **Testing**: Comprehensive pytest suite with >90% coverage
- **CI/CD**: GitHub Actions with multi-Python testing
- **Security**: Bandit + Safety dependency scanning

## ⚡ Performance & System Requirements

### Performance Benchmarks

**🧬 siRNA Design Performance:**
- **Small genes** (1-5 transcripts): ~2-5 seconds
- **Medium genes** (5-20 transcripts): ~10-30 seconds
- **Large genes** (20+ transcripts): ~1-2 minutes
- **Batch processing** (10 genes): ~5-15 minutes

**🔍 Off-target Analysis Performance:**
- **Per candidate** (single species): ~30-60 seconds
- **Multi-species** (3 genomes): ~2-5 minutes per candidate
- **Batch analysis** (50 candidates): ~1-3 hours (parallelized)

### System Requirements

**🔧 Minimum Requirements:**
- **CPU**: 2 cores, 2.0 GHz
- **RAM**: 4 GB (8 GB recommended for off-target analysis)
- **Storage**: 2 GB free space (+ 50 GB for genome indices)
- **Network**: Internet connection for gene searches and genome downloads

**⚡ Recommended Configuration:**
- **CPU**: 8+ cores, 3.0 GHz (for parallel Nextflow execution)
- **RAM**: 16-32 GB (for large-scale off-target analysis)
- **Storage**: SSD with 100+ GB (for genome indices and temporary files)
- **Network**: High-bandwidth connection for S3 genome downloads

**🐳 Docker Resource Allocation:**
```bash
# Recommended Docker settings
docker run --cpus="4" --memory="8g" \
  -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --genome-species human,rat,rhesus
```

## 🐳 Docker Usage

### Pre-built Images

```bash
# Pull latest stable release
docker pull ghcr.io/austin-s-h/sirnaforge:latest

# Run complete workflow
docker run --rm -v $(pwd):/data \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir /data/results

# Interactive development session
docker run -it --rm -v $(pwd):/data \
  ghcr.io/austin-s-h/sirnaforge:latest bash
```

### Building Custom Images

```bash
# Build production image
make docker

# Build with specific Python version
docker build --build-arg PYTHON_VERSION=3.11 \
  -f docker/Dockerfile -t sirnaforge:py311 .
```

The Docker image uses micromamba with `docker/environment-nextflow.yml` for consistent bioinformatics tool installations across all environments.

## 🧪 Testing & Quality Assurance

### Running Tests

| Command | Under the hood | When to use | Notes |
|---------|----------------|-------------|-------|
| `make test-local-python` | `uv run --group dev pytest -v -m "local_python"` | Fastest feedback loop during development | Python-only markers, no Docker/Nextflow required |
| `make test-unit` | `uv run --group dev pytest -v -m "unit"` | Validate core algorithms | Includes ~30 tests (~30s) |
| `make test-fast` | `uv run --group dev pytest -v -m "not slow"` | Pre-commit or PR checks | Skips slow/integration markers |
| `make test` | `uv run --group dev pytest -v` | Full Python suite | May include slow and docker-marked tests; expect >60s |
| `make test-ci` | `uv run --group dev pytest -m "ci" --junitxml=pytest-report.xml --cov=sirnaforge --cov-report=term-missing --cov-report=xml:coverage.xml -v` | CI pipelines needing artifacts | Produces coverage + JUnit reports |
| `make test-cov` | `uv run --group dev pytest --cov=sirnaforge --cov-report=html --cov-report=term-missing` | Local coverage runs | Outputs HTML coverage in `htmlcov/` |
| `make lint` | Ruff lint + Ruff format check + MyPy | Quick code-quality gate | No automatic fixes |
| `make check` | `make lint-fix` + `make test-fast` | Pre-commit parity | Applies Ruff fixes before running fast pytest subset |

Docker-powered tiers share the same pytest markers but execute inside the published image:

| Command | Container invocation | Resource profile | Purpose |
|---------|----------------------|------------------|---------|
| `make docker-test-smoke` | `docker run … python -m pytest -q -n 1 -m 'docker and smoke'` | 0.5 CPU / 256 MB | Minimal CI smoke (MUST PASS) |
| `make docker-test-fast` | `docker run … python -m pytest -q -n 1 -m 'docker and not slow'` | 1 CPU / 2 GB | Dev-friendly docker coverage |
| `make docker-test` | `docker run … python -m pytest -v -n 1 -m 'docker and (docker_integration or (not smoke))'` | 2 CPUs / 4 GB | Standard docker regression |
| `make docker-test-full` | `docker run … uv run --group dev pytest -v -n 2` | 4 CPUs / 8 GB | Release-grade validation |

> ℹ️ Run `make install-dev` once to install development dependencies and pre-commit hooks before using these targets. The full matrix of commands, filters, and expected runtimes lives in [`docs/testing_guide.md`](docs/testing_guide.md).

#### Docker smoke snapshot

For a quick environment sanity check, `make docker-test-smoke` exercises the published container image with toy data in ~40 seconds (0.5 CPU, 256 MB). A passing run prints **9 passed** with no failures; any remaining pytest collection warnings are tracked in the test suite and should disappear once the dataclass fix in this branch lands.

### Fast CI/CD with Toy Data ⚡

siRNAforge now includes an improved CI/CD workflow designed for quick feedback with minimal resources:

- **⚡ Ultra-fast execution**: < 15 minutes total
- **🪶 Minimal resources**: 256MB memory, 0.5 CPU cores
- **🧸 Toy data**: < 500 bytes of test sequences
- **🔥 Smoke tests**: Essential functionality validation

```bash
# Trigger fast CI/CD workflow locally
pytest -m "smoke" --tb=short

# Use toy data for quick validation
ls tests/unit/data/toy_*.fasta

# Fast workflow vs comprehensive workflow
# Fast:    15 min,  256MB RAM, toy data
# Full:    60 min,    8GB RAM, real datasets
```

See [`docs/ci-cd-fast.md`](docs/ci-cd-fast.md) for detailed documentation.

### Test Categories

- **Unit Tests** - Core algorithm validation
- **Integration Tests** - Component interaction testing
- **Pipeline Tests** - Nextflow workflow validation
- **Docker Tests** - Container functionality testing

## 📚 Documentation

### Local Documentation Building

```bash
# Install documentation dependencies
uv sync --group docs

# Build HTML documentation
make docs

# Generate CLI reference
make docs-cli

# Live-reload docs during editing
make docs-dev
```

### Generated Documentation

- `docs/_build/html/` - Complete Sphinx HTML documentation (via `make docs`)
- `docs/CLI_REFERENCE.md` - Auto-generated CLI help (via `make docs-cli`)
- `docs/api_reference.rst` - Python API reference source
- `docs/modification_annotation_spec.md` - Chemical modifications metadata specification

> 📖 See [docs/getting_started.md](docs/getting_started.md) for detailed tutorials and [docs/deployment.md](docs/deployment.md) for deployment guides.

### Chemical Modifications Metadata

siRNAforge supports structured annotation of chemical modifications, overhangs, and provenance information for siRNA sequences. This enables systematic tracking of modifications like 2'-O-methyl, 2'-fluoro, and phosphorothioate linkages.

**Quick Example:**
```bash
# Create metadata JSON file
cat > metadata.json << 'EOF'
{
  "patisiran_ttr_guide": {
    "id": "patisiran_ttr_guide",
    "sequence": "AUGGAAUACUCUUGGUUAC",
    "target_gene": "TTR",
    "strand_role": "guide",
    "overhang": "dTdT",
    "chem_mods": [
      {
        "type": "2OMe",
        "positions": [1, 4, 6, 11, 13, 16, 19]
      }
    ],
    "provenance": {
      "source_type": "patent",
      "identifier": "US10060921B2",
      "url": "https://patents.google.com/patent/US10060921B2"
    },
    "confirmation_status": "confirmed"
  }
}
EOF

# Annotate FASTA with metadata
sirnaforge sequences annotate sequences.fasta metadata.json -o annotated.fasta

# View sequences with metadata
sirnaforge sequences show annotated.fasta
sirnaforge sequences show annotated.fasta --format json
```

**Features:**
- 🧪 **Chemical Modifications** - Annotate 2'-O-methyl, 2'-fluoro, PS linkages, LNA, etc.
- 📍 **Position Tracking** - 1-based position numbering for each modification
- 🔗 **Overhang Support** - DNA (dTdT) or RNA (UU) overhangs
- 📚 **Provenance** - Track sources (patents, publications, clinical trials)
- ✅ **Confirmation Status** - Mark validated vs. predicted sequences
- 🗂️ **FASTA Headers** - Standardized key-value encoding in headers
- 📄 **JSON Sidecars** - Separate metadata files for easy curation

**Common Modification Types:**
- `2OMe` - 2'-O-methyl (nuclease resistance)
- `2F` - 2'-fluoro (enhanced stability)
- `PS` - Phosphorothioate (nuclease resistance)
- `LNA` - Locked Nucleic Acid (enhanced binding)
- `MOE` - 2'-O-methoxyethyl (improved pharmacokinetics)

**Python API:**
```python
from sirnaforge.models.modifications import (
    StrandMetadata,
    ChemicalModification,
    Provenance,
    SourceType
)

# Create metadata
metadata = StrandMetadata(
    id="my_sirna_guide",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 4, 6, 11])
    ],
    provenance=Provenance(
        source_type=SourceType.PUBLICATION,
        identifier="PMID12345678"
    )
)

# Generate FASTA with metadata
from sirnaforge.models.modifications import SequenceRecord, StrandRole
record = SequenceRecord(
    target_gene="BRCA1",
    strand_role=StrandRole.GUIDE,
    metadata=metadata
)
print(record.to_fasta())
```

📖 See [docs/modification_annotation_spec.md](docs/modification_annotation_spec.md) for complete specification, API reference, and examples.

## 🤝 Contributing

We welcome contributions to siRNAforge! Here's how to get started:

### Development Setup

1. **Fork** the repository on GitHub
2. **Clone** your fork: `git clone https://github.com/yourusername/sirnaforge`
3. **Setup** development environment: `make install-dev`
4. **Create** a feature branch: `git checkout -b feature/amazing-feature`

### Development Workflow

```bash
# Make your changes
# ...

# Ensure code quality
make lint           # Check code style and types
make format         # Auto-format code
make test-local-python  # Fast sanity suite
make check              # Auto-fix lint + fast pytest

# Commit and push
git add .
git commit -m 'Add amazing feature'
git push origin feature/amazing-feature
```

### Contribution Guidelines

- **Code Style**: Follow Black formatting and Ruff linting rules
- **Type Hints**: All new code must include type annotations
- **Tests**: Add tests for new functionality
- **Documentation**: Update docstrings and documentation
- **Commit Messages**: Use conventional commit format

### Pull Request Process

1. Ensure all tests pass and code is properly formatted
2. Update documentation for any API changes
3. Add entries to `CHANGELOG.md` for user-facing changes
4. Create a pull request with a clear description

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for detailed guidelines.

## 📄 License

This project is licensed under the MIT License - see the [`LICENSE`](LICENSE) file for details.

## 🙏 Acknowledgments

siRNAforge builds upon excellent open-source tools and libraries:

- **ViennaRNA Package** - RNA secondary structure prediction
- **BWA-MEM2** - Fast and accurate sequence alignment
- **Nextflow** - Workflow management and containerization
- **BioPython** - Python bioinformatics toolkit
- **Pydantic** - Data validation and type safety
- **Modern Python Stack** - uv, Typer, Rich for developer experience

> **Note**: Much of the code in this repository was developed with assistance from AI agents, but all code has been reviewed, tested, and validated by human developers.
