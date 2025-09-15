# 🧬 siRNAforge — Comprehensive siRNA Design Tool

<div align="center">
  <img src="docs/branding/sirnaforge_logo_3.svg" alt="siRNAforge Logo" width="200"/>

  **Multi-species gene to siRNA design, off-target prediction, and ranking**

  [![Python 3.9–3.12](https://img.shields.io/badge/python-3.9--3.12-blue.svg)](https://www.python.org/downloads/)
  [![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
  [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
  [![Docker](https://img.shields.io/badge/docker-available-blue?logo=docker)](https://github.com/users/Austin-s-h/packages/container/package/sirnaforge)
  [![Nextflow](https://img.shields.io/badge/nextflow-pipeline-brightgreen?logo=nextflow)](https://github.com/Austin-s-h/sirnaforge/tree/main/nextflow_pipeline)
  [![Tests](https://img.shields.io/badge/tests-passing-brightgreen?logo=pytest)](https://github.com/Austin-s-h/sirnaforge/actions)
  [![Coverage](https://img.shields.io/badge/coverage-%3E90%25-brightgreen?logo=codecov)](https://github.com/Austin-s-h/sirnaforge)
  [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Austin-s-h/sirnaforge/blob/main/LICENSE)
</div>


siRNAforge is a modern, comprehensive toolkit for designing high-quality siRNAs with integrated off-target analysis. Built with Python 3.9-3.12, it combines cutting-edge bioinformatics algorithms with robust software engineering practices to provide a complete gene silencing solution for researchers and biotechnology applications.

## ✨ Key Features

- 🎯 **Algorithm-driven design** - Comprehensive siRNA design with multi-component thermodynamic scoring
- 🔍 **Multi-species off-target analysis** - BWA-MEM2 and Bowtie alignment across human, rat, rhesus genomes
- 📊 **Advanced scoring system** - Composite scoring with seed-region specificity and secondary structure prediction
- 🧪 **ViennaRNA integration** - Secondary structure prediction for enhanced design accuracy
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
make test
```

**� Local Development Installation:**
```bash
# Install uv (lightning-fast Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup with development dependencies
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge
make install-dev

# Run tests to verify installation
make test
```

### Essential Dependencies for Off-target Analysis

The Docker image includes all bioinformatics dependencies via conda environment (`docker/environment-nextflow.yml`):

- ✅ **Nextflow** (≥25.04.0) - Workflow orchestration and parallelization
- ✅ **BWA-MEM2** (≥2.2.1) - High-performance genome alignment
- ✅ **Bowtie** (≥1.3.1) - miRNA seed-region analysis
- ✅ **SAMtools** (≥1.19.2) - SAM/BAM file processing and indexing
- ✅ **ViennaRNA** (≥2.7.0) - RNA secondary structure prediction
- ✅ **AWS CLI** (≥2.0) - Automated genome reference downloads
- ✅ **Java 17** - Nextflow runtime environment

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

# Workflow from pre-existing FASTA file
uv run sirnaforge workflow --input-fasta transcripts.fasta \
  --output-dir custom_analysis \
  --offtarget-n 25 \
  custom_gene_name
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

## 🏗️ Architecture & Workflow

### Complete Pipeline Overview

```
Gene Query → Transcript Search → ORF Validation → siRNA Design → Off-target Analysis → Ranked Results
     ↓              ↓                ↓               ↓               ↓                    ↓
Multi-database   Canonical       Coding Frame   Thermodynamic   Multi-species BWA    Scored & Filtered
Gene Search      Isoform         Validation     + Structure     Alignment +          siRNA Candidates
(Ensembl/        Selection                      Scoring         Bowtie miRNA         with Off-target
RefSeq/GENCODE)                                                Analysis             Predictions
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

**🔍 Off-target Analysis** (`sirnaforge.core.off_target`)
- **Dual-engine approach**:
  - **BWA-MEM2**: Sensitive genome-wide alignment for transcriptome off-targets
  - **Bowtie**: High-specificity miRNA seed-region analysis (positions 2-8)
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
│   │   ├── off_target.py          # BWA-MEM2/Bowtie off-target analysis
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
```

### Repository Structure

```
sirnaforge/
├── 📦 src/sirnaforge/              # Main package (modern src-layout)
│   ├── 🎯 core/                   # Core algorithms
│   │   ├── design.py              # siRNA design and scoring
│   │   ├── off_target.py          # Off-target analysis
│   │   └── scoring.py             # Scoring algorithms
│   ├── 📊 models/                 # Pydantic data models
│   │   ├── sirna.py              # siRNA candidate models
│   │   └── transcript.py         # Transcript models
│   ├── 💾 data/                   # Data access layer
│   │   ├── gene_search.py        # Multi-database gene search
│   │   ├── orf_analysis.py       # ORF validation
│   │   └── base.py               # Common data utilities
│   ├── 🔧 pipeline/               # Nextflow integration
│   ├── 🛠️ utils/                  # Helper utilities
│   ├── 📟 cli.py                  # Rich CLI interface
│   └── workflow.py               # Orchestration logic
├── 🧪 tests/                      # Comprehensive test suite
│   ├── unit/                     # Unit tests
│   ├── integration/              # Integration tests
│   └── pipeline/                 # Pipeline tests
├── 🌊 nextflow_pipeline/          # Nextflow workflow
│   ├── main.nf                   # Main workflow
│   ├── nextflow.config           # Configuration
│   ├── genomes.yaml             # Genome specifications
│   └── modules/local/           # Process modules
├── 🐳 docker/                     # Containerization
│   ├── Dockerfile               # Multi-stage build
│   └── environment-nextflow.yml # Conda environment
├── 📚 docs/                       # Documentation
├── 📋 examples/                   # Usage examples
├── ⚙️ pyproject.toml              # Modern Python packaging
├── 🔧 Makefile                    # Development commands
    └── uv.lock                   # Reproducible dependency resolution
```

## 📊 Output Formats & Results
```

## � Output Formats & Results

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

## �🔬 Nextflow Pipeline Integration

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
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge
make install-dev  # Installs all dev dependencies

# Core development commands
make test         # Full test suite with coverage
make lint         # Code quality checks (mypy, ruff, black)
make test-fast    # Quick tests (unit tests only)
make docs-build   # Generate documentation
make docker-build # Build container locally

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
make docker-build

# Build with specific Python version
docker build --build-arg PYTHON_VERSION=3.11 \
  -f docker/Dockerfile -t sirnaforge:py311 .
```

The Docker image uses micromamba with `docker/environment-nextflow.yml` for consistent bioinformatics tool installations across all environments.

## 🧪 Testing & Quality Assurance

### Running Tests

```bash
# Run complete test suite
make test

# Specific test categories
uv run pytest tests/unit/           # Unit tests
uv run pytest tests/integration/    # Integration tests
uv run pytest tests/pipeline/       # Pipeline tests

# With coverage reporting
uv run pytest --cov=sirnaforge --cov-report=html

# Run tests in parallel
uv run pytest -n auto

# Ultra-fast CI/CD smoke tests (NEW!)
make docker-test-smoke

# Validate fast CI/CD setup
python scripts/validate_fast_ci.py
```

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

# Generate usage examples
make docs-examples

# Build complete documentation set
make docs-full
```

### Generated Documentation

- `docs/_build/html/` - Complete Sphinx HTML documentation
- `docs/CLI_REFERENCE.md` - Auto-generated CLI help
- `docs/examples/USAGE_EXAMPLES.md` - Usage examples
- `docs/api_reference.rst` - Python API documentation

> 📖 See [docs/getting_started.md](docs/getting_started.md) for detailed tutorials and [docs/deployment.md](docs/deployment.md) for deployment guides.

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
make test          # Run test suite

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

---

<div align="center">
  <strong>siRNAforge — Comprehensive siRNA design toolkit for gene silencing</strong><br>
  Professional siRNA design for the modern researcher
</div>
