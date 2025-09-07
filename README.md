# 🧬 siRNAforge — Design. Verify. Deliver.

<div align="center">
  <img src="docs/branding/sirnaforge_logo_3.svg" alt="siRNAforge Logo" width="200"/>
  
  **siRNAforge - Multi-species gene to siRNA design, off-target prediction, and ranking. Comprehensive siRNA design toolkit for gene silencing.**
  
  [![Python 3.9–3.12](https://img.shields.io/badge/python-3.9--3.12-blue.svg)](https://www.python.org/downloads/)
  [![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
  [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
</div>

## ✨ Key Features

- 🎯 **Algorithm-driven design** - Comprehensive siRNA design with thermodynamic analysis
- 🔍 **Off-target analysis** - Integrated genome-wide off-target prediction  
- 📊 **Rich scoring system** - Multi-component scoring with customizable weights
- 🧪 **Secondary structure** - RNA folding prediction using ViennaRNA
- 🐍 **Modern Python** - Type-safe code with Pydantic models and rich CLI
- ⚡ **High performance** - Built with `uv` for lightning-fast dependency management
- 🐳 **Containerized** - Docker support for reproducible environments
- 🔬 **Pipeline ready** - Nextflow integration for scalable analysis

> **Note:** Supports Python 3.9-3.12. Python 3.13+ not supported due to ViennaRNA compatibility.

## 🚀 Quick Start

### Installation

**🐳 Docker (Recommended - Everything Included):**
```bash
# Build comprehensive image with all dependencies
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge
docker build -f docker/Dockerfile -t sirnaforge:latest .

# Quick start - complete workflow
docker run -v $(pwd):/workspace -w /workspace sirnaforge:latest \
    sirnaforge workflow TP53 --output-dir results
```

The Docker build uses the `docker/environment.yml` conda environment specification to install all bioinformatics tools and dependencies via micromamba, ensuring consistent versions across all environments.

**📦 Local Development:**
```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup Python environment
git clone https://github.com/Austin-s-h/sirnaforge
cd sirnaforge
make install-dev
```

### Additional Dependencies (Off-target Analysis)

**🐳 Docker (Recommended):**
The comprehensive Docker image includes everything via the `docker/environment.yml` conda environment:
- ✅ Nextflow (≥25.04.0) - Workflow management
- ✅ BWA-MEM2 (≥2.2.1) - Sequence alignment  
- ✅ SAMtools (≥1.19.2) - BAM/SAM processing
- ✅ ViennaRNA (≥2.7.0) - RNA folding prediction
- ✅ AWS CLI (≥2.0) - S3 genome downloads
- ✅ Java 17 - Nextflow runtime

**📦 Alternative Local Setup:**
For local development without Docker:
```bash
# Install bioinformatics tools via micromamba
micromamba env create -f environment.yml
micromamba activate nextflow
```

**Deployment Options:**
- 🐳 **Docker** - Complete environment in one image
- ☁️ **AWS Batch** - Scalable cloud computing  
- 🖥️ **HPC/SLURM** - High-performance computing
- ⚙️ **Kubernetes** - Container orchestration

> 📖 **See [docs/deployment.md](docs/deployment.md) for comprehensive deployment guides.**

### Basic Usage

```bash
# Complete workflow: gene → transcripts → siRNAs → off-target analysis
uv run sirnaforge workflow TP53 --output-dir results

# Design siRNAs from transcripts
uv run sirnaforge design transcripts.fasta -o results.tsv -n 20

# Search for gene transcripts
uv run sirnaforge search TP53 -o transcripts.fasta

# Validate input files
uv run sirnaforge validate input.fasta

# Show help
uv run sirnaforge --help
```

### Python API

```python
from sirnaforge.workflow import run_sirna_workflow
from sirnaforge import SiRNADesigner, DesignParameters
from sirnaforge.models import FilterCriteria
import asyncio

# Complete workflow
async def design_sirnas():
    results = await run_sirna_workflow(
        gene="TP53",
        output_dir="results",
        genome_species=["human", "rat"]
    )
    return results

# Individual components
params = DesignParameters(
    sirna_length=21,
    filters=FilterCriteria(gc_min=30, gc_max=60)
)

designer = SiRNADesigner(params)
results = designer.design_from_file("transcripts.fasta")

for candidate in results.top_candidates:
    print(f"{candidate.id}: {candidate.guide_sequence} (score: {candidate.composite_score:.1f})")
```

## 🏗️ Architecture

### Repository Structure

```
sirnaforge/
├── 📦 src/sirnaforge/           # Main package (modern src-layout)
│   ├── 🎯 core/                # Core algorithms (design, scoring, filtering)
│   ├── 📊 models/              # Pydantic data models  
│   ├── 💾 data/                # Data utilities (Ensembl, file I/O)
│   ├── 🔧 pipeline/            # Nextflow integration
│   ├── 🛠️ utils/               # Helper utilities
│   └── 📟 cli.py               # Rich CLI interface
├── 🧪 tests/                   # Comprehensive test suite
│   ├── unit/                   # Unit tests
│   ├── integration/            # Integration tests
│   └── pipeline/               # Pipeline tests
├── 🐳 docker/                  # Multi-stage Dockerfiles
├── 📚 docs/                    # Documentation
├── 📋 examples/                # Usage examples
├── ⚙️ pyproject.toml           # Modern Python packaging with uv
└── 🔧 Makefile                 # Development commands
```

### Algorithm Pipeline

1. **Gene Search** → Retrieve transcripts from Ensembl/RefSeq
2. **ORF Analysis** → Validate open reading frames
3. **siRNA Design** → Generate candidates with thermodynamic scoring
4. **Quality Control** → Filter based on GC content, structure, etc.
5. **Off-target Analysis** → Multi-species genome scanning (Nextflow)
6. **Ranking & Selection** → Composite scoring and top candidate selection

## 🛠️ Development

### Modern Dependency Management with uv

siRNAforge uses `uv`'s dependency groups for flexible, fast development:

```bash
# Standard development (installs dev group automatically)
uv sync

# Specific workflows
uv sync --group analysis    # Jupyter, plotting, data science
uv sync --group pipeline    # Nextflow, Docker integration
uv sync --group docs        # Sphinx documentation building

# Install everything
uv sync --all-groups

# Production only (no dev dependencies)
uv sync --no-dev
```

### Development Commands

```bash
# Setup
make install-dev          # Full development environment
make install-analysis     # Add Jupyter/plotting tools
make install-pipeline     # Add Nextflow/Docker tools

# Quality & Testing
make test                 # Run test suite
make lint                 # Code quality checks
make format               # Auto-format code
make type-check           # MyPy type checking

# Building
make build                # Build package
make docker-build         # Build Docker image
```

### Available Dependency Groups

| Group | Purpose | Key Tools |
|-------|---------|-----------|
| `dev` | Core development (default) | pytest, black, ruff |
| `test` | Testing frameworks | pytest-cov, pytest-xdist |
| `lint` | Code quality | mypy, ruff, black |
| `analysis` | Data science | jupyter, matplotlib, pandas |
| `pipeline` | Workflow integration | nextflow, docker |
| `docs` | Documentation | sphinx, sphinx-rtd-theme |

### Code Quality

- **Type Safety**: Full mypy coverage with Pydantic models
- **Formatting**: Black + Ruff for consistent style
- **Testing**: Comprehensive pytest suite with coverage
- **CI/CD**: GitHub Actions with multi-Python testing
- **Security**: Bandit + Safety scanning

## 🐳 Docker Usage

The Docker image uses micromamba with the `docker/environment.yml` environment specification for consistent, reproducible bioinformatics tool installations.

```bash
# Build optimized image
make docker-build

# Run design workflow
docker run --rm -v $(pwd):/data sirnaforge \
  sirnaforge design /data/transcripts.fasta -o /data/results.tsv

# Interactive development
docker run -it --rm -v $(pwd):/data sirnaforge bash
```

## 🔬 Nextflow Pipeline

The integrated Nextflow pipeline handles computationally intensive off-target analysis:

```bash
# Standalone pipeline
nextflow run nextflow_pipeline/main.nf \
  --input candidates.fasta \
  --download_indexes true

# Multi-species analysis
nextflow run nextflow_pipeline/main.nf \
  --input candidates.fasta \
  --genomes "human,rat,rhesus" \
  --profile docker
```

### Pipeline Features

- **Multi-Species Support** - Human, rat, rhesus macaque genomes
- **Dual Alignment** - Bowtie (fast) + BWA-MEM2 (sensitive)
- **Seed-Aware Scoring** - Prioritizes mismatches in seed regions
- **Auto Index Building** - Downloads and builds genome indexes on demand
- **Cloud Ready** - AWS Batch, Kubernetes support

## 📊 Output Formats

- **TSV/CSV** - Machine-readable results tables
- **FASTA** - Sequence outputs for downstream analysis
- **JSON** - Structured metadata and configuration
- **HTML** - Interactive reports with visualizations

## 🧪 Testing

```bash
# Run all tests
make test

# Specific test suites
uv run pytest tests/unit/           # Unit tests
uv run pytest tests/integration/    # Integration tests
uv run pytest tests/pipeline/       # Pipeline tests

# With coverage
uv run pytest --cov=sirnaforge --cov-report=html
```

## 🤝 Contributing

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature/amazing-feature`
3. **Setup** development environment: `make install-dev`
4. **Make** your changes with tests
5. **Verify** quality: `make lint && make test`
6. **Commit** changes: `git commit -m 'Add amazing feature'`
7. **Push** and **create** a Pull Request

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for detailed guidelines.

## 📄 License

MIT License - see [`LICENSE`](LICENSE) file for details.

## 🙏 Acknowledgments

- **ViennaRNA** - RNA secondary structure prediction
- **Bowtie/BWA** - Sequence alignment capabilities  
- **Python Scientific Stack** - NumPy, Pandas, BioPython
- **Modern Python Tooling** - uv, Pydantic, Typer, Rich
- **AI Assistance** - Much of the code in this repository was generated with the assistance of AI agents. However, all of it was cobbled together and reviewed by a human.

---

<div align="center">
  <strong>siRNAforge — Design. Verify. Deliver.</strong><br>
  Professional siRNA design for the modern researcher.
</div>
