# siRNAforge

<div align="center">
  <img src="docs/branding/sirnaforge_logo_3.svg" alt="siRNAforge Logo" width="200"/>

  **Computational platform for siRNA design and off-target analysis**

  [![Release](https://github.com/austin-s-h/sirnaforge/actions/workflows/release.yml/badge.svg?branch=master)](https://github.com/austin-s-h/sirnaforge/actions/workflows/release.yml)
  [![Python 3.9–3.12](https://img.shields.io/badge/python-3.9--3.12-blue.svg)](https://www.python.org/downloads/)
  [![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
  [![Docker](https://img.shields.io/badge/docker-available-blue?logo=docker)](https://github.com/users/austin-s-h/packages/container/package/sirnaforge)
  [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/austin-s-h/sirnaforge/blob/main/LICENSE)
</div>

---

siRNAforge is a computational toolkit for designing small interfering RNAs with integrated multi-species off-target analysis. The platform combines thermodynamic modeling, secondary structure prediction, and genome-wide alignment to identify high-specificity siRNA candidates for gene silencing applications.

## Overview

siRNAforge provides an end-to-end workflow from gene query to ranked siRNA candidates:

- **Multi-database transcript retrieval** - Ensembl, RefSeq, GENCODE integration
- **Thermodynamic and structural scoring** - ViennaRNA-based secondary structure prediction
- **Multi-species off-target analysis** - BWA-MEM2 alignment across human, rat, and rhesus macaque genomes
- **Nextflow pipeline integration** - Containerized, scalable workflow execution
- **Chemical modification tracking** - Metadata support for 2'-O-methyl, 2'-fluoro, phosphorothioate linkages

**Supported Python versions:** 3.9–3.12 (Python 3.13+ pending ViennaRNA compatibility)

## Documentation

**📚 [View Full Documentation](docs/)**

### Getting Started
- [Installation Guide](docs/getting_started.md#installation) - Docker, conda, and local setup
- [Quick Start Tutorial](docs/getting_started.md#quick-start) - Your first analysis in minutes
- [Usage Examples](docs/usage_examples.md) - Common workflows and patterns

### Technical Reference
- [API Reference](docs/api_reference.rst) - Complete Python API documentation
- [CLI Reference](docs/cli_reference.md) - Command-line interface
- [Thermodynamic Scoring](docs/thermodynamic_guide.md) - Algorithm details
- [Chemical Modifications](docs/modification_annotation_spec.md) - Metadata specification

### Development
- [Testing Guide](docs/developer/testing_guide.md) - Test suite and execution
- [Contributing Guide](CONTRIBUTING.md) - Development workflow
- [Developer Documentation](docs/developer/) - Architecture and design

## Installation

### Docker (Recommended)

Pre-built images include all bioinformatics dependencies:

```bash
docker pull ghcr.io/austin-s-h/sirnaforge:latest

docker run -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir results
```

### Local Installation

Using the `uv` package manager:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make install-dev
```

**→ [Complete Installation Guide](docs/getting_started.md#installation)**

## Quick Start

### Basic Workflow

Design siRNAs for a target gene:

```bash
sirnaforge workflow TP53 --output-dir results
```

### Advanced Configuration

Multi-species analysis with custom parameters:

```bash
sirnaforge workflow BRCA1 \
  --genome-species "human,rat,rhesus" \
  --gc-min 40 --gc-max 60 \
  --sirna-length 21 \
  --top-n 50 \
  --output-dir brca1_analysis
```

### Component Usage

Individual pipeline components:

```bash
# Retrieve transcript sequences
sirnaforge search TP53 --output transcripts.fasta

# Design siRNAs from sequences
sirnaforge design transcripts.fasta --output results.csv

# Validate input files
sirnaforge validate candidates.fasta
```

**→ [Usage Examples](docs/usage_examples.md) • [CLI Reference](docs/cli_reference.md)**

## Output Structure

siRNAforge generates structured output directories with comprehensive analysis results:

```
output_directory/
├── transcripts/              # Retrieved sequences
├── orf_reports/             # Coding sequence validation
├── sirnaforge/              # siRNA design results
│   ├── {gene}_sirna_results.csv
│   ├── {gene}_top_candidates.fasta
│   └── {gene}_candidate_summary.txt
├── off_target/              # Off-target analysis
│   ├── basic_analysis.json
│   └── results/
└── workflow_summary.json
```

**→ [Output Format Guide](docs/getting_started.md#output-formats)**

## Development

### Setup

```bash
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make install-dev
```

### Development Commands

```bash
make test-local-python  # Fast unit tests
make lint               # Code quality checks
make check              # Lint + tests
make docs               # Build documentation
```

**→ [Testing Guide](docs/developer/testing_guide.md) • [Contributing Guide](CONTRIBUTING.md)**

### Code Quality

- **Type Safety** - Full mypy coverage with Pydantic models
- **Formatting** - Black and ruff automated code formatting
- **Testing** - >90% test coverage with pytest
- **CI/CD** - Automated quality checks and multi-version testing



## System Requirements

### Docker Environment (Recommended)

The Docker image includes all required bioinformatics tools:

- Nextflow (≥25.04.0)
- BWA-MEM2 (≥2.2.1)
- SAMtools (≥1.19.2)
- ViennaRNA (≥2.7.0)
- AWS CLI (≥2.0)
- Java 17

### Local Development

Bioinformatics dependencies can be installed via conda for local development.

**→ [Dependency Guide](docs/getting_started.md#dependencies)**

## Architecture

siRNAforge implements a modular pipeline architecture:

```
Gene Query → Transcript Retrieval → siRNA Design → Off-target Analysis → Ranked Results
```

**Core components:**

- **Gene Search** (`sirnaforge.data.gene_search`) - Multi-database transcript retrieval
- **siRNA Design** (`sirnaforge.core.design`) - Thermodynamic and structural scoring
- **Off-target Analysis** (`sirnaforge.core.off_target`) - BWA-MEM2-based genome alignment
- **Nextflow Pipeline** (`nextflow_pipeline/`) - Scalable workflow execution

**→ [Developer Documentation](docs/developer/)**

## Citation

If you use siRNAforge in your research, please cite:

```
[Citation information to be added upon publication]
```

## License

This project is licensed under the MIT License. See **[LICENSE](LICENSE)** for details.

## Support

- **Issues & Bug Reports** - [GitHub Issues](https://github.com/austin-s-h/sirnaforge/issues)
- **Documentation** - [docs/](docs/)
- **Contributing** - [CONTRIBUTING.md](CONTRIBUTING.md)

## Acknowledgments

siRNAforge integrates several open-source bioinformatics tools:

- **ViennaRNA Package** - RNA secondary structure prediction
- **BWA-MEM2** - High-performance sequence alignment
- **Nextflow** - Scalable workflow orchestration
- **BioPython** - Computational biology utilities

---

**Note:** This software is provided for research use. Portions of the codebase were developed with AI assistance and have been reviewed and validated by human developers.
