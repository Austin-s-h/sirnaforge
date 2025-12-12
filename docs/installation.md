# Installation

## Recommended: pip/uv

```bash
# Using pip
pip install sirnaforge

# Using uv (faster)
uv pip install sirnaforge

# Verify
sirnaforge version
```

**Requirements:** Python 3.9-3.12

## Development Setup

For contributing or running from source:

```bash
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make dev    # Installs deps + pre-commit hooks
```

(docker-full-bioinformatics-stack)=
## Docker (Full Bioinformatics Stack)

The Docker image includes Nextflow, BWA-MEM2, SAMtools, and ViennaRNA for complete off-target analysis.

```bash
# Pull pre-built image
docker pull ghcr.io/austin-s-h/sirnaforge:latest

# Run a workflow
docker run --rm -v $(pwd):/data -w /data \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir results/
```

:::{tip}
Use Docker when you need:
- Nextflow pipeline execution
- BWA-MEM2 off-target alignment
- Reproducible analysis environment
:::

## Verify Installation

```bash
# Check version
sirnaforge version

# Run help
sirnaforge --help

# Quick test with sample data
sirnaforge design examples/sample_transcripts.fasta -o test.csv
```
