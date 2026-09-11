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

**Requirements:** Python 3.10-3.12

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
SIRNAFORGE_CACHE_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/sirnaforge"
OUTPUT_DIR="$(pwd)/results"
mkdir -p "$SIRNAFORGE_CACHE_DIR" "$OUTPUT_DIR"
docker run --rm --userns=host --user "$(id -u):$(id -g)" \
  -v "$(pwd)":/data:ro -w /tmp \
  -v "$OUTPUT_DIR":/output \
  -v "$SIRNAFORGE_CACHE_DIR":/home/sirnauser/.cache/sirnaforge \
  -e SIRNAFORGE_CACHE_DIR=/home/sirnauser/.cache/sirnaforge \
  -e NXF_HOME=/home/sirnauser/.cache/sirnaforge/nextflow/home \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir /output
```

This reuses the global siRNAforge cache across container runs while keeping bind-mounted outputs owned by the invoking host user.

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
