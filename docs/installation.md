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

(execution-modes-and-index-preparation)=
## Execution modes and index preparation

Reference preparation runs in the Python orchestrator: it downloads/decompresses cDNA, resolves
species and transcript annotations, and reuses a complete, trusted cached index when available.
A reference without an index is handed to Nextflow as a FASTA. **Index construction and alignment
then use the same selected execution provider and image.** The Python host does not need BWA-MEM2
when these tasks execute in Docker or Singularity.

| Deployment | Nextflow task profile | Where BWA-MEM2 runs | Required setup |
| --- | --- | --- | --- |
| Full siRNAforge Docker image | `local`, selected automatically inside the container | Inside that image for both indexing and alignment | Docker on the host; mount inputs, output and cache. No Docker socket or host BWA installation is needed. |
| Host pip/uv installation with working Docker | `docker`, selected automatically | Task containers using the configured Nextflow image | Python package, Java and Nextflow on the host; Docker daemon access and bind-mountable input/cache/work paths. |
| Native installation | `local` | Host toolchain on `PATH` | Python package, Java, Nextflow, BWA-MEM2, SAMtools and ViennaRNA. Set `SIRNAFORGE_USE_LOCAL_EXECUTION=1` to select it explicitly. |
| Activated Conda environment | `local`, or `conda` with the prepared environment | The activated environment | Install the complete toolchain and siRNAforge first. The embedded `conda`/`test` profiles do not bootstrap an environment or install the application. uv alone is not a Conda runtime. |
| Host with Singularity | `singularity` fallback when Docker is unavailable | The configured task image | Java/Nextflow, Singularity, a compatible image and accessible bind paths. This deployment is not exercised by the Docker CI suite. |

Automatic selection inside the full image stays `local`, including when the host requested another
profile. On a host, working Docker is preferred, followed by Singularity, Conda and local execution.
`SIRNAFORGE_USE_LOCAL_EXECUTION=1` selects local execution. Python callers can supply a profile in
`WorkflowConfig(nextflow_config={"profile": "local"})`; a local profile requires its native tools.

For host Docker execution, select an image containing both siRNAforge and its bioinformatics tools:

```bash
sirnaforge workflow TP53 --species human \
  --transcriptome-fasta human:/data/human_chr21.cdna.fa \
  --nextflow-docker-image ghcr.io/austin-s-h/sirnaforge:dev \
  --output-dir results
```

`SIRNAFORGE_NEXTFLOW_IMAGE` supplies the same image setting through the environment. Pin an image
digest for reproducible deployments. Inside a full container, select the outer container's tag or
digest instead; the tasks execute locally inside it. Native Conda/Singularity deployments require
their own installation validation; the container test suite validates the full Docker image.

The automatic CPU ceiling is the available host/affinity/container-quota count, capped at 16.
Direct Nextflow runs use the JVM-detected processor count with the same cap. Process resource
requests are bounded by `max_cpus`, `max_memory` and `max_time`; retries stay within those ceilings.
Python callers can set explicit ceilings through `WorkflowConfig(nextflow_config={"max_cpus": 2})`.

Pipeline-built indices live in the Nextflow work cache and can be reused with resume. Keep that
cache if you want to avoid rebuilding; an incomplete or failed build is not reused as a valid index.
The low-level Python `build_bwa_index()` helper intentionally remains a native-tool API.

## Prebuilt screening indices

```bash
sirnaforge workflow TP53 \
  --transcriptome-indices human:/data/human_cdna,mouse:/data/mouse_cdna \
  --output-dir results
```

Explicit indices replace automatic transcriptome downloads. The species labels on the entries
define their screening scope; `--species` / `--mirna-species` still configure the miRNA layer.
If `--transcriptome-fasta` is also supplied, that explicit reference is included alongside the
indices. Each prebuilt prefix needs its complete BWA-MEM2 files and a readable neighboring cDNA
FASTA (`<prefix>.fa`, `.fasta`, `.fna`, or the prefix itself) for hit annotation. Index files are
staged as Nextflow inputs so they are visible in task containers.

`sirnaforge design` only designs candidates. Its unused `--genome-index` option has been removed;
use `workflow` or `offtarget` with `--transcriptome-indices` to request screening.

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
