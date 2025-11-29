# Getting Started & Quick Reference

Get from installation to your first siRNA analysis in minutes, with essential commands and parameters at your fingertips.

## Installation

### Development (Makefile + uv - Recommended)

`````{tab-set}

````{tab-item} Linux / macOS
```bash
# Install uv package manager
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make dev              # One-command setup: installs dependencies + pre-commit hooks
# OR install development dependencies only:
uv sync --dev
```
````

````{tab-item} Windows
```powershell
# Install uv via PowerShell
irm https://astral.sh/uv/install.ps1 | iex

# Clone and setup
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
uv sync --dev        # Windows: use uv directly (make requires WSL/MinGW)
```
````

`````

**Prerequisites:** Python 3.9-3.12, Git, Make (Linux/macOS)

### Production (Docker - Full Stack)

`````{tab-set}

````{tab-item} Pre-built Image
```bash
# Pull from GitHub Container Registry
docker pull ghcr.io/austin-s-h/sirnaforge:latest

# Verify installation
docker run --rm ghcr.io/austin-s-h/sirnaforge:latest sirnaforge version
```
````

````{tab-item} Build Locally
```bash
# Clone and build
git clone https://github.com/austin-s-h/sirnaforge
cd sirnaforge
make docker-build     # Builds complete image with Nextflow, BWA-MEM2, SAMtools, ViennaRNA

# Test the build
docker run --rm sirnaforge:latest sirnaforge version
```
````

`````

**Includes**: Python packages, Nextflow, BWA-MEM2, SAMtools, ViennaRNA, AWS CLI


## Your First Analysis

### Verify Installation

`````{tab-set}

````{tab-item} Development (uv/make)
```bash
# Direct uv commands
uv run sirnaforge --help
uv run sirnaforge version

# Or using make for quick checks
make version          # Show version via make
make check            # Run linting + fast tests to verify setup
```
````

````{tab-item} Docker
```bash
docker run --rm ghcr.io/austin-s-h/sirnaforge:latest sirnaforge --help
docker run --rm ghcr.io/austin-s-h/sirnaforge:latest sirnaforge version
```
````

`````

### Complete Workflow

`````{tab-set}

````{tab-item} uv
```bash
# End-to-end analysis for TP53
uv run sirnaforge workflow TP53 --output-dir my_first_analysis
```
````

````{tab-item} Docker
```bash
# End-to-end analysis for TP53
docker run --rm \
  -v $(pwd):/workspace \
  -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --output-dir my_first_analysis
```

:::{note}
**Container Execution Note**: When running inside the container, siRNAforge automatically detects the container environment and uses the `local` profile for Nextflow workflows (avoiding Docker-in-Docker). All bioinformatics tools (BWA-MEM2, SAMtools, ViennaRNA, Nextflow) are pre-installed and ready to use.
:::
````

`````

**What this does:**
1. Search TP53 transcripts from Ensembl
2. Design siRNA candidates
3. Score using thermodynamic asymmetry
4. Rank and filter best candidates
5. Generate results

### Results Structure
```
my_first_analysis/
├── sirnaforge/
│   ├── TP53_pass.csv        # Best candidates for lab use
│   └── TP53_all.csv         # All candidates with scores
├── transcripts/
│   ├── TP53_transcripts.fasta      # All retrieved transcripts
│   ├── TP53_canonical.fasta        # Canonical isoform only
│   └── temp_for_design.fasta       # Transcripts used for design
├── logs/
│   ├── workflow_stream.log         # Console output
│   └── workflow_summary.json       # High-level summary
└── orf_reports/
    └── TP53_orf_validation.txt     # ORF analysis results
```

### Examine Your Results

`````{tab-set}

````{tab-item} uv
```bash
# View top candidates
head -6 my_first_analysis/sirnaforge/TP53_pass.csv

# Check workflow summary
cat my_first_analysis/logs/workflow_summary.json
```
````

````{tab-item} Docker
```bash
# View top candidates
docker run --rm -v $(pwd):/workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  head -6 /workspace/my_first_analysis/sirnaforge/TP53_pass.csv

# Check workflow summary (use host tools)
cat my_first_analysis/logs/workflow_summary.json
```
````

`````

## Understanding Your Results

### Key Output Files

| File | Purpose | Use For |
|------|---------|---------|
| `*_pass.csv` | High-quality candidates that pass all filters | **Laboratory experiments** |
| `*_all.csv` | Complete candidate list with detailed scores | Analysis and custom filtering |
| `workflow_summary.json` | High-level analysis statistics | QC and reporting |


## Quality Indicators

**Key metrics:**
- `asymmetry_score` ≥0.65 (optimal)
- `gc_content` 35-60% (balanced)
- `melting_temp` 55-65°C (effective)

### Interpreting Scores

**Score Ranges:**
- **Composite Score (7-10)**: Higher is better; combines all factors
- **Thermodynamic Score (0-1)**: Measures RNA folding favorability and asymmetry
- **GC Content (35-60%)**: Optimal range 40-55% for stability and accessibility
- **Asymmetry Score (0.65-1.0)**: Higher values indicate better guide strand selection
- **MFE (-2 to -8 kcal/mol)**: Moderate stability preferred for effective processing

## Customizing Your Analysis

### High-Quality Design

`````{tab-set}

````{tab-item} uv
```bash
# Stricter quality parameters for research publications
uv run sirnaforge workflow TP53 \
  --output-dir tp53_high_quality \
  --gc-min 35 --gc-max 50 \
  --top-n 30 \
  --verbose
```
````

````{tab-item} Docker
```bash
# Stricter quality parameters for research publications
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 \
    --output-dir tp53_high_quality \
    --gc-min 35 --gc-max 50 \
    --top-n 30 \
    --verbose
```
````

`````

### Step-by-Step Workflow

`````{tab-set}

````{tab-item} uv
```bash
# 1. Search transcripts
uv run sirnaforge search TP53 -o transcripts.fasta --verbose

# 2. Validate sequences
uv run sirnaforge validate transcripts.fasta

# 3. Design siRNAs
uv run sirnaforge design transcripts.fasta \
  -o results.csv \
  --top-n 25 \
  --verbose
```
````

````{tab-item} Docker
```bash
# 1. Search transcripts
docker run --rm -v $(pwd):/workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search TP53 -o transcripts.fasta --verbose

# 2. Validate sequences
docker run --rm -v $(pwd):/workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge validate transcripts.fasta

# 3. Design siRNAs
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge design transcripts.fasta \
    -o results.csv \
    --top-n 25 \
    --verbose
```
````

`````

### Batch Processing Multiple Genes

`````{tab-set}

````{tab-item} uv
```bash
# Process multiple genes
genes=("TP53" "BRCA1" "EGFR")
for gene in "${genes[@]}"; do
    uv run sirnaforge workflow "$gene" \
        --output-dir "analysis_$gene" \
        --gc-min 35 --gc-max 55 \
        --top-n 20
done
```
````

````{tab-item} Docker
```bash
# Process multiple genes
genes=("TP53" "BRCA1" "EGFR")
for gene in "${genes[@]}"; do
    docker run --rm -v $(pwd):/workspace -w /workspace \
        ghcr.io/austin-s-h/sirnaforge:latest \
        sirnaforge workflow "$gene" \
          --output-dir "analysis_$gene" \
          --gc-min 35 --gc-max 55 \
          --top-n 20
done
```
````

`````

## Essential Commands & Parameters

### Core Commands

| Command | Purpose | Key Options |
|---------|---------|------------|
| `workflow` | Gene → siRNA pipeline | `--output-dir`, `--top-n` |
| `search` | Find transcripts | `--database`, `--output` |
| `design` | Generate siRNAs | `--top-n`, `--gc-min/max` |
| `validate` | Check FASTA | `--verbose` |
| `version` | Show version | - |

### Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--gc-min` / `--gc-max` | 30 % / 60 % | GC content window |
| `--length` | 21 | Candidate length in nt |
| `--top-n` | 10 (`design`), 20 (`workflow`) | Number of candidates retained |
| `--output-dir` | `sirna_workflow_output` | Output directory for workflows |
| `--output` | `sirna_results.tsv` | Output path for single-step commands |
| `--verbose` | `false` | Emit detailed progress messages |

### Common Usage Patterns

`````{tab-set}

````{tab-item} Development (uv)
```bash
# High-quality candidates
uv run sirnaforge workflow GENE --gc-min 35 --gc-max 50 --top-n 30

# Fast screening
uv run sirnaforge design input.fasta --skip-structure --top-n 10

# Multi-species analysis
uv run sirnaforge workflow GENE --genome-species "human,mouse,rat"

# Custom parameters
uv run sirnaforge workflow GENE --length 19 --gc-min 40 --max-poly-runs 2

# Quick example using sample data
make example          # Runs design on examples/sample_transcripts.fasta
```
````

````{tab-item} Docker
```bash
# High-quality candidates
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow GENE --gc-min 35 --gc-max 50 --top-n 30

# Fast screening
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge design input.fasta --skip-structure --top-n 10

# Multi-species analysis
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow GENE --genome-species "human,mouse,rat"

# Custom parameters
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow GENE --length 19 --gc-min 40 --max-poly-runs 2

# Using make with Docker
make docker-run GENE=TP53  # Simplified Docker workflow execution
```
````

`````

### Quality Thresholds

**Recommended for lab experiments:**
- GC Content: 35-60% (optimal: 40-55%)
- Asymmetry Score: ≥0.65
- Melting Temperature: 55-65°C
- End Stability (ΔΔG): +2 to +6 kcal/mol


## Development Quick Reference

For contributors and developers, siRNAforge uses **Make commands** for streamlined workflows powered by `uv`:

### Essential Make Commands

**Quick Setup**
```bash
make dev              # One-command setup: install deps + pre-commit hooks
make install          # Install production dependencies only (uv sync --no-dev)
```

**Testing (Organized by Speed)**
```bash
make test-dev         # Fast unit tests (~15s) - for rapid iteration
make test-ci          # Smoke tests with coverage - for CI/CD
make test-release     # Full validation - for releases
make test             # All tests (may have skips/failures)
```

**Testing (By Requirement)**
```bash
make test-requires-docker   # Tests requiring Docker daemon
make test-requires-network  # Tests requiring network access
make test-requires-nextflow # Tests requiring Nextflow
```

**Code Quality**
```bash
make lint             # Check code quality (ruff check + mypy)
make format           # Auto-format and autofix style issues (ruff)
make check            # format + test-dev (mutating quick validation)
make pre-commit       # Run all pre-commit hooks locally
make security         # Run bandit + safety scans
```

**Docker**
```bash
make docker-build     # Build Docker image
make docker-test      # Run tests inside container
make docker-build-test # Clean, rebuild, and validate Docker image
make docker-run       # Run workflow (usage: make docker-run GENE=TP53)
make docker-shell     # Interactive shell
```

**Documentation & Utilities**
```bash
make docs             # Build documentation
make docs-serve       # Serve at localhost:8000
make clean            # Clean build artifacts and caches
make version          # Show version
make example          # Run bundled sample transcripts workflow
make cache-info       # Inspect local transcript/miRNA cache mounts
```

**Why Make?** The Makefile wraps `uv` commands for consistency, provides organized test tiers, and includes convenient shortcuts for Docker workflows. Run `make help` for complete command list.

📖 **[Full development guide →](developer/development.md)**

## Next Steps

- **📖 [Usage Examples](usage_examples.md)** - Comprehensive real-world examples
- **⚙️ [CLI Reference](cli_reference.md)** - Complete parameter documentation
- **🧬 [Custom Scoring Guide](tutorials/custom_scoring.md)** - Advanced thermodynamic principles
