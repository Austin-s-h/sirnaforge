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
uv run sirnaforge --help
uv run sirnaforge version
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
├── off_target/
│   └── input_candidates.fasta      # FASTA passed to Nextflow (includes dirty controls)
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

```bash
# View top candidates
head -6 my_first_analysis/sirnaforge/TP53_pass.csv

# Check workflow summary
cat my_first_analysis/logs/workflow_summary.json
```

## Understanding Your Results

### Key Output Files

| File | Purpose | Use For |
|------|---------|---------|
| `*_pass.csv` | High-quality candidates that pass all filters | **Laboratory experiments** |
| `*_all.csv` | Complete candidate list with detailed scores | Analysis and custom filtering |
| `workflow_summary.json` | High-level analysis statistics | QC and reporting |

### Quality Metrics

**Key metrics to look for:**
- `asymmetry_score` ≥0.65 (optimal)
- `gc_content` 35-60% (balanced)
- `melting_temp` 55-65°C (effective)

**Score Ranges:**
- **Composite Score (7-10)**: Higher is better; combines all factors
- **Thermodynamic Score (0-1)**: Measures RNA folding favorability and asymmetry
- **GC Content (35-60%)**: Optimal range 40-55% for stability and accessibility
- **Asymmetry Score (0.65-1.0)**: Higher values indicate better guide strand selection

## Essential Commands

| Command | Purpose | Example |
|---------|---------|---------|
| `workflow` | Complete gene → siRNA pipeline | `sirnaforge workflow TP53 --output-dir results` |
| `search` | Find gene transcripts | `sirnaforge search BRCA1 -o transcripts.fasta` |
| `design` | Generate siRNAs from FASTA | `sirnaforge design input.fasta --top-n 20` |
| `validate` | Check FASTA format | `sirnaforge validate input.fasta` |
| `version` | Show version | `sirnaforge version` |

### Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--gc-min` / `--gc-max` | 30% / 60% | GC content window |
| `--length` | 21 | siRNA length in nucleotides |
| `--top-n` | 10 (`design`), 20 (`workflow`) | Number of candidates retained |
| `--output-dir` | `sirna_workflow_output` | Output directory for workflows |
| `--verbose` | `false` | Emit detailed progress messages |

## Next Steps

- **📖 [Usage Examples](usage_examples.md)** - Batch processing, custom parameters, chemical modifications
- **⚙️ [CLI Reference](cli_reference.md)** - Complete parameter documentation
- **🧬 [Thermodynamic Guide](thermodynamic_guide.md)** - Understanding scoring metrics
- **🔧 [Developer Guide](developer/development.md)** - For contributors and developers
