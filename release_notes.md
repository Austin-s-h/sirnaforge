# 🧬 siRNAforge v0.4.0

**Comprehensive siRNA design toolkit with multi-species off-target analysis**

## [0.4.0] - 2025-1-5

### Breaking

- Drop support for Python 3.9; minimum supported Python is now 3.10. This
  enables PEP 604 union syntax (e.g. `list[str] | None`) and other modern
  typing features. Update developer environments, CI, and pre-commit hooks to
  use Python 3.10 or later.


### New Features
- **Variant Targeting Implementation**: Complete Phase 1-5 implementation for targeting specific genetic variants
  - Core variant models and resolver infrastructure with Parquet-based caching
  - Population-specific AF filtering for geographic variant targeting
  - Phase 5 workflow integration with CLI flags for variant_mode parameter
  - Comprehensive variant feature implementation summary and documentation
- **Enhanced CLI and Workflow Integration**:
  - Off-target-only entry point for pre-designed siRNA candidates
  - Support for any sequence length in thermodynamic calculations
  - Improved CLI enum integration for variant_mode parameter
- **Docker and Container Improvements**:
  - Python 3.12 upgrade with optimized Dockerfile targeting
  - Login shell PATH preservation with `/etc/profile.d/conda-path.sh`
  - Enhanced container testing with dedicated test categories
  - Improved Docker entrypoint and health checks

### Improvements
- **Performance Optimizations**:
  - Parquet-based variant cache for improved performance
  - Cache-first index reuse with complete validation
  - Transcriptome filtering and major refactoring for simplicity
  - Unified cache management with Pythonic interface
- **CI/CD Pipeline Enhancements**:
  - Python 3.12 requirement across all workflows
  - Aligned pre-commit mypy with uv package manager
  - Enhanced release workflow with comprehensive testing
  - Improved Docker build and test categorization
- **Code Quality and Maintenance**:
  - Python 3.10+ syntax modernization throughout codebase
  - Comprehensive typing improvements and linting fixes
  - Enhanced error handling and validation middleware
  - Improved documentation with live CLI output examples

### Dependencies
- **Python 3.12 Support**: Full upgrade to Python 3.12 with modern syntax
- **Enhanced Dependencies**: Added `pyarrow>=18.0.0` for Parquet support
- **Updated Packages**: Modernized dependency versions with improved compatibility
- **uv Package Manager**: Full alignment with uv for faster dependency resolution

### Performance
- **Variant Caching**: Parquet-based storage for improved variant data performance
- **Memory Optimization**: Reduced memory requirements for Docker-constrained environments
- **Parallel Processing**: Enhanced concurrent execution for variant analysis
- **Index Reuse**: Cache-first approach for transcriptome indices


## [0.3.4] - 2025-12-31

### Added
- **Transcript Annotation Provider Layer**: New data provider interface for fetching genomic transcript annotations
  - Added `AbstractTranscriptAnnotationClient` interface in `src/sirnaforge/data/base.py`
  - Implemented `EnsemblTranscriptModelClient` using Ensembl REST API (lookup/id and overlap/region endpoints)
  - Added `VepConsequenceClient` stub for optional VEP enrichment (behind config flag, placeholder implementation)
  - New Pydantic models: `Interval`, `TranscriptAnnotation`, and `TranscriptAnnotationBundle` in `src/sirnaforge/models/transcript_annotation.py`
  - In-memory LRU cache with TTL for transcript annotations
  - Support for fetching by stable IDs or genomic regions
  - Comprehensive unit tests with mocked REST responses
  - Integration tests for real Ensembl REST API (gated by `@pytest.mark.requires_network`)

### Improvements
- **Extensible Architecture**: Transcript annotation provider follows the same layered pattern as existing data providers (gene search, ORF analysis, transcriptome management)
- **Reference Tracking**: Annotations include provenance metadata (provider, endpoint, reference choice) for reproducibility
- **Error Handling**: Robust handling of unresolved IDs and network errors with fallback to unresolved list

### Documentation
- Added comprehensive docstrings for all new classes and methods
- Unit and integration tests serve as usage examples


## Installation

### Docker (recommended)
```bash
# Pull the latest release
docker pull ghcr.io/austin-s-h/sirnaforge:0.4.0

# Quick test - should complete in ~2 seconds
docker run --rm ghcr.io/austin-s-h/sirnaforge:0.4.0 sirnaforge --help
```

### Python package
```bash
# Via pip
pip install sirnaforge==0.4.0

# Via uv (recommended for speed)
uv add sirnaforge==0.4.0

# Verify installation
sirnaforge --help
```

### Development
```bash
# Clone repository
git clone https://github.com/austin-s-h/sirnaforge.git
cd sirnaforge

# Setup using the Makefile (uses uv under the hood)
make dev

# Verify with quick tests
make test-dev  # ~15 seconds
```

## Usage

### Basic workflow (Docker)
```bash
# Complete gene-to-siRNA workflow
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:0.4.0 \
  sirnaforge workflow TP53 --output-dir results --genome-species human

# Custom transcript file
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:0.4.0 \
  sirnaforge design transcripts.fasta --output results.csv
```

### Python API
```python
from sirnaforge import SiRNADesigner, GeneSearcher

# Search for gene transcripts
searcher = GeneSearcher(species="human")
transcripts = searcher.search_gene("TP53")

# Design siRNAs
designer = SiRNADesigner()
candidates = designer.design_from_transcripts(transcripts)
```

### Commands
```bash
sirnaforge workflow   # Complete pipeline: gene → siRNAs
sirnaforge search     # Gene/transcript search
sirnaforge design     # siRNA candidate generation
sirnaforge validate   # Input file validation
sirnaforge config     # Show configuration
sirnaforge cache      # Manage miRNA databases
sirnaforge version    # Version information
```

## Key features in v0.4.0

- **Smart siRNA design** - Thermodynamic scoring with 90% duplex binding weight
- **Off-target analysis** - BWA-MEM2 alignment across human/rat/rhesus genomes
- **ViennaRNA integration** - Secondary structure prediction for accuracy
- **Pandera data schemas** - Type-safe output validation and formatting
- **uv package manager** - Fast dependency resolution
- **Production Docker** - Pre-built images with all bioinformatics tools
- **Nextflow pipeline** - Scalable execution with automatic parallelization

## Testing & quality

**This release passed comprehensive validation:**
- **Unit tests** - Core algorithm validation
- **Local Python tests** - Fastest development iteration
- **Docker smoke tests** - Critical functionality verification
- **Integration tests** - End-to-end workflow validation
- **Code quality** - Ruff formatting and MyPy type checking

## Resources

### Documentation
- [**Full documentation**](https://austin-s-h.github.io/sirnaforge) - Complete user guide
- [**Quick start**](https://github.com/austin-s-h/sirnaforge#-quick-start) - Get running in minutes
- [**API reference**](https://austin-s-h.github.io/sirnaforge/api) - Python API documentation
- [**Development guide**](https://github.com/austin-s-h/sirnaforge/blob/main/CONTRIBUTING.md) - Contributing instructions

### Container images
- **Versioned:** `ghcr.io/austin-s-h/sirnaforge:0.4.0`
- **Latest:** `ghcr.io/austin-s-h/sirnaforge:latest`
- **Registry:** [ghcr.io/austin-s-h/sirnaforge](ghcr.io/austin-s-h/sirnaforge)

### Support
- [**Source code**](https://github.com/austin-s-h/sirnaforge)
- [**Issues**](https://github.com/austin-s-h/sirnaforge/issues) - Bug reports & feature requests
- [**Discussions**](https://github.com/austin-s-h/sirnaforge/discussions) - Community support
- [**Changelog**](https://github.com/austin-s-h/sirnaforge/blob/main/CHANGELOG.md) - Version history
- [**All releases**](https://github.com/austin-s-h/sirnaforge/releases) - Previous versions

---

**Quick verification:**
```bash
# Test Docker image (should complete in ~2 seconds)
docker run --rm ghcr.io/austin-s-h/sirnaforge:0.4.0 sirnaforge version

# Expected output: siRNAforge v0.4.0
```


[Unreleased]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/austin-s-h/sirnaforge/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/austin-s-h/sirnaforge/releases/tag/v0.1.0
