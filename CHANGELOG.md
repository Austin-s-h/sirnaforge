# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.2] - 2025-10-26

### ✨ New Features
- **miRNA Design Mode**: New `--design-mode mirna` option for microRNA-specific siRNA design
  - Specialized `MiRNADesigner` subclass with miRNA-biogenesis-aware scoring
  - Enhanced CSV schema with miRNA-specific columns (strand_role, biogenesis_score)
  - CLI support via `--design-mode` flag with automatic parameter adjustment
- **miRNA Seed Match Analysis**: Integrated miRNA off-target screening in Nextflow pipeline
  - Lightweight seed region matching (positions 2-8) against miRNA databases
  - Automatic miRNA database download and caching from MirGeneDB
  - Per-candidate and aggregated miRNA hit reports in TSV/JSON formats
  - Configurable via `--mirna-db` and `--mirna-species` flags
- **Species Registry System**: Canonical species name mapping and normalization
  - Unified species identifiers across genome and miRNA databases
  - Automatic species alias resolution (e.g., "human" → "Homo sapiens" → mirgenedb slug)
  - Support for multi-species analysis with consistent naming

### 🔧 Improvements
- **Nextflow Pipeline Enhancements**:
  - Reduced memory requirements for Docker-constrained environments (2GB → 1GB for most processes)
  - Added miRNA seed analysis module with BWA-based matching
  - Improved error handling and progress reporting
  - Better resource allocation with memory/CPU buffers
- **Data Validation**: Extended Pandera schemas for miRNA-specific columns
- **CSV Output**: New columns `transcript_hit_count` and `transcript_hit_fraction` track guide specificity
- **miRNA Database Manager**: Enhanced with species normalization and canonical name mapping

### 🐛 Bug Fixes
- Fixed Nextflow Docker configuration for resource-constrained CI environments
- Resolved schema validation errors for miRNA columns in mixed-mode workflows
- Fixed typing issues in pipeline CLI functions

### 📚 Documentation
- **Major Documentation Consolidation**: Reorganized structure for improved user experience
  - Simplified navigation from 4 to 3 main sections (Getting Started, User Guide, Reference, Developer)
  - Consolidated `getting_started.md` and `quick_reference.md` into comprehensive guide
  - Streamlined tutorials to 2 focused guides (pipeline integration, custom scoring)
  - Created dedicated developer section for advanced documentation
- **Complete API Reference**: Added 18 previously missing modules
  - Comprehensive coverage of all 27 sirnaforge modules
  - Auto-generated Sphinx documentation with proper cross-references
- **Quality Improvements**: Configured ruff D rules for docstring validation
  - Fixed 116 docstring formatting issues automatically
  - Clean Sphinx builds with no warnings
- **Usage Examples**: Added miRNA seed analysis workflow documentation

### 🧪 Testing
- **New Test Coverage**: 232 new tests for miRNA design mode
  - Comprehensive unit tests for MiRNADesigner scoring
  - Schema validation tests for miRNA-specific columns
  - Integration tests for miRNA database functionality
- **Test Organization**: Normalized test markers for consistent CI/CD workflows
- **Documentation Tests**: Verified all doc builds and cross-references work correctly

### 📦 Dependencies
- No new runtime dependencies (leverages existing httpx, pydantic, pandera)
- Enhanced development dependencies for documentation generation

---

## [0.2.1] - 2025-10-24

```markdown
## [X.Y.Z] - YYYY-MM-DD

### ✨ New Features
- Brief description of new features

### 🔧 Improvements
- Improvements to existing functionality
- Performance enhancements

### 🐛 Bug Fixes
- Fixed specific issues
- Resolved edge cases

### 📊 Performance
- Performance improvements with metrics if available

### 🧪 Testing
- New tests added
- Test coverage improvements
```

---

## [0.2.1] - 2025-10-24

### ✨ New Features
- **Chemical Modification System**: Comprehensive infrastructure for siRNA chemical modifications
  - Default modification patterns automatically applied to designed siRNAs (standard_2ome, minimal_terminal, maximal_stability)
  - New `--modifications` and `--overhang` CLI flags for workflow and design commands
  - FDA-approved Patisiran (Onpattro) pattern included in example library
- **Modification Metadata Models**: Pydantic models for StrandMetadata, ChemicalModification, Provenance tracking
- **FASTA Annotation System**: Merge modification metadata into FASTA headers with full roundtrip support
- **Remote FASTA Inputs**: Workflow supports `--input-fasta` with automatic HTTP download and caching
- **Enhanced Pandera Schemas**: Runtime DataFrame validation with @pa.check_types decorators, automatic addition of modification columns

### 🔧 Improvements
- Modification columns (guide/passenger overhangs and modifications) now included in CSV outputs
- CLI `sequences show` command with JSON/FASTA/table output formats
- CLI `sequences annotate` command for merging metadata into FASTA files
- Standardized `+` separators in modification headers (backward compatible with `|`)
- Resource resolver for flexible input handling (local files, HTTP URLs)
- Improved type safety with Pandera schema validation on DesignResult.save_csv() and _generate_orf_report()

### 🐛 Bug Fixes
- Fixed JSON metadata loading regression with StrandMetadata subscripting
- Resolved mypy typing issues for optional FASTA descriptions
- Fixed CLI output handling for modification metadata

### 📚 Documentation
- **Chemical Modification Review** (551 lines): Comprehensive analysis and integration guide
- **Modification Integration Guide** (543 lines): Developer documentation with code examples
- **Modification Annotation Spec** (381 lines): Complete FASTA header specification
- **Example Patterns Library**: 4 production-ready modification patterns with usage guide
- Updated README with chemical modifications feature documentation
- Remote FASTA usage documented in CLI and gene search guides

### 🧪 Testing
- **18 new tests** for chemical modifications (100% passing):
  - 11 integration tests for workflow roundtrip validation
  - 7 tests validating example pattern files
- Added resource resolver unit tests (local paths, HTTP downloads, schemes)
- Extended modification metadata tests for delimiter compatibility
- All 164 tests passing with enhanced Pandera validation

### 📦 Dependencies
- No new runtime dependencies added (uses existing Pydantic, Pandera, httpx)

### ⚡ Performance
- Removed Bowtie indexing (standardized on BWA-MEM2)
- Streamlined off-target analysis pipeline configuration

---

## [0.2.0] - 2025-09-27

### ✨ New Features
- **miRNA Database Cache System** (`sirnaforge cache`) - Local caching and management of miRNA databases from multiple sources with automatic updates
- **Comprehensive Data Validation** - Pandera DataFrameSchemas for type-safe output validation ensuring consistent CSV/TSV report formatting
- **Enhanced Thermodynamic Scoring** - Modified composite score to heavily favor (90%) duplex binding energy for improved siRNA selection accuracy
- **Workflow Input Flexibility** - Added FASTA file input support for custom transcript analysis workflows
- **Embedded Nextflow Pipeline** - Integrated Nextflow execution directly within Python API for scalable processing

### 🔧 Improvements
- **Performance Optimization** - Parallelized off-target analysis and improved memory efficiency for large transcript sets
- **CLI Enhancement** - Better Unicode support, cleaner help text, and improved error reporting
- **Data Schema Validation** - Robust output validation with detailed error messages using modern Pandera 0.26.1 patterns
- **Documentation Overhaul** - Comprehensive testing guide, thermodynamic documentation, and improved API references
- **Development Workflow** - Enhanced Makefile with Docker testing categories, release validation, and conda environment support

### � Bug Fixes
- **Security Improvements** - Resolved security linting issues and improved dependency management
- **Off-target Analysis** - Fixed alignment indexing and improved multi-species database handling
- **CI/CD Pipeline** - Resolved build failures, improved test categorization, and enhanced release automation
- **Unicode Handling** - Fixed CLI display issues in various terminal environments

### 📊 Performance
- **10-100x Faster Dependencies** - Full migration to uv package manager for ultra-fast installs and environment management
- **Optimized Algorithms** - Improved thermodynamic calculation efficiency with better filtering strategies
- **Parallel Processing** - Enhanced concurrent execution for off-target analysis across multiple genomes

### 🧪 Testing & Infrastructure
- **Enhanced Test Categories** - Smoke tests (256MB), integration tests (2GB), and full CI validation
- **Docker Improvements** - Multi-stage builds, intelligent entrypoint, and resource-aware testing
- **Release Automation** - Comprehensive GitHub Actions workflow with quality gates and artifact management

### 📚 Documentation
- **Testing Guide** - Comprehensive documentation for all test categories and Docker workflows
- **Thermodynamic Guide** - Detailed explanation of scoring algorithms and parameter optimization
- **CLI Reference** - Auto-generated command documentation with examples
- **Development Setup** - Streamlined onboarding with conda environment and uv integration

### 📦 Dependencies & Architecture
- **Modern Python Support** - Maintained compatibility across Python 3.9-3.12 with improved type safety
- **Pydantic Integration** - Enhanced data models with validation middleware and error handling
- **Containerization** - Production-ready Docker images with conda bioinformatics stack
- **Package Management** - Full uv adoption for dependency resolution and virtual environment management

## [0.1.0] - 2025-09-06

### Added
- Initial release of siRNAforge toolkit
- Core siRNA design algorithms with thermodynamic scoring
- Multi-database gene search (Ensembl, RefSeq, GENCODE)
- Rich command-line interface with Typer and Rich
- Comprehensive siRNA candidate scoring system
- Off-target prediction framework
- Nextflow pipeline integration for scalable analysis
- Docker containerization for reproducible environments
- Python API with Pydantic data models
- Comprehensive test suite with unit and integration tests
- Modern development tooling with uv, black, ruff, mypy

### Core Features
- **Gene Search**: Multi-database transcript retrieval
- **siRNA Design**: Algorithm-driven candidate generation
- **Quality Control**: GC content, structure, and specificity filters
- **Scoring System**: Composite scoring with multiple components
- **Workflow Orchestration**: End-to-end gene-to-siRNA pipeline
- **CLI Interface**: Rich, user-friendly command-line tools
- **Python API**: Programmatic access for automation

### Supported Operations
- `sirnaforge workflow`: Complete gene-to-siRNA analysis
- `sirnaforge search`: Gene and transcript search
- `sirnaforge design`: siRNA candidate generation
- `sirnaforge validate`: Input file validation
- `sirnaforge config`: Configuration display
- `sirnaforge version`: Version information

### Technical Stack
- **Language**: Python 3.9-3.12
- **Package Management**: uv for fast dependency resolution
- **Data Models**: Pydantic for type-safe data handling
- **CLI Framework**: Typer with Rich for beautiful output
- **Testing**: pytest with comprehensive coverage
- **Code Quality**: black, ruff, mypy for consistency
- **Containerization**: Multi-stage Docker builds
- **Pipeline**: Nextflow integration for scalability
- **Documentation**: Sphinx with MyST parser, Read the Docs theme

[Unreleased]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/austin-s-h/sirnaforge/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/austin-s-h/sirnaforge/releases/tag/v0.1.0
