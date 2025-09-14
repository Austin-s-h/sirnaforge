# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2025-09-13
- Improved thermodynamic calculation reporting and modified composite score to heavily favor (90%) duplex binding energy
- Improved output data formatting and introduced pandera DataFrameSchemas to describe key reports
- Cleaning up docs and depreciated/placeholder behavior

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

[Unreleased]: https://github.com/Austin-s-h/sirnaforge/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Austin-s-h/sirnaforge/releases/tag/v0.1.0
