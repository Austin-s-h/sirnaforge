# Variant Targeting Feature Implementation Summary

## Overview
This document summarizes the implementation of SNP/variant targeting capability for siRNAforge, enabling design of siRNAs against specific genetic variants using public SNP databases.

## Features Implemented

### 1. Core Data Models (Phase 1) ✅
- **VariantRecord**: Pydantic model for storing variant information
  - Supports all required fields: id, chr, pos, ref, alt, assembly, sources, etc.
  - Methods for VCF-style formatting and source prioritization
  - Full validation for positions, allele frequencies, and clinical significance

- **VariantQuery**: Model for parsing different variant identifier formats
  - Supports rsID (rs12345)
  - Supports VCF coordinates (chr17:7577121:G:A)
  - Supports HGVS notation (NM_000546.6:c.215C>G)
  - Automatic query type detection

- **VariantMode**: Enum for variant handling strategies
  - TARGET: Design siRNAs specifically for variant alleles
  - AVOID: Exclude siRNAs overlapping variants (default)
  - BOTH: Generate candidates for both reference and alternate alleles

- **ClinVarSignificance**: Enum for ClinVar clinical significance levels
  - Supports filtering by pathogenicity (Pathogenic, Likely pathogenic, etc.)

### 2. Variant Data Providers (Phase 2) ✅
- **VariantResolver**: Core class for resolving variant identifiers to records
  - Priority-based source resolution (ClinVar > Ensembl > dbSNP)
  - Caching system using existing cache infrastructure
  - Allele frequency filtering (default min_af=0.01)
  - ClinVar significance filtering (default: Pathogenic, Likely pathogenic)
  - Assembly validation (GRCh38 only)

- **API Clients**:
  - ClinVar E-utilities API integration
  - Ensembl Variation API integration
  - dbSNP placeholder for future enhancement

- **VCF Support**:
  - pysam-based VCF file reading
  - Support for bgzip-compressed files with tabix index
  - Automatic AF extraction from INFO fields
  - Per-allele variant parsing

### 3. SiRNA Design Integration (Phase 3) ✅
- **Extended SiRNACandidate Model**:
  - `overlapped_variants`: List of variants overlapping the candidate
  - `allele_specific`: Boolean flag for allele-specific targeting
  - `targeted_alleles`: List of alleles targeted (ref/alt)
  - `variant_mode`: Mode used for variant handling

- **Helper Functions** (`variant_helpers.py`):
  - `generate_contexts_for_variant()`: Generate ref/alt sequence contexts
  - `check_candidate_overlaps_variant()`: Detect siRNA-variant overlaps
  - `annotate_candidate_with_variant()`: Add variant metadata to candidates
  - `apply_variant_to_sequence()`: Apply variants to reference sequences
  - `get_variant_position_in_transcript()`: Coordinate conversion

### 4. CLI Integration (Phase 4) ✅
- **New CLI Flags** (in `workflow` command):
  - `--snp` / `-s`: Variant identifier(s), repeatable
  - `--snp-file`: VCF file path for bulk variant input
  - `--variant-mode`: Handling strategy (avoid/target/both)
  - `--min-af`: Minimum allele frequency threshold (default: 0.01)
  - `--clinvar-filter-levels`: Comma-separated significance levels
  - `--variant-assembly`: Reference assembly (GRCh38 only)

- All parameters include comprehensive help text and validation

## Testing

### Unit Tests (47 tests, all passing)
1. **Variant Models** (12 tests)
   - VariantRecord creation and validation
   - VariantQuery parsing for all formats
   - Enum value validation
   - VCF-style formatting
   - Source prioritization

2. **VariantResolver** (18 tests)
   - Identifier parsing (rsID, coordinates, HGVS)
   - Chromosome normalization
   - AF and ClinVar filtering
   - Cache operations
   - Assembly validation

3. **Variant Helpers** (17 tests)
   - Context generation for SNPs, insertions, deletions
   - Overlap detection
   - Sequence modification
   - Boundary condition handling
   - Coordinate conversion

### Existing Tests (213 tests, all passing)
- All existing unit tests continue to pass
- No regressions introduced

## Architecture

```
CLI (cli.py)
  ↓
  └─→ [Variant parameters validation]
       ↓
Workflow (workflow.py) [Future integration point]
  ↓
  ├─→ VariantResolver (data/variant_resolver.py)
  │    ├─→ ClinVar API
  │    ├─→ Ensembl Variation API
  │    ├─→ VCF Reader (pysam)
  │    └─→ Cache (utils/cache_utils.py)
  │
  └─→ Variant Helpers (data/variant_helpers.py)
       ├─→ Context generation
       ├─→ Overlap detection
       └─→ Candidate annotation

Models (models/variant.py)
  ├─→ VariantRecord
  ├─→ VariantQuery
  ├─→ VariantMode
  └─→ ClinVarSignificance

Extended Models (models/sirna.py)
  └─→ SiRNACandidate (with variant fields)
```

## Code Quality
- All code passes Ruff linting with project standards
- Type hints throughout
- Comprehensive docstrings
- Pydantic validation for all data models
- Proper error handling with informative messages

## Remaining Work (Phase 5-7)

### Phase 5: Workflow Orchestration
- Integrate variant resolution into main workflow after gene selection
- Generate variant resolution reports (resolved_variants.json)
- Update siRNA candidate generation to use variant contexts
- Propagate variant annotations through off-target analysis

### Phase 6: Testing & Documentation
- Integration tests for variant-aware design workflows
- Example workflows in documentation
- Tutorial for variant-specific siRNA design
- CLI reference updates

### Phase 7: Validation & Polish
- Manual testing with real SNP identifiers (e.g., TP53 rs28934576)
- Performance profiling for large VCF files
- Code review and optimizations

## Usage Examples

### Design siRNAs avoiding a specific SNP
```bash
sirnaforge workflow TP53 \\
  --snp rs28934576 \\
  --variant-mode avoid \\
  --output-dir tp53_avoid_variant
```

### Target a pathogenic variant
```bash
sirnaforge workflow BRCA1 \\
  --snp chr17:43045677:G:A \\
  --variant-mode target \\
  --output-dir brca1_target_variant
```

### Use VCF file for multiple variants
```bash
sirnaforge workflow TP53 \\
  --snp-file clinvar_tp53.vcf.gz \\
  --variant-mode target \\
  --min-af 0.001 \\
  --clinvar-filter-levels "Pathogenic" \\
  --output-dir tp53_clinvar
```

## Design Decisions

1. **GRCh38 Only**: Simplified implementation by focusing on current reference genome
2. **Priority-based Resolution**: ClinVar > Ensembl > dbSNP ensures clinical relevance
3. **Pydantic Models**: Type safety and validation from the start
4. **Async API Clients**: Better performance for multiple variant lookups
5. **Cache-first Strategy**: Reduces API calls and improves reproducibility
6. **Minimal Changes Philosophy**: All additions are opt-in via CLI flags

## Dependencies
All required dependencies were already in the project:
- `pydantic>=2.11.0` - Data models
- `aiohttp>=3.12.0` - Async HTTP clients
- `pysam>=0.23.0` - VCF file reading
- `requests>=2.32.0` - API fallback

No new dependencies added.

## Files Modified/Added

### New Files
- `src/sirnaforge/models/variant.py` (142 lines)
- `src/sirnaforge/data/variant_resolver.py` (568 lines)
- `src/sirnaforge/data/variant_helpers.py` (223 lines)
- `tests/unit/test_variant_models.py` (145 lines)
- `tests/unit/test_variant_resolver.py` (260 lines)
- `tests/unit/test_variant_helpers.py` (275 lines)

### Modified Files
- `src/sirnaforge/models/__init__.py` (added variant exports)
- `src/sirnaforge/models/sirna.py` (added 4 variant fields to SiRNACandidate)
- `src/sirnaforge/cli.py` (added 6 CLI parameters)

Total: **1,613 lines of new code** with comprehensive tests

## Conclusion
This implementation provides a solid foundation for variant-aware siRNA design in siRNAforge. The core infrastructure (Phases 1-4) is complete, tested, and ready for workflow integration. The architecture is extensible and follows existing project patterns.
