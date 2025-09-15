# Test Data for siRNAforge

## smoke_test.fasta

Ultra-minimal test dataset designed for Docker smoke tests and CI/CD pipelines.

### Characteristics:
- **Size**: 281 bytes (optimal for fast CI)
- **Sequences**: 3 minimal RNA sequences  
- **Format**: Valid FASTA with RNA sequences (using U instead of T)
- **Purpose**: Quick validation that siRNAforge can process input without requiring extensive computation

### Sequences:
1. `smoke_test_seq_minimal`: 30 nucleotides - basic test sequence
2. `smoke_test_seq_short`: 27 nucleotides - shorter variant for edge case testing  
3. `smoke_test_seq_gc_balanced`: 31 nucleotides - balanced GC content

### Usage:
- Docker smoke tests: `tests/docker/test_smoke_ci.py`
- CLI validation: `sirnaforge design tests/data/smoke_test.fasta`
- Validation: `python tests/validate_smoke_data.py`

### Design Principles:
- **Minimal**: Small enough for fastest possible CI execution
- **Valid**: Real RNA sequences that can be processed by siRNA design algorithms
- **Diverse**: Different lengths and GC content for broader testing
- **Realistic**: Based on actual transcript patterns but simplified

This data is specifically designed for smoke testing - it's not meant to produce meaningful biological results, just to verify that the software stack works correctly.