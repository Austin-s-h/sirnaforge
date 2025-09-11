# Pandera Schema Modernization

This document describes the improvements made to siRNAforge's pandera schemas to align with current best practices for Pandera 0.26.1.

## Summary of Changes

### 1. Modern Import Patterns
- **Before**: Used deprecated `import pandera as pa` with warnings
- **After**: Using recommended `import pandera.pandas as pa` to eliminate deprecation warnings

### 2. Class-Based Schema Models
- **Before**: Used legacy `DataFrameSchema` with dictionary definitions
- **After**: Modern `DataFrameModel` classes with type annotations for better IDE support and maintainability

### 3. Improved Type Safety
- **Before**: Basic column type definitions
- **After**: Typed `Series[dtype]` annotations with proper nullable field handling

### 4. Enhanced Validation
- **Before**: Basic regex and range checks
- **After**: Custom `@pa.dataframe_check` validators for bioinformatics-specific logic:
  - Nucleotide sequence validation (supports both DNA and RNA)
  - siRNA length range validation (19-23 nucleotides)
  - Comprehensive documentation with field descriptions

### 5. Better Error Reporting
- **Before**: Generic validation errors
- **After**: Detailed schema metadata with descriptions and titles for clearer debugging

### 6. Flexible Nullable Field Handling
- **Before**: Inconsistent nullable type handling
- **After**: Proper `object` types for nullable fields with configuration options

## Schema Files Modified

### `src/sirnaforge/models/schemas.py`
- Complete rewrite using modern pandera patterns
- Three main schemas: `SiRNACandidateSchema`, `ORFValidationSchema`, `OffTargetHitsSchema`
- Custom validation functions for biological data integrity

### `src/sirnaforge/workflow.py`
- Updated empty DataFrame type handling to match new schema expectations
- Minimal changes to maintain backward compatibility

## Benefits

1. **Eliminated deprecation warnings** - Code is future-proof for Pandera updates
2. **Improved type safety** - Better IDE support and early error detection
3. **Enhanced documentation** - Clear field descriptions and validation rules
4. **Better error messages** - More informative validation failures
5. **Biological validation** - Domain-specific checks for siRNA and genomic data
6. **Maintainability** - Cleaner, more readable schema definitions

## Usage Examples

```python
import pandas as pd
from sirnaforge.models.schemas import SiRNACandidateSchema

# Validate siRNA candidate data
df = pd.DataFrame({...})  # Your data
validated_df = SiRNACandidateSchema.validate(df)
```

## Testing

New comprehensive test suite in `tests/unit/test_schemas.py` validates:
- Valid data acceptance
- Invalid data rejection 
- Nullable field handling
- Bioinformatics-specific constraints

All existing tests continue to pass, ensuring backward compatibility.