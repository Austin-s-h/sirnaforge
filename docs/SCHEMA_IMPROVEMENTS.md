# Pandera DataFrameModel Implementation Guide

A practical guide for implementing type-safe data validation in bioinformatics workflows using modern Pandera patterns.

## Why Use DataFrameModel?

### Key Benefits
- **Type Safety**: IDE autocomplete and early error detection
- **Self-Documenting**: Schema serves as living documentation
- **Validation**: Automatic data integrity checks with custom validators
- **Maintainability**: Class-based approach is cleaner than dictionary definitions

### Modern vs Legacy Patterns
```python
# ❌ Legacy (deprecated warnings)
import pandera as pa
schema = pa.DataFrameSchema({"col": pa.Column(str)})

# ✅ Modern (recommended)
import pandera.pandas as pa
class MySchema(pa.DataFrameModel):
  col: pa.Series[str]
```

## Implementation Example: Protein Analysis Workflow

### 1. Define Your Schema

```python
import pandera.pandas as pa
from pandera import Field, dataframe_check
import pandas as pd

class ProteinAnalysisSchema(pa.DataFrameModel):
  """Schema for protein sequence analysis results."""

  # Required fields with validation
  protein_id: pa.Series[str] = Field(description="Unique protein identifier")
  sequence: pa.Series[str] = Field(min_length=10, description="Amino acid sequence")
  molecular_weight: pa.Series[float] = Field(gt=0, description="Molecular weight in Da")

  # Optional fields (nullable)
  expression_level: pa.Series[object] = Field(nullable=True, description="Expression level")

  class Config:
    title = "Protein Analysis Results"
    description = "Validated protein sequence data with computed properties"

  @dataframe_check
  def valid_amino_acids(cls, df: pd.DataFrame) -> bool:
    """Ensure sequences contain only valid amino acids."""
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    return df["sequence"].str.upper().apply(
      lambda seq: all(aa in valid_aa for aa in seq)
    ).all()

  @dataframe_check
  def reasonable_mw_range(cls, df: pd.DataFrame) -> bool:
    """Molecular weight should be realistic for proteins."""
    return df["molecular_weight"].between(1000, 500000).all()
```

### 2. Use in Your Workflow

```python
def process_protein_data(input_file: str) -> pd.DataFrame:
  """Process and validate protein analysis data."""

  # Load raw data
  df = pd.read_csv(input_file)

  # Validate with schema
  try:
    validated_df = ProteinAnalysisSchema.validate(df)
    print(f"✅ Validated {len(validated_df)} protein records")
    return validated_df
  except pa.errors.SchemaError as e:
    print(f"❌ Validation failed: {e}")
    raise

# Usage
results = process_protein_data("protein_results.csv")
```

### 3. Handle Empty DataFrames

```python
def create_empty_result() -> pd.DataFrame:
  """Create empty DataFrame matching schema."""
  return pd.DataFrame({
    "protein_id": pd.Series([], dtype="object"),
    "sequence": pd.Series([], dtype="object"),
    "molecular_weight": pd.Series([], dtype="float64"),
    "expression_level": pd.Series([], dtype="object")
  })
```

## Best Practices

1. **Start Simple**: Begin with basic types, add validators incrementally
2. **Document Everything**: Use `Field(description=...)` for all columns
3. **Test Edge Cases**: Include tests for empty data, nulls, and boundary values
4. **Custom Validators**: Use `@dataframe_check` for domain-specific validation
5. **Nullable Handling**: Use `object` dtype for nullable string/mixed columns

## Quick Start Template

```python
import pandera.pandas as pa
from pandera import Field

class YourWorkflowSchema(pa.DataFrameModel):
  # Required fields
  id: pa.Series[str] = Field(description="Unique identifier")
  value: pa.Series[float] = Field(gt=0, description="Measured value")

  # Optional fields
  notes: pa.Series[object] = Field(nullable=True, description="Optional notes")

  class Config:
    title = "Your Workflow Results"
    description = "Description of your data"
```

Replace field names and types with your workflow's specific requirements.
