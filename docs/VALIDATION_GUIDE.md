# Enhanced Data Validation System

- What it does: Keeps your inputs, intermediate data, and outputs sane using Pydantic + Pandera.
- How to use: Pick a preset (development, production, performance) and pass `ValidationConfig` into your `WorkflowConfig`.
- When it runs: At INPUT, TRANSCRIPT_RETRIEVAL, DESIGN, and OUTPUT stages with configurable strictness.
- Why it helps: Fewer silent errors, clearer failures, and actionable reports.

siRNAforge includes a comprehensive validation system that integrates Pydantic models and Pandera schemas to ensure data integrity throughout the siRNA design pipeline.

## Overview

The validation system provides:

- **Multi-level validation** with configurable strictness (STRICT, WARNING, DISABLED)
- **Stage-specific controls** for different parts of the pipeline
- **Bioinformatics-specific validations** for genomic data
- **Cross-system consistency** between Pydantic and Pandera
- **Comprehensive reporting** with detailed error tracking

## Quick Start

### Basic Usage

```python
from sirnaforge.validation import ValidationConfig, ValidationMiddleware
from sirnaforge.workflow import WorkflowConfig, SiRNAWorkflow

# Create validation config (defaults to STRICT level)
validation_config = ValidationConfig()

# Create workflow with validation
workflow_config = WorkflowConfig(
    output_dir="results",
    gene_query="TP53",
    validation_config=validation_config
)

workflow = SiRNAWorkflow(workflow_config)
results = await workflow.run_complete_workflow()
```

### Configuration Presets

```python
from sirnaforge.validation import ValidationPresets

# Development: More lenient, collect all errors
dev_config = ValidationPresets.development()

# Production: Strict validation, fail fast
prod_config = ValidationPresets.production()

# Performance: Minimal validation for speed
perf_config = ValidationPresets.performance()
```

## Validation Levels

### STRICT (Default)
- Validation errors cause workflow failure
- Best for production and critical analyses
- Ensures highest data quality

### WARNING
- Validation errors logged as warnings
- Workflow continues execution
- Good for development and exploratory analysis

### DISABLED
- Skip validation entirely
- Maximum performance
- Use only when data quality is guaranteed

## Stage-Specific Validation

The system validates data at key pipeline stages:

### INPUT
- Design parameter consistency
- Filter criteria validation
- Scoring weight normalization
- Cross-system constraint checking

### TRANSCRIPT_RETRIEVAL
- Sequence composition validation
- Required field presence
- Sequence length checks
- Nucleotide character validation

### DESIGN
- siRNA candidate consistency
- Sequence length constraints (19-23 nt)
- GC content validation
- Biological constraint checking

### OUTPUT
- DataFrame schema validation
- Data type consistency
- Required column presence
- Value range validation

## Bioinformatics Validations

### Sequence Validation
```python
from sirnaforge.validation import ValidationUtils

# Validate nucleotide sequence
result = ValidationUtils.validate_nucleotide_sequence("ATCGATCGATCG")
print(f"Valid: {result.is_valid}")
print(f"GC Content: {result.metadata['gc_content']:.1f}%")

# Validate siRNA length
result = ValidationUtils.validate_sirna_length("ATCGATCGATCGATCGATCGA")
print(f"Length valid: {result.is_valid}")
```

### Parameter Validation
```python
from sirnaforge.models.sirna import DesignParameters

params = DesignParameters()
result = ValidationUtils.validate_parameter_consistency(params)
print(f"Parameters valid: {result.is_valid}")
print(f"Total weight: {result.metadata['total_weight']}")
```

### Candidate Validation
```python
result = ValidationUtils.validate_candidate_consistency(candidate)
print(f"Candidate valid: {result.is_valid}")
print(f"Calculated GC: {result.metadata['calculated_gc']:.1f}%")
```

## Advanced Configuration

### Custom Stage Levels
```python
from sirnaforge.validation import ValidationConfig, ValidationLevel, ValidationStage

config = ValidationConfig(
    default_level=ValidationLevel.WARNING,
    stage_levels={
        ValidationStage.INPUT: ValidationLevel.STRICT,    # Critical inputs
        ValidationStage.OUTPUT: ValidationLevel.STRICT,   # Critical outputs
        ValidationStage.DESIGN: ValidationLevel.WARNING,  # Allow warnings
    }
)
```

### Selective Validation
```python
config = ValidationConfig(
    validate_sequences=True,      # Check nucleotide sequences
    validate_biology=True,        # Apply bio-specific rules
    validate_consistency=False,   # Skip cross-checks (faster)
    validate_ranges=True          # Check numerical ranges
)
```

### Performance Tuning
```python
config = ValidationConfig(
    batch_size=5000,              # Larger batches for big datasets
    enable_caching=True,          # Cache validation results
    collect_all_errors=False,     # Fail fast for performance
    max_validation_errors=10      # Stop after 10 errors
)
```

## Validation Reports

The system generates comprehensive validation reports:

```json
{
  "validation_config": {
    "default_level": "strict",
    "validate_sequences": true,
    "validate_biology": true
  },
  "stage_reports": [
    {
      "stage": "input",
      "overall_valid": true,
      "summary": {
        "total_items": 1,
        "valid_items": 1,
        "total_errors": 0,
        "total_warnings": 0,
        "success_rate": 100.0
      }
    }
  ],
  "summary": {
    "total_stages": 4,
    "stages_passed": 4,
    "total_errors": 0,
    "total_warnings": 2,
    "overall_success": true
  }
}
```

## Integration with Existing Code

### Workflow Integration
The validation system integrates seamlessly with existing workflows:

```python
# Validation is automatically applied at checkpoints
workflow = SiRNAWorkflow(config)
results = await workflow.run_complete_workflow()

# Access validation results
validation_summary = workflow.validation._generate_summary()
print(f"Validation success: {validation_summary['overall_success']}")
```

### Manual Validation
You can also use validation utilities independently:

```python
from sirnaforge.validation import ValidationUtils
import pandas as pd

# Validate a DataFrame against schema
df = pd.read_csv("sirna_results.csv")
result = ValidationUtils.validate_dataframe_schema(df, "sirna_candidates")
print(f"Schema valid: {result.is_valid}")
```

## Error Handling

### Graceful Degradation
```python
try:
    workflow = SiRNAWorkflow(config)
    results = await workflow.run_complete_workflow()
except ValueError as e:
    if "Validation failed" in str(e):
        print("Validation errors detected - check validation report")
        # Handle validation failure gracefully
```

### Warning Collection
```python
# Set config to collect warnings without failing
config = ValidationConfig(default_level=ValidationLevel.WARNING)
workflow = SiRNAWorkflow(WorkflowConfig(validation_config=config))

results = await workflow.run_complete_workflow()

# Check for warnings in validation report
validation_report = workflow.validation._generate_summary()
if validation_report['total_warnings'] > 0:
    print(f"Found {validation_report['total_warnings']} validation warnings")
```

## Best Practices

### Development Workflow
1. Start with **development** preset for exploratory work
2. Use **warning** level to identify potential issues
3. Switch to **strict** level before production

### Production Deployment
1. Always use **strict** validation in production
2. Enable comprehensive logging
3. Monitor validation reports for data quality trends
4. Set up alerts for validation failures

### Performance Optimization
1. Use **performance** preset for large datasets
2. Disable unnecessary validation types
3. Increase batch sizes for bulk processing
4. Enable caching for repeated validations

### Data Quality Monitoring
1. Save validation reports for trend analysis
2. Track validation success rates over time
3. Monitor warning patterns for early issue detection
4. Use validation metadata for quality metrics

## Migration Guide

### From Previous Versions
The enhanced validation system is backward compatible:

```python
# Old way (still works)
design_result.save_csv("results.csv")  # Uses built-in validation

# New way (enhanced validation)
workflow_config = WorkflowConfig(
    validation_config=ValidationConfig()  # Add validation config
)
```

### Integrating Custom Validations
Add your own validation rules:

```python
from sirnaforge.validation.utils import ValidationResult

def custom_validation(data):
    result = ValidationResult()
    # Add your validation logic
    if some_condition:
        result.add_error("Custom validation failed")
    return result
```

This enhanced validation system ensures robust data quality throughout your siRNA design pipeline while providing the flexibility to adapt to different use cases and environments.
