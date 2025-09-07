# Nextflow 25+ Integration Notes

## Key Changes in Nextflow 25+

### Command Line Changes
- ✅ **`-preview`** replaces the old `-dry-run` option
- ✅ **Enhanced reporting** with better HTML output
- ✅ **Improved Docker integration** with better error handling
- ✅ **Better resource management** and process monitoring

### New Features We're Using
1. **Enhanced Preview Mode**: `-preview` shows execution plan without running
2. **Improved Linting**: `nextflow lint` has better error detection
3. **Better Docker Support**: More reliable container execution
4. **Enhanced Reporting**: Richer HTML reports with more details

### Development Workflow
```bash
# Check syntax and preview execution
make nextflow-preview

# Lint pipeline scripts  
make lint-nextflow

# Run with test data
make nextflow-run

# Generate comprehensive reports
make nextflow-report
```

### Integration Test Results
Based on the test output:
- ✅ Nextflow 25.04.6 installed and working
- ✅ Linting passes with no errors
- ✅ Basic functionality works
- ✅ Docker integration functional
- ✅ Preview mode now correctly implemented

### Next Steps for Pipeline Development
1. **Add more comprehensive tests** for different genome species
2. **Implement stub-run mode** for faster testing
3. **Add workflow-level error handling** using Nextflow 25+ features
4. **Enhance reporting** with custom HTML templates
5. **Add parameter validation** using Nextflow schema features
