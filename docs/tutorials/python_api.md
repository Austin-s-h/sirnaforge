# Python API Tutorial

This tutorial demonstrates how to use siRNAforge programmatically through its Python API.

## Learning Objectives

- Use siRNAforge classes and functions directly in Python
- Automate siRNA design workflows
- Customize scoring and filtering algorithms
- Integrate with existing Python data analysis pipelines

## Prerequisites

- siRNAforge installed in development mode
- Basic Python programming experience
- Familiarity with basic siRNAforge concepts

For complete class and method documentation, see the {doc}`API Reference <../api_reference>`.


## Getting Started

### Import Core Components

```python
# Core design functionality
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters, FilterCriteria

# Data access
from sirnaforge.data.gene_search import GeneSearcher, DatabaseType
from sirnaforge.workflow import run_sirna_workflow

# Utilities
from pathlib import Path
import asyncio
import pandas as pd
```

```{tip}
Hover over class names in your IDE to see docstrings, or check the {doc}`API Reference <../api_reference>` for full documentation.
```


### Basic API Usage

#### 1. Simple siRNA Design

```python
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters

# Configure design parameters
# See API docs for all DesignParameters options:
# https://sirnaforge.readthedocs.io/en/latest/api_reference.html#sirnaforge.models.sirna.DesignParameters
params = DesignParameters(
    sirna_length=21,
    top_candidates=20,
    gc_content_range=(30.0, 60.0)
)

# Create designer instance
designer = SiRNADesigner(params)

# Design from sequences
sequences = [
    "AUGCGCGAUCUCGAUGCAUGUGCGAUCGAUGCGUAUCGAU",
    "GCCAUGCGAUCGAUGCGUAUCAUCGAUGCGCGAUGCUGUGAU"
]

results = designer.design_from_sequences(sequences)

# Access results (SiRNADesignResult object)
print(f"Generated {len(results.candidates)} siRNA candidates")
for candidate in results.top_candidates:
    print(f"{candidate.sirna_id}: {candidate.guide_sequence} (score: {candidate.composite_score:.1f})")
```

See the {doc}`API Reference <../api_reference>` for complete documentation of `DesignParameters`, `SiRNADesigner`, and `SiRNACandidate`.


#### 2. Design from FASTA File

```python
from pathlib import Path

# Design from file
input_file = Path("transcripts.fasta")
results = designer.design_from_file(input_file)

# Save results
output_file = Path("sirna_results.tsv")
results.save_tsv(output_file)

print(f"Results saved to {output_file}")
```

## Gene Search API

The gene search system retrieves transcripts from genomic databases. See the {doc}`API Reference <../api_reference>` for complete `GeneSearcher` documentation.


### Synchronous Gene Search

```python
from sirnaforge.data.gene_search import search_gene_sync, DatabaseType

# Simple gene search
result = search_gene_sync("TP53", database=DatabaseType.ENSEMBL)

if result.success:
    print(f"Found {len(result.transcripts)} transcripts for {result.gene_info.gene_name}")

    # Access transcript information
    for transcript in result.transcripts:
        print(f"  {transcript.transcript_id}: {len(transcript.sequence)} bp")
        if transcript.is_canonical:
            print("    (Canonical transcript)")
else:
    print(f"Search failed: {result.error}")
```

See the {doc}`API Reference <../api_reference>` for details on `GeneSearchResult`, `TranscriptInfo`, and `GeneInfo`.


### Asynchronous Gene Search

```python
import asyncio
from sirnaforge.data.gene_search import GeneSearcher

async def search_multiple_genes():
    searcher = GeneSearcher()

    genes = ["TP53", "BRCA1", "EGFR", "MYC"]
    results = {}

    for gene in genes:
        print(f"Searching for {gene}...")
        result = await searcher.search_gene(gene, DatabaseType.ENSEMBL)
        results[gene] = result

    return results

# Run the async search
results = asyncio.run(search_multiple_genes())

# Process results
for gene, result in results.items():
    if result.success:
        canonical_count = sum(1 for t in result.transcripts if t.is_canonical)
        print(f"{gene}: {len(result.transcripts)} transcripts ({canonical_count} canonical)")
```

See the {doc}`API Reference <../api_reference>` for details on async methods like `GeneSearcher.search_gene` and `GeneSearcher.search_multiple_databases`.


## Advanced API Usage

### Custom Scoring

siRNAforge allows you to implement custom scoring algorithms by extending the base scorer class. This is useful for research-specific requirements or novel siRNA design principles.

See the {doc}`API Reference <../api_reference>` for complete documentation of the `BaseScorer` interface and built-in scorer implementations.


```python
from sirnaforge.core.design import BaseScorer, SiRNADesigner
from sirnaforge.models.sirna import SiRNACandidate, DesignParameters

class CustomScorer(BaseScorer):
    """Custom scoring algorithm for specific research needs.

    This example prioritizes GC content and 5' position over structure.
    Extend BaseScorer and implement calculate_score() to create your own.
    """

    def __init__(self, weight_gc=0.3, weight_position=0.2, weight_structure=0.5):
        self.weight_gc = weight_gc
        self.weight_position = weight_position
        self.weight_structure = weight_structure

    def calculate_score(self, candidate: SiRNACandidate) -> float:
        """Calculate custom composite score.

        Args:
            candidate: SiRNACandidate with thermodynamic properties

        Returns:
            Composite score between 0-10
        """
        # GC content scoring (optimal around 45%)
        gc_score = 1.0 - abs(candidate.gc_content - 45.0) / 45.0

        # Position scoring (prefer 5' region)
        position_score = max(0, 1.0 - candidate.position / 1000)

        # Structure scoring (from thermodynamics)
        structure_score = candidate.thermodynamic_score

        # Weighted combination
        composite = (
            self.weight_gc * gc_score +
            self.weight_position * position_score +
            self.weight_structure * structure_score
        )

        return composite * 10.0  # Scale to 0-10

# Use custom scorer
custom_scorer = CustomScorer(weight_gc=0.4, weight_structure=0.6)
designer = SiRNADesigner(DesignParameters(), scorer=custom_scorer)
```

For more complex scoring, see the [Custom Scoring Tutorial](custom_scoring.md). The {doc}`API Reference <../api_reference>` provides details on `ThermodynamicAnalyzer` for thermodynamic calculations.

## Configuration and Customization

Configuration in siRNAforge is handled through Pydantic models for type safety and validation. See the {doc}`API Reference <../api_reference>` for complete documentation of `DesignParameters` and `FilterCriteria` with all available options and validation rules.


### Using Configuration Models

```python
from sirnaforge.models.sirna import DesignParameters, FilterCriteria

# Create parameters with validation
params = DesignParameters(
    sirna_length=21,
    top_candidates=20,
    filters=FilterCriteria(
        gc_min=30.0,
        gc_max=60.0,
        min_composite_score=5.0,
        max_poly_runs=4
    )
)

# Parameters are validated automatically
designer = SiRNADesigner(params)
```

### Environment-Based Configuration

```python
import os
from sirnaforge.models.sirna import DesignParameters

# Read from environment variables
gc_min = float(os.getenv("SIRNAFORGE_GC_MIN", "30.0"))
gc_max = float(os.getenv("SIRNAFORGE_GC_MAX", "60.0"))
top_n = int(os.getenv("SIRNAFORGE_TOP_N", "20"))

params = DesignParameters(
    sirna_length=21,
    top_candidates=top_n,
    gc_content_range=(gc_min, gc_max)
)
```

> **Validation**: All parameters are validated by Pydantic. Invalid values raise `ValidationError` with clear messages.

## Testing Your Code

### Unit Testing API Usage

```python
import pytest
from unittest.mock import patch, MagicMock
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters

def test_custom_api_usage():
    """Test custom API usage patterns."""

    # Test parameters
    params = DesignParameters(
        sirna_length=21,
        top_candidates=10,
        gc_content_range=(35.0, 50.0)
    )

    designer = SiRNADesigner(params)

    # Test with known sequences
    test_sequences = [
        "AUGCGCGAUCUCGAUGCAUGUGCGAUCGAUGCGUAUCGAU" * 2  # 80 nt
    ]

    results = designer.design_from_sequences(test_sequences)

    # Assertions
    assert len(results.candidates) > 0
    assert all(35.0 <= c.gc_content <= 50.0 for c in results.candidates)
    assert all(len(c.guide_sequence) == 21 for c in results.candidates)

@patch('sirnaforge.data.gene_search.httpx.AsyncClient')
async def test_mock_api_integration(mock_client):
    """Test API integration with mocked responses."""

    # Mock successful response
    mock_response = MagicMock()
    mock_response.json.return_value = {
        "id": "ENSG00000141510",
        "display_name": "TP53"
    }
    mock_response.raise_for_status.return_value = None
    mock_client.return_value.__aenter__.return_value.get.return_value = mock_response

    # Test search
    from sirnaforge.data.gene_search import GeneSearcher
    searcher = GeneSearcher()
    result = await searcher.search_gene("TP53")

    # Verify mock was called
    assert mock_client.called
```

## Best Practices

### 1. Resource Management

```python
import asyncio
from contextlib import asynccontextmanager

@asynccontextmanager
async def sirna_session():
    """Managed session for siRNAforge operations."""
    searcher = GeneSearcher()
    try:
        yield searcher
    finally:
        await searcher.close()  # Clean up resources

# Usage
async def main():
    async with sirna_session() as searcher:
        result = await searcher.search_gene("TP53")
        # Process result
```

### 2. Progress Monitoring

```python
from tqdm import tqdm
import time

def batch_design_with_progress(genes, output_dir):
    """Batch design with progress monitoring."""

    results = {}

    for gene in tqdm(genes, desc="Processing genes"):
        try:
            # Add small delay to show progress
            time.sleep(0.1)

            result = robust_sirna_design(gene)
            results[gene] = result

            # Update progress description
            tqdm.write(f"✅ Completed {gene}")

        except Exception as e:
            tqdm.write(f"❌ Failed {gene}: {e}")
            results[gene] = None

    return results
```

### 3. Memory-Efficient Processing

```python
def stream_large_fasta(file_path, chunk_size=1000):
    """Process large FASTA files in chunks."""

    from Bio import SeqIO

    with open(file_path, 'r') as handle:
        sequences = []

        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))

            if len(sequences) >= chunk_size:
                yield sequences
                sequences = []

        # Yield remaining sequences
        if sequences:
            yield sequences

# Usage
designer = SiRNADesigner(DesignParameters())
all_results = []

for chunk in stream_large_fasta("large_transcripts.fasta", chunk_size=100):
    chunk_results = designer.design_from_sequences(chunk)
    all_results.extend(chunk_results.candidates)

    # Optional: Save intermediate results
    print(f"Processed chunk, total candidates: {len(all_results)}")
```

## Next Steps

After mastering the Python API:

1. **[Usage Examples](../usage_examples.md)** - Complex multi-step analyses and automation
2. **[Pipeline Integration](pipeline_integration.md)** - Nextflow pipeline development
3. **[Custom Scoring](custom_scoring.md)** - Advanced algorithm development
4. **{doc}`API Reference <../api_reference>`** - Complete class and method documentation
5. **[Developer Guide](../developer/development.md)** - Contributing to siRNAforge

## Key API Modules

For complete module documentation, see the {doc}`API Reference <../api_reference>`:

- **Core Design**: `sirnaforge.core.design` - Main design engine
- **Models**: `sirnaforge.models.sirna` - Data models and schemas
- **Gene Search**: `sirnaforge.data.gene_search` - Database access
- **Thermodynamics**: `sirnaforge.core.thermodynamics` - RNA folding analysis
- **Workflow**: `sirnaforge.workflow` - High-level orchestration


The Python API provides the flexibility to create sophisticated, automated siRNA design workflows tailored to your specific research needs.
