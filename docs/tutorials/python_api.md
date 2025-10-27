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

### Basic API Usage

#### 1. Simple siRNA Design

```python
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters

# Configure design parameters
params = DesignParameters(
    sirna_length=21,
    top_candidates=20,
    gc_content_range=(30.0, 52.0)
)

# Create designer instance
designer = SiRNADesigner(params)

# Design from sequences
sequences = [
    "AUGCGCGAUCUCGAUGCAUGUGCGAUCGAUGCGUAUCGAU",
    "GCCAUGCGAUCGAUGCGUAUCAUCGAUGCGCGAUGCUGUGAU"
]

results = designer.design_from_sequences(sequences)

# Access results
print(f"Generated {len(results.candidates)} siRNA candidates")
for candidate in results.top_candidates:
    print(f"{candidate.sirna_id}: {candidate.guide_sequence} (score: {candidate.composite_score:.1f})")
```

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

## Advanced API Usage

### Custom Scoring

```python
from sirnaforge.core.design import BaseScorer
from sirnaforge.models.sirna import SiRNACandidate

class CustomScorer(BaseScorer):
    """Custom scoring algorithm for specific research needs."""

    def __init__(self, weight_gc=0.3, weight_position=0.2, weight_structure=0.5):
        self.weight_gc = weight_gc
        self.weight_position = weight_position
        self.weight_structure = weight_structure

    def calculate_score(self, candidate: SiRNACandidate) -> float:
        """Calculate custom composite score."""

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

        return composite

# Use custom scorer
custom_scorer = CustomScorer(weight_gc=0.4, weight_structure=0.6)
designer = SiRNADesigner(params, scorer=custom_scorer)
```

### Batch Processing

```python
from concurrent.futures import ProcessPoolExecutor
from functools import partial

def process_gene(gene_name, output_dir):
    """Process a single gene."""
    try:
        # Run workflow for gene
        result = asyncio.run(run_sirna_workflow(
            gene=gene_name,
            output_dir=f"{output_dir}/{gene_name}",
            top_candidates=20
        ))
        return {"gene": gene_name, "status": "success", "candidates": len(result.candidates)}
    except Exception as e:
        return {"gene": gene_name, "status": "error", "error": str(e)}

# Batch process multiple genes
genes = ["TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "PIK3CA", "AKT1", "KRAS"]
output_directory = "batch_analysis"

# Process in parallel
with ProcessPoolExecutor(max_workers=4) as executor:
    process_func = partial(process_gene, output_dir=output_directory)
    results = list(executor.map(process_func, genes))

# Summary of results
for result in results:
    if result["status"] == "success":
        print(f"✅ {result['gene']}: {result['candidates']} candidates")
    else:
        print(f"❌ {result['gene']}: {result['error']}")
```

## Data Analysis Integration

### Integration with Pandas

```python
import pandas as pd
from sirnaforge.core.design import SiRNADesigner

# Load existing data
gene_list = pd.read_csv("target_genes.csv")

results_data = []

for _, row in gene_list.iterrows():
    gene_name = row['gene_name']

    try:
        # Search for gene
        search_result = search_gene_sync(gene_name)

        if search_result.success:
            # Design siRNAs
            designer = SiRNADesigner(DesignParameters(top_candidates=10))
            design_result = designer.design_from_sequences([
                t.sequence for t in search_result.transcripts if t.sequence
            ])

            # Add to results
            for candidate in design_result.candidates:
                results_data.append({
                    'gene_name': gene_name,
                    'sirna_id': candidate.sirna_id,
                    'guide_sequence': candidate.guide_sequence,
                    'gc_content': candidate.gc_content,
                    'composite_score': candidate.composite_score,
                    'transcript_count': len(search_result.transcripts)
                })

    except Exception as e:
        print(f"Error processing {gene_name}: {e}")

# Create DataFrame
results_df = pd.DataFrame(results_data)

# Analysis and visualization
print(f"Processed {results_df['gene_name'].nunique()} genes")
print(f"Generated {len(results_df)} total siRNA candidates")

# Top candidates by gene
top_per_gene = results_df.groupby('gene_name').apply(
    lambda x: x.nlargest(3, 'composite_score')
).reset_index(drop=True)

print("\nTop 3 candidates per gene:")
print(top_per_gene[['gene_name', 'sirna_id', 'composite_score']])
```

### Integration with Jupyter Notebooks

```python
# Jupyter notebook cell
%matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

# Design siRNAs for analysis
results = designer.design_from_file("transcripts.fasta")

# Convert to DataFrame for analysis
df = pd.DataFrame([
    {
        'sirna_id': c.sirna_id,
        'gc_content': c.gc_content,
        'composite_score': c.composite_score,
        'thermodynamic_score': c.thermodynamic_score,
        'position': c.position
    }
    for c in results.candidates
])

# Visualize results
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# GC content distribution
sns.histplot(df['gc_content'], bins=20, ax=axes[0,0])
axes[0,0].set_title('GC Content Distribution')

# Score correlation
sns.scatterplot(data=df, x='thermodynamic_score', y='composite_score', ax=axes[0,1])
axes[0,1].set_title('Score Correlation')

# Position vs Score
sns.scatterplot(data=df, x='position', y='composite_score', ax=axes[1,0])
axes[1,0].set_title('Position vs Composite Score')

# Score distribution
sns.boxplot(data=df[['composite_score', 'thermodynamic_score']], ax=axes[1,1])
axes[1,1].set_title('Score Distributions')

plt.tight_layout()
plt.show()
```

## Error Handling and Validation

### Comprehensive Error Handling

```python
from sirnaforge.core.exceptions import (
    SiRNAForgeException, DesignException, ValidationException
)

def robust_sirna_design(gene_name, max_retries=3):
    """Robust siRNA design with error handling and retries."""

    for attempt in range(max_retries):
        try:
            # Search for gene
            print(f"Attempt {attempt + 1}: Searching for {gene_name}")
            search_result = search_gene_sync(gene_name)

            if not search_result.success:
                raise DesignException(f"Gene search failed: {search_result.error}")

            # Validate sequences
            valid_transcripts = []
            for transcript in search_result.transcripts:
                if transcript.sequence and len(transcript.sequence) > 100:
                    valid_transcripts.append(transcript)

            if not valid_transcripts:
                raise ValidationException("No valid transcripts found")

            # Design siRNAs
            designer = SiRNADesigner(DesignParameters(top_candidates=20))
            results = designer.design_from_sequences([
                t.sequence for t in valid_transcripts
            ])

            if not results.candidates:
                raise DesignException("No siRNA candidates generated")

            print(f"✅ Successfully designed {len(results.candidates)} candidates for {gene_name}")
            return results

        except ValidationException as e:
            print(f"❌ Validation error for {gene_name}: {e}")
            break  # Don't retry validation errors

        except (DesignException, SiRNAForgeException) as e:
            print(f"⚠️  Attempt {attempt + 1} failed for {gene_name}: {e}")
            if attempt == max_retries - 1:
                print(f"❌ All attempts failed for {gene_name}")
                return None

        except Exception as e:
            print(f"💥 Unexpected error for {gene_name}: {e}")
            if attempt == max_retries - 1:
                print(f"❌ All attempts failed for {gene_name}")
                return None

    return None

# Use robust design
genes = ["TP53", "BRCA1", "INVALID_GENE", "EGFR"]
for gene in genes:
    result = robust_sirna_design(gene)
    if result:
        print(f"Best candidate for {gene}: {result.top_candidates[0].guide_sequence}")
```

## Configuration and Customization

### Custom Configuration

```python
from sirnaforge.models.sirna import DesignParameters, FilterCriteria
from pydantic import BaseSettings

class CustomConfig(BaseSettings):
    """Custom configuration with environment variable support."""

    # Default parameters
    default_sirna_length: int = 21
    default_candidates: int = 20
    gc_min_default: float = 30.0
    gc_max_default: float = 52.0

    # API settings
    ensembl_base_url: str = "https://rest.ensembl.org"
    request_timeout: int = 30

    # Quality thresholds
    min_composite_score: float = 5.0
    max_off_targets: int = 5

    class Config:
        env_prefix = "SIRNAFORGE_"

# Load configuration
config = CustomConfig()

# Use in design
params = DesignParameters(
    sirna_length=config.default_sirna_length,
    top_candidates=config.default_candidates,
    filters=FilterCriteria(
        gc_min=config.gc_min_default,
        gc_max=config.gc_max_default,
        min_composite_score=config.min_composite_score
    )
)
```

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

1. **{doc}`../usage_examples`** - Complex multi-step analyses
2. **{doc}`pipeline_integration`** - Nextflow pipeline development
3. **{doc}`custom_scoring`** - Advanced algorithm development
4. **{doc}`../developer/development`** - Contributing to siRNAforge development

The Python API provides the flexibility to create sophisticated, automated siRNA design workflows tailored to your specific research needs.
