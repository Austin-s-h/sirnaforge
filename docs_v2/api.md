# Python API

Use siRNAforge programmatically for custom workflows and automation.

## Basic Usage

```python
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters

# Configure parameters
params = DesignParameters(
    sirna_length=21,
    top_candidates=20,
    gc_content_range=(35.0, 55.0)
)

# Design from sequences
designer = SiRNADesigner(params)
results = designer.design_from_sequences([
    "AUGCGCGAUCUCGAUGCAUGUGCGAUCGAUGCGUAUCGAU"
])

# Access results
for candidate in results.top_candidates:
    print(f"{candidate.sirna_id}: {candidate.guide_sequence}")
    print(f"  Score: {candidate.composite_score:.1f}")
    print(f"  GC: {candidate.gc_content:.1f}%")
```

## Gene Search

```python
from sirnaforge.data.gene_search import search_gene_sync, DatabaseType

# Search for gene
result = search_gene_sync("TP53", database=DatabaseType.ENSEMBL)

if result.success:
    print(f"Found {len(result.transcripts)} transcripts")
    for t in result.transcripts:
        print(f"  {t.transcript_id}: {len(t.sequence)} bp")
```

### Async Search

```python
import asyncio
from sirnaforge.data.gene_search import GeneSearcher

async def search_genes():
    searcher = GeneSearcher()
    results = {}
    for gene in ["TP53", "BRCA1", "EGFR"]:
        result = await searcher.search_gene(gene)
        results[gene] = result
    return results

results = asyncio.run(search_genes())
```

## Complete Workflow

```python
from sirnaforge.workflow import run_sirna_workflow
import asyncio

async def analyze_gene(gene: str):
    results = await run_sirna_workflow(
        gene_query=gene,
        output_dir=f"results/{gene}",
        top_n_candidates=20,
        gc_min=35.0,
        gc_max=55.0
    )
    return results

results = asyncio.run(analyze_gene("TP53"))
print(f"Generated {len(results['candidates'])} candidates")
```

## Key Classes

### DesignParameters

```python
from sirnaforge.models.sirna import DesignParameters

params = DesignParameters(
    sirna_length=21,           # 19-23 nt
    top_candidates=20,         # Number to return
    gc_content_range=(30, 60), # Min/max GC %
)
```

### SiRNACandidate

Each candidate includes:

```python
candidate.sirna_id          # Unique identifier
candidate.guide_sequence    # 21nt guide strand
candidate.passenger_sequence # Passenger strand
candidate.position          # Position in transcript
candidate.composite_score   # Overall score (0-10)
candidate.asymmetry_score   # Thermodynamic asymmetry
candidate.gc_content        # GC percentage
candidate.melting_temp_c    # Melting temperature
```

### Chemical Modifications

```python
from sirnaforge.models.modifications import (
    StrandMetadata,
    ChemicalModification,
    Provenance,
    SourceType
)

metadata = StrandMetadata(
    id="sirna_001",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 3, 5])
    ],
    provenance=Provenance(
        source_type=SourceType.DESIGNED,
        identifier="my_experiment_001"
    )
)
```

## File I/O

### Save Results

```python
# Save to CSV
results.save_csv("candidates.csv")

# Save to FASTA
results.save_fasta("candidates.fasta")
```

### Load Sequences

```python
from pathlib import Path
from sirnaforge.core.design import SiRNADesigner

designer = SiRNADesigner(params)
results = designer.design_from_file(Path("transcripts.fasta"))
```

## Error Handling

```python
from sirnaforge.core.design import DesignException
from sirnaforge.data.gene_search import GeneNotFoundError

try:
    result = search_gene_sync("INVALID_GENE")
except GeneNotFoundError as e:
    print(f"Gene not found: {e}")

try:
    results = designer.design_from_sequences([])
except DesignException as e:
    print(f"Design failed: {e}")
```

## See Also

- [Workflows](workflows.md) - CLI workflows
- [Scoring](scoring.md) - Understanding metrics
