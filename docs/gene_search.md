# Gene Search Functionality

The siRNAforge toolkit includes comprehensive gene search functionality that allows you to retrieve transcript sequences from multiple genomic databases.ene Search Functionality 🧬

The siRNA Design Toolkit includes comprehensive gene search functionality that allows you to retrieve transcript sequences from multiple genomic databases.

## Features

- 🔍 **Multi-database support**: Search Ensembl, RefSeq, and GENCODE databases
- 🧬 **Comprehensive metadata**: Gene info, transcript types, canonical annotations
- 📋 **FASTA output**: Save sequences in standard FASTA format
- ⚡ **Async/sync APIs**: Both asynchronous and synchronous interfaces
- 🎯 **CLI integration**: Rich command-line interface with progress indicators

## Quick Start

### Command Line Usage

```bash
# Search for a gene and save sequences
uv run sirna search TP53 -o tp53_transcripts.fasta

# Search specific database
uv run sirna search BRCA1 -d ensembl -o brca1_ensembl.fasta

# Search all databases
uv run sirna search MYC --all -o myc_all_databases.fasta

# Get metadata only (no sequences)
uv run sirna search EGFR --no-sequence
```

### Python API Usage

#### Synchronous API

```python
from sirna_design.data.gene_search import search_gene_sync, DatabaseType

# Simple gene search
result = search_gene_sync("TP53", database=DatabaseType.ENSEMBL)

if result.success:
    print(f"Found {len(result.transcripts)} transcripts for {result.gene_info.gene_name}")
    for transcript in result.transcripts:
        print(f"  {transcript.transcript_id}: {len(transcript.sequence)} bp")
```

#### Asynchronous API

```python
import asyncio
from sirna_design.data.gene_search import GeneSearcher, DatabaseType

async def search_genes():
    searcher = GeneSearcher()
    
    # Single database search
    result = await searcher.search_gene("BRCA1", DatabaseType.ENSEMBL)
    
    # Multi-database search
    results = await searcher.search_multiple_databases(
        "MYC", 
        databases=[DatabaseType.ENSEMBL, DatabaseType.REFSEQ]
    )
    
    return results

# Run the search
results = asyncio.run(search_genes())
```

## Supported Databases

### Ensembl ✅ (Fully Implemented)
- **Base URL**: https://rest.ensembl.org
- **Features**: Gene lookup, transcript retrieval, sequence fetching
- **Data**: Human genome (homo_sapiens)

### RefSeq 🚧 (Planned)
- **Base URL**: https://eutils.ncbi.nlm.nih.gov/entrez/eutils
- **Status**: Placeholder implementation
- **Features**: Will support NCBI E-utilities integration

### GENCODE 🚧 (Planned)
- **Base URL**: https://www.gencodegenes.org
- **Status**: Placeholder implementation
- **Features**: Will support GENCODE release data

## Data Models

### GeneInfo
```python
class GeneInfo(BaseModel):
    gene_id: str                    # e.g., "ENSG00000141510"
    gene_name: Optional[str]        # e.g., "TP53"
    gene_type: Optional[str]        # e.g., "protein_coding"
    chromosome: Optional[str]       # e.g., "17"
    start: Optional[int]            # Genomic start position
    end: Optional[int]              # Genomic end position
    strand: Optional[int]           # 1 or -1
    description: Optional[str]      # Gene description
    database: DatabaseType          # Source database
```

### TranscriptInfo
```python
class TranscriptInfo(BaseModel):
    transcript_id: str              # e.g., "ENST00000269305"
    transcript_name: Optional[str]  # Display name
    transcript_type: Optional[str]  # e.g., "protein_coding"
    gene_id: str                    # Parent gene ID
    gene_name: Optional[str]        # Parent gene name
    sequence: Optional[str]         # cDNA sequence
    length: Optional[int]           # Sequence length
    database: DatabaseType          # Source database
    is_canonical: bool              # Canonical transcript flag
```

## Example Workflows

### 1. Design siRNAs from Gene Search

```python
from sirna_design.data.gene_search import search_gene_sync
from sirna_design.core.design import SiRNADesigner
from sirna_design.models.sirna import DesignParameters

# Search for gene
result = search_gene_sync("TP53")

if result.success:
    # Save transcripts to FASTA
    searcher = GeneSearcher()
    searcher.save_transcripts_fasta(result.transcripts, "tp53.fasta")
    
    # Design siRNAs
    designer = SiRNADesigner(DesignParameters())
    sirna_results = designer.design_from_file("tp53.fasta")
    
    print(f"Designed {len(sirna_results.candidates)} siRNA candidates")
```

### 2. Compare Databases

```python
from sirna_design.data.gene_search import search_multiple_databases_sync

# Compare gene across databases
results = search_multiple_databases_sync("BRCA1")

for result in results:
    if result.success:
        canonical_count = sum(1 for t in result.transcripts if t.is_canonical)
        print(f"{result.database.value}: {len(result.transcripts)} transcripts "
              f"({canonical_count} canonical)")
```

## Error Handling

The gene search functionality includes robust error handling:

```python
result = search_gene_sync("NONEXISTENT_GENE")

if not result.success:
    print(f"Search failed: {result.error}")
    # Handle the error appropriately
else:
    # Process successful result
    print(f"Found {len(result.transcripts)} transcripts")
```

## Configuration

### Searcher Options

```python
searcher = GeneSearcher(
    preferred_db=DatabaseType.ENSEMBL,  # Default database
    timeout=60,                         # Request timeout (seconds)
    max_retries=3                       # Maximum retry attempts
)
```

### FASTA Output Options

```python
searcher.save_transcripts_fasta(
    transcripts=transcripts,
    output_path="sequences.fasta",
    include_metadata=True  # Include gene info in headers
)
```

## Integration with siRNA Design

The gene search functionality integrates seamlessly with the siRNA design pipeline:

```bash
# Complete pipeline: search → design → analyze
sirna search TP53 -o tp53.fasta
sirnaforge design tp53.fasta -o tp53_sirnas.tsv --top-n 20
```

## Testing

Run the gene search tests:

```bash
# Run all gene search tests
uv run pytest tests/unit/test_gene_search.py -v

# Run with coverage
uv run pytest tests/unit/test_gene_search.py --cov=sirnaforge.data.gene_search
```

## Contributing

To add support for additional databases:

1. Implement the `_search_<database>` method in `GeneSearcher`
2. Add database configuration to `db_configs`
3. Add the new database type to `DatabaseType` enum
4. Write comprehensive tests

See the existing Ensembl implementation as a reference.
