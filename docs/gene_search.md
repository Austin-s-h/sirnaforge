# Gene Search Functionality

TL;DR
-----

- Retrieve transcript sequences (FASTA + metadata) from Ensembl today; RefSeq and GENCODE are planned.
- Use `sirnaforge search <QUERY>` from the CLI or the Python API for programmatic use.
- Output plays nicely with the `sirnaforge design` step and the full `sirnaforge workflow`.

The siRNAforge toolkit includes comprehensive gene search functionality to retrieve transcript sequences and metadata from multiple genomic databases.

## Features

- 🔍 **Multi-database support**: Search Ensembl, RefSeq, and GENCODE databases
- 🧬 **Comprehensive metadata**: Gene info, transcript types, canonical annotations
- 📋 **FASTA output**: Save sequences in standard FASTA format
- ⚡ **Async/sync APIs**: Both asynchronous and synchronous interfaces
- 🎯 **CLI integration**: Rich command-line interface with progress indicators

## Quick Start

### Command Line Usage

`````{tab-set}

````{tab-item} uv
```bash
# Search for a gene and save sequences
uv run sirnaforge search TP53 -o tp53_transcripts.fasta

# Search specific database
uv run sirnaforge search BRCA1 -d ensembl -o brca1_ensembl.fasta

# Search all databases
uv run sirnaforge search MYC --all -o myc_all_databases.fasta

# Get metadata only (no sequences)
uv run sirnaforge search EGFR --no-sequence

# Use a FASTA from disk or remote URL as input (skip gene search)
uv run sirnaforge workflow TP53 --input-fasta /path/to/transcripts.fasta -o sirna_output
uv run sirnaforge workflow TP53 --input-fasta https://example.org/transcripts.fasta -o sirna_output
```
````

````{tab-item} Docker
```bash
# Search for a gene and save sequences
docker run --rm -v $(pwd):/workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search TP53 -o tp53_transcripts.fasta

# Search specific database
docker run --rm -v $(pwd):/workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search BRCA1 -d ensembl -o brca1_ensembl.fasta

# Search all databases
docker run --rm -v $(pwd):/workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search MYC --all -o myc_all_databases.fasta

# Get metadata only (no sequences)
docker run --rm ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge search EGFR --no-sequence

# Use a FASTA from disk or remote URL as input (skip gene search)
docker run --rm -v $(pwd):/workspace -w /workspace \
  ghcr.io/austin-s-h/sirnaforge:latest \
  sirnaforge workflow TP53 --input-fasta /workspace/transcripts.fasta -o sirna_output
```
````

`````

### Python API Usage

#### Synchronous API

```python
from sirnaforge.data.gene_search import search_gene_sync, DatabaseType

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
from sirnaforge.data.gene_search import GeneSearcher, DatabaseType

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

The gene search system uses Pydantic models for type-safe data handling and validation.

> **API Reference**: See complete model documentation at:
> - [`GeneInfo`](api_reference.rst) - Gene metadata
> - [`TranscriptInfo`](api_reference.rst) - Transcript data and sequences
> - [`GeneSearchResult`](api_reference.rst) - Search result container
> - [`DatabaseType`](api_reference.rst) - Database enumerations

**Key Fields:**
- **GeneInfo**: gene_id, gene_name, gene_type, chromosome, coordinates, description
- **TranscriptInfo**: transcript_id, sequence, length, is_canonical, database
- **GeneSearchResult**: gene_info, transcripts, success, error

## Example Workflows

### 1. Design siRNAs from Gene Search

```python
from sirnaforge.data.gene_search import search_gene_sync, GeneSearcher
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters

# Search for gene
result = search_gene_sync("TP53")

if result.success:
    # Save transcripts to FASTA
    searcher = GeneSearcher()
    searcher.save_transcripts_fasta(result.transcripts, "tp53.fasta")

    # Design siRNAs with thermodynamic analysis
    designer = SiRNADesigner(DesignParameters())
    sirna_results = designer.design_from_file("tp53.fasta")

    print(f"Designed {len(sirna_results.candidates)} siRNA candidates")
```

> **See Also**: [`GeneSearcher.save_transcripts_fasta`](api_reference.rst) for FASTA output options

### 2. Compare Databases

```python
from sirnaforge.data.gene_search import search_multiple_databases_sync

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
    if result.is_access_error:
        print("Database access issue - check network/permissions")
    # Handle the error appropriately
else:
    # Process successful result
    print(f"Found {len(result.transcripts)} transcripts")
```

> **Exception Types**: See [`DatabaseAccessError`](api_reference.rst) and [`GeneNotFoundError`](api_reference.rst) in the API docs

## Configuration

Configuration is handled through constructor parameters. See [`GeneSearcher.__init__`](api_reference.rst) for all options.

```python
searcher = GeneSearcher(
    timeout=60,         # Request timeout (seconds)
    max_retries=3      # Maximum retry attempts
)
```

## See Also

- **[Getting Started](getting_started.md)** - Complete workflow including gene search
- **[Usage Examples](usage_examples.md)** - More search examples and batch processing
- **[API Reference](api_reference.rst)** - Full Python API documentation
