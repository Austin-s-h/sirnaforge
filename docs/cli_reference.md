# CLI Reference

## Commands Overview

```bash
sirnaforge [command] [options]
```

| Command | Description |
|---------|-------------|
| `workflow` | Complete gene → siRNA pipeline |
| `search` | Find gene transcripts |
| `design` | Generate siRNA candidates |
| `validate` | Check FASTA format |
| `version` | Show version |
| `config` | Show settings |

## Commands

### `search`
```bash

 Usage: sirnaforge search [OPTIONS] QUERY

 Search for gene transcripts and retrieve sequences.

┌─ Arguments ─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ *    query      TEXT  Gene ID, gene name, or transcript ID to search for [required]                                 │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
┌─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ --output             -o                            PATH  Output FASTA file for transcript sequences                 │
│                                                          [default: transcripts.fasta]                               │
│ --database           -d                            TEXT  Database to search (ensembl, refseq, gencode)              │
│                                                          [default: ensembl]                                         │
│ --all                -a                                  Search all databases                                       │
│ --no-sequence                                            Skip sequence retrieval (metadata only)                    │
│ --canonical-only                                         Extract only canonical isoforms                            │
│ --extract-canonical      --no-extract-canonical          Automatically extract canonical isoforms to separate file  │
│                                                          [default: extract-canonical]                               │
│ --types              -t                            TEXT  Comma-separated list of transcript types to include (e.g., │
│                                                          protein_coding,lncRNA)                                     │
│                                                          [default: protein_coding,lncRNA]                           │
│ --exclude-types                                    TEXT  Comma-separated list of transcript types to exclude        │
│                                                          [default: nonsense_mediated_decay,retained_intron]         │
│ --verbose            -v                                  Enable verbose output                                      │
│ --help                                                   Show this message and exit.                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### `workflow`
```bash

 Usage: sirnaforge workflow [OPTIONS] GENE_QUERY

 Run complete siRNA design workflow from gene query to off-target analysis.

┌─ Arguments ─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ *    gene_query      TEXT  Gene name or ID to analyze [required]                                                     │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
┌─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ --input-fasta                                      TEXT  Local FASTA path or remote URI (http/https/ftp) to use     │
│                                                          instead of performing a gene search                         │
│ --output-dir         -o                            PATH  Output directory for all workflow results                  │
│                                                          [default: sirna_workflow_output]                           │
│ --database           -d                            TEXT  Database to search (ensembl, refseq, gencode)              │
│                                                          [default: ensembl]                                         │
│ --top-n              -n                            INTEGER RANGE [1<=x]  Number of top siRNA candidates to generate │
│                                                          [default: 20]
│                                                          analysis [default: 10]                                     │
│ --genome-species                                   TEXT  Comma-separated list of genome species for off-target      │
│                                                          analysis [default: human,rat,rhesus]                       │
│ --gc-min                                           FLOAT RANGE [0.0<=x<=100.0]  Minimum GC content percentage      │
│                                                          [default: 30.0]                                            │
│ --gc-max                                           FLOAT RANGE [0.0<=x<=100.0]  Maximum GC content percentage      │
│                                                          [default: 52.0]                                            │
│ --length             -l                            INTEGER RANGE [19<=x<=23]  siRNA length in nucleotides          │
│                                                          [default: 21]                                              │
│ --verbose            -v                                  Enable verbose output                                      │
│ --log-file                                         PATH  Path to centralized log file (overrides                   │
│                                                          SIRNAFORGE_LOG_FILE env)                                   │
│ --help                                                   Show this message and exit.                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### `design`
```bash

 Usage: sirnaforge design [OPTIONS] INPUT_FILE

 Design siRNA candidates from transcript sequences.

┌─ Arguments ─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ *    input_file      PATH  Input FASTA file containing transcript sequences [required]                              │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
┌─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ --output             -o                            PATH  Output file for siRNA candidates                           │
│                                                          [default: sirna_results.tsv]                              │
│ --length             -l                            INTEGER RANGE [19<=x<=23]  siRNA length in nucleotides          │
│                                                          [default: 21]                                              │
│ --top-n              -n                            INTEGER RANGE [1<=x]  Number of top candidates to return        │
│                                                          [default: 10]                                              │
│ --gc-min                                           FLOAT RANGE [0.0<=x<=100.0]  Minimum GC content percentage      │
│                                                          [default: 30.0]                                            │
│ --gc-max                                           FLOAT RANGE [0.0<=x<=100.0]  Maximum GC content percentage      │
│                                                          [default: 52.0]                                            │
│ --max-poly-runs                                    INTEGER RANGE [1<=x]  Maximum consecutive identical nucleotides  │
│                                                          [default: 3]                                               │
│ --genome-index                                     PATH  Genome index for off-target analysis                      │
│ --snp-file                                         PATH  VCF file with SNPs to avoid                               │
│ --skip-structure                                         Skip secondary structure prediction (faster)               │
│ --skip-off-targets                                       Skip off-target analysis (faster)                         │
│ --verbose            -v                                  Enable verbose output                                      │
│ --help                                                   Show this message and exit.                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### `validate`
```bash

 Usage: sirnaforge validate [OPTIONS] INPUT_FILE

 Validate input FASTA file format and content.

┌─ Arguments ─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ *    input_file      PATH  FASTA file to validate [required]                                                        │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
┌─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ --help          Show this message and exit.                                                                         │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### `version`
```bash

 Usage: sirnaforge version [OPTIONS]

 Show version information.

┌─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ --help          Show this message and exit.                                                                         │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### `config`
```bash

 Usage: sirnaforge config [OPTIONS]

 Show default configuration parameters.

┌─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ --help          Show this message and exit.                                                                         │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

## Examples

### Basic Gene Search
```bash
# Search for TP53 gene transcripts
sirnaforge search TP53

# Search with custom output and database
sirnaforge search BRCA1 --output brca1_transcripts.fasta --database ensembl

# Search all databases
sirnaforge search MYC --all
```

### siRNA Design from Transcripts
```bash
# Design siRNAs from transcript file
sirnaforge design transcripts.fasta --output results.tsv

# Custom parameters
sirnaforge design input.fasta --length 19 --gc-min 40 --gc-max 60 --top-n 20
```

### Full Workflow
```bash
# Complete workflow from gene name
sirnaforge workflow TP53 --output-dir tp53_analysis

# Workflow with custom parameters
sirnaforge workflow BRCA1 --top-n 30 --gc-min 35 --gc-max 55
```

### File Validation
```bash
# Validate FASTA file
sirnaforge validate my_sequences.fasta
```
