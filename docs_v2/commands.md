# Commands

Complete reference for all siRNAforge CLI commands.

## workflow

Run complete siRNA design pipeline from gene query to scored candidates.

```bash
sirnaforge workflow GENE [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--output-dir` | `sirna_workflow_output` | Output directory |
| `--gc-min` | 30 | Minimum GC content (%) |
| `--gc-max` | 60 | Maximum GC content (%) |
| `--length` | 21 | siRNA length (nt) |
| `--top-n` | 20 | Number of candidates to retain |
| `--verbose` | off | Show detailed progress |

**Examples:**
```bash
sirnaforge workflow TP53
sirnaforge workflow BRCA1 --gc-min 35 --gc-max 55 --top-n 30
sirnaforge workflow EGFR --output-dir egfr_analysis --verbose
```

---

## search

Search gene databases and retrieve transcript sequences.

```bash
sirnaforge search QUERY [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | `transcripts.fasta` | Output FASTA file |
| `-d, --database` | `ensembl` | Database (ensembl/refseq/gencode) |
| `-a, --all` | off | Search all databases |
| `--canonical-only` | off | Only canonical isoforms |
| `-t, --types` | `protein_coding,lncRNA` | Transcript types to include |
| `--verbose` | off | Show detailed progress |

**Examples:**
```bash
sirnaforge search TP53 -o tp53.fasta
sirnaforge search HOTAIR --database ensembl --types lncRNA
sirnaforge search BRCA1 --all --canonical-only
```

---

## design

Design siRNA candidates from FASTA sequences.

```bash
sirnaforge design INPUT [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | `sirna_results.tsv` | Output file |
| `--gc-min` | 30 | Minimum GC content (%) |
| `--gc-max` | 60 | Maximum GC content (%) |
| `--length` | 21 | siRNA length (nt) |
| `--top-n` | 10 | Number of candidates to retain |
| `--skip-structure` | off | Skip secondary structure prediction |
| `--verbose` | off | Show detailed progress |

**Examples:**
```bash
sirnaforge design transcripts.fasta -o results.csv
sirnaforge design input.fasta --gc-min 40 --gc-max 55 --top-n 25
sirnaforge design large_file.fasta --skip-structure  # Faster
```

---

## validate

Check FASTA file format and content.

```bash
sirnaforge validate INPUT [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--verbose` | off | Show detailed validation results |

**Examples:**
```bash
sirnaforge validate input.fasta
sirnaforge validate transcripts.fasta --verbose
```

---

## sequences

Manage siRNA sequences and chemical modification metadata.

### sequences show

Display sequences with modification metadata.

```bash
sirnaforge sequences show FILE [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--id` | - | Show specific sequence by ID |
| `--format` | `table` | Output format (table/json/fasta) |

### sequences annotate

Merge modification metadata into FASTA headers.

```bash
sirnaforge sequences annotate FASTA METADATA -o OUTPUT
```

**Examples:**
```bash
sirnaforge sequences show candidates.fasta
sirnaforge sequences show candidates.fasta --format json
sirnaforge sequences annotate candidates.fasta mods.json -o annotated.fasta
```

---

## config

Show default configuration parameters.

```bash
sirnaforge config
```

---

## cache

Manage miRNA database cache for off-target analysis.

```bash
sirnaforge cache [COMMAND]
```

| Command | Description |
|---------|-------------|
| `info` | Show cache location and status |
| `clear` | Clear cached databases |

---

## version

Show version information.

```bash
sirnaforge version
```
