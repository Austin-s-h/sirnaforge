# `sirnaforge`

siRNAforge - siRNA design toolkit for gene silencing

**Usage**:

```console
$ sirnaforge [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `search`: Search for gene transcripts and retrieve...
* `workflow`: Run complete siRNA design workflow from...
* `design`: Design siRNA candidates from transcript...
* `validate`: Validate input FASTA file format and content.
* `version`: Show version information.
* `config`: Show default configuration parameters.
* `cache`: Manage miRNA database cache.
* `sequences`: Manage siRNA sequences and metadata

## `sirnaforge search`

Search for gene transcripts and retrieve sequences.

**Usage**:

```console
$ sirnaforge search [OPTIONS] QUERY
```

**Arguments**:

* `QUERY`: Gene ID, gene name, or transcript ID to search for  [required]

**Options**:

* `-o, --output PATH`: Output FASTA file for transcript sequences  [default: transcripts.fasta]
* `-d, --database TEXT`: Database to search (ensembl, refseq, gencode)  [default: ensembl]
* `-a, --all`: Search all databases
* `--fallback / --no-fallback`: Enable automatic fallback to other databases if access is blocked  [default: fallback]
* `--no-sequence`: Skip sequence retrieval (metadata only)
* `--canonical-only`: Extract only canonical isoforms
* `--extract-canonical / --no-extract-canonical`: Automatically extract canonical isoforms to separate file  [default: extract-canonical]
* `-t, --types TEXT`: Comma-separated list of transcript types to include (e.g., protein_coding,lncRNA)  [default: protein_coding,lncRNA]
* `--exclude-types TEXT`: Comma-separated list of transcript types to exclude  [default: nonsense_mediated_decay,retained_intron]
* `-v, --verbose`: Enable verbose output
* `--help`: Show this message and exit.

## `sirnaforge workflow`

Run complete siRNA design workflow from gene query to off-target analysis.

**Usage**:

```console
$ sirnaforge workflow [OPTIONS] GENE_QUERY
```

**Arguments**:

* `GENE_QUERY`: Gene name or ID to analyze  [required]

**Options**:

* `--input-fasta TEXT`: Local path or remote URI to an input FASTA file (http/https/ftp)
* `-o, --output-dir PATH`: Output directory for all workflow results  [default: sirna_workflow_output]
* `-d, --database TEXT`: Database to search (ensembl, refseq, gencode)  [default: ensembl]
* `--design-mode TEXT`: Design mode: sirna (default) or mirna (miRNA-biogenesis-aware)  [default: sirna]
* `-n, --top-n INTEGER RANGE`: Number of top siRNA candidates to select (also used for off-target analysis)  [default: 20; x&gt;=1]
* `--species, --genome-species TEXT`: Comma-separated canonical species identifiers (genome+miRNA). Supported values include human, mouse, rhesus, rat, chicken  [default: human,mouse,rhesus,rat,chicken]
* `--mirna-db TEXT`: miRNA reference database to use for seed analysis  [default: mirgenedb]
* `--mirna-species TEXT`: Optional comma-separated override for miRNA species identifiers. Defaults to mapping the --species selections
* `--transcriptome-fasta TEXT`: Path or URL to transcriptome FASTA (local file, URL, or pre-configured source such as `ensembl_human_cdna`). Cached and indexed automatically; defaults to `ensembl_human_cdna` when omitted
* `--offtarget-indices TEXT`: Comma-separated overrides for genome indices used in off-target analysis. Format: `human:/abs/path/GRCh38,mouse:/abs/path/GRCm39`. When provided, these replace cached/default genome references.
* `--gc-min FLOAT RANGE`: Minimum GC content percentage  [default: 30.0; 0.0&lt;=x&lt;=100.0]
* `--gc-max FLOAT RANGE`: Maximum GC content percentage  [default: 60.0; 0.0&lt;=x&lt;=100.0]
* `-l, --length INTEGER RANGE`: siRNA length in nucleotides  [default: 21; 19&lt;=x&lt;=23]
* `-m, --modifications TEXT`: Chemical modification pattern (standard_2ome, minimal_terminal, maximal_stability, none)  [default: standard_2ome]
* `--overhang TEXT`: Overhang sequence (dTdT for DNA, UU for RNA)  [default: dTdT]
* `-v, --verbose`: Enable verbose output
* `--log-file PATH`: Path to centralized log file (overrides SIRNAFORGE_LOG_FILE env)
* `--json-summary / --no-json-summary`: Write logs/workflow_summary.json (disable to skip JSON output)  [default: json-summary]
* `--help`: Show this message and exit.

### Input Sources & Transcriptome References

siRNAforge accepts two complementary inputs:

- `--input-fasta` optionally replaces the transcript retrieval step. You can point it at a local file, HTTP(S) URL, or FTP location. The gene query still acts as a logical label for outputs, while the workflow designs guides from the supplied sequences. **When using `--input-fasta` without `--transcriptome-fasta`, transcriptome off-target analysis is disabled** (design-only mode).
- `--transcriptome-fasta` selects the dataset used for transcriptome off-target analysis. It accepts local/remote files as well as shortcuts such as `ensembl_human_cdna`, `ensembl_mouse_cdna`, and other presets listed in `sirnaforge cache --info`. **Explicitly provide this flag to enable transcriptome off-target when using `--input-fasta`.**
- `--offtarget-indices` overrides the genome indices used for off-target analysis with explicit `species:/index_prefix` entries. When present, these override cached/default genome references and drive the genome species selection for Nextflow.

Passing both flags is common: the input FASTA feeds the design engine, while the transcriptome FASTA controls which reference is indexed for the Nextflow/BWA-MEM2 stage. When `--transcriptome-fasta` is omitted the CLI defaults to `ensembl_human_cdna`. Remote files are cached under `~/.cache/sirnaforge/` and reused across runs.

## `sirnaforge design`

Design siRNA candidates from transcript sequences.

**Usage**:

```console
$ sirnaforge design [OPTIONS] INPUT_FILE
```

**Arguments**:

* `INPUT_FILE`: Input FASTA file containing transcript sequences  [required]

**Options**:

* `-o, --output PATH`: Output file for siRNA candidates  [default: sirna_results.tsv]
* `--design-mode TEXT`: Design mode: sirna (default) or mirna (miRNA-biogenesis-aware)  [default: sirna]
* `-l, --length INTEGER RANGE`: siRNA length in nucleotides  [default: 21; 19&lt;=x&lt;=23]
* `-n, --top-n INTEGER RANGE`: Number of top candidates to return  [default: 10; x&gt;=1]
* `--gc-min FLOAT RANGE`: Minimum GC content percentage  [default: 30.0; 0.0&lt;=x&lt;=100.0]
* `--gc-max FLOAT RANGE`: Maximum GC content percentage  [default: 60.0; 0.0&lt;=x&lt;=100.0]
* `--max-poly-runs INTEGER RANGE`: Maximum consecutive identical nucleotides  [default: 3; x&gt;=1]
* `--genome-index PATH`: Genome index for off-target analysis
* `--snp-file PATH`: VCF file with SNPs to avoid
* `--skip-structure`: Skip secondary structure prediction (faster)
* `--skip-off-targets`: Skip off-target analysis (faster)
* `-m, --modifications TEXT`: Chemical modification pattern (standard_2ome, minimal_terminal, maximal_stability, none)  [default: standard_2ome]
* `--overhang TEXT`: Overhang sequence (dTdT for DNA, UU for RNA)  [default: dTdT]
* `-v, --verbose`: Enable verbose output
* `--help`: Show this message and exit.

## `sirnaforge validate`

Validate input FASTA file format and content.

**Usage**:

```console
$ sirnaforge validate [OPTIONS] INPUT_FILE
```

**Arguments**:

* `INPUT_FILE`: FASTA file to validate  [required]

**Options**:

* `--help`: Show this message and exit.

## `sirnaforge version`

Show version information.

**Usage**:

```console
$ sirnaforge version [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.

## `sirnaforge config`

Show default configuration parameters.

**Usage**:

```console
$ sirnaforge config [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.

## `sirnaforge cache`

Manage miRNA database cache.

**Usage**:

```console
$ sirnaforge cache [OPTIONS]
```

**Options**:

* `--clear`: Clear all cached miRNA databases
* `--dry-run`: Show what would be deleted without actually deleting
* `--info`: Show cache information
* `--help`: Show this message and exit.

## `sirnaforge sequences`

Manage siRNA sequences and metadata

**Usage**:

```console
$ sirnaforge sequences [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `show`: Show sequences with their metadata from...
* `annotate`: Merge metadata from JSON into FASTA headers.

### `sirnaforge sequences show`

Show sequences with their metadata from FASTA file.

**Usage**:

```console
$ sirnaforge sequences show [OPTIONS] INPUT_FILE
```

**Arguments**:

* `INPUT_FILE`: FASTA file to display  [required]

**Options**:

* `--id TEXT`: Show only this sequence ID
* `-f, --format TEXT`: Output format (table, json, fasta)  [default: table]
* `--help`: Show this message and exit.

### `sirnaforge sequences annotate`

Merge metadata from JSON into FASTA headers.

**Usage**:

```console
$ sirnaforge sequences annotate [OPTIONS] INPUT_FASTA METADATA_JSON
```

**Arguments**:

* `INPUT_FASTA`: Input FASTA file  [required]
* `METADATA_JSON`: JSON file with metadata  [required]

**Options**:

* `-o, --output PATH`: Output FASTA file (default: &lt;input&gt;_annotated.fasta)
* `-v, --verbose`: Enable verbose output
* `--help`: Show this message and exit.
