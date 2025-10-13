# Nextflow Pipeline Integration with SiRNAForge Off-Target Analysis

This document describes how the Nextflow pipeline utilizes the core entrypoint functions defined in `src/sirnaforge/core/off_target.py`.

## Core Entrypoint Functions Used

The pipeline uses the following functions from the `__all__` export list in `off_target.py`:

### Primary Analysis Functions

1. **`run_comprehensive_offtarget_analysis`**
   - **Used in**: `modules/local/offtarget_analysis.nf`
   - **Purpose**: Main off-target analysis using BWA-MEM2
   - **Output**: TSV, JSON, and summary files
   - **Parameters**: species, sequences_file, index_path, output_prefix, bwa_k, bwa_T, max_hits, seed_start, seed_end

2. **`run_mirna_analysis_for_nextflow`**
   - **Used in**: `modules/local/mirna_offtarget_analysis.nf` (optional module)
   - **Purpose**: Specialized miRNA off-target analysis using BWA-MEM2 seed-mode parameters
   - **Output**: miRNA-specific TSV, JSON, and summary files

3. **`run_transcriptome_analysis_for_nextflow`**
   - **Used in**: `modules/local/transcriptome_offtarget_analysis.nf` (optional module)
   - **Purpose**: Specialized transcriptome off-target analysis using BWA-MEM2
   - **Output**: Transcriptome-specific TSV, JSON, and summary files

### Index Building Functions

4. **`build_bwa_index`**
   - **Used in**: `modules/local/build_bwa_index.nf`
   - **Purpose**: Build BWA-MEM2 indices from FASTA files (shared by transcriptome and miRNA seed modes)
   - **Parameters**: fasta_file, index_prefix

### Utility Functions

6. **`validate_and_write_sequences`**
   - **Used in**: `modules/local/prepare_candidates.nf`
   - **Purpose**: Validate siRNA candidate sequences
   - **Parameters**: input_file, output_file, expected_length

7. **`parse_fasta_file`** and **`write_fasta_file`**
   - **Used in**: `modules/local/aggregate_results.nf`
   - **Purpose**: FASTA file I/O operations

8. **`check_tool_availability`** and **`validate_index_files`**
   - **Available for**: Quality control and validation in modules

## Pipeline Architecture

### Main Workflow: `main.nf`
- Coordinates the entire off-target analysis pipeline
- Calls the `SIRNA_OFFTARGET_ANALYSIS` subworkflow
- Manages parameter validation and output organization

### Subworkflow: `subworkflows/local/sirna_offtarget_analysis.nf`
- Orchestrates the parallel analysis of siRNA candidates against multiple genomes
- Builds indices when needed using the core entrypoint functions
- Distributes candidates across analysis processes

### Core Analysis Modules

#### `modules/local/offtarget_analysis.nf`
- **Primary module** using `run_comprehensive_offtarget_analysis`
- Runs BWA-MEM2-based comprehensive off-target analysis
- Processes each candidate-genome combination independently
- Outputs: `.tsv` (detailed results), `.json` (structured data), `_summary.txt` (summary stats)

#### `modules/local/mirna_offtarget_analysis.nf` (Optional)
- Specialized module using `run_mirna_analysis_for_nextflow`
- Uses BWA-MEM2 with ultra-short read parameters for miRNA seed-match analysis
- Outputs: `_mirna_hits.tsv`, `_mirna_hits.json`, `_mirna_summary.txt`

#### `modules/local/transcriptome_offtarget_analysis.nf` (Optional)
- Specialized module using `run_transcriptome_analysis_for_nextflow`
- Uses BWA-MEM2 for transcriptome analysis
- Outputs: `_transcriptome_hits.tsv`, `_transcriptome_hits.json`, `_transcriptome_summary.txt`

### Supporting Modules

#### `modules/local/prepare_candidates.nf`
- Uses `validate_and_write_sequences` to validate siRNA sequences
- Ensures all candidates meet length and composition requirements
- Generates validation reports

#### `modules/local/build_bwa_index.nf`
- Uses `build_bwa_index`
- Builds genome indices when FASTA files are provided instead of pre-built indices
- Supports human, mouse, and custom genome references for both transcriptome and miRNA analysis

#### `modules/local/aggregate_results.nf`
- Aggregates results from all candidate-genome combinations
- Uses utility functions for file I/O operations
- Generates combined analysis files and final summary reports

## Function Parameter Mapping

### From Nextflow Parameters to Function Arguments

| Nextflow Parameter | Function Parameter | Description |
|-------------------|-------------------|-------------|
| `params.bwa_k` | `bwa_k` | BWA seed length |
| `params.bwa_T` | `bwa_T` | BWA minimum score threshold |
| `params.max_hits` | `max_hits` | Maximum hits per candidate |
| `params.seed_start` | `seed_start` | Seed region start position |
| `params.seed_end` | `seed_end` | Seed region end position |
| `params.genome_species` | `species` | Target species name |

### Output File Patterns

The entrypoint functions generate consistent output patterns:

- **Comprehensive analysis**: `{output_prefix}.tsv`, `{output_prefix}.json`, `{output_prefix}_summary.txt`
- **miRNA analysis**: `{output_prefix}_mirna_hits.tsv`, `{output_prefix}_mirna_hits.json`, `{output_prefix}_mirna_summary.txt`
- **Transcriptome analysis**: `{output_prefix}_transcriptome_hits.tsv`, `{output_prefix}_transcriptome_hits.json`, `{output_prefix}_transcriptome_summary.txt`
