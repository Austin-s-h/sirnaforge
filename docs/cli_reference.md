# CLI Reference

> **Auto-generated**: All command output below is captured live during documentation build from the actual `sirnaforge` CLI.

This reference shows each command with its real `--help` output and working examples.

## Help & Version

### Main Help

```{program-output} uv run sirnaforge --help
```

### Version

```{program-output} uv run sirnaforge version
```

---

## workflow

Run complete siRNA design from gene query to scored candidates.

### Help

```{program-output} uv run sirnaforge workflow --help
```

:::{note}
The workflow command searches for gene transcripts, designs siRNA candidates, scores them using thermodynamic analysis, and outputs ranked results.
:::

#### Input Sources & Transcriptome References

siRNAforge accepts complementary inputs when you need to bypass gene search or control the reference used for transcriptome off-target analysis:

* `--input-fasta` replaces the transcript retrieval step. Point it at a local FASTA file, HTTP(S) URL, or FTP location. The positional argument (`GENE_QUERY`) still names the outputs, while the workflow designs guides from the supplied sequences. **When you pass `--input-fasta` without `--transcriptome-fasta`, transcriptome off-target analysis is disabled** (design-only mode).
* `--transcriptome-fasta` selects the dataset used for transcriptome off-target analysis. It accepts local paths, remote URLs, or presets such as `ensembl_human_cdna` and `ensembl_mouse_cdna` (see `sirnaforge cache --info`). **Provide this flag to re-enable transcriptome off-target analysis when running from a custom FASTA.**
* `--offtarget-indices` overrides the genome indices used for Nextflow/BWA-MEM2 with explicit `species:/path/to/index_prefix` entries. When present, these drive the set of species processed by the off-target pipeline.

Passing both flags is common: the input FASTA feeds the design engine, while the transcriptome FASTA controls which reference is indexed for the Nextflow stage. Remote resources are cached under `~/.cache/sirnaforge/` and reused automatically.

Rows inside `off_target/results/*/analysis.tsv` and the aggregated `combined_offtargets.tsv` include a `species` column so you can filter hits directly. Aggregated summaries collapse those values into `human` vs `other` buckets, exposing `hits_per_species`, `human_hits`, and `other_species_hits` in `combined_summary.json` plus the workflow console output. The workflow also records the resolved reference decision in `logs/workflow_summary.json` (`reference_summary.transcriptome`) so each run documents whether the transcriptome reference was disabled, defaulted, or explicitly provided.

---

## search

Search gene databases and retrieve transcript sequences.

### Help

```{program-output} uv run sirnaforge search --help
```

---

## design

Design siRNA candidates from FASTA sequences.

### Help

```{program-output} uv run sirnaforge design --help
```

### Example: Design from Sample Data

```{program-output} uv run sirnaforge design ../examples/sample_transcripts.fasta -o /tmp/sirna_example.csv --top-n 5
```

#### Output Preview

```{program-output} head -6 /tmp/sirna_example.csv
:shell:
```

---

## validate

Check FASTA file format and content.

### Help

```{program-output} uv run sirnaforge validate --help
```

### Example: Validate Sample Data

```{program-output} uv run sirnaforge validate ../examples/sample_transcripts.fasta
```

---

## config

Show default configuration parameters.

```{program-output} uv run sirnaforge config
```

---

## sequences

Manage siRNA sequences and chemical modification metadata.

### Help

```{program-output} uv run sirnaforge sequences --help
```

---

## cache

Manage miRNA database cache for off-target analysis.

### Help

```{program-output} uv run sirnaforge cache --help
```
