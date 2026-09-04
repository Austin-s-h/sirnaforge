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

#### ZFN Notes

ZFN activity/off-target evaluation now has a dedicated command: `sirnaforge zfn`.
Use the `workflow` command for transcript-centric siRNA/miRNA runs.

#### Input Sources & Transcriptome References

siRNAforge accepts complementary inputs when you need to bypass gene search or control the reference used for transcriptome off-target analysis:

* `--input-fasta` replaces the transcript retrieval step. Point it at a local FASTA file, HTTP(S) URL, or FTP location. The positional argument (`GENE_QUERY`) still names the outputs, while the workflow designs guides from the supplied sequences. **When you pass `--input-fasta` without `--transcriptome-fasta`, transcriptome off-target analysis is disabled** (design-only mode).
* `--transcriptome-fasta` selects the dataset used for transcriptome off-target analysis. It accepts local paths, remote URLs, or presets such as `ensembl_human_cdna` and `ensembl_mouse_cdna` (see `sirnaforge cache --info`). **Provide this flag to re-enable transcriptome off-target analysis when running from a custom FASTA.**
* `--offtarget-indices` overrides the genome indices used for Nextflow/BWA-MEM2 with explicit `species:/path/to/index_prefix` entries. When present, these drive the set of species processed by the off-target pipeline.

Passing both flags is common: the input FASTA feeds the design engine, while the transcriptome FASTA controls which reference is indexed for the Nextflow stage. Remote resources are cached under `~/.cache/sirnaforge/` and reused automatically.

Design-only mode is a deliberate cost guard, not an oversight: resolving the built-in defaults means downloading and indexing four multi-gigabyte Ensembl cDNA references (human, mouse, rat, macaque). Supplying your own sequences never triggers that implicitly. Library callers get the same policy — `run_sirna_workflow(input_fasta=...)` is design-only unless you pass `transcriptome_fasta=...` or opt in with `allow_transcriptome_with_input_fasta=True`.

`--skip-off-targets` disables **all** reference-based screening for the run: no transcriptome reference is resolved, downloaded or indexed, the Nextflow off-target stage does not run, **and repeat-element detection is skipped as well**. Repeat detection scans guides against the query species' cDNA reference, so it cannot run without the very download the flag exists to avoid; `logs/workflow_summary.json` reports it as `repeat_summary.status = "skipped"` with `reason = "user_disabled"`, and candidates keep `repeat_flagged = false`. Drop `--skip-off-targets` (optionally with `--transcriptome-fasta`) whenever you need repeat verdicts.

Rows inside `off_target/results/*/analysis.tsv` and the aggregated `combined_offtargets.tsv` include a `species` column so you can filter hits directly. Aggregated summaries collapse those values into `human` vs `other` buckets, exposing `hits_per_species`, `human_hits`, and `other_species_hits` in `combined_summary.json` plus the workflow console output. The workflow also records the resolved reference decision in `logs/workflow_summary.json` (`reference_summary.transcriptome`) so each run documents whether the transcriptome reference was disabled, defaulted, or explicitly provided.

---

## search

Search gene databases and retrieve transcript sequences.

### Help

```{program-output} uv run sirnaforge search --help
```

---

## design

Design siRNA/miRNA candidates from FASTA sequences.

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

## zfn

Evaluate a ZFN pair and run exhaustive genome-wide off-target search.

### Help

```{program-output} uv run sirnaforge zfn --help
```

#### Notes

* `--zfn-left-half-site` and `--zfn-right-half-site` are required.
* `--zfn-search-space` accepts either a local/remote FASTA or a configured reference key.
* `--zfn-search-backend` selects the half-site scan engine: `pyahocorasick` (default), `exhaustive_python` (baseline), or `fm_index` (experimental).
* `--zfn-search-space-index` accepts a persisted index-bundle directory for indexed backends. This is currently supported by `fm_index`.
* `--zfn-algorithm` supports `homology`, `conserved_g`, and `zfn_v2`.
* Outputs are written as `sirnaforge/zfn_candidate_summary.json` and `sirnaforge/zfn_offtarget_sites.csv`, with run metadata in `logs/workflow_summary.json`.

Operational guidance from the backend tuning work:

- prefer `pyahocorasick` for the first run on large references
- use `fm_index` only for repeated persisted-index workflows; treat it as experimental on large references
- keep `exhaustive_python` as the baseline comparator and fallback implementation

For reproducible `fm_index` runs, prebuild one search-space bundle once, then reuse it across runs:

```bash
uv run sirnaforge internal zfn-build-search-index \
	--genome-fasta /path/to/hg38.fa \
	--search-backend fm_index
```

The command prints a JSON summary including `bundle_dir`; pass that directory to `--zfn-search-space-index` on subsequent `sirnaforge zfn` runs.

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
