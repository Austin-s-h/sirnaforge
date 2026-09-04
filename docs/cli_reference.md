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

:::{warning}
**EXPERIMENTAL.** `--design-mode zfn` runs the experimental ZFN arm, which has known unfixed defects
(half-site orientation handling, FokI seed-region weighting, off-target region classification,
inverted `worst_site_score`/`best_offtarget_score` exports) tracked in
[#82](https://github.com/Austin-s-h/sirnaforge/issues/82). Do not use ZFN output for any decision
without independent validation. See [ZFN Module Guide](zfn_module.md).
:::

ZFN activity/off-target evaluation now has a dedicated command: `sirnaforge zfn`.
Use the `workflow` command for transcript-centric siRNA/miRNA runs.

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

:::{warning}
**EXPERIMENTAL — results are not decision-grade.** The ZFN arm ships experimental in 0.6.0 with known
unfixed defects in half-site orientation handling, FokI seed-region weighting and off-target region
classification, tracked in [#82](https://github.com/Austin-s-h/sirnaforge/issues/82). Do not use ZFN
output for any decision without independent validation. The published CCR5 half-site pair does not
match its own on-target site under the default strand-pairing rule — pass
`--zfn-right-half-site CTTTTGCAGTTT` rather than the published `AAACTGCAAAAG` — which also
invalidates the recorded ZFN validation runs. Two further defects change nothing visible in the
output: the exported `worst_site_score` and `best_offtarget_score` fields are inverted
(`worst_site_score` is the minimum site score, `best_offtarget_score` the maximum, whereas the
highest-scoring off-target is the most dangerous one), and a site inside a large containing gene can
be classified `intergenic`, which undercounts the exonic/promoter tallies the pass/fail filters read.
See [ZFN Module Guide](zfn_module.md).
:::

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

Operational guidance from the backend tuning work — measured before the half-site convention issue was
found, so read it as a runtime observation only, not as a validated correctness result:

- prefer `pyahocorasick` for the first run on large references, **but only for
  `--zfn-max-mismatches` of 3 or less**
- use `fm_index` only for repeated persisted-index workflows; treat it as experimental on large references
- keep `exhaustive_python` as the baseline comparator and fallback implementation

:::{warning}
**The default `pyahocorasick` backend aborts above 3 mismatches on a 12 bp half-site.** Both
pattern-enumerating backends (`pyahocorasick`, `fm_index`) expand the query over the full 15-letter
IUPAC alphabet rather than the four bases a genome contains, and reject the search when the expansion
exceeds 1,000,000 patterns. A 12 bp half-site at `--zfn-max-mismatches 4` expands to 5,498,165
patterns and an 18 bp half-site at 3 mismatches to 1,717,605, so both raise:

```text
ValueError: ZFN L half-site is too complex for the pyahocorasick backend: 5498165 candidate
patterns exceed the safety limit of 1000000.
```

`--zfn-max-mismatches 4` is the budget the CCR5 benchmark needs, so pass
`--zfn-search-backend exhaustive_python` for those runs. Tracked in
[#82](https://github.com/Austin-s-h/sirnaforge/issues/82).
:::

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
