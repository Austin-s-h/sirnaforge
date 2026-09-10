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

- `--input-fasta` replaces the transcript retrieval step. Point it at a local FASTA file, HTTP(S) URL, or FTP location. The positional argument (`GENE_QUERY`) still names the outputs, while the workflow designs guides from the supplied sequences. **When you pass `--input-fasta` without `--transcriptome-fasta`, transcriptome off-target analysis is disabled** (design-only mode).
- `--transcriptome-fasta` selects the dataset used for transcriptome off-target analysis. It accepts local paths, remote URLs, or presets such as `ensembl_human_cdna` and `ensembl_mouse_cdna` (see `sirnaforge cache --info`). **Provide this flag to re-enable transcriptome off-target analysis when running from a custom FASTA.**
- `--offtarget-indices` overrides the genome indices used for Nextflow/BWA-MEM2 with explicit `species:/path/to/index_prefix` entries. When present, these drive the set of species processed by the off-target pipeline.

Passing both flags is common: the input FASTA feeds the design engine, while the transcriptome FASTA controls which reference is indexed for the Nextflow stage. Remote resources are cached under `~/.cache/sirnaforge/` and reused automatically.

Design-only mode is a deliberate cost guard, not an oversight: resolving the built-in defaults means downloading and indexing four multi-gigabyte Ensembl cDNA references (human, mouse, rat, macaque). Supplying your own sequences never triggers that implicitly. Library callers get the same policy — `run_sirna_workflow(input_fasta=...)` is design-only unless you pass `transcriptome_fasta=...` or opt in with `allow_transcriptome_with_input_fasta=True`.

`--skip-off-targets` disables **all** reference-based screening for the run: no transcriptome reference is resolved, downloaded or indexed, the Nextflow off-target stage does not run, **and repeat-element detection is skipped as well**. Repeat detection scans guides against the query species' cDNA reference, so it cannot run without the very download the flag exists to avoid; `logs/workflow_summary.json` reports it as `repeat_summary.status = "skipped"` with `reason = "user_disabled"`, and candidates keep `repeat_flagged = false`. Drop `--skip-off-targets` (optionally with `--transcriptome-fasta`) whenever you need repeat verdicts.

Rows inside `off_target/results/*/analysis.tsv` and the aggregated `combined_offtargets.tsv` include a `species` column so you can filter hits directly. The aggregator writes seven classification columns onto `combined_offtargets.tsv`, holding the literal `not_classified` until something classifies each row; `sirnaforge workflow` and `sirnaforge offtarget` fill those verdicts in place, and also onto each `genome/*_analysis.tsv` the run read.

Once something has classified a row, `hit_class` is exactly one of `on_target`, `ortholog`, `repeat`, `off_target` or `undetermined`; before that it holds the literal `not_classified`, which is not a verdict and is what a bare `nextflow run` publishes. `undetermined` means the hit's species has no transcript index at all, so orthology and query-gene membership could not be checked — `species_index_missing` on the same row says so. Undetermined hits are still counted in `off_target_count` and still feed every gate: absent evidence must not loosen a screen, and `undetermined_hits` on the candidate row reports how much of the count is unqualified. A table produced with no index is therefore distinguishable from one where every alignment genuinely is a liability.

Two symbol columns, and they answer different questions. `matched_symbol` is the symbol that _established_ the class — the ortholog symbol on an `ortholog` row, the query symbol on an `on_target` row recognised by symbol — with `symbol_lookup_missing` flagging the rows where that check could not run. It is **not** a per-hit gene name: it is the literal `unknown` on every `off_target`, `repeat` and `undetermined` row, and on an `on_target` row matched by transcript ID, even where the index does resolve a symbol. `hit_symbol` **is** the per-row gene name: the symbol the transcript index resolves for `rname`, independent of class, with `hit_symbol_missing` flagging the rows it could not resolve. Group off-targets by gene on `hit_symbol`, never on `matched_symbol`. Both use the literal `unknown` rather than an empty cell, and `unknown` is common and real — roughly one in seven reference transcripts carries no symbol — so read it as "not annotated", never as "no gene". A stale table read back with any classification cell blank counts as unannotated and is re-annotated in full before it is republished, so no published row carries an empty classification cell (`test_a_row_annotated_in_part_is_repaired_rather_than_republished_with_blank_cells`).

`ortholog_evidence` says _how_ an `ortholog` verdict was reached, because the two ways are not the same claim. `gene_id` means the hit's gene ID was in the orthologue mapping resolved from Ensembl Compara for the query gene — the only validated tier. `symbol_heuristic` means only uppercased gene-symbol equality matched, which HGNC and MGI being separate nomenclature authorities makes both incomplete and unsafe: human `TP53`'s mouse orthologue is `Trp53`, so symbol equality misses it, while unrelated genes sharing a symbol across species satisfy it. Every non-ortholog row reads the literal `not_applicable`, never an empty cell. Treat `symbol_heuristic` as a lead, not as conservation evidence; a species whose Compara lookup could not be completed is named in `offtarget_summary.filtering_stats.orthology.unresolved_species` in `logs/workflow_summary.json`, and its hits are exactly the ones that fall back to the heuristic. Screening against a custom `--transcriptome-fasta` resolves orthology only when the species can be read from the reference's own Ensembl cDNA headers; otherwise the run warns and every cross-species hit stays an unqualified `off_target`.

The liability population is `hit_class ∈ {off_target, undetermined}`, and it is the same population the per-candidate `off_target_count` is derived from; the published row count is reconciled against the hits that fed those counters, and any shortfall is reported as a run warning. So no counted hit can go unpublished silently. That reconciliation is a **row total**, not per-candidate attribution — a mis-attribution that nets to zero across the table would not raise a warning.

**`combined_offtargets.tsv` has one column shape; `genome/*_analysis.tsv` has two.** The aggregator writes the classification columns itself, so the aggregated table carries all 19 columns whichever entry point produced it — a direct `nextflow run` and `sirnaforge workflow` differ in whether the verdicts are filled (`not_classified` versus a real class), not in which columns exist. The per-species files are a different matter: they gain the columns only on the runs whose aggregate came back header-only and the fallback read them, so that artifact has both shapes and `GenomeAlignmentSchema` (`strict=True`, 12 columns) will reject the 19-column version of it. Validate against `sirnaforge.models.schemas.AggregatedOffTargetSchema`, which accepts either width and admits `not_classified` in every classification column.

Candidate rows carry `screen_query_id`, the id the candidate was screened under. Guides are deduplicated before alignment, so a hit row's `qname` is the _representative's_ id — join row-level evidence on `screen_query_id = qname`, not on `id` (which attributes one guide's entire hit set to a single candidate) and not on `guide_sequence` (which breaks between a U-spelled guide and its T-spelled twin). It is `None` when the candidate was **not submitted to the aligner's input FASTA**. A candidate that was submitted but whose alignment never ran (`nextflow_unavailable`, `nextflow_failed`, or a step-5 exception) still carries the key, and it will join to nothing — `off_target_screened=False` is the signal for that case, not an absent key.

Aggregated summaries collapse those values into `human` vs `other` buckets, exposing `hits_per_species`, `human_hits`, and `other_species_hits` in `combined_summary.json` plus the workflow console output.

`combined_summary.json` also states which species were actually screened, and `hits_per_species` covers only those. `species_screened` is the positive evidence that a search ran; a species named in `unscreened_species` has **no** hit count rather than a hit count of zero, and `rejected_species_files` says per file why (an empty table, an unreadable one, a schema rejection, or an index prefix that resolved to nothing). `status` is `completed` only when every requested species was screened. A species dropped before Nextflow — requested with no resolvable reference, or one whose index build failed — is reported in `offtarget_summary.filtering_stats.species_screening_shortfalls` in `logs/workflow_summary.json` and included in `unscreened_species`, and the run reports `partial`. The workflow also records the resolved reference decision in `logs/workflow_summary.json` (`reference_summary.transcriptome`) so each run documents whether the transcriptome reference was disabled, defaulted, or explicitly provided.

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

- `--zfn-left-half-site` and `--zfn-right-half-site` are required.
- `--zfn-search-space` accepts either a local/remote FASTA or a configured reference key.
- `--zfn-search-backend` selects the half-site scan engine: `pyahocorasick` (default), `exhaustive_python` (baseline), or `fm_index` (experimental).
- `--zfn-search-space-index` accepts a persisted index-bundle directory for indexed backends. This is currently supported by `fm_index`.
- `--zfn-algorithm` supports `homology`, `conserved_g`, and `zfn_v2`.
- Outputs are written as `sirnaforge/candidate_summary.json` and `sirnaforge/offtarget_sites.csv`, with run metadata in `logs/workflow_summary.json`.

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
