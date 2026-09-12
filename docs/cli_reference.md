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
- `--transcriptome-indices` (also spelled `--offtarget-indices`) names transcriptome BWA-MEM2 indices you have already
  built, as `species:/path/to/index_prefix` entries. The species you declare is the label that reference carries and is
  added to the screen; it must be one the registry recognises (the same set `--species` accepts), and it does **not**
  suppress the references `--species` resolves — pass `--transcriptome-fasta` if you do not want the Ensembl defaults
  fetched alongside it. The cDNA FASTA the index was built from must be readable as plain text beside the prefix (the
  prefix itself, or `<prefix>.fa`/`.fasta`/`.fna`; a gzipped or binary neighbour does not count), or that species is
  reported unscreened rather than screened against a reference its hits cannot be classified with.

Passing both flags is common: the input FASTA feeds the design engine, while the transcriptome FASTA controls which reference is indexed for the Nextflow stage. Remote resources are cached under `~/.cache/sirnaforge/` and reused automatically.

Design-only mode is a deliberate cost guard, not an oversight: resolving the built-in defaults means downloading and indexing four multi-gigabyte Ensembl cDNA references (human, mouse, rat, macaque). Supplying your own sequences never triggers that implicitly. Library callers get the same policy — `run_sirna_workflow(input_fasta=...)` is design-only unless you pass `transcriptome_fasta=...` or opt in with `allow_transcriptome_with_input_fasta=True`.

#### Run Policy: Run Modes, Filter Actions and Where a Default Comes From

Every threshold, the design mode and the run mode are resolved **once**, by
`sirnaforge.config.run_policy.resolve_run_policy`, before anything is downloaded and before the
output directory is created. The CLI, `run_sirna_workflow`, `run_offtarget_only_workflow` and a
direct `WorkflowConfig` all receive the same `ResolvedRunPolicy`, and `logs/.../manifest.json`
records it under `run_policy`.

Precedence, applied exactly once and in this order:

1. the built-in versioned profile (`legacy` today, and **experimental**: its numbers are expert
   priors, and the two that were measured were measured on one target);
2. the design-mode preset, for `--design-mode mirna`, applied **only where you said nothing**;
3. `--policy-config FILE` (JSON or TOML: `{"settings": {"gc_max": 65}}`);
4. an option you gave on the command line.

**An omitted option and an option typed with the default value are different inputs.** The decision
is made from Click's parameter source, not by comparing your value against a default, which is why
`--design-mode mirna --gc-max 60` now keeps 60. It used to be silently rewritten to 52, and
`--overhang dTdT` to `UU`, because a value equal to the siRNA default read as "unset". `--gc-max 65`
is a supported setting in either mode.

`--run-mode` declares how much screening evidence a run claims to have:

| Run mode      | Meaning                                                                                    | Default for                                   |
| ------------- | ------------------------------------------------------------------------------------------ | --------------------------------------------- |
| `design_only` | No screening; every post-screen gate is **not evaluated**, which is not the same as passed | `sirnaforge design`                           |
| `exploratory` | Screening ran; incomplete evidence is kept and labelled                                    | —                                             |
| `qualified`   | The transcriptome channel in the query species must complete                               | `sirnaforge workflow`, `sirnaforge offtarget` |

`--skip-off-targets` maps to `--run-mode design_only`, and the manifest records the rule that did it
rather than presenting the mode as something you chose. The two are one run: both suppress reference
resolution, so `reference_summary.transcriptome` reads `disabled` either way. Asking for `--run-mode
qualified` together with `--skip-off-targets` is rejected. `--run-mode exploratory` **is** accepted
with `--skip-off-targets` (a run that keeps incomplete evidence and screens nothing is a coherent
thing to ask for), and every post-screen gate is `off` for it, because what decides that is whether
the run screened — not the mode's name. Run mode is **independent of `--design-mode`**: choosing miRNA
design says nothing about evidence completeness.

`--filter-action filter_id=off|warn|fail` (repeatable) sets one gate's action. A gate can be off for
three different reasons — you turned it off, it has no declared threshold, or the run holds no
evidence of the kind it reads — and all three mean **not evaluated**, never passed. What cannot be
waived while staying `qualified` is required evidence completeness (`--run-mode exploratory` is how
you say that out loud).

The three words are three different mechanisms, and that is why their limits differ:

- `off` works by **clearing the gate's threshold**, which the existing gate code already reads as "no
  gate". That is available for the ten gates whose threshold has an absent state; the help text of
  `sirnaforge workflow` lists them. The seven design-stage gates (`gc_content_min`, `gc_content_max`,
  `max_poly_runs`, `max_repeat_transcript_fraction`, `max_paired_fraction`, `min_asymmetry_score`,
  `min_empirical_score`) read a plain float with no absent value, so switching them off is **refused**
  rather than faked with an inert number that would be reported as a threshold you chose. Widen the
  threshold instead.
- `warn` works by **leaving `passes_filters` alone**: the gate still measures, still compares, and
  publishes its outcome in `<filter_id>_verdict` and `<filter_id>_observed`, but a `fail` verdict does
  not remove the candidate. So a `fail` verdict beside `passes_filters=PASS` is consistent — read the
  verdict columns, not the single label. This is what the per-filter verdict columns bought: before
  them a gate had exactly one way to express a failure, and `warn` sat in the vocabulary unapplied.
  It is honoured by `max_paired_fraction`, `min_asymmetry_score`, `min_empirical_score` and **every**
  post-screen gate, each of which records through `SiRNACandidate.record_filter_verdict`, the one
  place that reads the action. `min_asymmetry_score`, `max_mirna_perfect_seed` and
  `fail_on_high_risk_mirna` already **ship** as `warn`.
- **Four gates accept `warn` and do not honour it — a known defect, not a design choice.**
  `gc_content_min`, `gc_content_max` and `max_poly_runs` decide during candidate enumeration and
  assign the rejection label directly, dropping the candidate into `DesignResult.rejected_candidates`
  (kept only for dirty controls and auditing) before anything reads the action; the verdict is
  recorded correctly, but the candidate is gone from the scored output. `max_repeat_transcript_fraction`
  is stamped by a static method that stamps `REPEAT_ELEMENT` without being given the action at all.
  For these four, `warn` resolves and appears in the manifest and the gate still rejects. Fixing them
  means routing their rejection through `record_filter_verdict`, as `min_isoform_coverage` was in #105.
- **A gate can only be turned off, never on.** `warn` and `fail` are refused when the gate would be
  off anyway — it has no threshold to compare (`min_isoform_coverage`, `max_transcriptome_seed_perfect`
  and `max_total_offtarget_hits` ship that way, so set a threshold first), the boolean it reads is
  `False`, the run holds no screening evidence, or no 0.7.1 code reads it at all
  (`max_mirna_1mm_seed`). Accepting it would report an enforced limit that does not exist. This is why
  `sirnaforge design`, which defaults to `--run-mode design_only`, refuses `warn` and `fail` on every
  post-screen gate: they evaluate nothing on that command.

One more asymmetry between the commands, because `--filter-action` is offered on all three:
`sirnaforge design` resolves the actions into its manifest but constructs its designer **without**
them, so a `warn` set there is inert and each design-stage gate applies its declared default. Use
`sirnaforge workflow` when you need `warn` honoured at the design stage. On `sirnaforge offtarget` the
design-stage gates never run on pre-designed guides, so an action set on one is accepted and applies
to nothing.

Two honesty notes the manifest carries per gate, because the code earns them and prose would not:

- Six of the nine off-target gates read a **human-stratified** counter, not the all-species column of
  the same name, so their `scope.species` is `["human"]` and `evidence_exported` is `false` — the
  counter they compare is not in the candidate CSV, so a client cannot re-apply them and get the
  pipeline's answer. Exporting those counters is separate filter-scope work. The six do **not** all
  stratify the same way, and each gate's `definition` says which it is: `max_transcriptome_hits_0mm`,
  `_1mm` and `_2mm` count hits that are human **or unlabelled** (a blank species label is read as the
  query species); `max_mirna_perfect_seed` and `fail_on_high_risk_mirna` count hits labelled human
  **only**, so an unlabelled miRNA hit reaches neither gate; and `max_total_offtarget_hits` **sums the
  two conventions** in one number.
- `max_mirna_1mm_seed` is read by no gate in 0.7.1, so it resolves to `off`, and asking for `warn` or
  `fail` on it is refused. It no longer declares a threshold either: the default was `10`, which read
  as a limit the run enforced, and it is now `None` so the setting states nothing it cannot apply.

`--skip-off-targets` disables **all** reference-based screening for the run: no transcriptome reference is resolved, downloaded or indexed, the Nextflow off-target stage does not run, **and repeat-element detection is skipped as well**. Repeat detection scans guides against the query species' cDNA reference, so it cannot run without the very download the flag exists to avoid; `logs/workflow_summary.json` reports it as `repeat_summary.status = "skipped"` with `reason = "user_disabled"`, and candidates keep `repeat_flagged = false`. Drop `--skip-off-targets` (optionally with `--transcriptome-fasta`) whenever you need repeat verdicts.

Rows inside `off_target/results/*/analysis.tsv` and the aggregated `combined_offtargets.tsv` include a `species` column so you can filter hits directly. The aggregator writes seven classification columns onto `combined_offtargets.tsv`, holding the literal `not_classified` until something classifies each row; `sirnaforge workflow` and `sirnaforge offtarget` fill those verdicts in place, and also onto each `transcriptome/*_analysis.tsv` the run read.

Once something has classified a row, `hit_class` is exactly one of `on_target`, `ortholog`, `repeat`, `off_target` or `undetermined`; before that it holds the literal `not_classified`, which is not a verdict and is what a bare `nextflow run` publishes. `undetermined` means the hit's species has no transcript index at all, so orthology and query-gene membership could not be checked — `species_index_missing` on the same row says so. Undetermined hits are still counted in `off_target_count` and still feed every gate: absent evidence must not loosen a screen, and `undetermined_hits` on the candidate row reports how much of the count is unqualified. A table produced with no index is therefore distinguishable from one where every alignment genuinely is a liability.

Two symbol columns, and they answer different questions. `matched_symbol` is the symbol that _established_ the class — the ortholog symbol on an `ortholog` row, the query symbol on an `on_target` row recognised by symbol — with `symbol_lookup_missing` flagging the rows where that check could not run. It is **not** a per-hit gene name: it is the literal `unknown` on every `off_target`, `repeat` and `undetermined` row, and on an `on_target` row matched by transcript ID, even where the index does resolve a symbol. `hit_symbol` **is** the per-row gene name: the symbol the transcript index resolves for `rname`, independent of class, with `hit_symbol_missing` flagging the rows it could not resolve. Group off-targets by gene on `hit_symbol`, never on `matched_symbol`. Both use the literal `unknown` rather than an empty cell, and `unknown` is common and real — roughly one in seven reference transcripts carries no symbol — so read it as "not annotated", never as "no gene". A stale table read back with any classification cell blank counts as unannotated and is re-annotated in full before it is republished, so no published row carries an empty classification cell (`test_a_row_annotated_in_part_is_repaired_rather_than_republished_with_blank_cells`).

`ortholog_evidence` says _how_ an `ortholog` verdict was reached, because the two ways are not the same claim. `gene_id` means the hit's gene ID was in the orthologue mapping resolved from Ensembl Compara for the query gene — the only validated tier. `symbol_heuristic` means only uppercased gene-symbol equality matched, which HGNC and MGI being separate nomenclature authorities makes both incomplete and unsafe: human `TP53`'s mouse orthologue is `Trp53`, so symbol equality misses it, while unrelated genes sharing a symbol across species satisfy it. Every non-ortholog row reads the literal `not_applicable`, never an empty cell. Treat `symbol_heuristic` as a lead, not as conservation evidence; a species whose Compara lookup could not be completed is named in `offtarget_summary.filtering_stats.orthology.unresolved_species` in `logs/workflow_summary.json`, and its hits are exactly the ones that fall back to the heuristic. Screening against a custom `--transcriptome-fasta` resolves orthology when you declare the species (`mouse:/path/custom_cdna.fa`) or when it can be read from the reference's own Ensembl cDNA headers; otherwise the run warns and every cross-species hit stays an unqualified `off_target`.

The liability population is `hit_class ∈ {off_target, undetermined}`, and it is the same population the per-candidate `off_target_count` is derived from; the published row count is reconciled against the hits that fed those counters, and any shortfall is reported as a run warning. So no counted hit can go unpublished silently. That reconciliation is a **row total**, not per-candidate attribution — a mis-attribution that nets to zero across the table would not raise a warning.

**`combined_offtargets.tsv` has one column shape; `transcriptome/*_analysis.tsv` has two.** The aggregator writes the classification columns itself, so the aggregated table carries all 19 columns whichever entry point produced it — a direct `nextflow run` and `sirnaforge workflow` differ in whether the verdicts are filled (`not_classified` versus a real class), not in which columns exist. The per-species files are a different matter: they gain the columns only on the runs whose aggregate came back header-only and the fallback read them, so that artifact has both shapes and `GenomeAlignmentSchema` (`strict=True`, 12 columns) will reject the 19-column version of it. Validate against `sirnaforge.models.schemas.AggregatedOffTargetSchema`, which accepts either width and admits `not_classified` in every classification column.

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
