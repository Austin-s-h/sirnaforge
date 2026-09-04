# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2026-09-03

Correctness release for the off-target arm, from the audit in
[#80](https://github.com/Austin-s-h/sirnaforge/issues/80). **Composite scoring now runs after
off-target screening**, using a seven-term weight set; off-target counts are decomposed into
on-target/ortholog/repeat/genuine-off-target classes; and every distinct candidate sequence is
screened rather than only the top-ranked subset. Scores from this release are not comparable to
0.5.x, and ranking order changes for essentially every candidate — this is the intended outcome
of the redefinition, not a regression. Off-target counts for the same target will drop
substantially, for the same reason.

### Fixed

- **On-target self-hit exclusion only matched a candidate's exact source transcript (or, in a
  first pass, other isoforms present in the input FASTA), so uncapping `max_hits` turned the
  off-target filter into a near-total kill switch: 3,259 of 3,288 candidates failed
  "perfect transcriptome match" and zero survivors sat on the primary isoform.** A guide derived
  from one isoform necessarily also perfectly matches sibling isoforms (shared exons) — including
  novel/NMD-variant isoforms the gene search never returned, so they aren't in the 13-transcript
  input FASTA either. A hit is now recognised as on-target when it maps to the queried gene's ID
  or symbol, not only when its transcript appears in the input FASTA, so sibling isoforms stop
  being counted against the guide. The transcript→gene lookup that makes this possible is built
  per screened species — see the ortholog-recognition entry under Changed, which replaced a
  first-pass version of this fix that only ever worked for human.
- **Exhaustive off-target search is now the default, and `--max-hits` is finally exposed on the CLI.** The embedded Nextflow off-target pipeline capped retained hits per candidate/species at `params.max_hits = 10000` with no way to override it from `sirnaforge workflow` (the option didn't exist; `nextflow_config_overrides` only forwarded `docker_image` and ZFN params). Because `_filter_and_rank` sorts ascending by `offtarget_score`, the retained rows were the _most_ dangerous hits, not the least — but for genes with many perfect-match paralogs, human/rhesus hits filled all 10,000 slots, silently censoring the true per-species hit count and giving ~80% of candidates a false "zero off-target hits" result. `params.max_hits` now defaults to effectively unlimited (100,000,000); pass `--max-hits <n>` to cap it again for speed on large gene families `sirnaforge.workflow.run_sirna_workflow()` gained a matching `max_hits` parameter.
- **`--design-mode mirna` no longer silently designs with siRNA GC bounds.** `_resolve_design_mode` only applied the miRNA-specific GC range when the CLI's GC flags still equaled `(30.0, 52.0)`, but the actual `--gc-max` default is `60.0` — so the guard was never true and every miRNA-mode run used the siRNA GC ceiling under a "mirna" label. The guard now compares against the real siRNA defaults `(30.0, 60.0)`.
- **Ensembl/RefSeq/ClinVar HTTP clients now honor `HTTP_PROXY`/`HTTPS_PROXY` environment
  variables.** All `aiohttp.ClientSession` instances in `sirnaforge.data.base`,
  `sirnaforge.data.transcript_annotation`, and `sirnaforge.data.variant_resolver` were created
  without `trust_env=True`, so requests bypassed configured corporate proxies entirely (unlike
  `curl`/`requests`, `aiohttp` does not read proxy env vars by default). This caused misleading
  "HTTP 503: Access denied" errors on networks that block direct outbound access to Ensembl/NCBI
  but allow it via an authenticated proxy.
- **Embedded Nextflow off-target pipeline failed to start under Nextflow's strict config parser
  (`ConfigParserV2`, Nextflow 24+).** The bundled `nextflow.config` defined a `check_max()` Groovy
  method to cap per-process resource requests, which the strict parser rejects (and, once
  converted to a closure, still can't resolve from nested `withLabel:` scopes). Replaced the
  helper entirely with Nextflow's native `resourceLimits` process directive. Also cleaned up
  remaining `nextflow lint` warnings (deprecated `Channel.xxx` factory usage, implicit `it`
  closure params, unused closure/workflow parameters) across `main.nf` and the local subworkflows.
- **`--design-mode mirna` had no effect on the ranking of any run that reached off-target
  screening.** `MiRNADesigner` folds the ago-start, position-1 pairing and 3' supplementary
  bonuses into `composite_score` only, never into the composite term set. Moving scoring after
  screening meant the composite was rebuilt from that term set and the bonuses were dropped, so
  every screened miRNA run was ranked by the plain siRNA composite. The bonus and its normalising
  maximum are now recorded on the candidate and reapplied through one shared helper. Rows that
  never passed through `MiRNADesigner` — injected dirty controls, which are copies of _rejected_
  candidates — fall back to the mode's own maximum bonus rather than to zero, so they are divided
  by the same `1 + max_mirna_bonus` as everything else instead of sitting ~25% high in the same
  CSV. `docs/scoring.md` now states the exact relation between the `score_*` contribution columns
  and `composite_score` in miRNA mode; they do not sum to it, and not only by the bonus.
- **A screen that partly or wholly failed was scored as a clean, complete screen — and could
  outrank one.** A candidate with no results entry took the no-hits branch and received
  `off_target_sub_score(0) = 1.0`, maximal specificity, with `off_target_screened = True`. Three
  separate shapes of this: a species whose alignment stage produced no files (a BWA-MEM2 index
  OOM), a candidate that never reached the aligner at all, and — the common case — a run where
  _aggregation itself_ failed, where the aggregate's own `missing_species` list is empty because
  nothing wrote one. The guard now keys on positive evidence, the per-species alignment files the
  aggregator actually read, so its absence cannot read as success; miRNA-only runs, which align
  against no transcriptome at all, are covered by the same rule. Candidates in that state keep
  their design-time score, report `off_target_screened = False`, and are counted in
  `offtarget_summary.filtering_stats.candidates_not_scored_after_screening`. Hits that _were_ found are still
  classified, still reported, and still applied as filters — real evidence can only fail a
  candidate, never wrongly pass one — so `off_target_screened = False` now means "this screen was
  incomplete, the counts are a lower bound", not "the counts are zero". Where only some candidates
  could be scored after screening, the rest are sorted below every scored one and held out of
  `top_candidates` (a third exclusion beyond passing and non-repeat-flagged, so `top_candidates`
  can be shorter than `min(top_n, n_passing)`) rather than competing on an optimistic, differently
  scaled number; the run logs an ERROR and a console line naming the count.
- **Conservation was scored against the CLI species list rather than the species actually handed
  to the aligner, and a failed alignment inflated it.** `--genome-indices`, `--genome-fastas` and
  `--transcriptome-indices` all reach the same alignment channel in `main.nf`, but only
  `--genome-species` fed the conservation denominator, so an ortholog hit in one of the extra
  species could make the numerator exceed the denominator; `conservation_sub_score` raised and the
  `except` around the whole scoring block abandoned scoring silently, leaving the candidate with a
  design-time score while its neighbours carried post-screen ones. The denominator is now the set
  of non-query species handed to the aligner and the numerator is intersected with it. A species
  whose alignment produced nothing **stays** in that denominator: removing it emptied the term for
  single-ortholog-species runs, `compute_composite` redistributed its 0.10 weight to the surviving
  terms, and a partial screen scored 57.9 where the complete screen scored 51.1. A failed
  alignment can now only lower conservation. A failure to score after screening is loud rather
  than silent: ERROR level, a counter, and a console line.
  **Re-run to pick this up** — rankings and `top_candidates` membership change, so
  `candidates_all.csv`/`candidates_pass.csv` from an earlier build must be regenerated. Runs that
  supply species through `--transcriptome-indices` also get a new Nextflow work-dir cache key (the
  key includes the resolved species list), so the first run after upgrading re-screens instead of
  resuming; nothing needs to be cleared by hand.

### Added

- **Four new pure modules** the hit-classification/scoring overhaul is built on:
  `sirnaforge.data.transcript_index` (per-species `TranscriptGeneIndex`/`SpeciesTranscriptIndex`),
  `sirnaforge.core.hit_classification` (`classify_hit`, `HitClass`, `HitClassCounts` — the
  four-way on-target/ortholog/repeat/off-target classifier), `sirnaforge.core.repeat_detection`
  (`RepeatDetector`, a numpy-vectorized k-mer scanner), and `sirnaforge.core.scoring`
  (`compute_composite` and the per-term sub-score helpers).
- Design-time repeat-element flagging: every distinct guide is scanned against the query
  species' cDNA reference, and a guide occurring in more than 0.1%
  (`DEFAULT_REPEAT_TRANSCRIPT_FRACTION`) of reference transcripts is flagged. New
  `SiRNACandidate.repeat_flagged` / `repeat_transcript_fraction` fields record the verdict;
  flagged, currently-passing candidates get the new `FilterStatus.REPEAT_ELEMENT` verdict and
  are excluded from ranking. Costs roughly 47 s against the ~1 GB human cDNA reference (measured:
  0.47 s for a 10.1 MB reference, 1,082 needles across three lengths), once per run against the
  query species' reference only.
- `FilterStatus.EXCESS_OFF_TARGETS`, raised when the redefined genuine off-target count exceeds
  `OffTargetFilterCriteria.max_off_target_count`.
- New output columns in `SiRNACandidate`, `SiRNACandidateSchema` and the results CSVs:
  `on_target_hits`, `ortholog_hits`, `repeat_hits`, `ortholog_species`, `repeat_flagged`,
  `repeat_transcript_fraction`, `isoform_coverage`, `conservation_score`, `score_asymmetry`,
  `score_gc_content`, `score_accessibility`, `score_empirical`, `score_off_target`,
  `score_isoform_coverage`, `score_conservation`, `scored_after_screening`, `weight_set_version`.
- `design_summary` gains `repeat_excluded_count` and `repeat_threshold_fraction`.
  `offtarget_summary` gains `hit_classes` (`on_target`/`ortholog`/`repeat`/`off_target` totals),
  `query_gene_transcripts_recognised`, `ortholog_symbol_lookup_misses` and
  `species_index_misses`.
- `manifest.json` gains a `scoring` block: `{weight_set_version, weights, active_terms}`.

### Changed

- **`composite_score` is now computed _after_ off-target screening, from a seven-term weight
  set (`asymmetry`, `gc_content`, `accessibility`, `empirical`, `off_target`, `isoform_coverage`,
  `conservation`); prior scores are NOT COMPARABLE, and ranking order changes for essentially
  every candidate. This is the intended outcome of the redefinition, not a regression.** The
  weight-set version that produced a score is recorded in `manifest.json` as
  `scoring.weight_set_version` (currently `"2.0.0"`; `1.x` denotes the pre-#80 five-term set),
  so the two regimes cannot be silently compared. Weights, term by term:

  | Term               | Old weight | New weight | Note                                                                                                 |
  | ------------------ | ---------- | ---------- | ---------------------------------------------------------------------------------------------------- |
  | `asymmetry`        | 0.15       | 0.12       |                                                                                                      |
  | `gc_content`       | 0.15       | 0.10       |                                                                                                      |
  | `accessibility`    | 0.20       | 0.13       |                                                                                                      |
  | `empirical`        | 0.20       | 0.15       |                                                                                                      |
  | `off_target`       | 0.30       | 0.25       | redefined: post-screen genuine-off-target specificity, not the design-time self-repetitiveness proxy |
  | `isoform_coverage` | —          | 0.15       | new, post-screen                                                                                     |
  | `conservation`     | —          | 0.10       | new, post-screen                                                                                     |

- **`off_target_count` is redefined to count GENUINE off-targets only.** On-target, ortholog and
  repeat-mediated hits are now classified into their own columns (`on_target_hits`,
  `ortholog_hits`, `repeat_hits`) instead of inflating `off_target_count`; off-target counts for
  the same target will drop substantially. The mismatch-stratified
  `transcriptome_hits_{0,1,2}mm` counters are likewise redefined to count genuine off-targets
  only, and the `max_transcriptome_hits_*` thresholds in `OffTargetFilterCriteria` are read
  against those.
- **Ortholog recognition previously worked only when the screened species was human.** The
  transcript→gene index was built behind an `if transcriptome_species == "human"` gate into one
  shared dict, so a mouse/rat/macaque hit could never be recognised as an ortholog of the query
  gene. `TranscriptGeneIndex` now builds a separate index per species for every species
  screened, keyed on gene symbol compared case-insensitively. Where a species' annotation
  carries no symbol for a hit transcript, the hit stays `off_target` and the shortfall is
  counted (`ortholog_symbol_lookup_misses`) rather than guessed at.
- **Every distinct candidate sequence is now screened, not just the top-ranked passing
  subset.** `_prepare_offtarget_input` deduplicates candidates by normalized guide sequence
  before writing the off-target FASTA, and result integration fans the screening result for
  each distinct sequence back out to every candidate sharing it, so per-record CSV output is
  unchanged by the dedup. `top_n` no longer gates which candidates are screened — it keeps its
  reporting meaning only (how many candidates land in `top_candidates`). Candidates are
  re-ranked by the post-screen composite score once screening completes, and `top_candidates`
  is rebuilt from the viable (passing, non-repeat-flagged) subset rather than staying frozen at
  its design-time order.
- The design-time self-repetitiveness proxy (repeated 7-mers _within_ the guide) is retained as
  an unweighted diagnostic, `component_scores["design_off_target_proxy"]`; it no longer feeds
  `composite_score` under any name. For the standalone `sirnaforge design` path (no screening),
  the same scorer runs with the three post-screen terms inactive and renormalised, so
  design-only and post-screen runs share one score definition.
- `ScoringWeights` validation moved from a positional `field_validator` on `empirical` (which
  summed only the fields Pydantic had already validated at that point, so it silently stopped
  checking new weights as the term set grew) to a `model_validator(mode="after")` that sums all
  seven fields by name; an out-of-tolerance sum now always raises at configuration time.
- `SiRNACandidate.FilterStatus` gained the six off-target verdicts that previously existed only
  as raw strings in `workflow.py` (`TRANSCRIPTOME_PERFECT_MATCH`, `TRANSCRIPTOME_1MM`,
  `TRANSCRIPTOME_2MM`, `MIRNA_PERFECT_SEED`, `HIGH_RISK_MIRNA`, `TOTAL_OFFTARGETS`), plus
  `REPEAT_ELEMENT` and `EXCESS_OFF_TARGETS`. The Pandera allow-list in `SiRNACandidateSchema`
  now derives from the enum instead of maintaining its own copy, so the two cannot drift again.

### Removed

- **Eight `FilterCriteria` thermodynamic windows that were declared but never enforced** — no
  CLI flag exposed any of them: `mfe_min`, `mfe_max`, `duplex_stability_min`,
  `duplex_stability_max`, `melting_temp_min`, `melting_temp_max`, `delta_dg_end_min`,
  `delta_dg_end_max`. Re-deriving correct windows is deliberately out of scope for this change.
- `FilterCriteria.max_off_target_count` — design time cannot know the genuine off-target count,
  so it never gated anything there. Moved to `OffTargetFilterCriteria.max_off_target_count`,
  where it is now enforced post-screen (see `FilterStatus.EXCESS_OFF_TARGETS` above).
- `SiRNACandidate.on_target_transcriptome_hits` (and its CSV column), superseded by
  `on_target_hits`: two columns for the same quantity would drift, and the older one was only
  one release old. `on_target_confirmed` is unchanged — it still answers "was any hit recognised
  as the query gene?", a different question from the new hit count.

## [0.5.2] - 2026-09-02

Correctness release for the candidate selection/scoring arm, from the audit in
[#78](https://github.com/Austin-s-h/sirnaforge/issues/78). **Behavior note: every
thermodynamic column and the pass/fail composition change. Results are not comparable with
0.5.1 or earlier, and re-running an earlier design will select different candidates.** The
off-target arm (transcriptome + miRNA seed) is unaffected.

### Fixed

- **Duplex thermodynamics folded each strand against itself.** `passenger_sequence` is
  already the reverse complement of `guide_sequence`, but both ViennaRNA call sites
  reverse-complemented it again, so `calculate_duplex_stability` folded `guide & guide` and
  `_calculate_end_stability` folded `guide[:7] & guide[-7:]`. ViennaRNA's cofold MFE is
  symmetric in its two strands, so at the default `sirna_length=21` the 5' and 3' end folds
  were the _same_ pair of 7-mers: `dg_5p == dg_3p` exactly, which pinned `asymmetry_score`
  to the constant 0.5 for every candidate and `delta_dg_end` to 0.0 (making the documented
  +2 to +6 kcal/mol optimum unreachable). 23 nt designs desynchronised the two windows and
  so produced a spread of wrong values instead of one wrong value. Measured on the 505
  candidates of the packaged GAPDH transcript: duplex ΔG moves from -24.2..+0.2 (median
  -8.8) to -43.3..-32.5 (median -39.2), and the 199 rows with a non-physical positive
  `dg_5p`/`dg_3p` drop to zero. Both end windows now use the same width, which also fixes
  19-20 nt guides comparing a 7-mer 5' window against a 5-mer 3' window.
- **`LOW_ASYMMETRY` was assigned by testing the empirical score.**
  `_calculate_empirical_score` compared its own value against
  `filters.min_asymmetry_score`; `asymmetry_score` was computed, stored, exported, and never
  gated. Since the simplified Reynolds rule caps at 0.70, the default threshold of 0.65
  reduced eligibility to `guide[0] in {G,C} and guide[18] == "A"` — a two-nucleotide test
  that rejected 84.6% of candidates as "low asymmetry", including ones the same run ranked
  at the top. Each threshold now gates the quantity it is named after.
- **The position-19 A/U test never matched T.** Guides are stored as DNA, so 141 of 505
  GAPDH guides were denied a bonus they qualified for. The same defect affected three
  miRNA-mode bonuses: `ago_start_bonus` skipped all 144 T-initiated guides; every A:U pair
  at guide position 1 (257 of 505) was classified `mismatch` and collected an undeserved
  `pos1_mismatch_bonus`, while the G:U wobble branch could never fire at all; and
  `supp_13_16_score` counted no T, capping at 0.75 with 176 candidates stuck at 0.0.
- **Melting temperature was a scaled ΔG.** `37 + 2 * -ΔG` cannot yield a physical Tm (that
  needs ΔH and ΔS separately) and only looked plausible because it was reading the self-fold
  ΔG; on corrected input it returns 102.0-123.6 °C (median 115.4). Replaced with Biopython's `Tm_NN` on the RNA
  nearest-neighbour table of Xia et al. (1998) at 100 mM Na+ / 100 nM strands: 66.3-81.7 °C,
  median 75.3 °C. These are physical values for a 21 bp RNA duplex; a quarter of them sit
  above the advisory 60-78 °C reference range, which is not enforced.
- **The miRNA composite score was clamped at its 100-point ceiling.** Bonuses of up to 0.25
  were added to a 0-100 score and the sum truncated, flattening the top of the ranking —
  nine GAPDH candidates tied at exactly 100.0. Now rescaled by the maximum attainable total,
  which is monotone in `base + bonus` and so preserves the order the clamp was destroying.
- **`duplex_stability_score` saturated.** Its [-40, -5] kcal/mol window was calibrated
  against the self-fold ΔG; against real duplex ΔG it put 209 of 505 candidates at exactly
  1.0. Now normalised per nucleotide (-2.1 kcal/mol/nt -> 1.0, -1.4 -> 0.0), which also
  stops longer designs scoring higher for their length alone. Reported only; it does not
  feed `composite_score`.
- **The run manifest was not a record of the run.** `manifest.json` listed four design
  parameters, so `min_asymmetry_score`, the homopolymer gate, the scoring weights and the
  design mode appeared nowhere in it. Both `manifest.json` and
  `logs/workflow_summary.json` now dump `DesignParameters` wholesale, so a newly added
  parameter cannot go unrecorded; the manifest also carries `tool_version`.
- Corrected the `delta_dg_end` sign convention in `FilterCriteria`, `SiRNACandidateSchema`
  and the tutorials: it is `dg_5p - dg_3p` (positive = destabilised 5' end), not the
  reverse.

### Added

- `FilterCriteria.min_empirical_score` (default 0.5), bounded by the empirical rule's
  attainable range 0.4-0.7, so a threshold the rule can never satisfy raises at
  construction instead of silently rejecting every candidate.
- `FilterStatus.LOW_EMPIRICAL_SCORE`, so an empirical-rule rejection is no longer reported
  as an asymmetry failure. Registered with the candidate CSV schema.
- `SiRNACandidate.off_target_screened` (new CSV column). The hit-count columns default to 0
  and only the top-N selection is screened, so previously a candidate that was never
  screened was indistinguishable from one that came back clean.
- `ThermodynamicCalculator.meets_asymmetry_threshold`, which gates an already-computed score
  without re-folding both duplex ends. `is_thermodynamically_favorable` — previously dead
  code — now delegates to it, and its default threshold matches
  `FilterCriteria.min_asymmetry_score` (it was 0.5 against a documented 0.65).

### Changed

- Default filter composition on the GAPDH set: PASS 12.3% -> 26.1%, `LOW_ASYMMETRY` 84.6%
  (mislabelled) -> 66.1% (real thermodynamic asymmetry), plus `LOW_EMPIRICAL_SCORE` 4.6%.
  `min_empirical_score` defaults to 0.5 rather than 0.65 deliberately: at 0.65 the empirical
  gate alone passes 4.8% of candidates, which is the degenerate behaviour being fixed. The
  rule still influences ranking through its 0.20 composite weight.
- Removed three `except (AttributeError, ValueError, TypeError): pass` blocks around filter
  assignment rather than making them log. `FilterStatus` is always present and assigning a
  str-enum to a Pydantic field without `validate_assignment` cannot raise, so the handlers
  only ever hid the real defect — that the wrong field was being compared.
- Documented that the composite score's 0.30 off-target weight is a design-time proxy
  (repeated 7-mers _within_ the guide). Transcriptome/miRNA screening runs after design and
  gates `passes_filters`; it does not contribute to `composite_score`.
- Documented that `mfe_*`, `duplex_stability_*`, `melting_temp_*`, `delta_dg_end_*` and
  `max_off_target_count` in `FilterCriteria` have no reader in `src/` — they are reference
  ranges, not filters. Wiring them up needs calibration against truth data, and the
  `duplex_stability` window in particular is unreachable now that duplex ΔG is correct.
- Rewrote the `awk` filter examples in `docs/thermodynamic_guide.md` to resolve columns by
  header name, so they survive column-set changes such as `off_target_screened`.

## [0.5.1] - 2026-07-21

This release adds the Zinc Finger Nuclease (ZFN) design/off-target module, multi-species genome
references, and a correctness fix to miRNA seed-hit scoping. **Behavior note:** the miRNA seed-hit
fix changes the reported hit census (previously-counted non-seed perfect matches are no longer
counted as hits), so results are not identical to 0.4.x.

### Added

- **Zinc Finger Nuclease (ZFN) module** (`sirnaforge.zfn`): typed ZFN pair models, exhaustive and
  indexed half-site search backends, candidate ranking, annotation, and an optional Nextflow
  execution bridge. New `sirnaforge zfn` CLI command and `zfn_candidate_summary.json` /
  `zfn_offtarget_sites.csv` outputs.
- **Multi-species genome references.** Genome (DNA) references now cover human, mouse, rat, and
  macaque, matching the transcriptome/miRNA species set. Previously only
  `ensembl_human_hg38_primary` was available, so ZFN off-target search could only use a human
  search space. New keys: `ensembl_mouse_grcm39_primary`, `ensembl_rat_grcr8_toplevel`,
  `ensembl_macaque_mmul10_toplevel`. The human key is unchanged for backward compatibility.
- New `sirnaforge.data.ensembl_references` module: a single Ensembl assembly table generates both
  the transcriptome cDNA and genome DNA `SOURCES`, so adding a species is a one-line change.

### Fixed

- **miRNA seed hits are now scoped to the seed region.** `run_mirna_seed_analysis` previously
  copied all raw alignments into the filtered output (`df_filtered = df_raw.copy()`), so perfect
  guide-seed matches anywhere in a miRNA — including the 3′ region — were counted as hits,
  inflating `total_hits`/`mirna_hits_0mm_seed`. A raw alignment is now counted as a hit only when
  the guide seed lands on the miRNA's own seed region (0-based `coord == seed_start - 1`);
  non-seed matches are retained in the `*_raw` outputs but excluded from filtered hits.
- **Seed-mismatch coordinate frame.** Mismatch positions from the seed scanners were window-relative
  (`1..len(seed)`) but filtered with the guide-relative `seed_start..seed_end` threshold, so a
  mismatch at guide position 2 was dropped and the hit mislabeled a perfect seed. Window positions
  are now mapped onto guide coordinates in all three scan paths (exhaustive, pyahocorasick, BWA).

### Changed

- Raised the max analyzable guide/passenger length from 23 nt to 40 nt (`SiRNACandidate`), so longer
  observed species (e.g. Dicer 3′ read-through isoforms) can be scored. A validator warns above the
  23 nt recommended biological range; `DesignParameters.sirna_length` (what the tool _designs_) stays
  capped at 23.
- Off-target hit reporting is now exhaustive by default: `_mirna_max_hits()` returns `None` (no cap),
  and the default `max_hits` for `BwaAnalyzer`, `run_bwa_alignment_analysis`, and
  `run_comprehensive_offtarget_analysis` is `None` instead of `10000`. `SIRNAFORGE_MIRNA_MAX_HITS`
  can still impose a cap. A silent cap biased downstream hit counts.
- `--zfn-search-space` CLI help now lists the built-in genome keys; `docs/zfn_module.md` documents
  the reference table. Rat/macaque use Ensembl's `dna.toplevel` file (verified: those assemblies do
  not publish `dna.primary_assembly`). Transcriptome cDNA keys/URLs (and their cache keys) are
  unchanged, so existing caches remain valid.

### Housekeeping

- Cleared embedded execution outputs from `notebooks/zfn_backend_runtime_comparison.ipynb`
  (~1 MB → 68 KB); benchmark methodology preserved.

## [0.4.3] - 2026-03-22

### New Features

- Changed default ZFN search backend from EXHAUSTIVE_PYTHON to PYAHOCORASICK for improved performance.
- Added new fields to ZFNShardingConfig: memory_budget_gb, memory_reserve_gb, target_cpu_utilization, and max_cpu_workers to allow for better resource management.
- Enhanced memory and CPU utilization handling in ExhaustiveZFNOffTargetSearcher to dynamically adjust worker counts based on available resources and user-defined limits.
- Introduced warnings for using fm_index on large references to prevent high memory consumption.
- Added integration tests for ZFN annotation and reference resolution, ensuring compatibility with real genomic data.
- Updated unit tests to validate new sharding parameters and their impact on worker allocation and chunk sizing.

- **Pluggable ZFN search backends**:
  - Added `ZFNSearchBackend` support across workflow, CLI, and typed ZFN contracts.
  - Added optional backend engines for `pyahocorasick` and `fm_index` alongside `exhaustive_python`.

- **Persisted ZFN search-space index bundles**:
  - Added `fm_index` bundle build/load support with manifest + artifact validation.
  - Added internal CLI command `sirnaforge internal zfn-build-search-index` for reproducible bundle generation.

### Improvements

- **Cache-aligned reproducibility for `fm_index`**:
  - Index bundle creation now resolves FASTA inputs through the shared genome cache pipeline.
  - Default bundle output paths are stable and cache-derived for repeatable runs.
- **miRNA backend rollout readiness**:
  - Clarified the internal miRNA seed-scanning seam around `pyahocorasick` as the operational default.
  - Kept `exhaustive_python` as the correctness oracle and BWA as the semantic comparison baseline for parity validation.
  - Preserved the existing miRNA output contract and embedded Nextflow batch behavior without adding a public backend-selection flag.
- **Strict typing and enum-driven CLI contracts**:
  - Replaced string validation paths for backend/algorithm/dimer mode with imported enum types.
  - Updated Nextflow bridge and shard execution paths to accept typed enums directly.
- **Documentation and tests**:
  - Expanded CLI reference with backend/index usage guidance and reproducible prebuild examples.
  - Added/updated unit coverage for persisted index flow and internal index command delegation.

### Dependencies

- Added runtime dependencies for optional accelerated backends:
  - `pyahocorasick>=2.1.0`
  - `fm-index>=2.3.4`

## [0.4.2] - 2026-03-01

### New Features

- **ZFN workflow maturity for release**:
  - Added robust ZFN workflow mode coverage for constrained half-site design and off-target site discovery.
  - Documented and stabilized ZFN algorithm selection (`homology`, `conserved_g`, `zfn_v2`) and dimer controls.
  - Added explicit reporting expectations for `zfn_candidate_summary.json` and `zfn_offtarget_sites.csv` outputs.

### Improvements

- **Documentation completeness and command correctness**:
  - Reconciled stale CLI examples with current command options (`--species`, `--modifications`, ZFN-required half-site flags).
  - Added curated ZFN guidance to complement auto-generated CLI help output.
  - Updated installation docs to reflect Python support (`3.10-3.12`).
- **Release metadata consistency**:
  - Aligned release-facing documentation for the 0.4.2 train across changelog and release notes.

## [0.4.0] - 2025-1-5

### Breaking

- Drop support for Python 3.9; minimum supported Python is now 3.10. This
  enables PEP 604 union syntax (e.g. `list[str] | None`) and other modern
  typing features. Update developer environments, CI, and pre-commit hooks to
  use Python 3.10 or later.

### New Features

- **Variant Targeting Implementation**: Complete Phase 1-5 implementation for targeting specific genetic variants
  - Core variant models and resolver infrastructure with Parquet-based caching
  - Population-specific AF filtering for geographic variant targeting
  - Phase 5 workflow integration with CLI flags for variant_mode parameter
  - Comprehensive variant feature implementation summary and documentation
- **Enhanced CLI and Workflow Integration**:
  - Off-target-only entry point for pre-designed siRNA candidates
  - Support for any sequence length in thermodynamic calculations
  - Improved CLI enum integration for variant_mode parameter
- **Docker and Container Improvements**:
  - Python 3.12 upgrade with optimized Dockerfile targeting
  - Login shell PATH preservation with `/etc/profile.d/conda-path.sh`
  - Enhanced container testing with dedicated test categories
  - Improved Docker entrypoint and health checks

### Improvements

- **Performance Optimizations**:
  - Parquet-based variant cache for improved performance
  - Cache-first index reuse with complete validation
  - Transcriptome filtering and major refactoring for simplicity
  - Unified cache management with Pythonic interface
- **CI/CD Pipeline Enhancements**:
  - Python 3.12 requirement across all workflows
  - Aligned pre-commit mypy with uv package manager
  - Enhanced release workflow with comprehensive testing
  - Improved Docker build and test categorization
- **Code Quality and Maintenance**:
  - Python 3.10+ syntax modernization throughout codebase
  - Comprehensive typing improvements and linting fixes
  - Enhanced error handling and validation middleware
  - Improved documentation with live CLI output examples

### Dependencies

- **Python 3.12 Support**: Full upgrade to Python 3.12 with modern syntax
- **Enhanced Dependencies**: Added `pyarrow>=18.0.0` for Parquet support
- **Updated Packages**: Modernized dependency versions with improved compatibility
- **uv Package Manager**: Full alignment with uv for faster dependency resolution

### Performance

- **Variant Caching**: Parquet-based storage for improved variant data performance
- **Memory Optimization**: Reduced memory requirements for Docker-constrained environments
- **Parallel Processing**: Enhanced concurrent execution for variant analysis
- **Index Reuse**: Cache-first approach for transcriptome indices

## [0.3.4] - 2025-12-31

### Added

- **Transcript Annotation Provider Layer**: New data provider interface for fetching genomic transcript annotations
  - Added `AbstractTranscriptAnnotationClient` interface in `src/sirnaforge/data/base.py`
  - Implemented `EnsemblTranscriptModelClient` using Ensembl REST API (lookup/id and overlap/region endpoints)
  - Added `VepConsequenceClient` stub for optional VEP enrichment (behind config flag, placeholder implementation)
  - New Pydantic models: `Interval`, `TranscriptAnnotation`, and `TranscriptAnnotationBundle` in `src/sirnaforge/models/transcript_annotation.py`
  - In-memory LRU cache with TTL for transcript annotations
  - Support for fetching by stable IDs or genomic regions
  - Comprehensive unit tests with mocked REST responses
  - Integration tests for real Ensembl REST API (gated by `@pytest.mark.requires_network`)

### Improvements

- **Extensible Architecture**: Transcript annotation provider follows the same layered pattern as existing data providers (gene search, ORF analysis, transcriptome management)
- **Reference Tracking**: Annotations include provenance metadata (provider, endpoint, reference choice) for reproducibility
- **Error Handling**: Robust handling of unresolved IDs and network errors with fallback to unresolved list

### Documentation

- Added comprehensive docstrings for all new classes and methods
- Unit and integration tests serve as usage examples

## [0.3.3] - 2025-12-15

### Bug Fixes

- **Docker Login Shell PATH**: Fixed issue #37 where login shells (`/bin/bash -lc`) would reset PATH and drop `/opt/conda/bin`, making `sirnaforge` and `nextflow` unavailable
  - Added `/etc/profile.d/conda-path.sh` to preserve conda toolchain paths in login shells
  - Non-login shells continue to work as before
  - Added regression test `test_docker_login_shell_path()` to container test suite
  - Added standalone test script `scripts/test-docker-login-shell.sh` for manual verification
- **Nextflow Off-target Aggregation**: Fixed a Groovy/DSL2 runtime crash during final aggregation (`No signature of method ... call(LinkedList)`) by correcting channel collection/defaulting semantics in the embedded workflow
  - Replaced invalid `ifEmpty([])`/`ifEmpty('')` usage with `ifEmpty { [] }`/`ifEmpty { '' }`
  - Switched from `collect()` to `toList()` for explicit channel materialization before combining genome + miRNA result lists

## [0.3.1] - 2025-12-04

### Added

- **Dirty Control Injection**: `workflow.py` now carries the worst rejected guides forward as "dirty control" candidates (see `sirnaforge/utils/control_candidates.py`) so every Nextflow/off-target run includes known-failing sentinels for health checks.

### Improvements

- **Resilient Aggregation & Reporting**: Nextflow modules (`modules/local/aggregate_results.nf`, `mirna_offtarget_analysis.nf`, `split_candidates.nf`, etc.) plus `pipeline/nextflow_cli.py` were refactored to emit consolidated TSV/JSON artefacts even when some analyses are skipped, ensuring miRNA/genome summaries always arrive in `workflow_output/`.
- **Deterministic Caching**: New cache utilities expose `SIRNAFORGE_CACHE_DIR`/XDG-aware paths and tag Nextflow workdirs with metadata, dramatically reducing repeated downloads and making cleanup predictable.
- **Workflow Parameter Safety**: CLI defaults now enforce valid GC range/length boundaries and automatically fall back to the bundled Ensembl transcriptome set (`ensembl_human_cdna`, `ensembl_mouse_cdna`, `ensembl_rat_cdna`, `ensembl_macaque_cdna`) when no input is provided, preventing empty design runs.

### Bug Fixes

- **Pipeline Robustness**: Aggregation handles missing combined TSVs, gracefully copies per-species miRNA batches, and logs explicit workdir pointers so failed runs can be recovered without manual spelunking.
- **Nextflow Reliability**: All embedded DSL2 modules gained scoped retries, consistent BWA-MEM2 index prep, and container profile detection, eliminating the intermittent crashes seen in long off-target analyses.

### Documentation

- **Docs v2 Stack**: Added a parallel `docs_v2/` tree with autogenerated CLI/API references, live `sirnaforge` command output, and refreshed installation guides focused on Docker + Nextflow workflows.
- **Workflow & Tutorial Refresh**: `docs/getting_started.md`, `docs/usage_examples.md`, and the new Nextflow tutorial now describe dirty controls, cache locations, and off-target artefacts so users can reproduce the updated pipeline end-to-end.

## [0.3.0] - 2025-11-21

### Improvements

- **Documentation Standardization**: Unified tab-based execution examples across all documentation
  - Added sphinx-design tab-sets for uv/Docker execution in all usage examples
  - Standardized command patterns in `usage_examples.md`, `gene_search.md`, `getting_started.md`
  - Simplified usage examples from redundant variations to minimal/comprehensive pattern
  - Improved user experience with consistent execution context switching
- **GC Content Default Update**: Increased default `--gc-max` from 52% to 60%
  - Updated across CLI, documentation, and tutorials
  - Better alignment with siRNA design best practices
  - Maintains conservative gc-min default of 35%
- **CI/CD Pipeline Enhancements**:
  - Fixed release.yml to use correct make targets (`docker-test`, `test-dev`)
  - Added comprehensive `test-release` job with full coverage reporting
  - Improved test tier structure: `test-dev` (unit), `test-ci` (smoke), `test-release` (comprehensive)
  - Coverage artifacts now uploaded with 30-day retention
  - Added coverage summary to GitHub Actions workflow UI

### Bug Fixes

- **Off Target and miRNA seed match search now works!**
- **Docker Test Environment**: Fixed environment conflicts in `make docker-test`
  - Removed uv sync from Docker container execution (conflicts with conda)
  - Explicit pytest installation and conda path execution
  - Resolved parallel execution issues with pytest-xdist
- **Missing Dependencies**: Added `psutil>=6.0.0` to production dependencies
  - Required by workflow.py but was previously undeclared
  - Ensures consistent environment across installations
- **Nextflow Integration Tests**: Fixed two previously skipped tests
  - Corrected Nextflow version flag: `--version` → `-version`
  - Fixed workflow access to use NextflowRunner API instead of module import
  - All 20 container tests now passing (was 18 passed, 2 skipped)

### Testing

- **Comprehensive Test Coverage**: Enhanced `make test-release` to run all test tiers
  - Now runs 179 tests (dev + ci + release markers) with 55% coverage
  - Generates XML, HTML, and terminal coverage reports
  - Full JUnit XML for CI/CD integration
  - Execution time: ~22 seconds for complete validation
- **Test Organization**: All tests properly tagged with tier markers
  - 9 expected skips (requires Docker/Nextflow/BWA-MEM2 tools)
  - Consistent marker structure across test suite
  - Better CI/CD integration with proper test selection

### Build & Infrastructure

- **Makefile Improvements**: Enhanced test targets with better coverage support
  - `test-release` now comprehensive (dev + ci + release tests)
  - All test targets include appropriate coverage/junit reporting
  - Clear documentation of tier structure in help text
- **GitHub Actions**: Aligned workflows with current Makefile structure
  - CI workflow: lint → security → test-ci → build
  - Release workflow: validate → ci → test-release → build artifacts → docker
  - Proper dependency chains ensure quality gates

## [0.2.2] - 2025-10-26

### New Features

- **miRNA Design Mode**: New `--design-mode mirna` option for microRNA-specific siRNA design
  - Specialized `MiRNADesigner` subclass with miRNA-biogenesis-aware scoring
  - Enhanced CSV schema with miRNA-specific columns (strand_role, biogenesis_score)
  - CLI support via `--design-mode` flag with automatic parameter adjustment
- **miRNA Seed Match Analysis**: Integrated miRNA off-target screening in Nextflow pipeline
  - Lightweight seed region matching (positions 2-8) against miRNA databases
  - Automatic miRNA database download and caching from MirGeneDB
  - Per-candidate and aggregated miRNA hit reports in TSV/JSON formats
  - Configurable via `--mirna-db` and `--mirna-species` flags
- **Species Registry System**: Canonical species name mapping and normalization
  - Unified species identifiers across genome and miRNA databases
  - Automatic species alias resolution (e.g., "human" → "Homo sapiens" → mirgenedb slug)
  - Support for multi-species analysis with consistent naming

### Improvements

- **Nextflow Pipeline Enhancements**:
  - Reduced memory requirements for Docker-constrained environments (2GB → 1GB for most processes)
  - Added miRNA seed analysis module with BWA-based matching
  - Improved error handling and progress reporting
  - Better resource allocation with memory/CPU buffers
- **Data Validation**: Extended Pandera schemas for miRNA-specific columns
- **CSV Output**: New columns `transcript_hit_count` and `transcript_hit_fraction` track guide specificity
- **miRNA Database Manager**: Enhanced with species normalization and canonical name mapping

### Bug Fixes

- Fixed Nextflow Docker configuration for resource-constrained CI environments
- Resolved schema validation errors for miRNA columns in mixed-mode workflows
- Fixed typing issues in pipeline CLI functions

### Documentation

- **Major Documentation Consolidation**: Reorganized structure for improved user experience
  - Simplified navigation from 4 to 3 main sections (Getting Started, User Guide, Reference, Developer)
  - Consolidated `getting_started.md` and `quick_reference.md` into comprehensive guide
  - Streamlined tutorials to 2 focused guides (pipeline integration, custom scoring)
  - Created dedicated developer section for advanced documentation
- **Complete API Reference**: Added 18 previously missing modules
  - Comprehensive coverage of all 27 sirnaforge modules
  - Auto-generated Sphinx documentation with proper cross-references
- **Quality Improvements**: Configured ruff D rules for docstring validation
  - Fixed 116 docstring formatting issues automatically
  - Clean Sphinx builds with no warnings
- **Usage Examples**: Added miRNA seed analysis workflow documentation

### Testing

- **New Test Coverage**: 232 new tests for miRNA design mode
  - Comprehensive unit tests for MiRNADesigner scoring
  - Schema validation tests for miRNA-specific columns
  - Integration tests for miRNA database functionality
- **Test Organization**: Normalized test markers for consistent CI/CD workflows
- **Documentation Tests**: Verified all doc builds and cross-references work correctly

### Dependencies

- No new runtime dependencies (leverages existing httpx, pydantic, pandera)
- Enhanced development dependencies for documentation generation

---

## [0.2.1] - 2025-10-24

```markdown
## [X.Y.Z] - YYYY-MM-DD

### New Features

- Brief description of new features

### Improvements

- Improvements to existing functionality
- Performance enhancements

### Bug Fixes

- Fixed specific issues
- Resolved edge cases

### 📊 Performance

- Performance improvements with metrics if available

### Testing

- New tests added
- Test coverage improvements
```

---

## [0.2.1] - 2025-10-24

### New Features

- **Chemical Modification System**: Comprehensive infrastructure for siRNA chemical modifications
  - Default modification patterns automatically applied to designed siRNAs (standard_2ome, minimal_terminal, maximal_stability)
  - New `--modifications` and `--overhang` CLI flags for workflow and design commands
  - FDA-approved Patisiran (Onpattro) pattern included in example library
- **Modification Metadata Models**: Pydantic models for StrandMetadata, ChemicalModification, Provenance tracking
- **FASTA Annotation System**: Merge modification metadata into FASTA headers with full roundtrip support
- **Remote FASTA Inputs**: Workflow supports `--input-fasta` with automatic HTTP download and caching
- **Enhanced Pandera Schemas**: Runtime DataFrame validation with @pa.check_types decorators, automatic addition of modification columns

### Improvements

- Modification columns (guide/passenger overhangs and modifications) now included in CSV outputs
- CLI `sequences show` command with JSON/FASTA/table output formats
- CLI `sequences annotate` command for merging metadata into FASTA files
- Standardized `+` separators in modification headers (backward compatible with `|`)
- Resource resolver for flexible input handling (local files, HTTP URLs)
- Improved type safety with Pandera schema validation on DesignResult.save_csv() and _generate_orf_report()

### Bug Fixes

- Fixed JSON metadata loading regression with StrandMetadata subscripting
- Resolved mypy typing issues for optional FASTA descriptions
- Fixed CLI output handling for modification metadata

### Documentation

- **Chemical Modification Review** (551 lines): Comprehensive analysis and integration guide
- **Modification Integration Guide** (543 lines): Developer documentation with code examples
- **Modification Annotation Spec** (381 lines): Complete FASTA header specification
- **Example Patterns Library**: 4 production-ready modification patterns with usage guide
- Updated README with chemical modifications feature documentation
- Remote FASTA usage documented in CLI and gene search guides

### Testing

- **18 new tests** for chemical modifications (100% passing):
  - 11 integration tests for workflow roundtrip validation
  - 7 tests validating example pattern files
- Added resource resolver unit tests (local paths, HTTP downloads, schemes)
- Extended modification metadata tests for delimiter compatibility
- All 164 tests passing with enhanced Pandera validation

### Dependencies

- No new runtime dependencies added (uses existing Pydantic, Pandera, httpx)

### Performance

- Removed Bowtie indexing (standardized on BWA-MEM2)
- Streamlined off-target analysis pipeline configuration

---

## [0.2.0] - 2025-09-27

### New Features

- **miRNA Database Cache System** (`sirnaforge cache`) - Local caching and management of miRNA databases from multiple sources with automatic updates
- **Comprehensive Data Validation** - Pandera DataFrameSchemas for type-safe output validation ensuring consistent CSV/TSV report formatting
- **Enhanced Thermodynamic Scoring** - Modified composite score to heavily favor (90%) duplex binding energy for improved siRNA selection accuracy
- **Workflow Input Flexibility** - Added FASTA file input support for custom transcript analysis workflows
- **Embedded Nextflow Pipeline** - Integrated Nextflow execution directly within Python API for scalable processing

### Improvements

- **Performance Optimization** - Parallelized off-target analysis and improved memory efficiency for large transcript sets
- **CLI Enhancement** - Better Unicode support, cleaner help text, and improved error reporting
- **Data Schema Validation** - Robust output validation with detailed error messages using modern Pandera 0.26.1 patterns
- **Documentation Overhaul** - Comprehensive testing guide, thermodynamic documentation, and improved API references
- **Development Workflow** - Enhanced Makefile with Docker testing categories, release validation, and conda environment support

### � Bug Fixes

- **Security Improvements** - Resolved security linting issues and improved dependency management
- **Off-target Analysis** - Fixed alignment indexing and improved multi-species database handling
- **CI/CD Pipeline** - Resolved build failures, improved test categorization, and enhanced release automation
- **Unicode Handling** - Fixed CLI display issues in various terminal environments

### 📊 Performance

- **10-100x Faster Dependencies** - Full migration to uv package manager for ultra-fast installs and environment management
- **Optimized Algorithms** - Improved thermodynamic calculation efficiency with better filtering strategies
- **Parallel Processing** - Enhanced concurrent execution for off-target analysis across multiple genomes

### Testing & Infrastructure

- **Enhanced Test Categories** - Smoke tests (256MB), integration tests (2GB), and full CI validation
- **Docker Improvements** - Multi-stage builds, intelligent entrypoint, and resource-aware testing
- **Release Automation** - Comprehensive GitHub Actions workflow with quality gates and artifact management

### Documentation

- **Testing Guide** - Comprehensive documentation for all test categories and Docker workflows
- **Thermodynamic Guide** - Detailed explanation of scoring algorithms and parameter optimization
- **CLI Reference** - Auto-generated command documentation with examples
- **Development Setup** - Streamlined onboarding with conda environment and uv integration

### Dependencies & Architecture

- **Modern Python Support** - Maintained compatibility across Python 3.9-3.12 with improved type safety
- **Pydantic Integration** - Enhanced data models with validation middleware and error handling
- **Containerization** - Production-ready Docker images with conda bioinformatics stack
- **Package Management** - Full uv adoption for dependency resolution and virtual environment management

## [0.1.0] - 2025-09-06

### Added

- Initial release of siRNAforge toolkit
- Core siRNA design algorithms with thermodynamic scoring
- Multi-database gene search (Ensembl, RefSeq, GENCODE)
- Rich command-line interface with Typer and Rich
- Comprehensive siRNA candidate scoring system
- Off-target prediction framework
- Nextflow pipeline integration for scalable analysis
- Docker containerization for reproducible environments
- Python API with Pydantic data models
- Comprehensive test suite with unit and integration tests
- Modern development tooling with uv, black, ruff, mypy

### Core Features

- **Gene Search**: Multi-database transcript retrieval
- **siRNA Design**: Algorithm-driven candidate generation
- **Quality Control**: GC content, structure, and specificity filters
- **Scoring System**: Composite scoring with multiple components
- **Workflow Orchestration**: End-to-end gene-to-siRNA pipeline
- **CLI Interface**: Rich, user-friendly command-line tools
- **Python API**: Programmatic access for automation

### Supported Operations

- `sirnaforge workflow`: Complete gene-to-siRNA analysis
- `sirnaforge search`: Gene and transcript search
- `sirnaforge design`: siRNA candidate generation
- `sirnaforge validate`: Input file validation
- `sirnaforge config`: Configuration display
- `sirnaforge version`: Version information

### Technical Stack

- **Language**: Python 3.9-3.12
- **Package Management**: uv for fast dependency resolution
- **Data Models**: Pydantic for type-safe data handling
- **CLI Framework**: Typer with Rich for beautiful output
- **Testing**: pytest with comprehensive coverage
- **Code Quality**: black, ruff, mypy for consistency
- **Containerization**: Multi-stage Docker builds
- **Pipeline**: Nextflow integration for scalability
- **Documentation**: Sphinx with MyST parser, Read the Docs theme

[Unreleased]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/austin-s-h/sirnaforge/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/austin-s-h/sirnaforge/releases/tag/v0.1.0
