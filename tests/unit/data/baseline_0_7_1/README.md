# `baseline_0_7_1` — tracked slice of the frozen 0.7.1 canonical public reference run

## Provenance

| Item | Value |
| --- | --- |
| Produced by | `scripts/build_baseline_0_7_1_fixture.py` |
| Cut from | the frozen canonical run under `work/baseline_0.7.1/run` (gitignored; see its `MANIFEST.md`) |
| Code that produced the run | `0ca8b88607ad19242370ae7ba9e055f420976a25` (`integration/0.7.1`) |
| Target | **TP53**, 34 protein-coding transcripts fetched live from Ensembl REST during the run |
| Screening reference | Ensembl GRCh38 cDNA release of 2026-03-24 (`Homo_sapiens.GRCh38.cdna.all.fa.gz`, md5 `b48fe40e46cff8a68ec28462feb82bf4`), reduced to the primary assembly, one transcript per gene (longest), **plus all 40 TP53 transcripts**: 35,757 records / 102,612,737 bases |
| Aligner | `bwa-mem2 2.2.1` in the siRNAforge container; `-a -k 12 -T 15 -w 100`, `max_hits` effectively unlimited |
| miRNA reference | MirGeneDB human (`https://www.mirgenedb.org/fasta/hsa?mat=1`, md5 `60e3ee19ea3bb76b4dcdde32ed10311b`) |
| Weight set | `SCORING_WEIGHT_SET_VERSION 4.0.0`, vector `postscreen_sirna_v4` |
| Redistribution | Ensembl (Apache 2.0 / no restriction) and MirGeneDB are public; these files contain derived coordinates and 21-mers only, no third-party sequence file is vendored |

## Files

| File | Rows | Content |
| --- | --- | --- |
| `candidates_sample.csv` | 292 | candidate rows, column layout **identical** to `candidates_all.csv` |
| `offtarget_hits_sample.tsv` | 2,173 | per-alignment rows, layout identical to `combined_offtargets.tsv` |
| `mirna_hits_sample.tsv` | 268 | per-hit rows, layout identical to `combined_mirna_hits.tsv` |

28 distinct guide sequences over 34 transcripts. Rows per guide: min 1, median 3, max 33.
Hit rows are keyed on `qname` = the candidate `id` that deduplication submitted to the aligner.

## What it is good for

* **Guide-vs-`id` join key (#103).** One guide carries up to 33 candidate rows here, so a report
  that joins hits on `id` reads 32 of them as zero-off-target.
* **`mirna_hits_high_risk == mirna_hits_0mm_seed` (#101).** True on all 292 rows, with 73 rows
  non-zero, so `max_mirna_perfect_seed` and `fail_on_high_risk_mirna` fail exactly the same rows.
* **`transcriptome_hits_total == off_target_count` (#101).** True on all 292 rows.
* **The `EXCESS_OFF_TARGETS` population and its `nm` rescope (#101/D3).** 42 rows across 9 guides
  exceed `max_off_target_count = 15`; the `nm` distribution is 1,700 / 425 / 0 / 1 / 1 / 2 / 44 for
  `nm` = 0…6, so both the `nm` ≤ 2 majority and the clipped-partial tail are represented.
* **Gate masking (#100/#103).** `passes_filters` here is `TRANSCRIPTOME_PERFECT_MATCH` on 254 rows,
  `MIRNA_PERFECT_SEED` on 10 and `PASS` on 28 — while several of those rows independently fail
  `LOW_ASYMMETRY` too. The single status column cannot express that.
* **Post-screen scoring columns.** `composite_score` plus the seven `score_*` columns the shared
  writer emits. The **four** `postscreen_sirna_v4` terms sum to `composite_score` exactly; the three
  miRNA-mode terms (`score_ago_start`, `score_pos1_mismatch`, `score_supp_13_16`) are null on every
  row here, which is what makes that decomposition exact.

**These claims are asserted, not merely written down:** `tests/unit/test_baseline_0_7_1_fixture.py`
pins the column layouts and dtypes, both redundancy identities, the guide-vs-`id` multiplicity, the
gate-masking population, the `nm` tail and the repeat-column trap below. Thresholds are asserted
there as the literals **this run used** (`max_off_target_count = 15`, `min_asymmetry_score = 0.65`),
never as live model defaults, because 0.7.1 moves several of those defaults.

⚠️ **Scale.** 292 rows is **0.8%** of the run's 34,863. This is a correctness fixture. #103's payload
budget, and anything else that needs magnitude, cannot be exercised from it.

## Verification limits — read before asserting on these numbers

1. **`hit_class`, `matched_symbol` and `symbol_lookup_missing` are absent.** `0ca8b88` does not
   persist them; adding them is WP0c's job. Anything needing per-hit classes must recompute them
   with `core/hit_classification.classify_hit`.
2. **`ortholog_hits` is 0 and `conservation_score` is null on every row** — the run was human-only.
   These files cannot test ortholog recognition or the conservation denominator.
3. **`repeat_hits` is 0 and `repeat_flagged` is False on every row** — the run recorded
   `repeat_summary: skipped / reference_unavailable` because repeat detection needs the
   transcriptome-materialisation path that this configuration had to disable. Do **not** read a zero
   here as "no repeats": an offline `RepeatDetector` scan of the full run flags 424 of its 2,400
   guides. These files cannot test the `REPEAT` hit class. No column records whether the scan ran,
   so "scanned, found none" and "never scanned" are byte-indistinguishable; the test asserting that
   fails deliberately once a repeat-evidence status column is added.
4. **Off-target counts are lower bounds.** The screened reference is 8% of the Ensembl human cDNA
   release by bases and holds one transcript per gene, because the full reference cannot be
   bwa-mem2-indexed on the machine that produced this run. Absolute magnitudes are not
   full-transcriptome magnitudes; assert on structure and relations, not on totals.
5. **`isoform_coverage` is populated but never gated** (`min_isoform_coverage` defaults to `None`).
6. **The 6 TP53 transcripts present in the reference but absent from the 34 designed-from
   transcripts are counted as genuine off-targets** in these rows, because the classifier had no
   transcript→gene index. That is a faithful record of `0ca8b88` behaviour under this
   configuration, not a target state — expect these counts to change once WP0c lands.
7. **`weight_set_version` is `4.0.0` and will break** when 0.7.1 bumps it. Any test asserting the
   literal value must be updated deliberately, not silently regenerated.
8. **The run's artifacts self-report `sirnaforge 0.6.0`.** The container image carries the 0.6.0
   wheel while the executing code was the mounted `0ca8b88` tree. The SHA above is authoritative.

## Regenerating

```bash
uv run python scripts/build_baseline_0_7_1_fixture.py \
    --run-dir work/baseline_0.7.1/run \
    --out-dir tests/unit/data/baseline_0_7_1
```

Deterministic given the same run: guides are chosen by sorted uniqueness within four buckets
(above the off-target cap, perfect miRNA seed, >5 rows per guide, passing) and rows are sorted.
Regenerating against a **different** run changes the numbers quoted above; update this README in
the same commit if that happens.
