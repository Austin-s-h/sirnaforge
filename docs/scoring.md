# Scoring Overview

siRNAforge uses research-backed thermodynamic metrics to rank siRNA candidates. Higher composite scores indicate better predicted efficacy.

## Quick Reference

| Metric                   | Optimal Range                   | What It Means                                               |
| ------------------------ | ------------------------------- | ----------------------------------------------------------- |
| `composite_score`        | 0-100 scale, higher is better   | Overall quality; **null before screening**                  |
| `design_score`           | 0-100 scale, higher is better   | Design-stage quality, not comparable with `composite_score` |
| `asymmetry_score`        | ≥0.65                           | Guide strand selection preference                           |
| `gc_content`             | 40-55%                          | Stability vs. accessibility balance                         |
| `melting_temp_c`         | 60-78°C                         | Duplex stability (nearest-neighbour Tm)                     |
| `mfe`                    | -4 to -7 kcal/mol               | Guide self-structure (0.0 = open chain)                     |
| `duplex_stability_dg`    | -32 to -43 kcal/mol for a 21mer | Guide:passenger duplex ΔG                                   |
| `target_accessibility_p` | higher is better, log-scaled    | P(mRNA site's seed-paired 8-mer open)                       |

## Two scores, three named weight vectors

Since issue #96 every weight vector is **hand-authored, named, sums to 1.0 and is written to the run
manifest**, and nothing rescales one at runtime. There is no renormalisation and no divisor. A vector
is _chosen_ by stage and design mode, never combined:

```
design_v4                        postscreen_sirna_v4          postscreen_mirna_v4
  target_accessibility  0.40       off_target           0.25    off_target            0.20
  asymmetry             0.35       target_accessibility 0.30    target_accessibility  0.24
  gc_content            0.25       asymmetry            0.25    asymmetry             0.20
                        ----       gc_content           0.20    gc_content            0.16
                        1.00                            ----    ago_start             0.10
                                                        1.00    supp_13_16            0.10
                                                                                      ----
                                                                                      1.00
```

`postscreen_mirna_v4`'s four shared terms are exactly `0.80 x postscreen_sirna_v4`; its two
biogenesis terms hold the remaining 0.20. It had a seventh term, `pos1_mismatch` at 0.05, until
issue #102 measured it **exactly constant at 0.0** on all 13,415 scored candidates of the public
baseline — constant by construction, because the passenger is the exact reverse complement of the
guide, so guide position 1 always pairs. A term that ranks nothing must not consume weight, so it
was removed from the vector. It is still computed, and the pairing state stays on the row as
`guide_pos1_base` and `pos1_pairing_state` — but `score_pos1_mismatch` is a *contribution* column, so
with the term in no vector that column is now **always null**.

`postscreen_sirna_v4` is exactly `design_v4`'s terms plus `off_target` — one extra term, cleanly
interpretable. Every row records the vector that produced it in `weight_vector`, and the manifest's
`scoring.vectors` block maps that name to the numbers, so any score is traceable to the weights that
made it.

| Field             | Vector                                        | Available                                                                                         |
| ----------------- | --------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| `design_score`    | `design_v4`                                   | at design time, from the three terms computable before screening                                  |
| `composite_score` | `postscreen_sirna_v4` / `postscreen_mirna_v4` | only after off-target screening produced usable evidence for that candidate; **null before that** |

**The two are not comparable.** They are different vectors over different term sets, and they also
weight their _shared_ terms differently, so **neither is systematically the larger**. Issue #102
deleted the claim that `design_score` is the more optimistic number: all features at 0.5 with
`off_target = 1.0` gives `design_v4` **50.0** against `postscreen_sirna_v4`'s **62.5**, because the
post-screen vector spends 0.25 on a term the design vector does not have and takes it from asymmetry
and accessibility. Do not rank a mixture of the two; the workflow does not (see
_What `top_candidates` excludes_).

The three declared terms:

- **Target accessibility** (0.40 design / 0.30 post-screen) — log-scaled RNAplfold probability that
  the 8 nt of the mRNA target site pairing guide positions 1-8 are unpaired. The guide seed pairs the
  target site's **3' end**, so that is the end scored; the 5'-end 8-mer is a measured near-null
  (ρ +0.07 vs +0.27 against knockdown). See `docs/models_and_scoring.md` §2.5.
- **Thermodynamic asymmetry** (0.35 / 0.25) — guide strand preferentially enters RISC.
- **GC content** (0.25 / 0.20) — balance between stability and accessibility. `exp(-((GC%-40)/10)^2)`,
  a Gaussian centred on 40% GC — note that 40 is not the midpoint of the default GC filter window
  (35-60), so the score's optimum and the gate's optimum are different numbers.

plus, post-screen only:

- **Off-target specificity** (0.25 siRNA / 0.20 miRNA) — decays with the _genuine_ off-target count
  (on-target, ortholog and repeat-mediated hits excluded), `exp(-count / 10)`.

and in `--design-mode mirna` only, two biogenesis terms that used to be an undeclared bonus:

- **`ago_start`** (0.10) — A/U at guide position 1, the Argonaute loading preference.
- **`supp_13_16`** (0.10) — low 3' supplementary pairing potential at guide positions 13-16. Its
  declared endpoint is _specificity_, and nothing has measured it against that endpoint.

Because both post-screen vectors sum to 1.0, the two modes are on **one scale**: a candidate whose
biogenesis sub-scores match its other sub-scores scores identically in either mode, at every level.
A miRNA candidate earning nothing on the two biogenesis terms scores 0.80 of the equivalent siRNA
candidate — and that 0.80 is four weights you can read in the manifest and attribute term by term on
the row, not a factor applied to the whole vector.

These are **declared expert priors, not fitted values, and every one of them is `experimental`.**
`src/sirnaforge/models/scoring_profile.py` is the term registry: per term it records the molecule,
strand and positions read, the endpoint claimed, the applicability condition, the formula and
transform, the missing-value policy, the declared and **attainable** feature ranges, the evidence
source and the evidence status. Nothing in siRNAforge is `validated` — no weight has been held out
and independently replicated — so a composite score ranks candidates and is neither a predicted
knockdown level nor a probability of any safety event.

Two terms carry benchmark evidence against measured knockdown (`target_accessibility` ρ +0.267,
`gc_content` cluster-robust β +0.152 on the Huesken panel, n = 2,816 / 41 transcripts) and
`off_target`'s weight rests on a screening-scope-dependent variance share rather than on any
efficacy or safety measurement. `design_v4`'s three numbers are round numbers awaiting sign-off.

### A/U content at guide positions 1-5

Computed on every candidate as `au_1_5` in `component_scores`, and **scored by no default vector**.
It is the strongest single feature on the Huesken panel after guide position 1 alone (ρ +0.378, against
ρ −0.017 / p = 0.37 for the same count at guide positions 17-21), and issue #97's decision D5 accepted
it as a scored term with the window **pre-declared** at 1-5 — no ρ-maximising sweep.

D5's pre-declared rule for whether it replaces `asymmetry` or joins it was executed in
`scripts/validate_scoring_profiles.py`: regress `au_1_5` out of `asymmetry_score`, then test the
residual against efficacy with by-transcript clustering. The residual survived (β +0.0533,
cluster-robust SE 0.0213, t = 2.50, p = 0.017, 41 clusters), so **both terms score** in the candidate
vector recorded as the `au_1_5_experimental` profile. It is not a default in 0.7.1 for two reasons,
both stated rather than implied: on predeclared held-out transcripts it buys 0.008 of ρ, and wiring
it into `postscreen_sirna_v4` needs the post-screen feature assembly to forward it, which is not part
of this change.

### Computed and reported, but not scored

Removing these from the composite is what makes "no renormalisation" achievable rather than merely
relocated: `conservation` is `None` on single-species runs and `isoform_coverage` is `None` on the
`design_from_sequence`/miRNA paths, so any vector containing them needs either variant vectors or
arithmetic. Every term that remains is universally computable.

| Quantity                        | Column               | What reads it                                                                         |
| ------------------------------- | -------------------- | ------------------------------------------------------------------------------------- |
| Empirical design rules          | `empirical_score`    | the `min_empirical_score` gate (LOW_EMPIRICAL_SCORE)                                  |
| Protein-coding isoform coverage | `isoform_coverage`   | the optional `--min-isoform-coverage` gate (LOW_ISOFORM_COVERAGE), **off by default** |
| Cross-species conservation      | `conservation_score` | nothing — reported for interpretation                                                 |
| Guide self-structure            | `paired_fraction`    | the `max_paired_fraction` gate (EXCESS_PAIRING)                                       |

The empirical rubric no longer judges guide position 1. It paid +0.1 for G/C there while the
biogenesis rule paid `ago_start` for A/U at the same base; measured over 29,605 candidates G/C gained
+1.6 empirical points and lost 7.9 to the biogenesis adjustment, so a declared 0.15-weight term was
overridden ~5× by an undeclared one and `empirical` ended up with a _negative_ variance share. A/U
wins. The rubric therefore attains only `{0.4, 0.5, 0.6}`, and `min_empirical_score` is bounded
accordingly.

### A term that cannot be computed yields no score

Weights are never redistributed, so there is nothing to fall back to. If a vector's term cannot be
computed for a candidate, that candidate carries no score on that vector — `design_score` or
`composite_score` stays null and the run logs it. In practice the only residual case is
`target_accessibility` for a site too close to the transcript 5' end for the scoring window, which is
per-candidate and numerically negligible (0 of 2,492 TP53 sites). A candidate scored with no
transcript context at all — possible only by calling the scorer directly — has no `design_score`.

`scored_after_screening` (bool) tells you whether a row's `composite_score` exists;
`weight_set_version` and `weight_vector` record which weights produced it. Scores are **not
comparable across major versions** — `1.x` is the pre-issue-#80 five-term set, `2.x` scored guide
self-structure in the slot `3.x` gave to real target-site accessibility (issue #95), and `4.x` is
this restructuring. Always compare candidates within one run, one weight-set version.

### What `top_candidates` excludes

`top_candidates` is rebuilt after screening from the candidates that clear **three** gates, so it
can be shorter than `min(top_n, number passing)`. (`candidates_all.csv` and `candidates_pass.csv`
are filtered on `passes_filters` alone, but both are written in the re-ranked order.) The gates:

1. `passes_filters` is `PASS` — a failed off-target filter is a rejection, not a low score;
2. `repeat_flagged` is `False` — a guide that saturates the query transcriptome is excluded even
   when screening never ran;
3. `scored_after_screening` is `True`, **whenever some but not all candidates were scored after
   screening**. A design-stage score is on a different vector and is systematically the more
   optimistic number, so letting it compete against post-screen neighbours would put exactly the
   candidates whose evidence is missing at the top. Those rows stay in `candidates_all.csv` with
   their `design_score`, and the run logs an ERROR naming the count. If _no_ candidate was scored
   after screening (a wholly failed or wholly pre-screen run) the list is internally consistent and
   this gate does not apply.

### Per-term contribution columns

Each candidate carries `score_off_target`, `score_target_accessibility`, `score_asymmetry`,
`score_gc_content` and — in miRNA mode — `score_ago_start` and `score_supp_13_16`: the declared
weight × sub-score × 100 contribution of each term. A column is empty for a term outside the vector
that scored that row. `score_pos1_mismatch` is still emitted for schema stability and is **always
empty** since #102 removed the term from the miRNA vector.

**They sum to the score, exactly, in both modes.** Nothing is added or divided afterwards. Before
issue #96, `--design-mode mirna` computed
`(Σ contributions + bonus × 100) / (1 + max_bonus)` with `max_bonus = 0.25`, so every declared weight
was silently scaled by 0.80 and a candidate earning no bonus reported a score 20% below the sum of
its own contribution columns.

### The design-time off-target proxy is a diagnostic only

The `off_target` term at design time was once a proxy for guide self-repetitiveness (repeated 7-mers
_within_ the guide, unrelated to alignment against a reference). That computation still runs, but it
feeds no score under any name — it survives only as the unweighted diagnostic
`component_scores["design_off_target_proxy"]`.

## Asymmetry Score

The most important single predictor of siRNA efficacy.

| Score    | Interpretation                           |
| -------- | ---------------------------------------- |
| 0.8-1.0  | Excellent - strong guide strand bias     |
| 0.65-0.8 | Good - likely correct strand selection   |
| 0.5-0.65 | Moderate - mixed strand loading possible |
| <0.5     | Poor - passenger strand may dominate     |

**Research basis:** Khvorova et al. (2003), Schwarz et al. (2003)

## GC Content

Affects duplex stability and target accessibility.

| Range      | Effect                               |
| ---------- | ------------------------------------ |
| <35%       | Unstable duplex, poor RISC loading   |
| 35-40%     | Acceptable, monitor stability        |
| **40-55%** | **Optimal range**                    |
| 55-60%     | Acceptable, may reduce accessibility |
| >60%       | Overly stable, poor target release   |

## Melting Temperature

Temperature at which 50% of duplexes dissociate.

- **<55°C**: Unstable, may dissociate prematurely
- **55-65°C**: Optimal for mammalian cells
- **65-75°C**: Acceptable, verify experimentally
- **>75°C**: May resist RISC processing

## Minimum Free Energy (MFE)

Predicts secondary structure stability of the guide strand.

- **>0 kcal/mol**: Unstable, poor structure
- **-2 to -4 kcal/mol**: Minimal structure (good)
- **-4 to -8 kcal/mol**: Moderate structure (optimal)
- **<-10 kcal/mol**: Strong self-structure (may reduce activity)

## Filtering Recommendations

### Standard (most applications)

```bash
sirnaforge workflow GENE --gc-min 35 --gc-max 60
```

### Stringent (publication-quality)

```bash
sirnaforge workflow GENE --gc-min 40 --gc-max 55 --top-n 30
```

### Relaxed (difficult targets)

```bash
sirnaforge workflow GENE --gc-min 30 --gc-max 65
```

## Output Columns

The `candidates_pass.csv` and `candidates_all.csv` files include:

| Column                                                            | Description                                                                                                            |
| ----------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `sirna_id`                                                        | Unique identifier                                                                                                      |
| `guide_sequence`                                                  | 21nt guide strand (5'→3')                                                                                              |
| `passenger_sequence`                                              | Passenger/sense strand                                                                                                 |
| `position`                                                        | Start position in transcript                                                                                           |
| `design_score`                                                    | Design-stage score on `design_v4` (3 terms), available without screening                                               |
| `composite_score`                                                 | Post-screen score on `postscreen_{sirna,mirna}_v4`; empty before screening                                             |
| `asymmetry_score`                                                 | Thermodynamic asymmetry                                                                                                |
| `gc_content`                                                      | GC percentage                                                                                                          |
| `melting_temp_c`                                                  | Melting temperature (°C)                                                                                               |
| `mfe`                                                             | Minimum free energy (kcal/mol)                                                                                         |
| `duplex_stability_dg`                                             | Guide:passenger duplex ΔG (kcal/mol)                                                                                   |
| `dg_5p` / `dg_3p`                                                 | Terminal 7 bp ΔG at each duplex end                                                                                    |
| `delta_dg_end`                                                    | `dg_5p - dg_3p`; positive favours guide loading                                                                        |
| `off_target_screened`                                             | `False` means the screen was incomplete, so the hit counts are a lower bound, not a total                              |
| `off_target_count`                                                | The liability population: `off_target` **plus** `undetermined` (on-target, ortholog and repeat-mediated hits excluded) |
| `on_target_hits` / `ortholog_hits` / `repeat_hits`                | Three of the other four classes from the same five-way split                                                           |
| `undetermined_hits`                                               | The part of `off_target_count` whose class could not be decided for want of a transcript index                         |
| `ortholog_species`                                                | Comma-joined canonical species with an ortholog hit; `ortholog_evidence` says which tier decided each                  |
| `repeat_flagged` / `repeat_transcript_fraction`                   | Design-time k-mer repeat verdict and the frequency it was based on                                                     |
| `isoform_coverage` / `conservation_score`                         | Reported, unscored (empty when not computable); isoform coverage feeds an optional gate                                |
| `empirical_score`                                                 | Reported, unscored; the `min_empirical_score` gate input                                                               |
| `score_*`                                                         | Per-term contribution, summing exactly to the score (see above)                                                        |
| `scored_after_screening` / `weight_set_version` / `weight_vector` | Which stage, weight set and named vector produced this row's score                                                     |
| `passes_filters`                                                  | `PASS` or the first failed filter                                                                                      |

## References

1. Khvorova A et al. (2003) - Thermodynamic asymmetry and RISC loading
2. Schwarz DS et al. (2003) - Asymmetry rule for siRNA strand selection
3. Reynolds A et al. (2004) - Rational siRNA design guidelines
4. Ui-Tei K et al. (2004) - Guidelines for effective siRNAs
