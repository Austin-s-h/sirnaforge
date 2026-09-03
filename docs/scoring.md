# Scoring Overview

siRNAforge uses research-backed thermodynamic metrics to rank siRNA candidates. Higher composite scores indicate better predicted efficacy.

## Quick Reference

| Metric                | Optimal Range                   | What It Means                           |
| --------------------- | ------------------------------- | --------------------------------------- |
| `composite_score`     | 0-100 scale, higher is better   | Overall quality                         |
| `asymmetry_score`     | ≥0.65                           | Guide strand selection preference       |
| `gc_content`          | 40-55%                          | Stability vs. accessibility balance     |
| `melting_temp_c`      | 60-78°C                         | Duplex stability (nearest-neighbour Tm) |
| `mfe`                 | -4 to -7 kcal/mol               | Secondary structure stability           |
| `duplex_stability_dg` | -32 to -43 kcal/mol for a 21mer | Guide:passenger duplex ΔG               |

## Composite Score

`composite_score` is a weighted sum of seven sub-scores, each normalised to `[0, 1]`, computed
once per candidate by `sirnaforge.core.scoring.compute_composite`:

- **Thermodynamic asymmetry** (weight 0.12) - Guide strand preferentially enters RISC
- **GC content** (weight 0.10) - Balance between stability and accessibility
- **Target accessibility** (weight 0.13) - mRNA region accessibility
- **Empirical rules** (weight 0.15) - Position-specific sequence features
- **Off-target specificity** (weight 0.25) - Post-screen: decays with the _genuine_ off-target
  count (on-target, ortholog and repeat-mediated hits excluded), `exp(-count / 10)`
- **Isoform coverage** (weight 0.15) - Post-screen: fraction of the query gene's protein-coding
  isoforms the guide hits
- **Conservation** (weight 0.10) - Post-screen: fraction of requested non-query species with an
  ortholog hit

These weights (`ScoringWeights`, `weight_set_version = "2.0.0"`) sum to 1.00 and are the
_post-screen_ set. Off-target, isoform coverage and conservation cannot be evaluated until
transcriptome/miRNA screening has run, so **the score is computed once, after screening, not at
design time.** A candidate produced by the standalone `sirnaforge design` path (no screening)
still gets a `composite_score` from the same scorer, with those three terms simply inactive (see
below) — it is not a second, cheaper score.

### Weights are renormalised over the active term set

A term is _active_ for a candidate only when its sub-score could actually be computed:

- `off_target`, `isoform_coverage` and `conservation` are inactive before screening has run.
- `isoform_coverage` stays inactive if the query gene has no protein-coding transcript (an
  annotation gap, not a candidate defect).
- `conservation` stays inactive when the user requested no species beyond the query species —
  a single-species run has no evidence to compute it from.

`compute_composite` renormalises the _remaining_ weights to sum to 1 before combining them, so a
missing term is dropped rather than scored zero. This is why a single-species run is **neither
rewarded nor penalised** for lacking a conservation term: the weight that would have gone to
`conservation` is redistributed proportionally across the terms that did run, not left on the
table as an implicit zero.

`scored_after_screening` (bool) tells you which regime produced a given row's score, and
`weight_set_version` records which weight set. Scores from `weight_set_version` `1.x` (the
pre-issue-#80 five-term set) are **not comparable** to `2.x` scores — always compare candidates
within one run, one weight-set version.

### Per-term contribution columns

Each candidate carries `score_asymmetry`, `score_gc_content`, `score_accessibility`,
`score_empirical`, `score_off_target`, `score_isoform_coverage` and `score_conservation`: the
renormalised-weight × sub-score × 100 contribution of each active term. These sum to
`composite_score` (to floating-point tolerance) and are `None`/empty for any term that was
inactive for that candidate, so you can see exactly which terms carried a given score rather than
inferring it from the total alone.

### The design-time off-target proxy is now a diagnostic only

Before this weight set, the `off_target` term at design time was a proxy for guide
self-repetitiveness (repeated 7-mers _within_ the guide, unrelated to alignment against a
reference). That computation still runs, but it no longer feeds `composite_score` under any
name — it survives only as the unweighted diagnostic
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

The `*_pass.csv` and `*_all.csv` files include:

| Column                                             | Description                                                                                     |
| -------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `sirna_id`                                         | Unique identifier                                                                               |
| `guide_sequence`                                   | 21nt guide strand (5'→3')                                                                       |
| `passenger_sequence`                               | Passenger/sense strand                                                                          |
| `position`                                         | Start position in transcript                                                                    |
| `composite_score`                                  | Overall quality score                                                                           |
| `asymmetry_score`                                  | Thermodynamic asymmetry                                                                         |
| `gc_content`                                       | GC percentage                                                                                   |
| `melting_temp_c`                                   | Melting temperature (°C)                                                                        |
| `mfe`                                              | Minimum free energy (kcal/mol)                                                                  |
| `duplex_stability_dg`                              | Guide:passenger duplex ΔG (kcal/mol)                                                            |
| `dg_5p` / `dg_3p`                                  | Terminal 7 bp ΔG at each duplex end                                                             |
| `delta_dg_end`                                     | `dg_5p - dg_3p`; positive favours guide loading                                                 |
| `off_target_screened`                              | `False` means this candidate was never screened, so the hit counts are unknown rather than zero |
| `off_target_count`                                 | Genuine off-target hits only (on-target, ortholog and repeat-mediated hits excluded)            |
| `on_target_hits` / `ortholog_hits` / `repeat_hits` | The other three classes from the same four-way split                                            |
| `ortholog_species`                                 | Comma-joined canonical species names with at least one ortholog hit                             |
| `repeat_flagged` / `repeat_transcript_fraction`    | Design-time k-mer repeat verdict and the frequency it was based on                              |
| `isoform_coverage` / `conservation_score`          | Post-screen sub-scores (empty when inactive for that candidate)                                 |
| `score_*`                                          | Per-term contribution to `composite_score` (see Composite Score above)                          |
| `scored_after_screening` / `weight_set_version`    | Which scoring regime produced this row's `composite_score`                                      |
| `passes_filters`                                   | `PASS` or the first failed filter                                                               |

## References

1. Khvorova A et al. (2003) - Thermodynamic asymmetry and RISC loading
2. Schwarz DS et al. (2003) - Asymmetry rule for siRNA strand selection
3. Reynolds A et al. (2004) - Rational siRNA design guidelines
4. Ui-Tei K et al. (2004) - Guidelines for effective siRNAs
