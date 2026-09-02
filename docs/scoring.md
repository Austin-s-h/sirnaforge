# Scoring Overview

siRNAforge uses research-backed thermodynamic metrics to rank siRNA candidates. Higher composite scores indicate better predicted efficacy.

## Quick Reference

| Metric | Optimal Range | What It Means |
|--------|---------------|---------------|
| `composite_score` | 0-100 scale, higher is better | Overall quality |
| `asymmetry_score` | ≥0.65 | Guide strand selection preference |
| `gc_content` | 40-55% | Stability vs. accessibility balance |
| `melting_temp_c` | 60-78°C | Duplex stability (nearest-neighbour Tm) |
| `mfe` | -4 to -7 kcal/mol | Secondary structure stability |
| `duplex_stability_dg` | -32 to -43 kcal/mol for a 21mer | Guide:passenger duplex ΔG |

## Composite Score

The composite score combines multiple factors with research-validated weights:

- **Thermodynamic asymmetry** (15%) - Guide strand preferentially enters RISC
- **GC content** (15%) - Balance between stability and accessibility
- **Target accessibility** (20%) - mRNA region accessibility
- **Off-target potential** (30%) - Design-time proxy: repeated 7-mers *within* the guide.
  Transcriptome and miRNA screening run after design and gate `passes_filters`; they do not
  contribute to `composite_score`.
- **Empirical rules** (20%) - Position-specific sequence features

The score is reported on a 0-100 scale. Compare candidates within one run rather than
against a fixed cut-off: the attainable maximum depends on the configured weights, and
scores are not comparable across versions (0.5.2 corrected the thermodynamics).

## Asymmetry Score

The most important single predictor of siRNA efficacy.

| Score | Interpretation |
|-------|----------------|
| 0.8-1.0 | Excellent - strong guide strand bias |
| 0.65-0.8 | Good - likely correct strand selection |
| 0.5-0.65 | Moderate - mixed strand loading possible |
| <0.5 | Poor - passenger strand may dominate |

**Research basis:** Khvorova et al. (2003), Schwarz et al. (2003)

## GC Content

Affects duplex stability and target accessibility.

| Range | Effect |
|-------|--------|
| <35% | Unstable duplex, poor RISC loading |
| 35-40% | Acceptable, monitor stability |
| **40-55%** | **Optimal range** |
| 55-60% | Acceptable, may reduce accessibility |
| >60% | Overly stable, poor target release |

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

| Column | Description |
|--------|-------------|
| `sirna_id` | Unique identifier |
| `guide_sequence` | 21nt guide strand (5'→3') |
| `passenger_sequence` | Passenger/sense strand |
| `position` | Start position in transcript |
| `composite_score` | Overall quality score |
| `asymmetry_score` | Thermodynamic asymmetry |
| `gc_content` | GC percentage |
| `melting_temp_c` | Melting temperature (°C) |
| `mfe` | Minimum free energy (kcal/mol) |
| `duplex_stability_dg` | Guide:passenger duplex ΔG (kcal/mol) |
| `dg_5p` / `dg_3p` | Terminal 7 bp ΔG at each duplex end |
| `delta_dg_end` | `dg_5p - dg_3p`; positive favours guide loading |
| `off_target_screened` | `False` means this candidate was never screened, so the hit counts are unknown rather than zero |
| `passes_filters` | `PASS` or the first failed filter |

## References

1. Khvorova A et al. (2003) - Thermodynamic asymmetry and RISC loading
2. Schwarz DS et al. (2003) - Asymmetry rule for siRNA strand selection
3. Reynolds A et al. (2004) - Rational siRNA design guidelines
4. Ui-Tei K et al. (2004) - Guidelines for effective siRNAs
