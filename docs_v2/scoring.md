# Scoring

siRNAforge uses research-backed thermodynamic metrics to rank siRNA candidates. Higher composite scores indicate better predicted efficacy.

## Quick Reference

| Metric | Optimal Range | What It Means |
|--------|---------------|---------------|
| `composite_score` | 7-10 | Overall quality (higher = better) |
| `asymmetry_score` | ≥0.65 | Guide strand selection preference |
| `gc_content` | 40-55% | Stability vs. accessibility balance |
| `melting_temp` | 55-65°C | Duplex stability |
| `mfe` | -4 to -7 kcal/mol | Secondary structure stability |

## Composite Score

The composite score combines multiple factors with research-validated weights:

- **Thermodynamic asymmetry** (25%) - Guide strand preferentially enters RISC
- **GC content** (20%) - Balance between stability and accessibility
- **Target accessibility** (25%) - mRNA region accessibility
- **Off-target potential** (20%) - Specificity prediction
- **Empirical rules** (10%) - Position-specific sequence features

**Interpreting scores:**
- **8-10**: Excellent candidates for experiments
- **6-8**: Good candidates, may need validation
- **<6**: Consider alternatives

## Asymmetry Score

The most important single predictor of siRNA efficacy.

**Why it matters:** RISC preferentially loads the strand with the less stable 5' end. Higher asymmetry = correct strand selection = better knockdown.

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
| `quality_flags` | Any warnings or notes |

## References

1. Khvorova A et al. (2003) - Thermodynamic asymmetry and RISC loading
2. Schwarz DS et al. (2003) - Asymmetry rule for siRNA strand selection
3. Reynolds A et al. (2004) - Rational siRNA design guidelines
4. Ui-Tei K et al. (2004) - Guidelines for effective siRNAs
