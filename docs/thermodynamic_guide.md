# Thermodynamic Metrics Interpretation Guide

This guide helps you interpret the thermodynamic metrics in siRNAforge output files and optimize siRNA selection based on specific experimental needs.

## Quick Reference Table

| Metric | Optimal Range | Good Range | Poor Range | Units |
|--------|---------------|------------|------------|-------|
| **GC Content** | 40-55% | 35-60% | <35% or >65% | % |
| **Asymmetry Score** | 0.7-1.0 | 0.6-0.8 | <0.5 | 0-1 scale |
| **Paired Fraction** | 0.5-0.7 | 0.4-0.8 | <0.3 or >0.9 | 0-1 scale |
| **MFE** | -4 to -7 | -2 to -8 | <-10 or >0 | kcal/mol |
| **Duplex Stability ΔG** | see note | see note | — | kcal/mol |
| **Delta ΔG End** (`dg_5p - dg_3p`) | +2 to +5 | +1 to +6 | <0 | kcal/mol |
| **Melting Temp** | 65-75°C | 60-78°C | <50°C or >80°C | °C |
| **Off-target Count** | 0-1 | 0-3 | >5 | count |

**Note on `duplex_stability_dg`**: since 0.5.2 this is the ΔG of the full guide:passenger
duplex, which scales with length — a fully paired 21mer measures **-32 to -43 kcal/mol**
(median -39). There is no efficacy-validated window for it in siRNAforge: it is not
enforced as a filter and does not feed `composite_score`. Compare candidates by the
length-normalised `duplex_stability_score` within a run instead. Pre-0.5.2 values (roughly
-5 to -25) came from folding the guide against itself and are not comparable.

## Analyzing Your Results

### Step 1: Load and Examine Data

```bash
# View top candidates sorted by composite score
head -10 your_results/sirna_design/GENE_pass.csv

# Check distribution of key metrics
csvcut -c gc_content,asymmetry_score,mfe,duplex_stability_dg,delta_dg_end,melting_temp_c \
  your_results/sirna_design/GENE_pass.csv | head -20
```

The awk filters below resolve columns by header name rather than by position, so they keep
working when the column set changes (0.5.2 added `off_target_screened`).

### Step 2: Identify High-Quality Candidates

Look for siRNAs that meet multiple criteria:

```bash
# Filter for high-quality candidates (example thresholds)
# GC 40-55%, asymmetry >=0.7, MFE -7..-4, less stable guide 5-prime end, <=3 off-targets
awk -F',' '
NR==1 { for (i=1; i<=NF; i++) c[$i]=i; print; next }
$c["gc_content"]>=40 && $c["gc_content"]<=55 &&
$c["asymmetry_score"]>=0.7 &&
$c["mfe"]>=-7 && $c["mfe"]<=-4 &&
$c["delta_dg_end"]>=2 &&
$c["off_target_count"]<=3 { print }
' your_results/sirna_design/GENE_pass.csv
```

### Step 3: Troubleshoot Poor Performance

If few candidates meet optimal criteria:

#### **Low GC Content Issues**
- **Problem**: GC content consistently <35%
- **Solution**: Consider relaxing GC minimum to 30% or target different transcript regions
- **Alternative**: Focus on asymmetry and MFE scores instead

#### **Poor Asymmetry Scores**
- **Problem**: Most candidates have asymmetry_score <0.6
- **Solution**: Prioritize candidates with highest available asymmetry scores
- **Check**: Verify end stability differences (delta_dg_end should be positive)

#### **Overly Stable Duplexes**
- **Problem**: MFE values <-10 kcal/mol (guide secondary structure)
- **Solution**: Consider candidates with less negative (higher) MFE values
- **Alternative**: Test experimentally as some cell types handle stable duplexes better

## Application-Specific Guidelines

### For High-Efficiency Applications

When maximum knockdown is critical:

```bash
# Prioritize asymmetry and low off-targets
# High asymmetry, actually screened, very few off-targets, good end asymmetry
awk -F',' '
NR==1 { for (i=1; i<=NF; i++) c[$i]=i; print; next }
$c["asymmetry_score"]>=0.8 &&
$c["off_target_screened"]=="True" &&
$c["off_target_count"]<=1 &&
$c["delta_dg_end"]>=2 { print }
' results.csv
```

### For Broad Target Coverage

When targeting multiple isoforms:

```bash
# Balance efficiency with transcript coverage
# Moderate asymmetry, hits at least half the input transcripts, few off-targets
awk -F',' '
NR==1 { for (i=1; i<=NF; i++) c[$i]=i; print; next }
$c["asymmetry_score"]>=0.6 &&
$c["transcript_hit_fraction"]>=0.5 &&
$c["off_target_count"]<=5 { print }
' results.csv
```

### For Sensitive Cell Types

When working with difficult-to-transfect cells:

```bash
# Favor stability and moderate parameters
# Higher GC for stability, strong length-normalised duplex, mid-range Tm
awk -F',' '
NR==1 { for (i=1; i<=NF; i++) c[$i]=i; print; next }
$c["gc_content"]>=45 && $c["gc_content"]<=60 &&
$c["duplex_stability_score"]>=0.7 &&
$c["melting_temp_c"]>=65 && $c["melting_temp_c"]<=75 { print }
' results.csv
```

## Experimental Validation Tips

### Testing Multiple Candidates

1. **Select 3-5 candidates** with diverse metric profiles
2. **Include positive controls** with known effective siRNAs
3. **Test dose response** to optimize concentration

### Metric Correlation Analysis

Monitor which metrics correlate with experimental success:

```python
import pandas as pd
import seaborn as sns

# Load results and experimental data
results = pd.read_csv('sirna_results.csv')
experimental = pd.read_csv('knockdown_efficiency.csv')

# Merge and analyze correlations
combined = results.merge(experimental, on='id')
correlations = combined[['asymmetry_score', 'gc_content', 'mfe',
                        'delta_dg_end', 'knockdown_percent']].corr()

sns.heatmap(correlations, annot=True)
```

### Iterative Optimization

1. **Baseline Test**: Use default siRNAforge parameters
2. **Analyze Results**: Identify which metrics correlate with success
3. **Refine Parameters**: Adjust thresholds based on your system
4. **Validate**: Test refined predictions experimentally

## Common Troubleshooting

### No High-Quality Candidates

**Possible Causes:**
- Target sequence has unfavorable composition
- Overly strict filtering parameters
- Transcript region lacks optimal sites

**Solutions:**
1. Relax one parameter at a time (start with GC content)
2. Increase candidate pool size (`--top-n 50`)
3. Try different transcript isoforms
4. Consider alternative target regions

### All Candidates Have High Off-targets

**Approach:**
1. Prioritize candidates with lowest off-target counts
2. Use experimental validation to test specificity
3. Consider tissue-specific expression of off-targets
4. Implement additional experimental controls

### System-Specific Optimization

**Different organisms/cell types may require adjusted thresholds:**

- **Plant cells**: Often tolerate higher GC content (45-65%)
- **Primary cells**: May need more stable duplexes
- **Cancer cell lines**: Often more permissive of various parameters
- **In vivo applications**: Require stricter off-target criteria

## Advanced Analysis

### Composite Score Interpretation

The composite score integrates multiple factors. Understanding its components helps optimization:

```python
# Estimate component contributions (example weights)
def estimate_composite_components(row):
    """Estimate how each metric contributes to composite score"""
    gc_component = score_gc_content(row['gc_content']) * 0.2
    asymmetry_component = row['asymmetry_score'] * 0.3
    structure_component = score_mfe(row['mfe']) * 0.2
    offtarget_component = score_offtargets(row['off_target_count']) * 0.3

    return {
        'gc': gc_component,
        'asymmetry': asymmetry_component,
        'structure': structure_component,
        'offtarget': offtarget_component
    }
```

### Custom Scoring Functions

For specialized applications, implement custom scoring (see {doc}`tutorials/custom_scoring`):

```python
def custom_therapeutic_score(candidate):
    """Scoring optimized for therapeutic applications"""
    # Heavily weight safety (low off-targets)
    safety_score = 1.0 / (1.0 + candidate.off_target_count) * 0.5

    # Moderate weight on efficiency
    efficiency_score = candidate.asymmetry_score * 0.3

    # Stability for in vivo delivery
    stability_score = score_stability(candidate.gc_content, candidate.mfe) * 0.2

    return safety_score + efficiency_score + stability_score
```

---

**Remember**: These guidelines provide starting points. Experimental validation remains essential for confirming siRNA effectiveness in your specific system.
