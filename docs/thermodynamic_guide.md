# Thermodynamic Metrics Interpretation Guide

This guide helps you interpret the thermodynamic metrics in siRNAforge output files and optimize siRNA selection based on specific experimental needs.

## Quick Reference Table

| Metric | Optimal Range | Good Range | Poor Range | Units |
|--------|---------------|------------|------------|-------|
| **GC Content** | 40-55% | 35-60% | <35% or >65% | % |
| **Asymmetry Score** | 0.7-1.0 | 0.6-0.8 | <0.5 | 0-1 scale |
| **Paired Fraction** | 0.5-0.7 | 0.4-0.8 | <0.3 or >0.9 | 0-1 scale |
| **MFE** | -4 to -7 | -2 to -8 | <-10 or >0 | kcal/mol |
| **Duplex Stability ΔG** | -18 to -22 | -15 to -25 | <-30 or >-10 | kcal/mol |
| **Delta ΔG End** | +2 to +5 | +1 to +6 | <0 | kcal/mol |
| **Melting Temp** | 65-75°C | 60-78°C | <50°C or >80°C | °C |
| **Off-target Count** | 0-1 | 0-3 | >5 | count |

## Analyzing Your Results

### Step 1: Load and Examine Data

```bash
# View top candidates sorted by composite score
head -10 your_results/sirna_design/GENE_pass.csv

# Check distribution of key metrics
cut -d',' -f6,7,9,11,15,16,17 your_results/sirna_design/GENE_pass.csv | head -20
```

### Step 2: Identify High-Quality Candidates

Look for siRNAs that meet multiple criteria:

```bash
# Filter for high-quality candidates (example thresholds)
awk -F',' '
NR==1 {print; next}  # Print header
$6>=40 && $6<=55 &&    # GC content 40-55%
$7>=0.7 &&             # Asymmetry score ≥0.7
$9>=-7 && $9<=-4 &&    # MFE between -7 and -4
$15>=+2 &&             # Positive delta_dg_end
$17<=3                 # Low off-target count
{print}' your_results/sirna_design/GENE_pass.csv
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
- **Problem**: MFE values <-10 kcal/mol or very negative duplex_stability_dg
- **Solution**: Consider candidates with less negative (higher) MFE values
- **Alternative**: Test experimentally as some cell types handle stable duplexes better

## Application-Specific Guidelines

### For High-Efficiency Applications

When maximum knockdown is critical:

```bash
# Prioritize asymmetry and low off-targets
awk -F',' '
NR==1 {print; next}
$7>=0.8 &&             # High asymmetry score
$17<=1 &&              # Very low off-targets
$15>=+2                # Good end asymmetry
{print}' results.csv
```

### For Broad Target Coverage

When targeting multiple isoforms:

```bash
# Balance efficiency with transcript coverage
awk -F',' '
NR==1 {print; next}
$7>=0.6 &&             # Moderate asymmetry acceptable
$18>=0.5 &&            # Good transcript hit fraction
$17<=5                 # Moderate off-target tolerance
{print}' results.csv
```

### For Sensitive Cell Types

When working with difficult-to-transfect cells:

```bash
# Favor stability and moderate parameters
awk -F',' '
NR==1 {print; next}
$6>=45 && $6<=60 &&    # Higher GC for stability
$11>=-20 &&            # Moderate duplex stability
$16>=65 && $16<=75     # Optimal melting temp
{print}' results.csv
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
