# Data Models, Filtering & Scoring: Technical Reference

> **Academic rigor note**: This document provides comprehensive documentation of siRNAforge's data models, filtering criteria, and scoring algorithms with citations and justification for critical thresholds. Sections marked with `[REVIEW NEEDED]` indicate areas requiring additional expert review or validation.

## Overview

siRNAforge implements a multi-stage computational pipeline for siRNA design that relies on validated data models and research-backed scoring algorithms. This document describes:

1. **Data Models** - Pydantic-validated structures for siRNA candidates and analysis results
2. **Filter Criteria** - Evidence-based thresholds for candidate quality control
3. **Scoring Algorithms** - Composite scoring with thermodynamic and empirical components
4. **Threshold Justification** - Literature citations and rationale for default parameters

---

## 1. Core Data Models

### 1.1 SiRNACandidate

The `SiRNACandidate` model represents a complete siRNA duplex with all computed properties.

```python
class SiRNACandidate(BaseModel):
    """Individual siRNA candidate with computed thermodynamic and efficacy properties."""
    
    # Identity (unique identifier and source)
    id: str                    # Format: SIRNAF_{transcript}_{start}_{end}
    transcript_id: str         # Source transcript (e.g., ENST00000269305)
    position: int              # 1-based start position in transcript
    
    # Duplex sequences
    guide_sequence: str        # Antisense strand (loaded into RISC)
    passenger_sequence: str    # Sense strand (typically degraded)
    
    # Basic properties
    gc_content: float          # GC percentage (optimal: 35-60%)
    length: int                # Duplex length (typically 21 nt)
    
    # Thermodynamic properties
    asymmetry_score: float     # RISC loading preference (optimal: ≥0.65)
    duplex_stability: float    # ΔG in kcal/mol (optimal: -15 to -25)
    
    # Secondary structure
    structure: str             # Dot-bracket notation
    mfe: float                 # Minimum free energy (optimal: -2 to -8 kcal/mol)
    paired_fraction: float     # Fraction paired bases (optimal: 0.4-0.6)
    
    # Off-target metrics
    off_target_count: int      # Potential off-target sites (goal: ≤3)
    transcriptome_hits_0mm: int   # Perfect match hits
    transcriptome_hits_1mm: int   # 1-mismatch hits
    transcriptome_hits_2mm: int   # 2-mismatch hits
    mirna_hits_total: int         # Total miRNA seed matches
    mirna_hits_0mm_seed: int      # Perfect seed matches
    
    # Scoring
    composite_score: float     # Overall quality (0-100 scale)
    component_scores: dict     # Individual scoring components
    passes_filters: bool|FilterStatus  # Quality control status
```

#### Sequence Validation

All sequences undergo validation:
- **Allowed nucleotides**: A, T, C, G, U
- **Length constraints**: 19-23 nucleotides (siRNA length)
- **Strand matching**: Guide and passenger must be same length

### 1.2 DesignParameters

Configuration model for the design workflow:

```python
class DesignParameters(BaseModel):
    """Complete configuration for siRNA design workflow."""
    
    # Design mode
    design_mode: DesignMode     # "sirna" or "mirna"
    
    # Sequence parameters
    sirna_length: int = 21     # Duplex length (19-23 nt)
    top_n: int = 50            # Number of candidates to return
    
    # Quality control
    filters: FilterCriteria    # Threshold parameters
    scoring: ScoringWeights    # Component weights
    
    # Chemical modifications
    apply_modifications: bool = True
    modification_pattern: str = "standard_2ome"
    default_overhang: str = "dTdT"
```

### 1.3 FilterCriteria

Threshold parameters for quality filtering:

```python
class FilterCriteria(BaseModel):
    """Quality filters based on thermodynamic and empirical criteria."""
    
    # GC content (literature: 30-60%, optimal: 40-55%)
    gc_min: float = 35.0
    gc_max: float = 60.0
    
    # Sequence composition
    max_poly_runs: int = 3     # Max consecutive identical nucleotides
    
    # Secondary structure
    max_paired_fraction: float = 0.6  # Prevent rigid structures
    
    # Thermodynamic asymmetry
    min_asymmetry_score: float = 0.65  # Guide strand selection
    
    # MFE thresholds (kcal/mol)
    mfe_min: float = -8.0      # Too stable (more negative)
    mfe_max: float = -2.0      # Too unstable (less negative)
    
    # Duplex stability (kcal/mol)
    duplex_stability_min: float = -25.0
    duplex_stability_max: float = -15.0
    
    # Melting temperature (°C, for mammalian cells)
    melting_temp_min: float = 60.0
    melting_temp_max: float = 78.0
    
    # End asymmetry ΔΔG (kcal/mol)
    delta_dg_end_min: float = 2.0
    delta_dg_end_max: float = 6.0
    
    # Off-target limits
    max_off_target_count: int = 3
```

### 1.4 OffTargetFilterCriteria

Specialized filtering for off-target analysis results:

```python
class OffTargetFilterCriteria(BaseModel):
    """Off-target analysis filtering criteria."""
    
    # Transcriptome off-targets (mismatch tolerance)
    max_transcriptome_hits_0mm: int = 0    # Perfect matches
    max_transcriptome_hits_1mm: int = 5    # 1-mismatch hits
    max_transcriptome_hits_2mm: int = 20   # 2-mismatch hits
    
    # miRNA seed matches (positions 2-8)
    max_mirna_perfect_seed: int = 3
    max_mirna_1mm_seed: int = 10
    fail_on_high_risk_mirna: bool = True
```

### 1.5 ScoringWeights

Relative weights for composite scoring:

```python
class ScoringWeights(BaseModel):
    """Component weights for composite scoring (must sum to 1.0)."""
    
    asymmetry: float = 0.25      # Thermodynamic asymmetry
    gc_content: float = 0.20     # GC optimization
    accessibility: float = 0.25  # Target accessibility
    off_target: float = 0.20     # Specificity
    empirical: float = 0.10      # Position-specific rules
```

---

## 2. Scoring Algorithms

### 2.1 Composite Score Calculation

The composite score integrates multiple evidence-based components:

$$\text{Composite} = \sum_{i} w_i \times S_i \times 100$$

Where $w_i$ are configurable weights and $S_i$ are normalized component scores (0-1).

### 2.2 Thermodynamic Asymmetry Score

**Research basis**: Khvorova et al. (2003), Schwarz et al. (2003)

RISC preferentially loads the strand with the less thermodynamically stable 5' end. The asymmetry score measures this preference:

**Algorithm**:
1. Calculate 5' end stability (positions 1-7): $\Delta G_{5'}$
2. Calculate 3' end stability (positions 15-21): $\Delta G_{3'}$  
3. Compute asymmetry: $\text{raw} = \Delta G_{5'} - \Delta G_{3'}$
4. Normalize: $\text{score} = \max(0, \min(1, (\text{raw} + 5) / 10))$

**Implementation** (ViennaRNA):
```python
def calculate_asymmetry_score(candidate) -> tuple[float, float, float]:
    """Returns (dg_5p, dg_3p, asymmetry_score)"""
    dg_5p = calculate_end_stability(guide[:7], passenger[:7])
    dg_3p = calculate_end_stability(guide[14:21], passenger[14:21])
    asymmetry_raw = dg_5p - dg_3p
    asymmetry_score = max(0.0, min(1.0, (asymmetry_raw + 5.0) / 10.0))
    return dg_5p, dg_3p, asymmetry_score
```

**Interpretation**:
| Score | Interpretation |
|-------|----------------|
| 0.8-1.0 | Excellent - strong guide strand bias |
| 0.65-0.8 | Good - likely correct strand selection |
| 0.5-0.65 | Moderate - mixed strand loading possible |
| <0.5 | Poor - passenger strand may dominate |

### 2.3 GC Content Score

**Research basis**: Reynolds et al. (2004), Ui-Tei et al. (2004)

GC content affects duplex stability and target accessibility. The scoring uses a Gaussian penalty centered at optimal GC (40%):

$$\text{GC\_score} = \exp\left(-\left(\frac{\text{GC} - 40}{10}\right)^2\right)$$

**Implementation**:
```python
def _calculate_gc_score(gc_content: float) -> float:
    """Gaussian penalty around 40% GC."""
    return math.exp(-(((gc_content - 40) / 10) ** 2))
```

**Interpretation**:
| GC Range | Effect |
|----------|--------|
| <35% | Unstable duplex, poor RISC loading |
| 35-40% | Acceptable, monitor stability |
| **40-55%** | **Optimal range** |
| 55-60% | Acceptable, may reduce accessibility |
| >60% | Overly stable, poor target release |

### 2.4 Duplex Stability Score

**Research basis**: Naito et al. (2009), Ichihara et al. (2017)

Duplex formation ΔG affects RISC loading efficiency. Score normalized from ΔG range [-40, -5] kcal/mol:

$$\text{score} = \frac{-\Delta G - 5}{40 - 5}$$

**Implementation**:
```python
def _calculate_duplex_score(candidate) -> tuple[float, float]:
    """Returns (normalized_score, dg_value)"""
    dg = calculate_duplex_stability(guide, passenger)
    dg_clamped = max(-40.0, min(-5.0, dg))
    score = (-(dg_clamped) - 5.0) / (40.0 - 5.0)
    return max(0.0, min(1.0, score)), dg
```

**Optimal range**: -15 to -25 kcal/mol

### 2.5 Target Accessibility Score

**Research basis**: Tafer et al. (2008)

Target site accessibility affects siRNA efficacy. Score based on guide strand secondary structure:

$$\text{Accessibility} = 1 - \text{paired\_fraction}$$

**Implementation** (ViennaRNA):
```python
def _calculate_accessibility_score(candidate) -> float:
    """Accessibility inversely related to secondary structure."""
    structure, mfe, paired_fraction = calculate_secondary_structure(guide)
    return 1.0 - paired_fraction
```

**Optimal**: paired_fraction 0.4-0.6 (moderate structure)

### 2.6 Off-Target Score

Specificity prediction based on internal repetitive sequences:

$$\text{OT\_score} = \exp\left(-\frac{\text{penalty}}{50}\right)$$

**Implementation**:
```python
def _calculate_off_target_score(candidate) -> float:
    """Penalty for repetitive 7-mer sequences."""
    penalty = 0
    for i in range(len(guide) - 6):
        seed = guide[i:i+7]
        if guide.count(seed) > 1:
            penalty += 10
    return math.exp(-penalty / 50)
```

**`[REVIEW NEEDED]`**: Current implementation is simplified. Full off-target analysis uses BWA-MEM2 alignment against reference genomes in the Nextflow pipeline.

### 2.7 Empirical Score (Reynolds Rules)

**Research basis**: Reynolds et al. (2004)

Position-specific sequence preferences:

```python
def _calculate_empirical_score(candidate) -> float:
    """Simplified Reynolds rules."""
    score = 0.5  # Base score
    
    # Prefer A/U at position 19 (3' end)
    if guide[18] in ["A", "U"]:
        score += 0.1
    
    # Prefer G/C at position 1
    if guide[0] in ["G", "C"]:
        score += 0.1
    
    # Avoid C at position 19
    if guide[18] == "C":
        score -= 0.1
    
    return max(0.0, min(1.0, score))
```

**`[REVIEW NEEDED]`**: Additional Reynolds criteria could be implemented:
- Position 10 preferences
- A/U content in positions 15-19
- Avoid GGG stretches

---

## 3. Filter Implementation

### 3.1 Early Filtering (Enumeration Stage)

During candidate enumeration, fast filters are applied to reduce computational load:

```python
def _enumerate_candidates(sequence, transcript_id):
    for i in range(len(sequence) - sirna_length + 1):
        target_seq = sequence[i:i+sirna_length]
        guide_seq = reverse_complement(target_seq)
        gc_content = calculate_gc_content(guide_seq)
        
        # Fast rejection
        fail_reason = None
        if not (gc_min <= gc_content <= gc_max):
            fail_reason = FilterStatus.GC_OUT_OF_RANGE
        elif has_poly_runs(guide_seq, max_poly_runs):
            fail_reason = FilterStatus.POLY_RUNS
        
        if fail_reason:
            # Store in rejected pool for "dirty control" analysis
            rejected.append(candidate)
        else:
            candidates.append(candidate)
```

### 3.2 Quality Filters (Scoring Stage)

Additional filters applied during scoring:

| Filter | Condition | Rationale |
|--------|-----------|-----------|
| `EXCESS_PAIRING` | paired_fraction > 0.6 | Prevents rigid structures |
| `LOW_ASYMMETRY` | asymmetry_score < min_asymmetry_score | Ensures guide strand selection |

### 3.3 Filter Status Codes

```python
class FilterStatus(str, Enum):
    PASS = "PASS"                    # All criteria met
    GC_OUT_OF_RANGE = "GC_OUT_OF_RANGE"  # GC content outside bounds
    POLY_RUNS = "POLY_RUNS"          # Homopolymer runs exceed limit
    EXCESS_PAIRING = "EXCESS_PAIRING"    # Too much secondary structure
    LOW_ASYMMETRY = "LOW_ASYMMETRY"  # Poor thermodynamic asymmetry
    DIRTY_CONTROL = "DIRTY_CONTROL"  # Reserved for controls
```

---

## 4. Threshold Justification

### 4.1 GC Content: 35-60%

**Literature support**:
- Reynolds et al. (2004): Optimal 30-52% for maximum silencing
- Ui-Tei et al. (2004): Functional siRNAs have 35-65% GC
- Jackson et al. (2006): Higher GC correlates with off-targets

**Rationale**: Balance between:
- **Lower bound (35%)**: Minimum duplex stability for RISC loading
- **Upper bound (60%)**: Maximum to prevent over-stabilization and off-targeting

### 4.2 Asymmetry Score: ≥0.65

**Literature support**:
- Khvorova et al. (2003): Thermodynamic asymmetry determines strand selection
- Schwarz et al. (2003): ΔΔG of 2+ kcal/mol ensures correct loading

**Rationale**: Score of 0.65 corresponds to approximately ΔΔG = 1.5 kcal/mol, providing >80% probability of correct strand selection.

### 4.3 Poly-runs: ≤3 consecutive

**Literature support**:
- Jackson et al. (2003): AAAA runs associated with off-targets
- Synthesis considerations: Long homopolymers cause synthesis issues

**Rationale**: Practical limit balancing efficacy and manufacturability.

### 4.4 MFE: -2 to -8 kcal/mol

**Literature support**:
- Tafer et al. (2008): Moderate structure optimal for target binding
- Too stable (<-10): Impaired target access
- Too unstable (>0): Poor duplex integrity

### 4.5 Melting Temperature: 60-78°C

**Literature support**:
- Standard for mammalian cell culture at 37°C
- Allows duplex stability while permitting RISC-mediated unwinding

**`[REVIEW NEEDED]`**: Temperature thresholds may need adjustment for:
- Plant cells (different optimal ranges)
- In vivo applications (serum stability requirements)

---

## 5. miRNA-Biogenesis Mode

### 5.1 MiRNADesignConfig

Specialized parameters for miRNA-like siRNA design:

```python
class MiRNADesignConfig(BaseModel):
    """miRNA-biogenesis-aware configuration."""
    
    # Conservative thresholds
    gc_min: float = 30.0       # Relaxed lower bound
    gc_max: float = 52.0       # Stricter upper bound
    asymmetry_min: float = 0.65
    
    # Argonaute loading preferences
    scoring_weights: dict = {
        "ago_start_bonus": 0.1,      # A/U at position 1
        "pos1_mismatch_bonus": 0.05, # G:U wobble preferred
        "seed_clean_bonus": 0.15,    # Clean seed region
        "supp_13_16_bonus": 0.1,     # 3' supplementary pairing
    }
```

### 5.2 miRNA-Specific Scoring

**Position 1 analysis**:
- Argonaute preferentially loads strands with A/U at position 1
- G:U wobble or mismatch at position 1 improves loading

**Seed region (positions 2-8)**:
- Critical for target recognition
- Clean seed = lower off-target potential

**3' Supplementary pairing (positions 13-16)**:
- Contributes to target specificity
- Lower stability preferred (more specific)

**`[REVIEW NEEDED]`**: miRNA-specific scoring weights are based on general principles. Experimental validation recommended for therapeutic applications.

---

## 6. Off-Target Analysis Models

### 6.1 OffTargetHit

```python
class OffTargetHit(BaseModel):
    """Single off-target alignment from BWA analysis."""
    
    qname: str           # siRNA identifier
    qseq: str            # Query sequence
    rname: str           # Reference (chromosome/transcript)
    coord: int           # Alignment position
    strand: str          # + or -
    nm: int              # Edit distance
    seed_mismatches: int # Mismatches in seed (pos 2-8)
    offtarget_score: float
```

### 6.2 MiRNAHit

```python
class MiRNAHit(BaseModel):
    """miRNA seed match from database alignment."""
    
    species: str         # e.g., "hsa" (human)
    database: str        # mirgenedb, mirbase, etc.
    mirna_id: str        # e.g., hsa-miR-21-5p
    seed_mismatches: int # Seed region mismatches
```

### 6.3 Supported miRNA Databases

| Database | Description |
|----------|-------------|
| `mirgenedb` | High-confidence, manually curated |
| `mirbase` | Comprehensive, all mature miRNAs |
| `mirbase_high_conf` | miRBase high-confidence subset |
| `targetscan` | miRNA family conservation data |

---

## 7. Chemical Modification Models

### 7.1 StrandMetadata

```python
class StrandMetadata(BaseModel):
    """Complete strand metadata with modifications."""
    
    id: str
    sequence: str
    overhang: str          # e.g., "dTdT", "UU"
    chem_mods: list[ChemicalModification]
    provenance: Provenance
```

### 7.2 ChemicalModification

```python
class ChemicalModification(BaseModel):
    """Position-specific chemical modification."""
    
    type: str              # 2OMe, 2F, PS, LNA
    positions: list[int]   # 1-based positions
```

### 7.3 Supported Modification Patterns

| Pattern | Description |
|---------|-------------|
| `standard_2ome` | 2'-O-methyl at alternating positions |
| `minimal_terminal` | Terminal modifications only |
| `maximal_stability` | Full backbone modifications |
| `none` | No modifications |

---

## 8. Workflows Requiring Documentation

**`[DOCUMENTATION NEEDED]`**: The following workflows exist but require detailed documentation:

### 8.1 Nextflow Pipeline
- Multi-genome off-target analysis
- BWA-MEM2 alignment parameters
- Species-specific reference handling

### 8.2 ORF Validation
- Start/stop codon detection
- Frame shift analysis
- Kozak sequence scoring

### 8.3 Transcript Retrieval
- Ensembl/RefSeq/GENCODE integration
- Isoform selection criteria
- Sequence validation

---

## 9. References

1. **Khvorova A, Reynolds A, Jayasena SD** (2003). Functional siRNAs and miRNAs exhibit strand bias. *Cell* 115(2):209-216.

2. **Schwarz DS, Hutvágner G, Du T, Xu Z, Aronin N, Bhatt DP** (2003). Asymmetry in the assembly of the RNAi enzyme complex. *Cell* 115(2):199-208.

3. **Reynolds A, Leake D, Boese Q, Scaringe S, Marshall WS, Khvorova A** (2004). Rational siRNA design for RNA interference. *Nature Biotechnology* 22(3):326-330.

4. **Ui-Tei K, Naito Y, Takahashi F, Haraguchi T, Ohki-Hamazaki H, Juni A, Ueda R, Saigo K** (2004). Guidelines for the selection of highly effective siRNA sequences for mammalian and chick RNA interference. *Nucleic Acids Research* 32(3):936-948.

5. **Naito Y, Yoshimura J, Morishita S, Ui-Tei K** (2009). siDirect 2.0: updated software for designing functional siRNA with reduced seed-dependent off-target effect. *BMC Bioinformatics* 10:392.

6. **Ichihara M, Murakumo Y, Masuda A, Matsuura T, Asai N, Jijiwa M, Ishida M, Shinmi J, Yatsuya H, Qiao S, Takahashi M, Ohno K** (2007). Thermodynamic instability of siRNA duplex is a prerequisite for dependable prediction of siRNA activities. *Nucleic Acids Research* 35(18):e123.

7. **Tafer H, Ameres SL, Obernosterer G, Gebeshuber CA, Schroeder R, Martinez J, Hofacker IL** (2008). The impact of target site accessibility on the design of effective siRNAs. *Nature Biotechnology* 26(5):578-583.

8. **Jackson AL, Bartz SR, Schelter J, Kobayashi SV, Burchard J, Mao M, Li B, Cavet G, Linsley PS** (2003). Expression profiling reveals off-target gene regulation by RNAi. *Nature Biotechnology* 21(6):635-637.

---

## Appendix A: Default Parameter Summary

| Parameter | Default | Range | Justification |
|-----------|---------|-------|---------------|
| `sirna_length` | 21 | 19-23 | Standard duplex length |
| `gc_min` | 35.0 | 0-100 | Minimum stability |
| `gc_max` | 60.0 | 0-100 | Maximum stability |
| `max_poly_runs` | 3 | 1+ | Synthesis/specificity |
| `max_paired_fraction` | 0.6 | 0-1 | Accessibility |
| `min_asymmetry_score` | 0.65 | 0.3-1 | Strand selection |
| `mfe_min` | -8.0 | kcal/mol | Structure stability |
| `mfe_max` | -2.0 | kcal/mol | Structure stability |
| `duplex_stability_min` | -25.0 | kcal/mol | Duplex formation |
| `duplex_stability_max` | -15.0 | kcal/mol | Duplex formation |
| `melting_temp_min` | 60.0 | °C | Mammalian cells |
| `melting_temp_max` | 78.0 | °C | Mammalian cells |
| `max_off_target_count` | 3 | 0+ | Specificity |

## Appendix B: Scoring Weight Defaults

| Component | Weight | Rationale |
|-----------|--------|-----------|
| Asymmetry | 0.25 | Most predictive single factor |
| GC Content | 0.20 | Stability/accessibility balance |
| Accessibility | 0.25 | Target site availability |
| Off-target | 0.20 | Specificity importance |
| Empirical | 0.10 | Position-specific fine-tuning |

---

*Document version: 1.0*
*Last updated: Auto-generated from source code*
*Review status: Initial draft - Expert review recommended*
