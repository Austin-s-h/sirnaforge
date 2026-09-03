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
    duplex_stability: float    # ΔG in kcal/mol (fully paired 21mer: -32 to -43)

    # Secondary structure
    structure: str             # Dot-bracket notation
    mfe: float                 # Minimum free energy (optimal: -2 to -8 kcal/mol)
    paired_fraction: float     # Fraction paired bases (optimal: 0.4-0.6)

    # Off-target metrics
    off_target_screened: bool  # False = never screened, so the counts below are unknown
    off_target_count: int      # Genuine off-target sites only (goal: ≤3); on-target,
                                # ortholog and repeat-mediated hits are excluded
    transcriptome_hits_0mm: int   # Perfect match GENUINE off-target hits
    transcriptome_hits_1mm: int   # 1-mismatch GENUINE off-target hits
    transcriptome_hits_2mm: int   # 2-mismatch GENUINE off-target hits
    mirna_hits_total: int         # Total miRNA seed matches
    mirna_hits_0mm_seed: int      # Perfect seed matches

    # Hit classification (four-way: on-target, ortholog, repeat, off-target)
    on_target_hits: int        # Hits on the query gene in the query species
    ortholog_hits: int         # Hits on the query gene's ortholog in another screened species
    repeat_hits: int           # Hits attributable to a flagged repeat element
    ortholog_species: str      # Comma-separated canonical species with an ortholog hit

    # Repeat detection (design-time k-mer frequency check)
    repeat_flagged: bool             # Guide exceeds the repeat transcript-fraction threshold
    repeat_transcript_fraction: float  # Fraction of reference transcripts containing the guide

    # Post-screen sub-scores (inactive/None until screening has run, or when there is no
    # evidence to compute them from -- see ScoringWeights below)
    isoform_coverage: float | None     # Protein-coding isoform coverage sub-score
    conservation_score: float | None   # Cross-species ortholog conservation sub-score

    # Scoring
    composite_score: float     # Overall quality (0-100 scale), computed once, after screening
    component_scores: dict     # Individual scoring components (design-time diagnostics)
    score_asymmetry: float | None       # Per-term composite contributions
    score_gc_content: float | None
    score_accessibility: float | None
    score_empirical: float | None
    score_off_target: float | None
    score_isoform_coverage: float | None
    score_conservation: float | None
    scored_after_screening: bool  # True once off_target/isoform_coverage/conservation are active
    weight_set_version: str       # Weight set that produced composite_score ("" = not yet scored)
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

    # Thermodynamic asymmetry -> gates asymmetry_score
    min_asymmetry_score: float = 0.65  # Guide strand selection

    # Empirical design rules -> gates the empirical component score (range 0.4-0.7)
    min_empirical_score: float = 0.5
```

Each threshold gates the quantity it is named after: `min_asymmetry_score` is compared
against `asymmetry_score`, `min_empirical_score` against the empirical component score.
`min_empirical_score` is bounded by the empirical rule's attainable range (0.4-0.7), so a
value the rule can never reach is rejected at construction instead of silently failing
every candidate.

Issue #80 removed eight thermodynamic windows that used to live here (`mfe_min`, `mfe_max`,
`duplex_stability_min`, `duplex_stability_max`, `melting_temp_min`, `melting_temp_max`,
`delta_dg_end_min`, `delta_dg_end_max`): they were declared but never enforced by
`SiRNADesigner` (no CLI flag exposed any of them), and re-deriving correct windows against
truth data is deliberately out of scope. `max_off_target_count` also used to live here; it
moved to `OffTargetFilterCriteria` below, since design time cannot know a genuine off-target
count that only exists after screening.

### 1.4 OffTargetFilterCriteria

Specialized filtering for off-target analysis results, applied after screening:

```python
class OffTargetFilterCriteria(BaseModel):
    """Off-target analysis filtering criteria."""

    # Genuine off-target count (on-target, ortholog and repeat hits excluded)
    max_off_target_count: int = 3

    # Transcriptome GENUINE off-targets (mismatch tolerance)
    max_transcriptome_hits_0mm: int = 1    # Perfect matches
    max_transcriptome_hits_1mm: int = 10   # 1-mismatch hits
    max_transcriptome_hits_2mm: int = 50   # 2-mismatch hits

    # miRNA seed matches (positions 2-8)
    max_mirna_perfect_seed: int = 0
    max_mirna_1mm_seed: int = 10
    fail_on_high_risk_mirna: bool = True
```

### 1.5 ScoringWeights

Relative weights for composite scoring (seven terms as of issue #80; must sum to 1.0):

```python
class ScoringWeights(BaseModel):
    """Component weights for composite scoring (must sum to 1.0)."""

    asymmetry: float = 0.12         # Thermodynamic asymmetry
    gc_content: float = 0.10        # GC optimization
    accessibility: float = 0.13     # Target accessibility
    empirical: float = 0.15         # Position-specific rules
    off_target: float = 0.25        # Post-screen genuine off-target specificity (see 2.6)
    isoform_coverage: float = 0.15  # Post-screen protein-coding isoform coverage (new)
    conservation: float = 0.10      # Post-screen cross-species ortholog conservation (new)
```

Validation is a `@model_validator(mode="after")` that sums all seven fields by name
(`COMPOSITE_TERM_NAMES`), not a `field_validator` that only sees fields declared earlier in the
class -- the previous positional check silently stopped seeing new weights as the term set grew.

`off_target`, `isoform_coverage` and `conservation` cannot be evaluated until off-target
screening has run, so they are inactive at design time. `compute_composite` (see 2.1) renormalises
the weights of whichever terms _are_ active, so a design-time-only score and a post-screen score
are both legitimate 0-100 values from the same weight set, not two incompatible scales.

---

## 2. Scoring Algorithms

### 2.1 Composite Score Calculation

The composite score integrates seven evidence-based components, computed once by
`compute_composite` over whichever terms are _active_ for that candidate:

$$\text{Composite} = \sum_{i \in \text{active}} w_i' \times S_i \times 100, \quad w_i' = \frac{w_i}{\sum_{j \in \text{active}} w_j}$$

Where $w_i$ are the configured `ScoringWeights`, $w_i'$ is the weight renormalised over the
active term set, and $S_i$ are normalized component scores (0-1). A term is active when its
sub-score could be computed for that candidate (see 1.5); inactive terms are omitted, not
scored zero, so $\sum_i w_i'$ over the active set is always 1.

### 2.2 Thermodynamic Asymmetry Score

**Research basis**: Khvorova et al. (2003), Schwarz et al. (2003)

RISC preferentially loads the strand with the less thermodynamically stable 5' end. The asymmetry score measures this preference:

**Algorithm**:

1. Calculate 5' end stability (guide 5' 7-mer against the passenger 3' end): $\Delta G_{5'}$
2. Calculate 3' end stability (guide 3' 7-mer against the passenger 5' end): $\Delta G_{3'}$
3. Compute asymmetry: $\text{raw} = \Delta G_{5'} - \Delta G_{3'}$
4. Normalize: $\text{score} = \max(0, \min(1, (\text{raw} + 5) / 10))$

The duplex is antiparallel, so each end window pairs with the _opposite_ end of the other
strand. Both windows use the same width, so the two ΔG values stay comparable.

**Implementation** (ViennaRNA):

```python
def calculate_asymmetry_score(candidate) -> tuple[float, float, float]:
    """Returns (dg_5p, dg_3p, asymmetry_score)"""
    window = min(END_WINDOW_NT, len(guide), len(passenger))
    dg_5p = calculate_end_stability(guide[:window], passenger[-window:])
    dg_3p = calculate_end_stability(guide[-window:], passenger[:window])
    asymmetry_raw = dg_5p - dg_3p
    asymmetry_score = max(0.0, min(1.0, (asymmetry_raw + 5.0) / 10.0))
    return dg_5p, dg_3p, asymmetry_score
```

`passenger` is already the guide's complement, so neither strand is reverse-complemented
before folding — ViennaRNA's `&` cofold notation pairs two 5'->3' strands antiparallel by
itself. Reverse-complementing the passenger folded the guide against itself, which pinned
`asymmetry_score` to exactly 0.5 for every 21 nt candidate (fixed in 0.5.2).

**Interpretation**:

| Score    | Interpretation                           |
| -------- | ---------------------------------------- |
| 0.8-1.0  | Excellent - strong guide strand bias     |
| 0.65-0.8 | Good - likely correct strand selection   |
| 0.5-0.65 | Moderate - mixed strand loading possible |
| <0.5     | Poor - passenger strand may dominate     |

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

| GC Range   | Effect                               |
| ---------- | ------------------------------------ |
| <35%       | Unstable duplex, poor RISC loading   |
| 35-40%     | Acceptable, monitor stability        |
| **40-55%** | **Optimal range**                    |
| 55-60%     | Acceptable, may reduce accessibility |
| >60%       | Overly stable, poor target release   |

### 2.4 Duplex Stability Score

**Research basis**: Naito et al. (2009), Ichihara et al. (2017)

Duplex formation ΔG affects RISC loading efficiency. ΔG scales with duplex length, so it
is normalised per nucleotide before scoring: -2.1 kcal/mol/nt maps to 1.0, -1.4 to 0.0.

$$\text{score} = \frac{-1.4 - \Delta G / L}{-1.4 - (-2.1)}$$

**Implementation**:

```python
def _calculate_duplex_score(candidate) -> tuple[float, float]:
    """Returns (normalized_score, dg_value)"""
    dg = calculate_duplex_stability(guide, passenger)
    dg_per_nt = dg / len(candidate.guide_sequence)
    span = DUPLEX_DG_PER_NT_WEAK - DUPLEX_DG_PER_NT_STRONG
    score = (DUPLEX_DG_PER_NT_WEAK - dg_per_nt) / span
    return max(0.0, min(1.0, score)), dg
```

**Measured range** (ViennaRNA, 37 °C): a fully paired 21mer duplex is -32 to -43 kcal/mol,
i.e. -1.55 to -2.06 kcal/mol/nt. The fixed [-40, -5] kcal/mol window used before 0.5.2 was
calibrated against the pre-fix self-fold ΔG; against real duplex ΔG it put 41% of
candidates at exactly 1.0, and it rewarded longer designs for their length alone.

`duplex_stability_score` is reported but does not feed `composite_score`.

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

As of issue #80, the `off_target` term contributing to `composite_score` is the **post-screen
genuine-off-target specificity** sub-score, computed by
`sirnaforge.core.scoring.off_target_sub_score` from the redefined `off_target_count` (on-target,
ortholog and repeat-mediated hits excluded):

$$\text{off\_target} = \exp\left(-\frac{\text{genuine\_off\_target\_count}}{10}\right)$$

Zero genuine off-targets scores 1.0; the score decays toward 0 as the count grows. This term is
inactive until screening has run (see 1.5), so it is absent from a design-time-only score.

**Design-time diagnostic (not part of `composite_score`)**: `SiRNADesigner` still computes a
self-repetitiveness proxy over internal repeated 7-mers, based on internal repetitive sequences:

$$\text{design\_off\_target\_proxy} = \exp\left(-\frac{\text{penalty}}{50}\right)$$

```python
def _calculate_off_target_score(candidate) -> float:
    """Diagnostic only: penalty for repetitive 7-mer sequences within the guide."""
    penalty = 0
    for i in range(len(guide) - 6):
        seed = guide[i:i+7]
        if guide.count(seed) > 1:
            penalty += 10
    return math.exp(-penalty / 50)
```

This value is stored in `component_scores["design_off_target_proxy"]` and reported for
diagnostic purposes only; it does not feed `composite_score` under any name. Full off-target
analysis uses BWA-MEM2 alignment against reference transcriptomes in the Nextflow pipeline,
followed by the four-way hit classifier (`sirnaforge.core.hit_classification`); that
decomposition, not the internal-repeat proxy, is what `off_target_count` and the `off_target`
scoring term are based on.

### 2.7 Empirical Score (Reynolds Rules)

**Research basis**: Reynolds et al. (2004)

Position-specific sequence preferences:

```python
def _calculate_empirical_score(candidate) -> float:
    """Simplified Reynolds rules."""
    guide = candidate.guide_sequence.upper().replace("T", "U")  # guides are stored as DNA
    score = 0.5  # Base score

    # Prefer A/U at position 19 (3' end)
    if guide[18] in ("A", "U"):
        score += 0.1

    # Prefer G/C at position 1
    if guide[0] in ("G", "C"):
        score += 0.1

    # Avoid C at position 19
    if guide[18] == "C":
        score -= 0.1

    return max(EMPIRICAL_SCORE_MIN, min(EMPIRICAL_SCORE_MAX, score))
```

The attainable range is **0.4-0.7**, not 0-1: three ±0.1 adjustments on a 0.5 base, and the
A/U and C tests at position 19 are mutually exclusive. `min_empirical_score` is bounded by
that range. Reading the guide as RNA matters — guides are stored as DNA, so before 0.5.2 a
T at position 19 never earned the A/U bonus.

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

| Filter                | Condition                             | Rationale                              |
| --------------------- | ------------------------------------- | -------------------------------------- |
| `EXCESS_PAIRING`      | paired_fraction > 0.6                 | Prevents rigid structures              |
| `LOW_ASYMMETRY`       | asymmetry_score < min_asymmetry_score | Ensures guide strand selection         |
| `LOW_EMPIRICAL_SCORE` | empirical score < min_empirical_score | Position-specific sequence preferences |

Only the first failure is recorded. Before 0.5.2, `LOW_ASYMMETRY` was assigned by comparing
the _empirical_ score against `min_asymmetry_score`: it reported the wrong reason, and
`asymmetry_score` was never gated at all.

Between design-time filtering and post-screen filtering, issue #80 added one more design-time
gate, applied to every distinct guide against the query species' cDNA reference:

| Filter           | Condition                                       | Rationale                                 |
| ---------------- | ----------------------------------------------- | ----------------------------------------- |
| `REPEAT_ELEMENT` | guide occurs in > 0.1% of reference transcripts | Flags guides overlapping a repeat element |

`REPEAT_ELEMENT` is applied only if the candidate is currently passing (existing failures are not
overwritten), and it excludes the candidate from ranking.

### 3.3 Post-Screen Off-Target Filters

Applied once transcriptome/miRNA screening has run, against `OffTargetFilterCriteria`
(see 1.4). These now gate the **redefined, genuine-off-target-only** counts:

| Filter                        | Condition                                                        | Rationale                          |
| ----------------------------- | ---------------------------------------------------------------- | ---------------------------------- |
| `TRANSCRIPTOME_PERFECT_MATCH` | 0-mismatch genuine off-targets > `max_transcriptome_hits_0mm`    | Perfect-match specificity          |
| `TRANSCRIPTOME_1MM`           | 1-mismatch genuine off-targets > `max_transcriptome_hits_1mm`    | Near-miss specificity              |
| `TRANSCRIPTOME_2MM`           | 2-mismatch genuine off-targets > `max_transcriptome_hits_2mm`    | Broader specificity                |
| `MIRNA_PERFECT_SEED`          | perfect miRNA seed hits > `max_mirna_perfect_seed`               | miRNA-mimicry risk                 |
| `HIGH_RISK_MIRNA`             | perfect seed + `offtarget_score` < 5.0                           | Strong-binding miRNA mimicry       |
| `TOTAL_OFFTARGETS`            | combined transcriptome + miRNA hits > `max_total_offtarget_hits` | Aggregate specificity              |
| `EXCESS_OFF_TARGETS`          | genuine off-target count > `max_off_target_count`                | Overall genuine-off-target ceiling |

These six verdicts previously existed only as raw strings assembled in `workflow.py`; they are
now members of `FilterStatus` (see 3.4), and the Pandera allow-list for the `passes_filters`
column derives from the enum instead of maintaining an independent copy that could drift from it.

### 3.4 Filter Status Codes

```python
class FilterStatus(str, Enum):
    PASS = "PASS"                    # All criteria met
    GC_OUT_OF_RANGE = "GC_OUT_OF_RANGE"  # GC content outside bounds
    POLY_RUNS = "POLY_RUNS"          # Homopolymer runs exceed limit
    EXCESS_PAIRING = "EXCESS_PAIRING"    # Too much secondary structure
    LOW_ASYMMETRY = "LOW_ASYMMETRY"  # Poor thermodynamic asymmetry
    LOW_EMPIRICAL_SCORE = "LOW_EMPIRICAL_SCORE"  # Fails the empirical design rules
    DIRTY_CONTROL = "DIRTY_CONTROL"  # Reserved for controls
    REPEAT_ELEMENT = "REPEAT_ELEMENT"  # Design-time repeat k-mer verdict
    EXCESS_OFF_TARGETS = "EXCESS_OFF_TARGETS"  # Post-screen genuine off-target ceiling
    TRANSCRIPTOME_PERFECT_MATCH = "TRANSCRIPTOME_PERFECT_MATCH"
    TRANSCRIPTOME_1MM = "TRANSCRIPTOME_1MM"
    TRANSCRIPTOME_2MM = "TRANSCRIPTOME_2MM"
    MIRNA_PERFECT_SEED = "MIRNA_PERFECT_SEED"
    HIGH_RISK_MIRNA = "HIGH_RISK_MIRNA"
    TOTAL_OFFTARGETS = "TOTAL_OFFTARGETS"
```

### 3.5 Hit Classes (`sirnaforge.core.hit_classification.HitClass`)

Every transcriptome hit is classified into exactly one of four mutually exclusive classes,
checked in this precedence order (on-target and ortholog deliberately outrank repeat, so a
repeat-flagged guide's hits on its own gene or orthologs are still counted as such):

| Class        | Meaning                                                                                                                                         |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `ON_TARGET`  | Hit is in the query species and matches the query gene (transcript ID, gene ID or symbol)                                                       |
| `ORTHOLOG`   | Hit is in a different screened species, and that species' gene symbol for the hit transcript matches the query gene's symbol (case-insensitive) |
| `REPEAT`     | The guide is in the design-time repeat-flagged set                                                                                              |
| `OFF_TARGET` | Everything else -- this is what `off_target_count` counts                                                                                       |

A hit whose species label is blank/missing is treated as the query species. When an ortholog
check cannot run because the hit species' annotation carries no gene symbol for that transcript,
the hit stays `OFF_TARGET` (never guessed at) and the shortfall is counted separately
(`ortholog_symbol_lookup_misses` in `offtarget_summary`).

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

| Database            | Description                       |
| ------------------- | --------------------------------- |
| `mirgenedb`         | High-confidence, manually curated |
| `mirbase`           | Comprehensive, all mature miRNAs  |
| `mirbase_high_conf` | miRBase high-confidence subset    |
| `targetscan`        | miRNA family conservation data    |

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

| Pattern             | Description                          |
| ------------------- | ------------------------------------ |
| `standard_2ome`     | 2'-O-methyl at alternating positions |
| `minimal_terminal`  | Terminal modifications only          |
| `maximal_stability` | Full backbone modifications          |
| `none`              | No modifications                     |

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

1. **Khvorova A, Reynolds A, Jayasena SD** (2003). Functional siRNAs and miRNAs exhibit strand bias. _Cell_ 115(2):209-216.

2. **Schwarz DS, Hutvágner G, Du T, Xu Z, Aronin N, Bhatt DP** (2003). Asymmetry in the assembly of the RNAi enzyme complex. _Cell_ 115(2):199-208.

3. **Reynolds A, Leake D, Boese Q, Scaringe S, Marshall WS, Khvorova A** (2004). Rational siRNA design for RNA interference. _Nature Biotechnology_ 22(3):326-330.

4. **Ui-Tei K, Naito Y, Takahashi F, Haraguchi T, Ohki-Hamazaki H, Juni A, Ueda R, Saigo K** (2004). Guidelines for the selection of highly effective siRNA sequences for mammalian and chick RNA interference. _Nucleic Acids Research_ 32(3):936-948.

5. **Naito Y, Yoshimura J, Morishita S, Ui-Tei K** (2009). siDirect 2.0: updated software for designing functional siRNA with reduced seed-dependent off-target effect. _BMC Bioinformatics_ 10:392.

6. **Ichihara M, Murakumo Y, Masuda A, Matsuura T, Asai N, Jijiwa M, Ishida M, Shinmi J, Yatsuya H, Qiao S, Takahashi M, Ohno K** (2007). Thermodynamic instability of siRNA duplex is a prerequisite for dependable prediction of siRNA activities. _Nucleic Acids Research_ 35(18):e123.

7. **Tafer H, Ameres SL, Obernosterer G, Gebeshuber CA, Schroeder R, Martinez J, Hofacker IL** (2008). The impact of target site accessibility on the design of effective siRNAs. _Nature Biotechnology_ 26(5):578-583.

8. **Jackson AL, Bartz SR, Schelter J, Kobayashi SV, Burchard J, Mao M, Li B, Cavet G, Linsley PS** (2003). Expression profiling reveals off-target gene regulation by RNAi. _Nature Biotechnology_ 21(6):635-637.

---

## Appendix A: Default Parameter Summary

| Parameter             | Default | Range   | Justification          |
| --------------------- | ------- | ------- | ---------------------- |
| `sirna_length`        | 21      | 19-23   | Standard duplex length |
| `gc_min`              | 35.0    | 0-100   | Minimum stability      |
| `gc_max`              | 60.0    | 0-100   | Maximum stability      |
| `max_poly_runs`       | 3       | 1+      | Synthesis/specificity  |
| `max_paired_fraction` | 0.6     | 0-1     | Accessibility          |
| `min_asymmetry_score` | 0.65    | 0.3-1   | Strand selection       |
| `min_empirical_score` | 0.5     | 0.4-0.7 | Position preferences   |

`mfe_min`, `mfe_max`, `duplex_stability_min`, `duplex_stability_max`, `melting_temp_min`,
`melting_temp_max`, `delta_dg_end_min` and `delta_dg_end_max` were removed from `FilterCriteria`
in issue #80: they were declared but never enforced by `SiRNADesigner`, and re-deriving correct
windows against truth data is deliberately out of scope. `max_off_target_count` (default `3`,
range `0+`, justification: specificity) moved to `OffTargetFilterCriteria`, where it is enforced
post-screen against the genuine off-target count (see 1.4 and 3.3).

## Appendix B: Scoring Weight Defaults

| Component        | Weight | Rationale                                           |
| ---------------- | ------ | --------------------------------------------------- |
| Asymmetry        | 0.12   | Most predictive single factor                       |
| GC Content       | 0.10   | Stability/accessibility balance                     |
| Accessibility    | 0.13   | Target site availability                            |
| Empirical        | 0.15   | Position-specific fine-tuning                       |
| Off-target       | 0.25   | Post-screen genuine off-target specificity          |
| Isoform coverage | 0.15   | Protein-coding isoform targeting completeness (new) |
| Conservation     | 0.10   | Cross-species ortholog specificity check (new)      |

Weight-set version `2.0.0`; `1.x` denotes the pre-issue-#80 five-term set and is not comparable.

---

_Document version: 1.0_
_Last updated: Auto-generated from source code_
_Review status: Initial draft - Expert review recommended_
