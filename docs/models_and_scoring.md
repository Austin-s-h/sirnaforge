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
    off_target_screened: bool  # False = the screen was INCOMPLETE (never screened, or a run whose
                                # query-species alignment produced nothing), so the counts below
                                # are a lower bound, not a total — hits found in the species that
                                # did align are still reported, and still applied as filters
    off_target_count: int      # Genuine off-target sites only (gated at 15); on-target,
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

    # Target-site accessibility (RNAplfold on the transcript). None when it could not be
    # computed, in which case the term is inactive -- never substituted with a default.
    target_accessibility_p: float | None        # P(seed-end 8-mer unpaired); the scored input
    target_accessibility_p_17mer: float | None  # Reported only, never scored
    target_accessibility_p_site: float | None   # Reported only, never scored

    # Scoring
    composite_score: float     # Overall quality (0-100 scale), computed once, after screening
    component_scores: dict     # Individual scoring components (design-time diagnostics)
    score_asymmetry: float | None       # Per-term composite contributions
    score_gc_content: float | None
    score_target_accessibility: float | None
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
    top_n: int | None = None   # Candidates to report (None = all, the default)

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

    # Empirical design rules -> gates the empirical component score (range 0.4-0.6)
    min_empirical_score: float = 0.4
```

Each threshold gates the quantity it is named after: `min_asymmetry_score` is compared
against `asymmetry_score`, `min_empirical_score` against the empirical component score.
`min_empirical_score` is bounded by the empirical rule's attainable range (0.4-0.6), so a
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

    # Genuine off-target count (on-target, ortholog and repeat hits excluded).
    # None means "no gate", which is what --filter-action <id>=off sets.
    max_off_target_count: int | None = 15

    # Transcriptome GENUINE off-targets (mismatch tolerance)
    max_transcriptome_hits_0mm: int | None = 1    # Perfect matches
    max_transcriptome_hits_1mm: int | None = 10   # 1-mismatch hits
    max_transcriptome_hits_2mm: int | None = 50   # 2-mismatch hits
    max_transcriptome_seed_perfect: int | None = None  # off by default; see the warning below

    # miRNA seed matches (positions 2-8)
    max_mirna_perfect_seed: int | None = 0
    max_mirna_1mm_seed: int | None = 10   # read by no gate in 0.7.1; resolves to off
    fail_on_high_risk_mirna: bool = True
```

> **A seed-perfect partial hit falls in no mismatch stratum.** The three
> `max_transcriptome_hits_{0,1,2}mm` thresholds are read against counters stratified by the
> **guide-level** `nm` (see 6.1), so a clipped or gapped hit is stratified by how many guide bases
> failed to pair, not by the aligner's `NM` tag. A minus-strand `6S15M` / `MD:Z:15` / `NM:i:0`
> record puts the clip on guide positions 16-21 and leaves guide positions 2-8 pairing perfectly:
> it carries `nm = 6`, so it is counted in `transcriptome_hits_seed_0mm`,
> `transcriptome_hits_total` and `off_target_count`, but in **none** of `_0mm`/`_1mm`/`_2mm`. Two
> such hits report `0mm=0 1mm=0 2mm=0 seed_0mm=2 total=2 off_target_count=2`. With stock defaults
> the only thing gating them is `max_off_target_count` (15), so a run of them can pass. Set
> `max_transcriptome_seed_perfect` to gate them directly — it is enforced (verdict
> `TRANSCRIPTOME_SEED_PERFECT`) but ships as `None` because no ceiling has been calibrated against
> truth data. Unlike the three mismatch thresholds it is **not** species-split: it is compared
> against the reported `transcriptome_hits_seed_0mm` column across all screened species.

### 1.5 ScoringWeights and the named weight vectors

`ScoringWeights` is a **container of three hand-authored vectors**, not a flat weight list. Each
vector subclasses `WeightVector`, declares its own `VECTOR_NAME` and `TERM_NAMES`, and refuses to
construct unless it is named and already sums to 1.0 (tolerance 1e-9 — enough for float
representation of two-decimal literals, not enough to be approximately normalised):

```python
class DesignWeights(WeightVector):            # design_v4 -> SiRNACandidate.design_score
    target_accessibility: float = 0.35
    asymmetry: float = 0.40
    gc_content: float = 0.25

class PostScreenSiRNAWeights(WeightVector):   # postscreen_sirna_v4 -> composite_score
    off_target: float = 0.25
    target_accessibility: float = 0.30
    asymmetry: float = 0.25
    gc_content: float = 0.20

class PostScreenMiRNAWeights(WeightVector):   # postscreen_mirna_v4 -> composite_score
    off_target: float = 0.20      # exactly 0.80 x postscreen_sirna_v4, term by term
    target_accessibility: float = 0.24
    asymmetry: float = 0.20
    gc_content: float = 0.16
    ago_start: float = 0.10       # A/U at guide position 1
    supp_13_16: float = 0.10      # low 3' supplementary pairing, guide positions 13-16
```

Issue #102 removed a seventh term, `pos1_mismatch` at 0.05: it is **exactly constant at 0.0** for the
exact-reverse-complement passenger every design uses, so it ranked nothing while consuming weight. Its
0.05 was not reassigned by judgement — the four shared terms were restored to exactly
`0.80 x postscreen_sirna_v4`, the proportional-scaling rule the vector already declared, which lands
on two decimal places without rounding. `pos1_mismatch` is still computed, and the pairing state it
derives from stays on the row as `guide_pos1_base` and `pos1_pairing_state`, so the removal is
auditable and reversible if mismatched passengers are ever designed. `score_pos1_mismatch` is a
*contribution* column, so with the term in no vector it is now **always null**.

Validation is a `@model_validator(mode="after")` reading the subclass's own `TERM_NAMES`.
`COMPOSITE_TERM_NAMES` is now only the ordered **union** of scored terms, used for column ordering;
nothing validates against it. It could not stay the validation input: the three vectors score 3, 4
and 6 different terms, and one global tuple made a missing term look like a licence to renormalise.

`ScoringWeights.vector_for(post_screen=..., design_mode=...)` returns exactly one vector; there is no
API for combining them and no runtime arithmetic on their values. `off_target` cannot be evaluated
until off-target screening has run, so the design stage scores `design_v4` into `design_score` and
`composite_score` stays `None`. The two fields are **different vectors over different term sets and
are not comparable**; `ranking_score(candidate)` is the single place that decides which number a
candidate currently has.

`empirical`, `isoform_coverage` and `conservation` appear in no vector. They are still computed and
reported (§2.7); removing them from scoring is what makes "no renormalisation" reachable, since the
latter two are legitimately `None` on some run shapes while every remaining term is universally
computable. `MiRNADesignConfig.scoring_weights` is gone with them: the miRNA weights live in
`postscreen_mirna_v4`, and the two entries that were declared there and read nowhere in `src/`
(`seed_clean_bonus`, `five_p_end_destabilization_bonus` — 0.25 of declared weight doing nothing) were
deleted.

---

## 2. Scoring Algorithms

### 2.1 Score Calculation

A score is a plain weighted sum over **exactly** the terms its vector declares:

$$\text{Score} = \sum_{i \in \text{vector}} w_i \times S_i \times 100$$

Where $w_i$ are the vector's declared weights, unchanged, and $S_i$ are normalized component scores
(0-1). Because $\sum_i w_i = 1$ by construction, the score spans [0, 100] with no clamping and no
rescaling, and the per-term contributions sum to it exactly.

`compute_composite(features, vector)` **raises** `ScoringError` if `features` is missing any of
`vector.terms`. That is the design: a caller that cannot compute a term has no score to report on
that vector, and must record the absence rather than a number derived from a smaller term set. Extra
keys in `features` are ignored, since `component_scores` also carries diagnostics.

> Before issue #96 this function divided the active weight vector by its own sum. Only four terms
> were computable before screening, so those four shared a 0.50 budget and every design-stage weight
> doubled: `target_accessibility`'s nominal 0.13 applied as 0.26. Two candidates scored under
> different active sets were not comparable, and nothing on the row recorded which vector had
> applied. `CompositeScore` now carries `vector_name`, and so does every candidate row
> (`weight_vector`) and the manifest (`scoring.vectors`).

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

**What it measures**: the equilibrium probability that the stretch of mRNA the guide seed has to
pair is already unpaired, computed with RNAplfold on the transcript.

Before issue #95 this term folded the **guide against itself** and reported the result as target
accessibility. Those are different physical quantities: against a real RNAplfold target-site
opening probability the old term explained r² = 0.01-0.04. Guide self-structure survives as
`paired_fraction` -- reported, and the `EXCESS_PAIRING` gate input -- but is no longer scored, and
the "optimal 0.4-0.6" note that used to appear beside it has been deleted (nothing scored it, and
the code it claimed to describe paid maximum at 0.0).

#### Geometry, which is the whole point

Guide and target are antiparallel. For a target site $T[1..L]$ 5'→3' and a guide $G[1..L]$ 5'→3',
$G[i]$ pairs $T[L+1-i]$. The guide seed $G[2..8]$ therefore pairs the **3' end of the target
site**, so the region RISC must open first is the 3'-most nucleotides.

Anchoring the window on the other end is a measured near-null. Against 2,779 siRNAs with measured
knockdown (Huesken et al. 2005 plus 385 further sequences), at W=150/L=100:

| Statistic                             | Spearman ρ vs measured inhibition | Dynamic range |
| ------------------------------------- | --------------------------------- | ------------- |
| **8-mer at the seed (3') end** ✅     | **+0.267**                        | 6.4 decades   |
| 17-mer, seed-anchored (reported only) | +0.244                            | 10.4 decades  |
| whole 21-nt site (reported only)      | +0.237                            | 13.2 decades  |
| mean P(base unpaired) over the site   | +0.167                            | none (linear) |
| 8-mer at the non-seed (5') end        | +0.067                            | 6.5 decades   |

The four-fold gap between the two 8-mers -- same length, same transcript, same fold, differing only
in which end of the site they cover -- is the strongest evidence the term reflects mechanism rather
than a composition artefact. It is also the regression guard: `tests/unit/test_target_accessibility.py`
asserts the control stays near null, which catches the site being indexed backwards.

**Honest limit**: ρ ≈ +0.27 is about 7% of rank variance. This is a genuine but weak term,
consistent with its 0.13 weight and not an argument for weighting it higher.

#### Statistic and normalisation

$$P = \Pr\left(\text{the 3'-most 8 nt of the target site are unpaired}\right)$$

$$\text{feature} = \frac{\mathrm{clamp}\left(\log_{10} P,\; \text{LOG\_FLOOR},\; 0\right) - \text{LOG\_FLOOR}}{-\text{LOG\_FLOOR}}$$

`LOG_FLOOR` defaults to −5.0, which captures ~99% of observed sites (benchmark distribution of
$\log_{10} P$: p1 = −4.93, p50 = −1.72, p99 = −0.18). It is a **fixed** floor, not a per-transcript
one: a self-calibrating scale would make `composite_score` incomparable between targets and between
runs.

**Implementation** (ViennaRNA `RNA.pfl_fold_up`, one fold per transcript, not per candidate):

```python
profile = TargetAccessibilityProfile.fold(transcript, window_size=150, max_bp_span=100, u_max=21)
site = profile.site_accessibility(start_0based, site_length)   # 3'-anchored windows
feature = target_accessibility_sub_score(site.seed_end_8mer)   # None -> term inactive
```

A global `fold_compound.pf()` is deliberately not used: it is O(n³) on a multi-kb transcript, and
the windowed form is the quantity local-opening models are defined on. Folding once per transcript
is cheaper than the per-candidate guide folds it replaced.

#### Missing values

When there is no transcript context, or a site sits too close to the transcript 5' end for the
window to fit, `target_accessibility_p` is `None`. Since no weight is ever redistributed, the
candidate then has **no score at all** on any vector declaring the term — `design_score` and
`composite_score` stay `None` and the run logs it. A missing input must never score as a good one,
and it must not be quietly rescaled into looking like a complete one either. In practice the 5'-end
case is per-candidate and numerically negligible: 0 of 2,492 TP53 sites, because the scored 8-mer
needs only `site_end >= 8`.

#### Configuration

| Field                                  | Default | Meaning                               |
| -------------------------------------- | ------- | ------------------------------------- |
| `target_accessibility.window_size` (W) | 150     | RNAplfold averaging window            |
| `target_accessibility.max_bp_span` (L) | 100     | RNAplfold maximum base-pair span      |
| `target_accessibility.log_floor`       | −5.0    | log₁₀ P treated as zero accessibility |

W and L are configuration rather than constants because they move a site's accessibility percentile
substantially. ρ rises only mildly with W (+0.249 at W=40 to +0.269 at W=240), and W=150-240 is the
plateau; 150 sits on it without W=240's cost. All three are recorded in the run manifest, and
changing any of them changes the numeric scale of `composite_score`.

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

### 2.7 Empirical Score (Reynolds Rules) — gate only, not a scoring term

**Research basis**: Reynolds et al. (2004)

Position-specific sequence preferences. Since issue #96 this appears in **no** weight vector: it is
computed on every candidate, written to `component_scores["empirical"]` and the `empirical_score`
column, and read only by the `min_empirical_score` gate.

```python
def _calculate_empirical_score(candidate) -> float:
    """Simplified Reynolds rules."""
    guide = candidate.guide_sequence.upper().replace("T", "U")  # guides are stored as DNA
    score = 0.5  # Base score

    # Prefer A/U at position 19 (3' end)
    if guide[18] in ("A", "U"):
        score += 0.1

    # Avoid C at position 19
    if guide[18] == "C":
        score -= 0.1

    return max(EMPIRICAL_SCORE_MIN, min(EMPIRICAL_SCORE_MAX, score))
```

The attainable range is **0.4-0.6**, not 0-1: two mutually exclusive ±0.1 adjustments on a 0.5 base.
`min_empirical_score` is bounded by that range, so `EMPIRICAL_SCORE_MAX` had to move 0.7 → 0.6 with
the clause below. Reading the guide as RNA matters — guides are stored as DNA, so before 0.5.2 a T at
position 19 never earned the A/U bonus.

#### The deleted G/C-at-position-1 clause

There used to be a third rule, `+0.1 if guide[0] in ("G", "C")`. It contradicted the biogenesis rule
rewarding **A/U** at the same base. Measured over 29,605 candidates:

| pos-1 base | n     | empirical | biogenesis adj. | composite |
| ---------- | ----- | --------- | --------------- | --------- |
| A          | 8,917 | 8.02      | **+2.83**       | **52.93** |
| T          | 9,558 | 8.06      | **+2.60**       | **53.24** |
| C          | 5,854 | **9.62**  | −5.30           | 45.21     |
| G          | 5,276 | **9.69**  | −5.00           | 44.56     |

G/C gained +1.6 on the empirical term and lost 7.9 on the biogenesis adjustment — net ~8 composite
points worse. A declared 0.15-weight term was overridden ~5× by an undeclared one, `empirical`
correlated r = −0.48 with the adjustment, and its variance share came out **negative**. A/U wins; the
clause is gone, and a test asserts the empirical score is now invariant to guide position 1.

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

| Filter                        | Condition                                                                        | Rationale                                          |
| ----------------------------- | -------------------------------------------------------------------------------- | -------------------------------------------------- |
| `TRANSCRIPTOME_PERFECT_MATCH` | 0-mismatch genuine off-targets > `max_transcriptome_hits_0mm`                    | Perfect-match specificity                          |
| `TRANSCRIPTOME_1MM`           | 1-mismatch genuine off-targets > `max_transcriptome_hits_1mm`                    | Near-miss specificity                              |
| `TRANSCRIPTOME_2MM`           | 2-mismatch genuine off-targets > `max_transcriptome_hits_2mm`                    | Broader specificity                                |
| `TRANSCRIPTOME_SEED_PERFECT`  | seed-perfect hits > `max_transcriptome_seed_perfect` (`None` = off, all species) | Catches partial hits that no mismatch stratum sees |
| `MIRNA_PERFECT_SEED`          | perfect miRNA seed hits > `max_mirna_perfect_seed`                               | miRNA-mimicry risk                                 |
| `HIGH_RISK_MIRNA`             | perfect seed + `offtarget_score` < 5.0                                           | Strong-binding miRNA mimicry                       |
| `TOTAL_OFFTARGETS`            | combined transcriptome + miRNA hits > `max_total_offtarget_hits`                 | Aggregate specificity                              |
| `EXCESS_OFF_TARGETS`          | genuine off-target count > `max_off_target_count`                                | Overall genuine-off-target ceiling                 |

Six of these seven verdicts previously existed only as raw strings assembled in `workflow.py`
(`TRANSCRIPTOME_SEED_PERFECT` is new); they are all now members of `FilterStatus` (see 3.4), and
the Pandera allow-list for the `passes_filters` column derives from the enum instead of
maintaining an independent copy that could drift from it.

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
    TRANSCRIPTOME_SEED_PERFECT = "TRANSCRIPTOME_SEED_PERFECT"  # opt-in, see 1.4
    MIRNA_PERFECT_SEED = "MIRNA_PERFECT_SEED"
    HIGH_RISK_MIRNA = "HIGH_RISK_MIRNA"
    TOTAL_OFFTARGETS = "TOTAL_OFFTARGETS"
```

### 3.5 Hit Classes (`sirnaforge.core.hit_classification.HitClass`)

Every transcriptome hit is classified into exactly one of five mutually exclusive classes,
checked in this precedence order (on-target and ortholog deliberately outrank repeat, so a
repeat-flagged guide's hits on its own gene or orthologs are still counted as such):

| Class          | Meaning                                                                                                                                                                       |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ON_TARGET`    | Hit is in the query species and matches the query gene (transcript ID, gene ID or symbol)                                                                                     |
| `ORTHOLOG`     | Hit is in a different screened species, on an orthologue of the query gene — by gene ID against the Compara mapping, else by symbol equality (`ortholog_evidence` says which) |
| `REPEAT`       | The guide is in the design-time repeat-flagged set                                                                                                                            |
| `OFF_TARGET`   | Everything else that could be decided                                                                                                                                         |
| `UNDETERMINED` | The hit's species has no transcript index, so orthology and the query gene could not be checked at all. Assigned by the annotation layer, never returned by `classify_hit`    |

`off_target_count` counts **both** `OFF_TARGET` and `UNDETERMINED` — the liability population — so
that removing a reference cannot loosen the screen. `undetermined_hits` reports how much of that
count is unqualified.

A hit whose species label is blank/missing is treated as the query species. The gene-ID tier needs no
symbol, so it is tried first; only when it finds nothing does a missing gene symbol stop the check, in
which case the hit stays `OFF_TARGET` (never guessed at) and the shortfall is counted separately
(`ortholog_symbol_lookup_misses` in `offtarget_summary`, which therefore counts symbol-tier shortfalls
only).

**Orthology evidence.** The two tiers are not the same claim, so which one fired is published per row
as `ortholog_evidence` (`gene_id` / `symbol_heuristic`, and `not_applicable` on every non-ortholog
row). `gene_id` is a lookup in the orthologue gene-ID set `sirnaforge.data.orthology` resolves from
Ensembl Compara (`ortholog_one2one`/`one2many`/`many2many` only — paralogues are liabilities, not
conservation) and handed to `ClassificationContext.ortholog_gene_ids`; `classify_hit` stays a pure
function and performs no lookup itself. `symbol_heuristic` is uppercased symbol equality, kept because
it is free offline and right for the many genes whose symbol is conserved, but wrong in both
directions: HGNC and MGI are separate authorities, so human `TP53`'s mouse orthologue is `Trp53`
(`TRP53` after ingest) and symbol equality cannot reach it, while unrelated same-symbol genes satisfy
it. Species whose lookup could not be completed are named in
`offtarget_summary.filtering_stats.orthology.unresolved_species` in `logs/workflow_summary.json` and
fall back to the heuristic; species with no hits at all are never looked up, so a single-species screen
makes no network request. `conservation_score` is derived from `ortholog_species`, so it inherits
whichever tiers those verdicts came from.

The **query species** is the organism the _target_ transcripts belong to. It is read from the
database the gene query was answered by (Ensembl/RefSeq/GENCODE are all human-only), never from
`--species`, which is an unordered set of genomes to screen _against_ and whose order carries no
meaning. Pass `--query-species` when designing against an input FASTA from another organism — it
also decides which species' alignment must have succeeded before candidates can be scored after
screening.

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

    # Thresholds and format defaults only -- the miRNA scoring weights live in
    # PostScreenMiRNAWeights (postscreen_mirna_v4), see 1.5.
```

### 5.2 miRNA-Specific Scoring

Two of the three biogenesis quantities are **ordinary declared terms** of `postscreen_mirna_v4`
(`ago_start` 0.10, `supp_13_16` 0.10), each a feature in [0, 1] like every other term. The third,
`pos1_mismatch`, is computed and reported and scored by nothing — issue #102 measured it exactly
constant. `biogenesis_features(guide, passenger)` derives all three from sequence alone, so they
are available for any candidate — including rows that never passed through `MiRNADesigner`, such as
the dirty controls cloned from rejected candidates.

The design stage scores `design_v4` in **both** modes, so a design-stage miRNA score is bit-identical
to the siRNA one; the biogenesis terms enter only once `off_target` exists.

> Before issue #96 the three were bonuses folded into `composite_score` and the result divided by
> `1 + max_bonus = 1.25`, which scaled every declared weight by 0.80 in miRNA mode. `off_target`'s
> nominal 0.25 applied as 0.20, a candidate earning no bonus kept only 80% of its score, and the
> adjustment appeared in neither `ScoringWeights` nor the manifest. Both post-screen vectors now sum
> to 1.0, so the two modes are on one scale: a candidate whose biogenesis sub-scores match its other
> sub-scores scores identically in either mode, at every level.

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
    nm: int              # Guide mismatch-equivalents (>= the aligner's NM tag)
    seed_mismatches: int # Mismatches in seed (pos 2-8)
    offtarget_score: float
```

`nm` is a **guide-level** distance, not a copy of the aligner's `NM` tag. `NM` counts only
differences inside the aligned block, so a `bwa mem -T 15` hit reported as `15M6S` with `NM:i:0`
looks like a perfect match even though 6 of the guide's 21 bases never paired with the target.
`nm` therefore counts the aligner's edit distance **plus** every guide base left unpaired by soft
or hard clipping and by insertions, plus the reference bases skipped by a deletion. Consequences
worth knowing:

- `mismatch_positions` and `seed_mismatches` are in **guide coordinates**. BWA stores minus-strand
  records reverse-complemented, and `design.py` builds every guide as
  `reverse_complement(target_seq)`, so minus-strand is the common case; read position `r` of a
  length-`L` record is guide position `L + 1 - r`. In `mirna_seed` mode the aligner only sees the
  extracted seed window, so positions are additionally shifted by `seed_start - 1`.
- `mismatch_positions` lists only positions that exist on the guide. A deletion has no guide
  position of its own, so it raises `nm` and the score without adding an entry — that is the one
  case where `len(mismatch_positions) < nm`.
- `offtarget_score` is a penalty: `0.0` means highest risk, and it is reserved for full-length
  exact matches. `_filter_and_rank` sorts ascending on it and `max_hits` keeps the head of that
  list, so a hit with `nm > 0` never scores `0.0`.
- The mismatch-stratified `transcriptome_hits_{0,1,2}mm` counters are read off `nm`, so a clipped
  partial hit lands in the stratum matching its guide-level distance rather than in `0mm` — and if
  that distance exceeds 2 it lands in **no** stratum, which is why a seed-perfect partial hit needs
  `max_transcriptome_seed_perfect` (see the warning in 1.4) rather than
  `max_transcriptome_hits_0mm`.
- **`AnalysisSummary.mean_mismatches`** (the `mean_mismatches` key of every
  `*_summary.json`) is the mean of this `nm`, so it also changed meaning: it is the mean
  guide mismatch-equivalent count, **not** the mean aligner edit distance, and it rose for any run
  containing clipped or gapped hits. Two `15M6S`/`NM:i:0` hits now report `mean_mismatches = 6.0`
  where they used to report `0.0`. Values are not comparable across the 0.6.0 boundary.
  `mean_seed_mismatches` and `mean_mapq` are unaffected in definition (though
  `mean_seed_mismatches` changes numerically, since the frame fix corrected which positions are
  in the seed).
- **`SiRNACandidate.off_target_penalty` is reporting only** and its direction depends on which
  stage wrote it last, so treat it as a diagnostic, not a risk metric. At design time
  `SiRNADesigner._calculate_off_target_score` writes an internal-repeat 7-mer penalty where higher
  is worse; after screening `_integrate_offtarget_results` overwrites it with the **maximum**
  `offtarget_score` over the candidate's hits, where higher is _safer_ and `0.0` is reserved for a
  full-length exact match. Taking the **max** means the field reports a candidate's _least_
  worrying hit: a candidate with one full-length perfect off-target reports
  `off_target_penalty = 0.0`, while a candidate whose only hit is a clipped partial reports 76
  (`6S15M`) or 98 (`15M6S`). Widening `nm` widened these numbers too. Judge risk from
  `off_target_count` and the hit strata.

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

7. **Tafer H, Ameres SL, Obernosterer G, Gebeshuber CA, Schroeder R, Martinez J, Hofacker IL** (2008). The impact of target site accessibility on the design of effective siRNAs. _Nature Biotechnology_ 26(5):578-583. — Cited for the concept of scoring mRNA local opening, not as the source of the exact statistic used here. Before issue #95 this reference was attached to a term that folded the guide, which the paper does not describe.

8. **Huesken D, Lange J, Mickanin C, Weiler J, Asselbergs F, Warner J, Meloon B, Engel S, Rosenberg A, Cohen D, Labow M, Reinhardt M, Natt F, Hall J** (2005). Design of a genome-wide siRNA library using an artificial neural network. _Nature Biotechnology_ 23(8):995-1001. — The measured-knockdown benchmark the accessibility statistic was selected against (via a third-party redistribution; the values were not verified against the primary paper, which is not open access).

9. **Jackson AL, Bartz SR, Schelter J, Kobayashi SV, Burchard J, Mao M, Li B, Cavet G, Linsley PS** (2003). Expression profiling reveals off-target gene regulation by RNAi. _Nature Biotechnology_ 21(6):635-637.

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

Target-site accessibility settings (see 2.5), on `DesignParameters.target_accessibility`:

| Parameter     | Default | Range   | Justification                                       |
| ------------- | ------- | ------- | --------------------------------------------------- |
| `window_size` | 150     | 20-1000 | RNAplfold W; on the benchmark ρ plateau at 150-240  |
| `max_bp_span` | 100     | 10-1000 | RNAplfold L; must not exceed `window_size`          |
| `log_floor`   | −5.0    | < 0     | Captures ~99% of observed sites; fixed, not per-run |

`mfe_min`, `mfe_max`, `duplex_stability_min`, `duplex_stability_max`, `melting_temp_min`,
`melting_temp_max`, `delta_dg_end_min` and `delta_dg_end_max` were removed from `FilterCriteria`
in issue #80: they were declared but never enforced by `SiRNADesigner`, and re-deriving correct
windows against truth data is deliberately out of scope. `max_off_target_count` (default `15`,
range `0+`, justification: specificity) moved to `OffTargetFilterCriteria`, where it is enforced
post-screen against the genuine off-target count (see 1.4 and 3.3).

## Appendix B: Scoring Weight Defaults

Weight-set version `4.0.0`. Three hand-authored vectors, each summing to exactly 1.0, never rescaled
at runtime. These are **declared expert priors**: only `target_accessibility` and `off_target` have
benchmark evidence behind them.

**`design_v4`** → `design_score` (design stage, both modes)

| Term                 | Weight | Rationale                                                                                              |
| -------------------- | ------ | ------------------------------------------------------------------------------------------------------ |
| Asymmetry            | 0.40   | Indistinguishable from accessibility on the benchmark (ρ +0.273 vs +0.267), and additionally ρ +0.53 with A/U at guide positions 1-5 |
| Target accessibility | 0.35   | RNAplfold local opening at the seed-paired end; the only term with a knockdown benchmark at this stage |
| GC content           | 0.25   | Stability/accessibility balance                                                                        |

⚠️ These three numbers are round numbers **awaiting sign-off** — they are the one part of the weight
set not chosen by the repo owner.

**`postscreen_sirna_v4`** → `composite_score` (post-screen, siRNA mode). `design_v4`'s terms plus one.

| Term                 | Weight | Rationale                                                                                          |
| -------------------- | ------ | -------------------------------------------------------------------------------------------------- |
| Off-target           | 0.25   | Post-screen genuine off-target specificity; measured 2.24× its nominal share of composite variance |
| Target accessibility | 0.30   | as above                                                                                           |
| Asymmetry            | 0.25   | as above                                                                                           |
| GC content           | 0.20   | as above                                                                                           |

> Holding `off_target` at 0.25 while the scored budget shrank from six terms to four **reduces** its
> relative influence, from 0.25/0.60 of the old scored budget to 0.25/1.00. Given it measured at 56%
> of composite variance that is probably the right direction, but it arrives as a side effect of the
> restructuring rather than as an explicit choice, and should be defensible deliberately.

**`postscreen_mirna_v4`** → `composite_score` (post-screen, `--design-mode mirna`)

| Term                 | Weight | Rationale                                                             |
| -------------------- | ------ | --------------------------------------------------------------------- |
| Off-target           | 0.20   | exactly 0.80 × the siRNA vector's 0.25                                |
| Target accessibility | 0.24   | exactly 0.80 × 0.30                                                   |
| Asymmetry            | 0.20   | exactly 0.80 × 0.25; ties `off_target`, as the siRNA vector also does  |
| GC content           | 0.16   | exactly 0.80 × 0.20                                                   |
| `ago_start`          | 0.10   | A/U at guide position 1 (Argonaute loading)                           |
| `supp_13_16`         | 0.10   | low 3' supplementary pairing potential; endpoint claimed, not measured |

The two biogenesis terms hold 0.20 and the four shared terms hold exactly 0.80 × the siRNA values.
Before issue #102 they were 0.22 / 0.18 / 0.15 with `off_target` at 0.20 — a rounded 0.75× scaling
with the resulting `off_target`/`asymmetry` tie broken in favour of `off_target`. Removing the
constant `pos1_mismatch` released 0.05 and made the exact scaling reachable at two decimal places,
so the rounding and its tie-break are both gone. This is still a **declared** scaling of a
hand-authored vector, not a fit; the numbers moved because a term left, not because anything was
tuned.

**Computed and reported, in no vector:** `empirical` (the `min_empirical_score` gate),
`isoform_coverage` (the optional `min_isoform_coverage` gate, default off), `conservation` (reporting
only) and `paired_fraction` (the `max_paired_fraction` gate).

Version comparability: `4.x` removed both hidden normalisations and restructured the term sets, so no
`3.x` score is comparable with it. `3.x` replaced guide self-structure with real target-site
accessibility in the 0.13 slot (issue #95), and `1.x` denotes the pre-issue-#80 five-term set.

---

_Document version: 1.0_
_Last updated: Auto-generated from source code_
_Review status: Initial draft - Expert review recommended_
