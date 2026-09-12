"""Composite siRNA scoring against hand-authored, named weight vectors.

This module provides the canonical scorer. It applies **no arithmetic to weights at runtime**:
there is no renormalisation, no rescaling and no division anywhere below. A caller hands in one
`WeightVector` -- already named, already summing to 1.0, validated at construction -- together with
a feature for **every** term that vector declares. A missing term is an error, not an invitation to
redistribute its weight, because redistributing it silently doubled every remaining weight and made
two candidates scored under different term sets incomparable while nothing on the row recorded which
weights had applied.

Which vector applies is decided by stage and design mode, never combined (`ScoringWeights.vector_for`):

    design_v4             design stage, both modes -> SiRNACandidate.design_score
    postscreen_sirna_v4   post-screen, siRNA mode  -> SiRNACandidate.composite_score
    postscreen_mirna_v4   post-screen, miRNA mode  -> SiRNACandidate.composite_score

`design_score` and `composite_score` are different vectors over different term sets and are not
comparable with each other. Every vector sums to 1.0, so every score spans the same [0, 100] and its
name is written to both the candidate row and the run manifest.

Terms that are computed but deliberately **not** scored: `empirical` (gate only, via
`min_empirical_score`), `isoform_coverage` (reported, plus the optional `min_isoform_coverage` gate)
and `conservation` (reported only). Their sub-score helpers still live here because they still have
to be computed and validated -- they simply no longer appear in any vector. Removing them is what
makes "no renormalisation" reachable rather than merely relocated: both are legitimately None on some
run shapes, and every term that remains is universally computable.

Version history:
    - 1.x: pre-issue-#80 five-term set (asymmetry, gc_content, accessibility,
      empirical, off_target proxy). Not comparable to 2.x.
    - 2.0.0: seven-term set with post-screen off-target redefined, isoform_coverage
      and conservation added. Weights retuned.
    - 3.0.0: issue #95. `accessibility` (which folded the guide against itself and
      reported it as target accessibility) is replaced by `target_accessibility`, a
      real RNAplfold local-opening probability on the transcript. Same 0.13 weight,
      different quantity, so 2.x scores are not comparable. Guide self-structure
      survives as `paired_fraction` -- reported, and the EXCESS_PAIRING gate input --
      but is no longer a scoring term.
    - 4.0.0: issues #96 and #102, one comparability break between them. Both hidden
      normalisations removed: the active-set renormalisation here, and the
      `1 + max_mirna_bonus` divisor that scaled every miRNA score by 0.80. One flat
      weight vector becomes three named ones; `empirical`, `conservation` and
      `isoform_coverage` leave the composite; the miRNA biogenesis bonuses become
      declared terms. `pos1_mismatch` is **not** among them -- #102's audit measured it
      exactly constant at 0.0 on all 13,415 scored baseline candidates, so it holds no
      weight and `postscreen_mirna_v4` has six terms, its shared four at exactly
      0.80 x `postscreen_sirna_v4`. No 3.x score is comparable with a 4.x one.
      Bump this version whenever any DEFAULT weight or any vector's term set changes.

Every weight and threshold in this release is `experimental`. `models/scoring_profile.py`
records, per term, what it reads and how far its evidence goes; nothing has been held out and
independently replicated, so no term is `validated` and these scores rank candidates rather
than predicting a knockdown level or a safety probability.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from math import exp, log10

from sirnaforge.models.sirna import (
    COMPOSITE_TERM_NAMES,
    DEFAULT_ACCESSIBILITY_LOG_FLOOR,
    WeightVector,
)

# Scoring weight set version. Bump when DEFAULT weights or any vector's term set change.
SCORING_WEIGHT_SET_VERSION = "4.0.0"

# Every scored term name in reporting order -- the union over the three vectors, not a vector
# itself. Nothing validates against it; each vector validates against its own TERM_NAMES.
COMPOSITE_TERMS = COMPOSITE_TERM_NAMES

# Off-target sub-score decay constant.
OFF_TARGET_DECAY = 10.0


class ScoringError(ValueError):
    """Raised when scoring inputs are invalid or unnormalizable."""

    pass


@dataclass(frozen=True)
class CompositeScore:
    """Result of scoring one candidate against one named weight vector.

    Attributes:
        score: Score in [0, 100].
        weight_set_version: Version string identifying the weight set used.
        vector_name: Name of the weight vector that produced the score, recorded on the
            candidate row so it resolves to the exact weights in the run manifest.
        terms: The vector's terms, in declaration order. Every one contributed -- there is no
            such thing as an inactive term any more.
        contributions: Per-term contributions (weight * feature * 100), summing to score.
    """

    score: float
    weight_set_version: str
    vector_name: str
    terms: tuple[str, ...]
    contributions: dict[str, float]


def off_target_sub_score(genuine_off_target_count: int) -> float:
    """Compute the off-target specificity sub-score from genuine off-target count.

    The score decays exponentially with off-target count: exp(-count / OFF_TARGET_DECAY).
    Zero off-targets yield 1.0 (perfect specificity); the score approaches 0 as count grows.

    Args:
        genuine_off_target_count: Number of genuine off-target sites (on-target,
            ortholog and repeat hits excluded).

    Returns:
        Off-target sub-score in [0, 1].

    Raises:
        ValueError: If count is negative.
    """
    if genuine_off_target_count < 0:
        raise ValueError(f"genuine_off_target_count must be non-negative, got {genuine_off_target_count}")
    score = exp(-genuine_off_target_count / OFF_TARGET_DECAY)
    return max(0.0, min(1.0, score))  # Clamp to [0, 1] for numerical safety.


def target_accessibility_sub_score(
    p_unpaired: float | None, log_floor: float = DEFAULT_ACCESSIBILITY_LOG_FLOOR
) -> float | None:
    """Normalise a target-site opening probability into a [0, 1] sub-score.

    The probability spans ~6.4 decades across real target sites, so it is log-scaled before
    being mapped linearly onto [0, 1] against a fixed floor:

        feature = (clamp(log10 P, log_floor, 0) - log_floor) / -log_floor

    The floor is fixed rather than per-transcript on purpose: a self-calibrating scale would make
    composite scores incomparable between targets and between runs.

    Args:
        p_unpaired: P(the scored window of the target site is unpaired), or None when it could
            not be computed (no transcript context, or the site is too close to the 5' end).
        log_floor: log10 probability treated as zero accessibility. Must be negative.

    Returns:
        Sub-score in [0, 1], or None when p_unpaired is None (term inactive). A probability at
        or below the floor -- including exactly 0, which underflow can produce -- returns 0.0:
        that is a computed "inaccessible", not a missing value.

    Raises:
        ValueError: If log_floor is not negative, or p_unpaired is outside [0, 1].
    """
    if log_floor >= 0.0:
        raise ValueError(f"log_floor must be negative, got {log_floor}")

    if p_unpaired is None:
        return None  # Term inactive.

    if not (0.0 <= p_unpaired <= 1.0):
        raise ValueError(f"p_unpaired must be a probability in [0, 1], got {p_unpaired}")

    if p_unpaired <= 0.0:
        return 0.0  # log10(0) is -inf, which clamps to the floor anyway.

    clamped = max(log_floor, min(0.0, log10(p_unpaired)))
    return (clamped - log_floor) / -log_floor


def isoform_coverage_sub_score(protein_coding_hit: int, protein_coding_total: int) -> float | None:
    """Compute protein-coding isoform coverage. Reported and gate-only; not a scoring term.

    Left the composite in issue #96 (it is None on `design_from_sequence` and miRNA paths, so any
    vector containing it needed either variant vectors or arithmetic), and is still computed on
    every candidate: it is written to `SiRNACandidate.isoform_coverage` and read by the optional
    `FilterCriteria.min_isoform_coverage` gate.

    Returns hit / total when total > 0, or None when total == 0 (an annotation gap or a gene with
    no protein-coding transcripts must not read as zero coverage; that would penalise the candidate
    for the annotation's shortcoming, and the gate skips a None).

    Args:
        protein_coding_hit: Number of protein-coding transcripts hit.
        protein_coding_total: Total number of protein-coding transcripts in the gene.

    Returns:
        Isoform coverage sub-score in [0, 1], or None when total == 0 (inactive).

    Raises:
        ValueError: If either count is negative, or hit > total.
    """
    if protein_coding_hit < 0:
        raise ValueError(f"protein_coding_hit must be non-negative, got {protein_coding_hit}")
    if protein_coding_total < 0:
        raise ValueError(f"protein_coding_total must be non-negative, got {protein_coding_total}")

    if protein_coding_total == 0:
        return None  # Term inactive.

    if protein_coding_hit > protein_coding_total:
        raise ValueError(
            f"protein_coding_hit ({protein_coding_hit}) cannot exceed protein_coding_total ({protein_coding_total})"
        )

    return protein_coding_hit / protein_coding_total


def conservation_sub_score(ortholog_species_hit: int, requested_non_query_species: int) -> float | None:
    """Compute cross-species conservation. Reported only since issue #96; not a scoring term.

    It is still computed on every candidate and written to `SiRNACandidate.conservation_score`; no
    weight reads it. It left the composite because it is None on single-species runs, and its
    returning None used to be the other half of what forced the deleted renormalisation.

    Returns hit / requested when requested > 0, or None when requested == 0 (a run that screened
    only the query species has no ortholog evidence, which is not the same as no conservation).

    Args:
        ortholog_species_hit: Number of non-query species with at least one ortholog hit.
        requested_non_query_species: Number of non-query species the user requested.

    Returns:
        Conservation sub-score in [0, 1], or None when requested == 0 (inactive).

    Raises:
        ValueError: If either count is negative, or hit > requested.
    """
    if ortholog_species_hit < 0:
        raise ValueError(f"ortholog_species_hit must be non-negative, got {ortholog_species_hit}")
    if requested_non_query_species < 0:
        raise ValueError(f"requested_non_query_species must be non-negative, got {requested_non_query_species}")

    if requested_non_query_species == 0:
        return None  # Term inactive.

    if ortholog_species_hit > requested_non_query_species:
        raise ValueError(
            f"ortholog_species_hit ({ortholog_species_hit}) cannot exceed "
            f"requested_non_query_species ({requested_non_query_species})"
        )

    return ortholog_species_hit / requested_non_query_species


def compute_composite(features: Mapping[str, float], vector: WeightVector) -> CompositeScore:
    """Score one candidate against one named weight vector.

    A pure function: same inputs always yield the same output, with no I/O, no global state and no
    logging side effects. The vector's weights are used **exactly as declared** -- this function
    contains no arithmetic on them beyond multiplying each by its feature.

    Every term the vector declares must be present in `features`. That is the whole point: scoring
    whatever happened to be available used to renormalise the weights over it -- now deleted -- so
    the same nominal 0.13 was worth 0.13 or 0.26 depending on the run stage. A caller that cannot
    compute a term has no score to report, and must record that instead of a rescaled one.

    Args:
        features: Mapping from term name to sub-score in [0, 1]. Must cover every term in
            `vector.terms`; extra keys are ignored (component_scores carries diagnostics too).
            Sub-scores outside [0, 1] indicate an upstream bug and raise ScoringError.
        vector: The named, hand-authored weight vector to score against.

    Returns:
        CompositeScore with score in [0, 100], version, vector name, terms, and per-term
        contributions (each in [0, 100], summing to score).

    Raises:
        ScoringError: If the vector declares a term `features` does not supply, or if any
            sub-score is outside [0, 1].
    """
    terms = vector.terms

    missing = [term for term in terms if term not in features]
    if missing:
        raise ScoringError(
            f"Weight vector '{vector.name}' requires {list(terms)} but {missing} were not supplied. "
            "Weights are never renormalised, so a candidate missing a term has no score on this "
            "vector; record that rather than a rescaled one."
        )

    for term in terms:
        value = features[term]
        if not (0.0 <= value <= 1.0):
            raise ScoringError(
                f"Sub-score for term '{term}' is {value:.4f}, outside [0, 1]. "
                "This indicates an upstream bug in the sub-score computation."
            )

    weights = vector.as_mapping()
    contributions = {term: weights[term] * features[term] * 100.0 for term in terms}
    # An all-1.0 feature vector sums to 100 only up to float error, and SiRNACandidate declares its
    # score fields as le=100, so pin the endpoints rather than let rounding fail validation.
    score = min(100.0, max(0.0, sum(contributions.values())))

    return CompositeScore(
        score=score,
        weight_set_version=SCORING_WEIGHT_SET_VERSION,
        vector_name=vector.name,
        terms=terms,
        contributions=contributions,
    )
