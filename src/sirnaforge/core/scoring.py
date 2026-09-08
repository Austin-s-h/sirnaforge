"""Composite siRNA scoring with renormalized sub-score contributions.

This module provides the canonical composite scorer for siRNA candidates, computing
a single 0-100 quality score from weighted sub-scores (asymmetry, GC content,
accessibility, empirical rules, off-target specificity, isoform coverage, and
ortholog conservation).

The scorer renormalizes weights over the active term set, so a run that lacks
evidence for some terms (e.g., no conservation term in a single-species run, or no
isoform_coverage when the gene has no protein-coding transcripts) is neither
rewarded nor penalised for the missing terms. This design lets the same scorer
produce a legitimate single-scale composite for both design-only and post-screen
contexts without special-casing.

Version history:
    - 1.x: pre-issue-#80 five-term set (asymmetry, gc_content, accessibility,
      empirical, off_target proxy). Not comparable to 2.x.
    - 2.0.0: seven-term set with post-screen off-target redefined, isoform_coverage
      and conservation added. Weights retuned. Bump this version whenever a DEFAULT
      weight in ScoringWeights changes.
"""

from collections.abc import Collection, Mapping
from dataclasses import dataclass
from math import exp

from sirnaforge.models.sirna import COMPOSITE_TERM_NAMES, ScoringWeights

# Scoring weight set version. Bump when DEFAULT weights change.
SCORING_WEIGHT_SET_VERSION = "2.0.0"

# Canonical composite-score term names (imported from models for consistency).
COMPOSITE_TERMS = COMPOSITE_TERM_NAMES

# Off-target sub-score decay constant.
OFF_TARGET_DECAY = 10.0


class ScoringError(ValueError):
    """Raised when scoring inputs are invalid or unnormalizable."""

    pass


@dataclass(frozen=True)
class CompositeScore:
    """Result of composite scoring with per-term contribution breakdown.

    Attributes:
        score: Composite score in [0, 100].
        weight_set_version: Version string identifying the weight set used.
        active_terms: Terms that contributed to the score, in declaration order.
        contributions: Per-term contributions (renormalized_weight * feature * 100),
            summing to score. Only active terms appear; a missing term was inactive,
            not zero-valued.
    """

    score: float
    weight_set_version: str
    active_terms: tuple[str, ...]
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


def isoform_coverage_sub_score(protein_coding_hit: int, protein_coding_total: int) -> float | None:
    """Compute the isoform coverage sub-score for protein-coding transcripts.

    Returns hit / total when total > 0, or None when total == 0 (the term is
    inactive — an annotation gap or a gene with no protein-coding transcripts must
    not be scored as zero; that would penalise the candidate for the annotation's
    shortcoming).

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
    """Compute the cross-species conservation sub-score.

    Returns hit / requested when requested > 0, or None when requested == 0 (the
    term is inactive for a single-species run — a run that screened only the query
    species must not be penalised for lacking ortholog evidence).

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


def compute_composite(
    features: Mapping[str, float],
    weights: ScoringWeights,
    active_terms: Collection[str] | None = None,
) -> CompositeScore:
    """Compute the composite siRNA quality score from weighted sub-scores.

    This is a pure function: same inputs always yield the same output, with no I/O,
    no global state, and no logging side effects. Weights are renormalized over the
    active term set so that missing terms (e.g., no conservation term in a
    single-species run) do not bias the composite.

    Args:
        features: Mapping from term name to sub-score in [0, 1]. Terms not present
            in this mapping are considered inactive. Sub-scores outside [0, 1]
            indicate an upstream bug and raise ScoringError.
        weights: ScoringWeights instance with weights for all known terms. The
            ScoringWeights validator ensures the full weight vector sums to ~1.0.
        active_terms: Optional explicit active term set. Defaults to the terms
            actually present in `features`. Only active terms contribute; their
            weights are renormalized to sum to 1.

    Returns:
        CompositeScore with score in [0, 100], version, active terms, and per-term
        contributions (each in [0, 100], summing to score).

    Raises:
        ScoringError: If any sub-score is outside [0, 1], the active term set is
            empty or contains only unknown terms, or the active weight vector is
            unnormalizable (all weights zero).
    """
    # Default active_terms to the keys actually present in features.
    if active_terms is None:
        resolved_active = tuple(term for term in COMPOSITE_TERMS if term in features)
    else:
        # Keep only known terms that are present in features, in declaration order.
        resolved_active = tuple(term for term in COMPOSITE_TERMS if term in active_terms and term in features)

    if not resolved_active:
        raise ScoringError("No active terms: either features is empty or active_terms names no known term in features")

    # Validate that all active sub-scores are in [0, 1].
    for term in resolved_active:
        value = features[term]
        if not (0.0 <= value <= 1.0):
            raise ScoringError(
                f"Sub-score for term '{term}' is {value:.4f}, outside [0, 1]. "
                "This indicates an upstream bug in the sub-score computation."
            )

    # Extract and renormalize the active weight vector.
    active_weights = {term: getattr(weights, term) for term in resolved_active}
    weight_sum = sum(active_weights.values())
    if weight_sum <= 0.0:
        raise ScoringError(
            f"Active weight sum is {weight_sum:.6f} (non-positive, cannot normalize). Active terms: {resolved_active}"
        )

    renorm_weights = {term: w / weight_sum for term, w in active_weights.items()}

    # Compute per-term contributions and composite score.
    contributions = {term: renorm_weights[term] * features[term] * 100.0 for term in resolved_active}
    # An all-1.0 feature vector sums to 100 only up to float error, and SiRNACandidate declares
    # composite_score as le=100, so pin the endpoints rather than let rounding fail validation.
    score = min(100.0, max(0.0, sum(contributions.values())))

    return CompositeScore(
        score=score,
        weight_set_version=SCORING_WEIGHT_SET_VERSION,
        active_terms=resolved_active,
        contributions=contributions,
    )
