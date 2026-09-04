"""Unit tests for composite siRNA scoring with renormalized sub-score contributions."""

import math

import pytest

from sirnaforge.core.scoring import (
    COMPOSITE_TERMS,
    SCORING_WEIGHT_SET_VERSION,
    ScoringError,
    compute_composite,
    conservation_sub_score,
    isoform_coverage_sub_score,
    off_target_sub_score,
)
from sirnaforge.models.sirna import ScoringWeights


@pytest.mark.unit
class TestOffTargetSubScore:
    """Tests for the off-target specificity sub-score helper."""

    def test_zero_off_targets_yields_perfect_score(self) -> None:
        """Zero off-targets should yield 1.0 (perfect specificity)."""
        assert off_target_sub_score(0) == 1.0

    def test_known_off_target_counts_yield_known_scores(self) -> None:
        """Known off-target counts should produce known exponential decay scores."""
        # Hand-computed: exp(-count / 10.0)
        assert off_target_sub_score(1) == pytest.approx(math.exp(-1 / 10.0), abs=1e-9)
        assert off_target_sub_score(5) == pytest.approx(math.exp(-5 / 10.0), abs=1e-9)
        assert off_target_sub_score(10) == pytest.approx(math.exp(-1), abs=1e-9)
        assert off_target_sub_score(100) == pytest.approx(math.exp(-10), abs=1e-9)

    def test_large_count_approaches_zero(self) -> None:
        """Large off-target counts should approach zero."""
        # At count = 10 * OFF_TARGET_DECAY = 100, exp(-10) ≈ 0.000045
        score = off_target_sub_score(100)
        assert 0.0 <= score < 0.001

    def test_negative_count_raises(self) -> None:
        """Negative off-target count should raise ValueError."""
        with pytest.raises(ValueError, match="must be non-negative"):
            off_target_sub_score(-1)


@pytest.mark.unit
class TestIsoformCoverageSubScore:
    """Tests for the isoform coverage sub-score helper."""

    def test_full_coverage_yields_one(self) -> None:
        """Hitting all protein-coding isoforms should yield 1.0."""
        assert isoform_coverage_sub_score(5, 5) == 1.0

    def test_partial_coverage_yields_fraction(self) -> None:
        """Partial isoform coverage should yield the exact hit/total fraction."""
        assert isoform_coverage_sub_score(3, 5) == 0.6
        assert isoform_coverage_sub_score(1, 4) == 0.25
        assert isoform_coverage_sub_score(7, 10) == 0.7

    def test_zero_hits_yields_zero(self) -> None:
        """Zero isoform hits with nonzero total should yield 0.0."""
        assert isoform_coverage_sub_score(0, 5) == 0.0

    def test_zero_total_returns_none(self) -> None:
        """Zero total protein-coding isoforms should return None (inactive term)."""
        # An annotation gap or gene with no protein-coding transcripts must not be
        # scored as zero; that would penalise the candidate for the annotation.
        assert isoform_coverage_sub_score(0, 0) is None

    def test_hit_exceeds_total_raises(self) -> None:
        """Hit count exceeding total should raise ValueError."""
        with pytest.raises(ValueError, match="cannot exceed"):
            isoform_coverage_sub_score(6, 5)

    def test_negative_hit_raises(self) -> None:
        """Negative hit count should raise ValueError."""
        with pytest.raises(ValueError, match="must be non-negative"):
            isoform_coverage_sub_score(-1, 5)

    def test_negative_total_raises(self) -> None:
        """Negative total count should raise ValueError."""
        with pytest.raises(ValueError, match="must be non-negative"):
            isoform_coverage_sub_score(0, -1)


@pytest.mark.unit
class TestConservationSubScore:
    """Tests for the cross-species conservation sub-score helper."""

    def test_full_conservation_yields_one(self) -> None:
        """Ortholog hits in all requested species should yield 1.0."""
        assert conservation_sub_score(3, 3) == 1.0

    def test_partial_conservation_yields_fraction(self) -> None:
        """Partial conservation should yield the exact hit/requested fraction."""
        assert conservation_sub_score(2, 3) == pytest.approx(2 / 3)
        assert conservation_sub_score(1, 4) == 0.25
        assert conservation_sub_score(5, 10) == 0.5

    def test_zero_hits_yields_zero(self) -> None:
        """Zero ortholog hits with nonzero requested should yield 0.0."""
        assert conservation_sub_score(0, 3) == 0.0

    def test_zero_requested_returns_none(self) -> None:
        """Zero requested non-query species should return None (inactive term)."""
        # A single-species run must not be penalised for lacking ortholog evidence.
        assert conservation_sub_score(0, 0) is None

    def test_hit_exceeds_requested_raises(self) -> None:
        """Ortholog hit count exceeding requested species should raise ValueError."""
        with pytest.raises(ValueError, match="cannot exceed"):
            conservation_sub_score(4, 3)

    def test_negative_hit_raises(self) -> None:
        """Negative ortholog hit count should raise ValueError."""
        with pytest.raises(ValueError, match="must be non-negative"):
            conservation_sub_score(-1, 3)

    def test_negative_requested_raises(self) -> None:
        """Negative requested count should raise ValueError."""
        with pytest.raises(ValueError, match="must be non-negative"):
            conservation_sub_score(0, -1)


@pytest.mark.unit
class TestComputeComposite:
    """Tests for the composite scorer with renormalization."""

    def test_known_features_and_weights_yield_known_score(self) -> None:
        """Known feature mapping and known weights should yield a hand-computed score.

        This is the anchor test: the arithmetic is written out by hand, not by
        calling the function twice. If the function changes to return a constant,
        this test will fail.
        """
        # Use default weights from ScoringWeights:
        # asymmetry=0.12, gc_content=0.10, accessibility=0.13, empirical=0.15,
        # off_target=0.25, isoform_coverage=0.15, conservation=0.10
        # Sum = 1.00

        features = {
            "asymmetry": 0.8,
            "gc_content": 0.9,
            "accessibility": 0.7,
            "empirical": 0.6,
            "off_target": 1.0,
            "isoform_coverage": 0.5,
            "conservation": 0.4,
        }
        weights = ScoringWeights()

        # Hand-computed score: sum(weight * feature * 100) for all terms.
        # asymmetry:        0.12 * 0.8 * 100 = 9.6
        # gc_content:       0.10 * 0.9 * 100 = 9.0
        # accessibility:    0.13 * 0.7 * 100 = 9.1
        # empirical:        0.15 * 0.6 * 100 = 9.0
        # off_target:       0.25 * 1.0 * 100 = 25.0
        # isoform_coverage: 0.15 * 0.5 * 100 = 7.5
        # conservation:     0.10 * 0.4 * 100 = 4.0
        # Total = 73.2
        expected_score = 73.2

        result = compute_composite(features, weights)

        assert result.score == pytest.approx(expected_score, abs=1e-9)
        assert result.weight_set_version == SCORING_WEIGHT_SET_VERSION
        assert result.active_terms == COMPOSITE_TERMS
        assert len(result.contributions) == 7

    def test_contributions_sum_to_score_full_active_set(self) -> None:
        """Contributions should sum to the composite score for a full active set."""
        features = {
            "asymmetry": 0.8,
            "gc_content": 0.9,
            "accessibility": 0.7,
            "empirical": 0.6,
            "off_target": 1.0,
            "isoform_coverage": 0.5,
            "conservation": 0.4,
        }
        weights = ScoringWeights()
        result = compute_composite(features, weights)

        contribution_sum = sum(result.contributions.values())
        assert contribution_sum == pytest.approx(result.score, abs=1e-9)

    def test_contributions_sum_to_score_partial_active_set(self) -> None:
        """Contributions should sum to the composite score for a partial active set."""
        # Only five of seven terms active (no isoform_coverage, no conservation).
        features = {
            "asymmetry": 0.8,
            "gc_content": 0.9,
            "accessibility": 0.7,
            "empirical": 0.6,
            "off_target": 1.0,
        }
        weights = ScoringWeights()
        result = compute_composite(features, weights)

        contribution_sum = sum(result.contributions.values())
        assert contribution_sum == pytest.approx(result.score, abs=1e-9)

        # Only active terms should appear in contributions.
        assert len(result.contributions) == 5
        assert "isoform_coverage" not in result.contributions
        assert "conservation" not in result.contributions

    def test_inactive_terms_renormalize_to_preserve_scale(self) -> None:
        """Dropping terms should not change the composite of a candidate with all equal sub-scores.

        When all active sub-scores are equal, renormalization ensures the composite
        stays the same regardless of which terms are active. This is the "neither
        rewarded nor penalised" property.
        """
        # All sub-scores = 1.0, full active set -> should score 100.0.
        features_full = dict.fromkeys(COMPOSITE_TERMS, 1.0)
        weights = ScoringWeights()
        result_full = compute_composite(features_full, weights)
        assert result_full.score == pytest.approx(100.0, abs=1e-9)

        # All sub-scores = 1.0, partial active set (drop two terms) -> should still score 100.0.
        features_partial = {
            "asymmetry": 1.0,
            "gc_content": 1.0,
            "accessibility": 1.0,
            "empirical": 1.0,
            "off_target": 1.0,
            # isoform_coverage and conservation absent
        }
        result_partial = compute_composite(features_partial, weights)
        assert result_partial.score == pytest.approx(100.0, abs=1e-9)

        # All sub-scores = 0.5, full active set -> should score 50.0.
        features_half_full = dict.fromkeys(COMPOSITE_TERMS, 0.5)
        result_half_full = compute_composite(features_half_full, weights)
        assert result_half_full.score == pytest.approx(50.0, abs=1e-9)

        # All sub-scores = 0.5, partial active set -> should still score 50.0.
        features_half_partial = {
            "asymmetry": 0.5,
            "gc_content": 0.5,
            "accessibility": 0.5,
            "empirical": 0.5,
            "off_target": 0.5,
        }
        result_half_partial = compute_composite(features_half_partial, weights)
        assert result_half_partial.score == pytest.approx(50.0, abs=1e-9)

    def test_empty_features_raises(self) -> None:
        """Empty features mapping should raise ScoringError (no active terms)."""
        features = {}
        weights = ScoringWeights()
        with pytest.raises(ScoringError, match="No active terms"):
            compute_composite(features, weights)

    def test_all_zero_active_weights_raises(self) -> None:
        """All-zero active weight vector should raise ScoringError at call time."""
        # Construct a weight set with one nonzero term to pass the ScoringWeights
        # validator, then call compute_composite with features that exclude it.
        weights = ScoringWeights(
            asymmetry=0.0,
            gc_content=0.0,
            accessibility=0.0,
            empirical=0.0,
            off_target=1.0,  # Only this is nonzero, so the full vector sums to 1.0.
            isoform_coverage=0.0,
            conservation=0.0,
        )
        # Features that exclude off_target -> active weight sum is 0.
        features = {
            "asymmetry": 0.5,
            "gc_content": 0.5,
        }
        with pytest.raises(ScoringError, match="non-positive, cannot normalize"):
            compute_composite(features, weights)

    def test_feature_below_zero_raises(self) -> None:
        """Feature value below 0.0 should raise ScoringError."""
        features = {
            "asymmetry": -0.1,  # Invalid
            "gc_content": 0.9,
        }
        weights = ScoringWeights()
        with pytest.raises(ScoringError, match="outside \\[0, 1\\]"):
            compute_composite(features, weights)

    def test_feature_above_one_raises(self) -> None:
        """Feature value above 1.0 should raise ScoringError."""
        features = {
            "asymmetry": 0.8,
            "gc_content": 1.1,  # Invalid
        }
        weights = ScoringWeights()
        with pytest.raises(ScoringError, match="outside \\[0, 1\\]"):
            compute_composite(features, weights)

    def test_explicit_active_terms_subset_is_respected(self) -> None:
        """Explicitly passing active_terms should restrict the active set."""
        features = {
            "asymmetry": 0.8,
            "gc_content": 0.9,
            "accessibility": 0.7,
            "empirical": 0.6,
            "off_target": 1.0,
        }
        weights = ScoringWeights()

        # Request only two terms be active.
        result = compute_composite(features, weights, active_terms=["asymmetry", "gc_content"])

        assert len(result.active_terms) == 2
        assert result.active_terms == ("asymmetry", "gc_content")
        assert len(result.contributions) == 2

        # Renormalized weights: asymmetry=0.12, gc_content=0.10 -> sum=0.22
        # Renorm: asymmetry=0.12/0.22, gc_content=0.10/0.22
        # Score = (0.12/0.22)*0.8*100 + (0.10/0.22)*0.9*100
        renorm_asym = 0.12 / 0.22
        renorm_gc = 0.10 / 0.22
        expected_score = renorm_asym * 0.8 * 100 + renorm_gc * 0.9 * 100
        assert result.score == pytest.approx(expected_score, abs=1e-9)

    def test_unknown_term_in_active_terms_is_silently_ignored(self) -> None:
        """Unknown term names in active_terms should be silently ignored."""
        features = {
            "asymmetry": 0.8,
            "gc_content": 0.9,
        }
        weights = ScoringWeights()

        # Request a nonexistent term alongside valid ones.
        result = compute_composite(features, weights, active_terms=["asymmetry", "unknown_term", "gc_content"])

        # Only the known terms present in features should be active.
        assert result.active_terms == ("asymmetry", "gc_content")
        assert len(result.contributions) == 2


@pytest.mark.unit
class TestScoringWeightsValidator:
    """Tests for the ScoringWeights model validator."""

    def test_default_weights_sum_to_one(self) -> None:
        """Default weights should sum to exactly 1.0 as specified in the contract.

        This test asserts the actual default values term by term so a silent retune
        fails the suite. If you change the defaults, update this test AND bump
        SCORING_WEIGHT_SET_VERSION.
        """
        weights = ScoringWeights()
        assert weights.asymmetry == 0.12
        assert weights.gc_content == 0.10
        assert weights.accessibility == 0.13
        assert weights.empirical == 0.15
        assert weights.off_target == 0.25
        assert weights.isoform_coverage == 0.15
        assert weights.conservation == 0.10

        # Sum to 1.00 exactly.
        total = sum(getattr(weights, term) for term in COMPOSITE_TERMS)
        assert total == pytest.approx(1.0, abs=1e-9)

    def test_unnormalizable_weight_set_raises_at_construction(self) -> None:
        """A weight set summing far from 1.0 should raise ValueError at construction."""
        with pytest.raises(ValueError, match="must sum to 1.0"):
            ScoringWeights(
                asymmetry=0.9,
                gc_content=0.0,
                accessibility=0.0,
                empirical=0.0,
                off_target=0.0,
                isoform_coverage=0.0,
                conservation=0.0,
            )

    def test_custom_weights_summing_to_one_are_accepted(self) -> None:
        """Custom weights summing to 1.0 should pass validation."""
        weights = ScoringWeights(
            asymmetry=0.2,
            gc_content=0.2,
            accessibility=0.2,
            empirical=0.2,
            off_target=0.1,
            isoform_coverage=0.05,
            conservation=0.05,
        )
        assert weights.asymmetry == 0.2


@pytest.mark.unit
class TestVersionConstant:
    """Tests for the SCORING_WEIGHT_SET_VERSION constant."""

    def test_version_is_2_0_0(self) -> None:
        """SCORING_WEIGHT_SET_VERSION should be "2.0.0".

        If you change the default weights in ScoringWeights, you MUST bump this
        version. This test exists to catch silent retuning without version bumps.
        """
        assert SCORING_WEIGHT_SET_VERSION == "2.0.0"

    def test_composite_score_records_version(self) -> None:
        """CompositeScore should record the weight set version."""
        features = {"asymmetry": 0.8, "gc_content": 0.9}
        weights = ScoringWeights()
        result = compute_composite(features, weights)
        assert result.weight_set_version == "2.0.0"
