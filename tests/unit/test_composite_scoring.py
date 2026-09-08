"""Unit tests for scoring against named, hand-authored weight vectors (issue #96).

The scorer applies no arithmetic to weights: no renormalisation, no rescaling, no division. The
tests below are written so that reintroducing any of those fails, rather than merely exercising
the current code path.
"""

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
from sirnaforge.models.sirna import (
    DesignMode,
    DesignWeights,
    PostScreenMiRNAWeights,
    PostScreenSiRNAWeights,
    ScoringWeights,
    WeightVector,
)


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
    """Tests for the scorer: one named vector, its exact term set, no arithmetic on weights."""

    def test_known_features_and_weights_yield_known_score(self) -> None:
        """A hand-computed score, term by term, against postscreen_sirna_v4.

        The anchor test: the arithmetic is written out by hand rather than by calling the function
        twice, so a scorer that returned a constant would fail here.
        """
        features = {
            "off_target": 1.0,
            "target_accessibility": 0.7,
            "asymmetry": 0.8,
            "gc_content": 0.9,
        }
        vector = PostScreenSiRNAWeights()

        # off_target:           0.25 * 1.0 * 100 = 25.0
        # target_accessibility: 0.30 * 0.7 * 100 = 21.0
        # asymmetry:            0.25 * 0.8 * 100 = 20.0
        # gc_content:           0.20 * 0.9 * 100 = 18.0
        # Total = 84.0
        result = compute_composite(features, vector)

        assert result.score == pytest.approx(84.0, abs=1e-9)
        assert result.weight_set_version == SCORING_WEIGHT_SET_VERSION
        assert result.vector_name == "postscreen_sirna_v4"
        assert result.terms == ("off_target", "target_accessibility", "asymmetry", "gc_content")
        assert len(result.contributions) == 4

    def test_design_vector_scores_its_own_three_terms(self) -> None:
        """design_v4 is a different vector over a different term set, hand-computed too."""
        features = {"target_accessibility": 0.5, "asymmetry": 1.0, "gc_content": 0.0}

        # 0.40 * 0.5 * 100 + 0.35 * 1.0 * 100 + 0.25 * 0.0 * 100 = 20.0 + 35.0 + 0.0
        result = compute_composite(features, DesignWeights())

        assert result.score == pytest.approx(55.0, abs=1e-9)
        assert result.vector_name == "design_v4"
        assert result.terms == ("target_accessibility", "asymmetry", "gc_content")

    def test_extra_features_are_ignored_not_scored(self) -> None:
        """component_scores carries diagnostics; only the vector's own terms may contribute."""
        vector = DesignWeights()
        lean = {"target_accessibility": 0.5, "asymmetry": 0.5, "gc_content": 0.5}
        fat = {**lean, "empirical": 1.0, "conservation": 1.0, "isoform_coverage": 1.0, "off_target": 1.0}

        assert compute_composite(fat, vector).score == pytest.approx(compute_composite(lean, vector).score)
        assert set(compute_composite(fat, vector).contributions) == set(vector.terms)

    def test_contributions_sum_to_score(self) -> None:
        """The score is exactly the sum of its declared contributions -- nothing is applied after."""
        features = {"off_target": 1.0, "target_accessibility": 0.7, "asymmetry": 0.8, "gc_content": 0.9}
        result = compute_composite(features, PostScreenSiRNAWeights())

        assert sum(result.contributions.values()) == pytest.approx(result.score, abs=1e-9)

    def test_a_missing_term_raises_instead_of_renormalising(self) -> None:
        """The defect this issue removes: a partial term set must not be scored at all.

        Renormalising over whichever terms happened to be populated doubled every design-stage
        weight, so the same nominal 0.30 was worth 0.30 or 0.60 depending on the run stage and
        nothing on the row recorded which. There is now no score to report instead.
        """
        vector = PostScreenSiRNAWeights()
        features = {"target_accessibility": 0.7, "asymmetry": 0.8, "gc_content": 0.9}

        with pytest.raises(ScoringError, match="off_target"):
            compute_composite(features, vector)

    def test_dropping_a_term_cannot_leave_the_score_unchanged(self) -> None:
        """The positive form of the same guarantee, stated as a value.

        Under renormalisation an all-equal candidate scored identically whichever terms were
        active -- that property is exactly what made two candidates incomparable. Here the smaller
        vector is a genuinely different number, and the larger one refuses the partial input.
        """
        equal = dict.fromkeys(("off_target", "target_accessibility", "asymmetry", "gc_content"), 0.5)
        post_screen = compute_composite(equal, PostScreenSiRNAWeights())
        design = compute_composite(equal, DesignWeights())

        # Both vectors sum to 1.0, so an all-0.5 candidate scores 50 on either -- by construction,
        # not by rescaling. What is refused is scoring one vector's terms with the other's weights.
        assert post_screen.score == pytest.approx(50.0, abs=1e-9)
        assert design.score == pytest.approx(50.0, abs=1e-9)
        with pytest.raises(ScoringError):
            compute_composite({"target_accessibility": 0.5}, DesignWeights())

    def test_empty_features_raises(self) -> None:
        """An empty feature mapping supplies no term, so every vector refuses it."""
        with pytest.raises(ScoringError, match="requires"):
            compute_composite({}, PostScreenSiRNAWeights())

    def test_feature_below_zero_raises(self) -> None:
        """Feature value below 0.0 should raise ScoringError."""
        features = {"target_accessibility": 0.5, "asymmetry": -0.1, "gc_content": 0.9}
        with pytest.raises(ScoringError, match="outside \\[0, 1\\]"):
            compute_composite(features, DesignWeights())

    def test_feature_above_one_raises(self) -> None:
        """Feature value above 1.0 should raise ScoringError."""
        features = {"target_accessibility": 0.5, "asymmetry": 0.8, "gc_content": 1.1}
        with pytest.raises(ScoringError, match="outside \\[0, 1\\]"):
            compute_composite(features, DesignWeights())

    def test_scorer_takes_no_active_term_argument(self) -> None:
        """There is no way to ask for a subset: that argument was the renormalisation hook."""
        with pytest.raises(TypeError):
            compute_composite(  # type: ignore[call-arg]
                {"target_accessibility": 0.5, "asymmetry": 0.5, "gc_content": 0.5},
                DesignWeights(),
                active_terms=["asymmetry"],
            )


@pytest.mark.unit
class TestWeightVectors:
    """Every vector is named, hand-authored, and sums to 1.0 -- enforced at construction."""

    def test_declared_default_weights(self) -> None:
        """Assert the defaults term by term so a silent retune fails the suite.

        If you change any of these, bump SCORING_WEIGHT_SET_VERSION: the numbers are declared
        expert priors, and a run scored under different ones is not comparable.
        """
        assert DesignWeights().as_mapping() == {
            "target_accessibility": 0.40,
            "asymmetry": 0.35,
            "gc_content": 0.25,
        }
        assert PostScreenSiRNAWeights().as_mapping() == {
            "off_target": 0.25,
            "target_accessibility": 0.30,
            "asymmetry": 0.25,
            "gc_content": 0.20,
        }
        assert PostScreenMiRNAWeights().as_mapping() == {
            "off_target": 0.20,
            "target_accessibility": 0.22,
            "asymmetry": 0.18,
            "gc_content": 0.15,
            "ago_start": 0.10,
            "pos1_mismatch": 0.05,
            "supp_13_16": 0.10,
        }

    def test_every_vector_sums_to_one_over_its_own_terms(self) -> None:
        """Each vector validates against its own term set; there is no global term tuple."""
        for vector in ScoringWeights().all_vectors():
            assert sum(vector.as_mapping().values()) == pytest.approx(1.0, abs=1e-9)
            assert set(vector.as_mapping()) == set(vector.terms)
            assert set(vector.terms) <= set(COMPOSITE_TERMS)

    def test_the_three_vectors_score_different_term_sets(self) -> None:
        """The reason the single global COMPOSITE_TERM_NAMES tuple could no longer validate."""
        sizes = {vector.name: len(vector.terms) for vector in ScoringWeights().all_vectors()}
        assert sizes == {"design_v4": 3, "postscreen_sirna_v4": 4, "postscreen_mirna_v4": 7}

    def test_a_mis_summed_vector_is_a_construction_error(self) -> None:
        """Not a silent pass, and not renormalised at use: refused outright."""
        with pytest.raises(ValueError, match="must sum to exactly 1.0"):
            DesignWeights(target_accessibility=0.9, asymmetry=0.35, gc_content=0.25)
        with pytest.raises(ValueError, match="must sum to exactly 1.0"):
            PostScreenSiRNAWeights(off_target=0.25, target_accessibility=0.30, asymmetry=0.25, gc_content=0.10)

    def test_the_old_wide_tolerance_is_gone(self) -> None:
        """0.95-1.05 used to pass, i.e. up to 5% of undeclared rescaling per run."""
        with pytest.raises(ValueError, match="must sum to exactly 1.0"):
            DesignWeights(target_accessibility=0.42, asymmetry=0.35, gc_content=0.25)

    def test_an_unnamed_vector_is_a_construction_error(self) -> None:
        """A vector whose name cannot be recorded cannot be traced from a row, so it is refused."""

        class UnnamedWeights(WeightVector):
            """A subclass that forgot VECTOR_NAME."""

            TERM_NAMES = ("asymmetry",)
            asymmetry: float = 1.0

        with pytest.raises(ValueError, match="no VECTOR_NAME"):
            UnnamedWeights()

    def test_a_vector_with_no_terms_is_a_construction_error(self) -> None:
        """A vector must name what it scores."""

        class TermlessWeights(WeightVector):
            """A subclass that forgot TERM_NAMES."""

            VECTOR_NAME = "termless"
            asymmetry: float = 1.0

        with pytest.raises(ValueError, match="no TERM_NAMES"):
            TermlessWeights()

    def test_custom_weights_summing_to_one_are_accepted(self) -> None:
        """Hand-authoring your own vector is supported -- it just has to sum to 1.0."""
        vector = DesignWeights(target_accessibility=0.5, asymmetry=0.3, gc_content=0.2)
        assert vector.as_mapping() == {"target_accessibility": 0.5, "asymmetry": 0.3, "gc_content": 0.2}
        assert vector.name == "design_v4"

    def test_vector_selection_is_by_stage_and_mode(self) -> None:
        """A vector is chosen, never combined."""
        weights = ScoringWeights()

        assert weights.vector_for(post_screen=False).name == "design_v4"
        assert weights.vector_for(post_screen=False, design_mode=DesignMode.MIRNA).name == "design_v4"
        assert weights.vector_for(post_screen=True, design_mode=DesignMode.SIRNA).name == "postscreen_sirna_v4"
        assert weights.vector_for(post_screen=True, design_mode=DesignMode.MIRNA).name == "postscreen_mirna_v4"

    def test_manifest_records_name_alongside_weights(self) -> None:
        """A row's weight_vector column must resolve to the numbers that produced it."""
        manifest = ScoringWeights().as_manifest()

        assert set(manifest) == {"design_v4", "postscreen_sirna_v4", "postscreen_mirna_v4"}
        assert manifest["postscreen_mirna_v4"]["ago_start"] == 0.10


@pytest.mark.unit
class TestVersionConstant:
    """Tests for the SCORING_WEIGHT_SET_VERSION constant."""

    def test_version_is_4_0_0(self) -> None:
        """SCORING_WEIGHT_SET_VERSION should be "4.0.0".

        Bump it whenever a default weight or a vector's term set changes. 4.0.0 marks issue #96:
        both hidden normalisations removed (the active-set renormalisation and the miRNA 1.25
        divisor), one flat vector replaced by three named ones, and empirical / conservation /
        isoform_coverage out of the composite. No 3.x score is comparable with a 4.x one.
        """
        assert SCORING_WEIGHT_SET_VERSION == "4.0.0"

    def test_composite_score_records_version_and_vector(self) -> None:
        """CompositeScore records both, so a row is traceable to the weights that made it."""
        features = {"target_accessibility": 0.5, "asymmetry": 0.8, "gc_content": 0.9}
        result = compute_composite(features, DesignWeights())

        assert result.weight_set_version == "4.0.0"
        assert result.vector_name == "design_v4"
