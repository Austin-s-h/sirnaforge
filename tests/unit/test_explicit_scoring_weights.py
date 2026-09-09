"""Issue #96: every weight is declared, and nothing rescales one at runtime.

The two defects removed here were both *hidden* normalisations -- arithmetic applied to a declared
weight vector between the manifest and the score. Each test below is written to fail if either
returns in any form, including a new one:

- the scorer renormalised over whichever terms were populated, doubling every design-stage weight;
- miRNA mode divided every score by ``1 + max_mirna_bonus`` (0.80x), which is absent from the
  manifest and which fought the empirical term over guide position 1.

The `no runtime arithmetic` test is a source-level guard rather than a behavioural one, because the
behavioural symptom of renormalisation is precisely that everything still looks plausible.
"""

import ast
import inspect
from pathlib import Path

import pytest

from sirnaforge.core import design as design_module
from sirnaforge.core import scoring as scoring_module
from sirnaforge.core.design import SiRNADesigner, biogenesis_features
from sirnaforge.core.scoring import compute_composite
from sirnaforge.models.sirna import (
    EMPIRICAL_SCORE_MAX,
    EMPIRICAL_SCORE_MIN,
    DesignMode,
    DesignParameters,
    FilterCriteria,
    ScoringWeights,
    SiRNACandidate,
)

SHARED_TERMS = ("off_target", "target_accessibility", "asymmetry", "gc_content")

# The old divisor, kept here as a number so the miRNA tests can state what they replaced.
OLD_MIRNA_DIVISOR = 1.25


def _candidate(guide: str) -> SiRNACandidate:
    """An unscored candidate whose passenger is the exact reverse complement of the guide."""
    passenger = guide.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    gc = (guide.count("G") + guide.count("C")) / len(guide) * 100
    return SiRNACandidate(
        id=f"SIRNAF_TP53_1_{len(guide)}",
        transcript_id="ENST00000269305",
        position=1,
        guide_sequence=guide,
        passenger_sequence=passenger,
        gc_content=gc,
        length=len(guide),
        asymmetry_score=0.7,
    )


@pytest.mark.unit
def test_mirna_and_sirna_composites_are_on_one_declared_scale() -> None:
    """A candidate as good at biogenesis as at everything else scores the same in both modes.

    This is the identity the 1.25 divisor broke. Both post-screen vectors sum to 1.0, so for a
    candidate whose sub-scores are all the same value `f`, both modes return exactly `100 * f` --
    for every `f`, not just at the endpoints. Under the divisor the same candidate returned
    `(100f + 0.25f * 100) / 1.25`, which only coincidentally agrees at f = 0.
    """
    weights = ScoringWeights()
    sirna_vector = weights.vector_for(post_screen=True, design_mode=DesignMode.SIRNA)
    mirna_vector = weights.vector_for(post_screen=True, design_mode=DesignMode.MIRNA)

    for level in (0.0, 0.25, 0.5, 0.75, 1.0):
        sirna = compute_composite(dict.fromkeys(sirna_vector.terms, level), sirna_vector)
        mirna = compute_composite(dict.fromkeys(mirna_vector.terms, level), mirna_vector)

        assert sirna.score == pytest.approx(100.0 * level, abs=1e-9)
        assert mirna.score == pytest.approx(sirna.score, abs=1e-9), (
            f"the two modes disagree at level {level}: {mirna.score} vs {sirna.score}"
        )


@pytest.mark.unit
def test_no_divisor_is_applied_after_the_miRNA_composite() -> None:
    """The miRNA score is exactly the sum of its declared contributions, undivided.

    The single most direct statement of the fix: with the divisor present the score was the sum
    over 1.25, so this equality failed by 20% on every miRNA candidate.
    """
    vector = ScoringWeights().vector_for(post_screen=True, design_mode=DesignMode.MIRNA)
    features = {
        "off_target": 1.0,
        "target_accessibility": 0.6,
        "asymmetry": 0.8,
        "gc_content": 0.5,
        "ago_start": 1.0,
        "pos1_mismatch": 0.0,
        "supp_13_16": 0.75,
    }

    result = compute_composite(features, vector)
    by_hand = sum(vector.as_mapping()[term] * features[term] * 100.0 for term in vector.terms)

    assert result.score == pytest.approx(by_hand, abs=1e-9)
    assert result.score == pytest.approx(sum(result.contributions.values()), abs=1e-9)
    assert result.score != pytest.approx(by_hand / OLD_MIRNA_DIVISOR)


@pytest.mark.unit
def test_a_zero_biogenesis_candidate_loses_only_declared_weight() -> None:
    """What the literal "same score" reading costs, stated as a measured number.

    A candidate scoring 0 on all three biogenesis terms cannot score identically in both modes:
    the three terms carry 0.25 of postscreen_mirna_v4, so it retains the shared four's 0.75. That
    is arithmetically similar to the old 0.80 factor and deliberately different in kind -- 0.75 is
    the sum of three weights written in the manifest, attributable term by term on the row, where
    0.80 was an undeclared divisor applied to the whole vector including off_target.
    """
    weights = ScoringWeights()
    sirna_vector = weights.vector_for(post_screen=True, design_mode=DesignMode.SIRNA)
    mirna_vector = weights.vector_for(post_screen=True, design_mode=DesignMode.MIRNA)

    shared = dict.fromkeys(SHARED_TERMS, 1.0)
    sirna = compute_composite(shared, sirna_vector)
    mirna = compute_composite({**shared, "ago_start": 0.0, "pos1_mismatch": 0.0, "supp_13_16": 0.0}, mirna_vector)

    retained = mirna.score / sirna.score
    assert retained == pytest.approx(0.75, abs=1e-9), f"zero-biogenesis candidates retain {retained:.4f}"

    # The shortfall is exactly the three declared biogenesis weights, not a factor on the vector.
    shortfall = (sirna.score - mirna.score) / 100.0
    assert shortfall == pytest.approx(
        sum(mirna_vector.as_mapping()[term] for term in ("ago_start", "pos1_mismatch", "supp_13_16")),
        abs=1e-9,
    )


@pytest.mark.unit
def test_no_two_guide_position_rules_disagree_on_the_same_base() -> None:
    """Guide position 1 must be judged by one rule, not two pulling opposite ways.

    The empirical rubric paid +0.1 for G/C at guide position 1 while the biogenesis rule paid
    ago_start for A/U at the same base. Measured over 29,605 candidates, G/C gained +1.6 empirical
    points and lost 7.9 to the biogenesis adjustment -- a declared 0.15-weight term overridden ~5x
    by an undeclared one, which is why `empirical` had a negative variance share. A/U wins; the
    empirical rubric must now be blind to position 1.
    """
    designer = SiRNADesigner(DesignParameters())
    middle = "CGTACGTACGTACGTAC"  # 17 nt between guide position 1 and position 19

    for last in "ACGT":
        scores = {designer._calculate_empirical_score(_candidate(f"{first}{middle}{last}")) for first in "ACGT"}
        assert len(scores) == 1, f"the empirical score still depends on guide position 1 (last={last}): {scores}"

    # And the surviving rule at that base does depend on it, in the A/U direction.
    au_first = biogenesis_features(f"A{middle}A", "T" * 19)["ago_start"]
    gc_first = biogenesis_features(f"G{middle}A", "T" * 19)["ago_start"]
    assert (au_first, gc_first) == (1.0, 0.0)


@pytest.mark.unit
def test_empirical_rubric_attains_exactly_the_declared_bounds() -> None:
    """With the position-1 clause gone the rubric attains {0.4, 0.5, 0.6}, so the max must be 0.6.

    EMPIRICAL_SCORE_MAX is the `le=` bound on min_empirical_score, so leaving it at 0.7 would let a
    user configure a threshold no candidate can reach.
    """
    designer = SiRNADesigner(DesignParameters())
    middle = "CGTACGTACGTACGTAC"
    attained = {
        designer._calculate_empirical_score(_candidate(f"{first}{middle}{last}")) for first in "ACGT" for last in "ACGT"
    }

    assert sorted(attained) == pytest.approx([0.4, 0.5, 0.6])
    assert EMPIRICAL_SCORE_MIN == 0.4
    assert EMPIRICAL_SCORE_MAX == 0.6
    assert FilterCriteria(min_empirical_score=EMPIRICAL_SCORE_MAX).min_empirical_score == 0.6


@pytest.mark.unit
def test_empirical_is_gate_only_and_not_a_scored_term() -> None:
    """It is computed and reported on every candidate; no vector reads it."""
    for vector in ScoringWeights().all_vectors():
        assert "empirical" not in vector.terms
        assert not hasattr(vector, "empirical")

    candidate = _candidate("GTAGATTACCACTGGAGTCAA")
    SiRNADesigner(DesignParameters())._score_candidates([candidate], "A" * 30 + candidate.passenger_sequence + "A" * 60)

    assert "empirical" in candidate.component_scores, "the empirical score must still be reported"
    assert candidate.design_score is not None
    assert not hasattr(candidate, "score_empirical"), "a gate has no contribution to the score"


@pytest.mark.unit
def test_conservation_and_isoform_coverage_are_reported_but_unscored() -> None:
    """Both left the composite; both are still fields on every candidate."""
    for vector in ScoringWeights().all_vectors():
        assert "conservation" not in vector.terms
        assert "isoform_coverage" not in vector.terms

    fields = set(SiRNACandidate.model_fields)
    assert {"conservation_score", "isoform_coverage"} <= fields
    assert not {"score_conservation", "score_isoform_coverage"} & fields


@pytest.mark.unit
def test_isoform_coverage_gate_defaults_off_and_bites_when_configured() -> None:
    """The gate is opt-in, and a coverage value that could not be computed never fails."""
    assert FilterCriteria().min_isoform_coverage is None
    assert FilterCriteria(min_isoform_coverage=0.5).min_isoform_coverage == 0.5

    # Bounds still apply because it is set through the constructor, not model_copy.
    with pytest.raises(ValueError, match="less than or equal to 1"):
        FilterCriteria(min_isoform_coverage=1.5)

    assert SiRNACandidate.FilterStatus.LOW_ISOFORM_COVERAGE.value == "LOW_ISOFORM_COVERAGE"


@pytest.mark.unit
def test_design_score_and_composite_score_are_separate_stages() -> None:
    """The design stage cannot produce a composite: off_target does not exist yet."""
    candidate = _candidate("GTAGATTACCACTGGAGTCAA")
    transcript = "A" * 30 + candidate.passenger_sequence + "A" * 60
    candidate.position = 31

    SiRNADesigner(DesignParameters())._score_candidates([candidate], transcript)

    assert candidate.design_score is not None
    assert candidate.composite_score is None, "composite_score must be null before screening"
    assert candidate.scored_after_screening is False
    assert candidate.weight_vector == "design_v4"
    assert candidate.score_off_target is None

    # design_score is exactly its three declared contributions -- nothing else is folded in.
    contributions = (
        candidate.score_target_accessibility,
        candidate.score_asymmetry,
        candidate.score_gc_content,
    )
    assert all(value is not None for value in contributions)
    assert candidate.design_score == pytest.approx(sum(contributions))  # type: ignore[arg-type]


def _weight_arithmetic_offences(module: object) -> list[str]:
    """Names of functions in `module` that divide or multiply a weight-shaped expression.

    Deliberately syntactic. Renormalisation is invisible in behaviour -- every score still lands in
    [0, 100] and still ranks plausibly -- so the guard has to look at the code. It flags any
    division whose operands mention a weight, and any `sum(` over something weight-shaped.
    """
    source = inspect.getsource(module)  # type: ignore[arg-type]
    tree = ast.parse(source)
    offences: list[str] = []

    def mentions_weight(node: ast.AST) -> bool:
        return any(
            "weight" in fragment.lower()
            for child in ast.walk(node)
            for fragment in (
                [child.id] if isinstance(child, ast.Name) else [child.attr] if isinstance(child, ast.Attribute) else []
            )
        )

    for node in ast.walk(tree):
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div) and mentions_weight(node):
            offences.append(ast.unparse(node))
        # A bare `sum(weights)` is the shape a normaliser takes; a generator expression is how
        # the validators legitimately total a vector to check it.
        is_weight_sum = (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "sum"
            and bool(node.args)
            and mentions_weight(node.args[0])
            and not isinstance(node.args[0], ast.GeneratorExp)
        )
        if is_weight_sum:
            offences.append(ast.unparse(node))
    return offences


@pytest.mark.unit
def test_no_module_divides_a_weight_vector_at_runtime() -> None:
    """Grep-verifiable form of the central rule, over the two modules that own scoring.

    The validators still sum a vector to check it -- that is the point, they refuse a vector that
    is not already normalised -- so this looks for division by a weight-shaped expression, which is
    what renormalisation and the miRNA divisor both were.
    """
    for module in (scoring_module, design_module):
        offences = _weight_arithmetic_offences(module)
        assert not offences, f"{module.__name__} performs runtime arithmetic on weights: {offences}"


@pytest.mark.unit
def test_the_deleted_normalisations_are_not_reachable_by_name() -> None:
    """The divisor helpers are gone, not deprecated: nothing can call them back into service."""
    assert not hasattr(design_module, "apply_mirna_biogenesis_bonus")
    assert not hasattr(design_module, "mirna_max_biogenesis_bonus")
    assert not hasattr(design_module, "MIRNA_BONUS_MAX_KEY")

    # And the two bonus weights that were declared but never read are gone from the config.
    from sirnaforge.models.sirna import MiRNADesignConfig  # noqa: PLC0415

    assert not hasattr(MiRNADesignConfig(), "scoring_weights")


@pytest.mark.unit
def test_source_mentions_no_renormalisation_of_the_active_term_set() -> None:
    """The words went with the behaviour: no live source claims weights are renormalised.

    scoring.py's module docstring used to document renormalisation as a deliberate feature, which
    would now contradict the code it introduces.
    """
    src = Path(scoring_module.__file__).resolve().parent.parent
    live = list(src.rglob("*.py"))
    assert live, "no source files found"

    # A mention is allowed only alongside a word that negates or historicises it. Prose wraps, so
    # the window is the line plus its two neighbours -- the sentence, not the line.
    negations = ("never", "not ", "no ", "used to", "would", "removed", "deleted", "gone", "rather than")

    for path in live:
        lines = path.read_text().splitlines()
        for index, line in enumerate(lines):
            if "renormalis" not in line.lower() and "renormaliz" not in line.lower():
                continue
            window = " ".join(lines[max(0, index - 1) : index + 2]).lower()
            assert any(marker in window for marker in negations), (
                f"{path.name}:{index + 1} describes renormalisation as current behaviour: {line.strip()}"
            )
