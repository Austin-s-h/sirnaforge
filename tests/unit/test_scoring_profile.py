"""Issue #102: the term registry, and the invariants that make it worth having.

The registry's value is not that it stores prose -- it is that three claims are now checkable by a
machine rather than by a careful reader:

- a term paid weight in any profile has a record describing it;
- a term whose attainable range is a single point cannot be paid weight at all (`pos1_mismatch`
  held 0.05 of `postscreen_mirna_v4` while being exactly constant at 0.0 on 13,415 scored rows);
- no profile claims `validated`, because nothing in siRNAforge has been held out and independently
  replicated (D1).

The panel measurements themselves are not tested here -- they need untracked third-party data and
live in `scripts/validate_scoring_profiles.py`. What is tested here is deterministic.
"""

from typing import ClassVar

import pytest
from pydantic import ValidationError

from sirnaforge.core.design import SiRNADesigner, au_content_5p_score, biogenesis_features
from sirnaforge.models.scoring_profile import (
    EXPERIMENTAL_AU_PROFILE,
    PROFILES,
    SHIPPED_PROFILE,
    TERM_REGISTRY,
    EvidenceStatus,
    ScoringProfile,
    TermRecord,
)
from sirnaforge.models.sirna import (
    COMPOSITE_TERM_NAMES,
    DesignParameters,
    PostScreenMiRNAWeights,
    ScoringWeights,
    SiRNACandidate,
)


def _record(**overrides: object) -> TermRecord:
    """A minimal valid TermRecord, so a test can vary exactly one field."""
    fields: dict[str, object] = {
        "term": "probe",
        "molecule": "guide",
        "strand": "guide 5'->3'",
        "positions": "1",
        "endpoint": "efficacy",
        "applicability": "any guide",
        "formula": "1.0 if the base is A else 0.0",
        "units": "indicator",
        "transform": "none",
        "missing_value_policy": "not reachable",
        "declared_range": (0.0, 1.0),
        "attainable_range": (0.0, 1.0),
        "evidence_source": "test fixture",
        "evidence_status": "experimental",
    }
    fields.update(overrides)
    return TermRecord(**fields)  # type: ignore[arg-type]


@pytest.mark.unit
def test_every_scored_term_of_every_shipped_vector_is_registered() -> None:
    """A weight with no record behind it is exactly what #102 was filed about."""
    for vector in ScoringWeights().all_vectors():
        for term in vector.terms:
            assert term in TERM_REGISTRY, f"{vector.name} scores '{term}' with no TermRecord"

    # And the union tuple used for column ordering agrees with the vectors it summarises.
    scored = {term for vector in ScoringWeights().all_vectors() for term in vector.terms}
    assert scored == set(COMPOSITE_TERM_NAMES)


@pytest.mark.unit
def test_pos1_mismatch_holds_no_weight_and_is_still_reported() -> None:
    """The audit's disposition, pinned: out of the vector, still computed.

    Its range is unattainable by construction rather than by accident -- the passenger is the exact
    reverse complement of the guide, so guide position 1 always forms a Watson-Crick pair.
    """
    assert "pos1_mismatch" not in PostScreenMiRNAWeights().TERM_NAMES
    assert not hasattr(PostScreenMiRNAWeights(), "pos1_mismatch")
    for vector in ScoringWeights().all_vectors():
        assert "pos1_mismatch" not in vector.terms

    assert TERM_REGISTRY["pos1_mismatch"].evidence_status is EvidenceStatus.DEPRECATED
    assert TERM_REGISTRY["pos1_mismatch"].is_constant

    # Still computed and still reported, so the removal is auditable and reversible.
    guide = "ATCGATCGATCGATCGATCGA"
    passenger = guide.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    assert biogenesis_features(guide, passenger)["pos1_mismatch"] == 0.0
    assert "score_pos1_mismatch" in SiRNACandidate.model_fields


@pytest.mark.unit
def test_an_exact_reverse_complement_passenger_never_earns_pos1_mismatch() -> None:
    """The structural argument, over every possible first base rather than one example."""
    for first in "ACGT":
        guide = f"{first}TCGATCGATCGATCGATCGA"
        passenger = guide.translate(str.maketrans("ATCG", "TAGC"))[::-1]
        assert biogenesis_features(guide, passenger)["pos1_mismatch"] == 0.0, (
            f"guide starting {first} earned pos1_mismatch against an exact reverse complement"
        )


@pytest.mark.unit
def test_a_profile_refuses_to_pay_weight_to_a_constant_term() -> None:
    """The machine-checked form of "a constant term must not consume weight"."""
    with pytest.raises(ValidationError, match="constant term ranks nothing"):
        ScoringProfile(
            profile_id="regression_probe",
            description="pays the constant term, as postscreen_mirna_v4 used to",
            status=EvidenceStatus.EXPERIMENTAL,
            ships_as_default=False,
            vectors=(_ConstantPayingWeights(),),
        )


@pytest.mark.unit
def test_a_profile_refuses_an_unregistered_term_and_the_validated_claim() -> None:
    """Two more things the registry will not let a profile say."""
    with pytest.raises(ValidationError, match="no TermRecord"):
        ScoringProfile(
            profile_id="unregistered_probe",
            description="scores a term nothing describes",
            status=EvidenceStatus.EXPERIMENTAL,
            ships_as_default=False,
            vectors=(_UnregisteredTermWeights(),),
        )

    with pytest.raises(ValidationError, match="claims VALIDATED"):
        ScoringProfile(
            profile_id="promotion_probe",
            description="claims a promotion that has not happened",
            status=EvidenceStatus.VALIDATED,
            ships_as_default=False,
            vectors=(),
        )


@pytest.mark.unit
def test_nothing_in_the_registry_is_promoted() -> None:
    """D1: 0.7.1 ships registry and audits, and promotes no weight or threshold."""
    assert all(record.evidence_status is not EvidenceStatus.VALIDATED for record in TERM_REGISTRY.values())
    assert all(profile.status is EvidenceStatus.EXPERIMENTAL for profile in PROFILES.values())
    assert SHIPPED_PROFILE.ships_as_default is True
    assert EXPERIMENTAL_AU_PROFILE.ships_as_default is False


@pytest.mark.unit
def test_the_shipped_profile_is_the_actual_default_vectors() -> None:
    """A registry that drifted from the code it documents would be worse than none."""
    declared = {vector.name: vector.as_mapping() for vector in SHIPPED_PROFILE.vectors}
    actual = ScoringWeights().as_manifest()
    assert declared == actual


@pytest.mark.unit
def test_au_1_5_is_registered_reported_and_scored_by_no_default_vector() -> None:
    """D5 settled the term and its window; D1 keeps it out of the shipped vectors."""
    record = TERM_REGISTRY["au_1_5"]
    assert record.positions == "1-5"
    assert record.evidence_status is EvidenceStatus.EXPERIMENTAL
    for vector in ScoringWeights().all_vectors():
        assert "au_1_5" not in vector.terms

    # It is in the experimental profile, alongside asymmetry -- D5's residual test kept both.
    experimental = {term for vector in EXPERIMENTAL_AU_PROFILE.vectors for term in vector.terms}
    assert {"au_1_5", "asymmetry"} <= experimental
    # And ago_start is absent there, because position 1 is inside 1-5 (#97 question 5).
    assert "ago_start" not in experimental


@pytest.mark.unit
def test_au_1_5_counts_the_declared_window_read_as_rna() -> None:
    """Guides are stored as DNA, so T must count as U or the term reads a different molecule."""
    assert au_content_5p_score("AAAAACGCGCGCGCGCGCGCG") == 1.0
    assert au_content_5p_score("TTTTTCGCGCGCGCGCGCGCG") == 1.0
    assert au_content_5p_score("GCGCGAAAAACGCGCGCGCGC") == 0.0
    assert au_content_5p_score("ATGCGCGCGCGCGCGCGCGCG") == pytest.approx(0.4)
    # Positions beyond 5 are outside the pre-declared window and must not move the score.
    assert au_content_5p_score("GGGGGAAAAAAAAAAAAAAAA") == 0.0
    # Too short to fill the window: None, not a substituted midpoint.
    assert au_content_5p_score("ATGC") is None


@pytest.mark.unit
def test_the_design_path_reports_au_1_5_without_scoring_it() -> None:
    """It has to be on the row to be auditable, and off every vector to stay unpromoted."""
    guide = "ATTGACTCCAGTGGTAATCT"
    passenger = guide.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    candidate = SiRNACandidate(
        id="au_probe",
        transcript_id="ENST00000269305",
        position=31,
        guide_sequence=guide,
        passenger_sequence=passenger,
        gc_content=(guide.count("G") + guide.count("C")) / len(guide) * 100,
        length=len(guide),
        asymmetry_score=0.7,
    )
    transcript = "A" * 30 + passenger + "A" * 60
    SiRNADesigner(DesignParameters())._score_candidates([candidate], transcript)

    assert candidate.component_scores["au_1_5"] == pytest.approx(au_content_5p_score(guide))
    assert candidate.design_score is not None
    # No contribution column, because no vector scores it.
    assert not hasattr(candidate, "score_au_1_5")
    assert "au_1_5" not in SiRNACandidate.model_fields


@pytest.mark.unit
def test_a_term_record_rejects_an_attainable_range_outside_its_declared_range() -> None:
    """The two ranges mean different things, so the narrower one must actually be narrower."""
    with pytest.raises(ValidationError, match="escapes declared_range"):
        _record(declared_range=(0.0, 1.0), attainable_range=(0.0, 1.5))
    with pytest.raises(ValidationError, match="is inverted"):
        _record(attainable_range=(1.0, 0.0))
    # A single point is legal to *record* -- it is only illegal to pay weight to.
    assert _record(attainable_range=(0.0, 0.0)).is_constant


@pytest.mark.unit
def test_every_registry_record_states_a_missing_value_policy_and_an_applicability() -> None:
    """The two fields that stop a term being quietly scored on a substituted value."""
    for term, record in TERM_REGISTRY.items():
        assert record.missing_value_policy.strip(), f"{term} declares no missing-value policy"
        assert record.applicability.strip(), f"{term} declares no applicability condition"
        assert record.evidence_source.strip(), f"{term} cites no evidence"
        assert record.term == term


# Probe vectors for the validator tests. They are deliberately invalid *as profiles*, not as
# vectors: each sums to 1.0, so it is the profile's registry check that rejects them.
class _ConstantPayingWeights(PostScreenMiRNAWeights):
    """Pays weight to `pos1_mismatch`, exactly as 4.0.0 did before #102."""

    VECTOR_NAME: ClassVar[str] = "probe_constant_paying"
    TERM_NAMES: ClassVar[tuple[str, ...]] = (
        "off_target",
        "target_accessibility",
        "asymmetry",
        "gc_content",
        "ago_start",
        "supp_13_16",
        "pos1_mismatch",
    )

    pos1_mismatch: float = 0.05
    target_accessibility: float = 0.19


class _UnregisteredTermWeights(PostScreenMiRNAWeights):
    """Pays weight to a term the registry has never heard of."""

    VECTOR_NAME: ClassVar[str] = "probe_unregistered"
    TERM_NAMES: ClassVar[tuple[str, ...]] = (
        "off_target",
        "target_accessibility",
        "asymmetry",
        "gc_content",
        "ago_start",
        "supp_13_16",
        "vibes",
    )

    vibes: float = 0.05
    target_accessibility: float = 0.19
