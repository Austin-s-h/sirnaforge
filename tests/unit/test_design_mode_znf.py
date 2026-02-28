"""Unit tests for ZNF design mode parsing and mutation-constraint wiring."""

import pytest

from sirnaforge.cli import _parse_znf_mutation_constraints, _resolve_design_mode
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    ZNFDefaultSubfingerMutationConstraint,
    ZNFMutationType,
    ZNFOverallMutationConstraint,
    ZNFSubfingerMutationConstraint,
)


def test_resolve_design_mode_accepts_znf() -> None:
    """ZNF should be accepted as a valid design mode."""
    mode, gc_min, gc_max, overhang, modification_pattern = _resolve_design_mode("znf", 30.0, 60.0, "dTdT", "none")
    assert mode == DesignMode.ZNF
    assert gc_min == 30.0
    assert gc_max == 60.0
    assert overhang == "dTdT"
    assert modification_pattern == "none"


def test_parse_znf_mutation_constraints_valid() -> None:
    """Valid CLI mutation constraints should parse into typed models."""
    constraints, default_constraint, overall_constraints = _parse_znf_mutation_constraints(
        ["1:2:substitution,transition", "2:0:deletion", "*:1:mismatch", "overall:3:substitution"]
    )
    assert len(constraints) == 2
    assert default_constraint == ZNFDefaultSubfingerMutationConstraint(
        max_mutations=1,
        mutation_types=[ZNFMutationType.SUBSTITUTION],
    )
    assert overall_constraints == [
        ZNFOverallMutationConstraint(
            max_mutations=3,
            mutation_types=[ZNFMutationType.SUBSTITUTION],
        )
    ]
    assert constraints[0] == ZNFSubfingerMutationConstraint(
        subfinger_index=1,
        max_mutations=2,
        mutation_types=[ZNFMutationType.SUBSTITUTION, ZNFMutationType.TRANSITION],
    )
    assert constraints[1].mutation_types == [ZNFMutationType.DELETION]


def test_parse_znf_mutation_constraints_invalid_format() -> None:
    """Invalid raw entries should raise ValueError."""
    with pytest.raises(ValueError):
        _parse_znf_mutation_constraints(["1:2"])


def test_design_parameters_store_znf_constraints() -> None:
    """DesignParameters should persist per-sub-finger mutation allowances."""
    params = DesignParameters(
        design_mode=DesignMode.ZNF,
        znf_subfinger_mutations=[
            ZNFSubfingerMutationConstraint(
                subfinger_index=3,
                max_mutations=1,
                mutation_types=[ZNFMutationType.TRANSVERSION],
            )
        ],
        znf_default_subfinger_mutation=ZNFDefaultSubfingerMutationConstraint(
            max_mutations=1,
            mutation_types=[ZNFMutationType.SUBSTITUTION],
        ),
        znf_overall_mutations=[
            ZNFOverallMutationConstraint(
                max_mutations=3,
                mutation_types=[ZNFMutationType.SUBSTITUTION],
            )
        ],
    )
    assert params.design_mode == DesignMode.ZNF
    assert params.znf_subfinger_mutations[0].subfinger_index == 3
    assert params.znf_default_subfinger_mutation is not None
    assert params.znf_overall_mutations[0].max_mutations == 3
