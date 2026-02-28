"""Unit tests for ZFN design mode parsing and mutation-constraint wiring."""

import pytest

from sirnaforge.cli import _parse_zfn_mutation_constraints, _resolve_design_mode
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    ZFNDefaultSubfingerMutationConstraint,
    ZFNMutationType,
    ZFNOverallMutationConstraint,
    ZFNSubfingerMutationConstraint,
)


def test_resolve_design_mode_accepts_zfn() -> None:
    """ZFN should be accepted as a valid design mode."""
    mode, gc_min, gc_max, overhang, modification_pattern = _resolve_design_mode("zfn", 30.0, 60.0, "dTdT", "none")
    assert mode == DesignMode.ZFN
    assert gc_min == 30.0
    assert gc_max == 60.0
    assert overhang == "dTdT"
    assert modification_pattern == "none"


def test_parse_zfn_mutation_constraints_valid() -> None:
    """Valid CLI mutation constraints should parse into typed models."""
    constraints, default_constraint, overall_constraints = _parse_zfn_mutation_constraints(
        ["1:2:substitution,transition", "2:0:deletion", "*:1:mismatch", "overall:3:substitution"]
    )
    assert len(constraints) == 2
    assert default_constraint == ZFNDefaultSubfingerMutationConstraint(
        max_mutations=1,
        mutation_types=[ZFNMutationType.SUBSTITUTION],
    )
    assert overall_constraints == [
        ZFNOverallMutationConstraint(
            max_mutations=3,
            mutation_types=[ZFNMutationType.SUBSTITUTION],
        )
    ]
    assert constraints[0] == ZFNSubfingerMutationConstraint(
        subfinger_index=1,
        max_mutations=2,
        mutation_types=[ZFNMutationType.SUBSTITUTION, ZFNMutationType.TRANSITION],
    )
    assert constraints[1].mutation_types == [ZFNMutationType.DELETION]


def test_parse_zfn_mutation_constraints_invalid_format() -> None:
    """Invalid raw entries should raise ValueError."""
    with pytest.raises(ValueError):
        _parse_zfn_mutation_constraints(["1:2"])


def test_design_parameters_store_zfn_constraints() -> None:
    """DesignParameters should persist per-sub-finger mutation allowances."""
    params = DesignParameters(
        design_mode=DesignMode.ZFN,
        zfn_subfinger_mutations=[
            ZFNSubfingerMutationConstraint(
                subfinger_index=3,
                max_mutations=1,
                mutation_types=[ZFNMutationType.TRANSVERSION],
            )
        ],
        zfn_default_subfinger_mutation=ZFNDefaultSubfingerMutationConstraint(
            max_mutations=1,
            mutation_types=[ZFNMutationType.SUBSTITUTION],
        ),
        zfn_overall_mutations=[
            ZFNOverallMutationConstraint(
                max_mutations=3,
                mutation_types=[ZFNMutationType.SUBSTITUTION],
            )
        ],
    )
    assert params.design_mode == DesignMode.ZFN
    assert params.zfn_subfinger_mutations[0].subfinger_index == 3
    assert params.zfn_default_subfinger_mutation is not None
    assert params.zfn_overall_mutations[0].max_mutations == 3
