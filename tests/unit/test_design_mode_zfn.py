"""Unit tests for ZFN design mode parsing and mutation-constraint wiring."""

import pytest

from sirnaforge.cli import _build_zfn_design_configuration, _parse_zfn_mutation_constraints, _resolve_design_mode
from sirnaforge.models.sirna import (
    DesignMode,
)
from sirnaforge.models.zfn import (
    ZFNDefaultSubfingerMutationConstraint,
    ZFNDesignParameters,
    ZFNMutationConstraints,
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
    """ZFNMutationConstraints should persist per-sub-finger mutation allowances on ZFNDesignParameters."""
    mc = ZFNMutationConstraints(
        subfinger_mutations=[
            ZFNSubfingerMutationConstraint(
                subfinger_index=3,
                max_mutations=1,
                mutation_types=[ZFNMutationType.TRANSVERSION],
            )
        ],
        default_subfinger_mutation=ZFNDefaultSubfingerMutationConstraint(
            max_mutations=1,
            mutation_types=[ZFNMutationType.SUBSTITUTION],
        ),
        overall_mutations=[
            ZFNOverallMutationConstraint(
                max_mutations=3,
                mutation_types=[ZFNMutationType.SUBSTITUTION],
            )
        ],
    )
    params = ZFNDesignParameters(
        left_half_site="GCCCCACTGTGGGGTG",
        right_half_site="ACCAGATGACTGATGA",
        mutation_constraints=mc,
    )
    assert params.mutation_constraints is not None
    assert params.mutation_constraints.subfinger_mutations[0].subfinger_index == 3
    assert params.mutation_constraints.default_subfinger_mutation is not None
    assert params.mutation_constraints.overall_mutations[0].max_mutations == 3


def test_design_parameters_canonical_contract_normalizes_spacers() -> None:
    """Canonical contract should expose a stable normalized internal representation."""
    params = ZFNDesignParameters(
        left_half_site="gcCCCACTG",
        right_half_site="accagatga",
        spacer_constraints={"allowed_spacer_lengths": [6, 5, 6]},
    )
    contract = params.canonical_search_contract()

    assert contract.left_half_site == "GCCCCACTG"
    assert contract.right_half_site == "ACCAGATGA"
    assert contract.allowed_spacer_lengths == [5, 6]
    assert contract.orientation_convention == "L...R genomic ordering"


def test_build_zfn_design_configuration_applies_cli_sharding_overrides() -> None:
    """CLI sharding options should propagate into typed ZFN sharding config."""
    params, _, _, _, _ = _build_zfn_design_configuration(
        zfn_subfinger_mutation=[],
        zfn_max_mismatches_per_subfinger=None,
        zfn_max_substitutions_overall=None,
        zfn_left_half_site="GCGTGGGCG",
        zfn_right_half_site="GCCCACGCG",
        zfn_search_space="ensembl_human_hg38_primary",
        zfn_algorithm="zfn_v2",
        zfn_dimer_mode="heterodimer_only",
        zfn_spacer_lengths="5,6,7",
        zfn_max_mismatches=2,
        zfn_annotation=None,
        zfn_shard_max_workers=6,
        zfn_shard_chunk_mb=40.0,
        zfn_shard_chromosomes="autosomes,sex",
    )

    assert params.sharding.max_workers == 6
    assert params.sharding.chunk_size_bp == 40_000_000
    assert params.sharding.chromosomes == ["autosomes", "sex"]
