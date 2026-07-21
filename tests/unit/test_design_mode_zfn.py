"""Unit tests for ZFN design mode parsing and mutation-constraint wiring."""

import json
from pathlib import Path

import pytest

import sirnaforge.cli as cli_module
from sirnaforge.cli import (
    _autotune_zfn_sharding,
    _build_zfn_design_configuration,
    _parse_zfn_mutation_constraints,
    _resolve_design_mode,
)
from sirnaforge.models.sirna import (
    DesignMode,
)
from sirnaforge.models.zfn import (
    DimerMode,
    ZFNAlgorithm,
    ZFNDefaultSubfingerMutationConstraint,
    ZFNDesignParameters,
    ZFNMutationConstraints,
    ZFNMutationType,
    ZFNOverallMutationConstraint,
    ZFNSearchBackend,
    ZFNShardingConfig,
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


def test_zfn_sharding_model_defaults_use_internal_parallel_profile() -> None:
    """Typed ZFN sharding defaults should request the internal two-worker profile."""
    sharding = ZFNShardingConfig()

    assert sharding.enabled is True
    assert sharding.chunk_size_bp == 12_000_000
    assert sharding.overlap_bp == 50
    assert sharding.chromosomes == []
    assert sharding.max_workers == 2


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


def test_build_zfn_design_configuration_autotunes_internal_sharding() -> None:
    """Sharding should be internalized with practical auto-tuned defaults."""
    params, _, _, _, _ = _build_zfn_design_configuration(
        zfn_subfinger_mutation=[],
        zfn_max_mismatches_per_subfinger=None,
        zfn_max_substitutions_overall=None,
        zfn_left_half_site="GCGTGGGCG",
        zfn_right_half_site="GCCCACGCG",
        zfn_search_space="ensembl_human_hg38_primary",
        zfn_search_space_index=None,
        zfn_search_backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON,
        zfn_algorithm=ZFNAlgorithm.ZFN_V2,
        zfn_dimer_mode=DimerMode.HETERODIMER_ONLY,
        zfn_spacer_lengths="5,6,7",
        zfn_max_mismatches=2,
        zfn_annotation=None,
    )

    assert params.sharding.enabled is True
    assert params.sharding.max_workers >= 1
    assert params.sharding.max_workers <= 8
    assert params.sharding.chunk_size_bp in {8_000_000, 12_000_000}
    assert params.sharding.chromosomes == []


def test_build_zfn_design_configuration_uses_workflow_cores_for_sharding() -> None:
    """Workflow core budget should cap/drive internal ZFN shard workers."""
    params, _, _, _, _ = _build_zfn_design_configuration(
        zfn_subfinger_mutation=[],
        zfn_max_mismatches_per_subfinger=None,
        zfn_max_substitutions_overall=None,
        zfn_left_half_site="GCGTGGGCG",
        zfn_right_half_site="GCCCACGCG",
        zfn_search_space="ensembl_human_hg38_primary",
        zfn_search_space_index=None,
        zfn_search_backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON,
        zfn_algorithm=ZFNAlgorithm.ZFN_V2,
        zfn_dimer_mode=DimerMode.HETERODIMER_ONLY,
        zfn_spacer_lengths="5,6,7",
        zfn_max_mismatches=2,
        workflow_cores=2,
        zfn_annotation=None,
    )

    assert params.sharding.max_workers == 2
    assert params.sharding.chunk_size_bp == 12_000_000


def test_build_zfn_design_configuration_caps_shard_workers_for_large_cores() -> None:
    """Internal sharding workers should remain bounded on high-core hosts."""
    params, _, _, _, _ = _build_zfn_design_configuration(
        zfn_subfinger_mutation=[],
        zfn_max_mismatches_per_subfinger=None,
        zfn_max_substitutions_overall=None,
        zfn_left_half_site="GCGTGGGCG",
        zfn_right_half_site="GCCCACGCG",
        zfn_search_space="ensembl_human_hg38_primary",
        zfn_search_space_index=None,
        zfn_search_backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON,
        zfn_algorithm=ZFNAlgorithm.ZFN_V2,
        zfn_dimer_mode=DimerMode.HETERODIMER_ONLY,
        zfn_spacer_lengths="5,6,7",
        zfn_max_mismatches=2,
        workflow_cores=16,
        zfn_annotation=None,
    )

    assert params.sharding.max_workers == 8
    assert params.sharding.chunk_size_bp == 8_000_000


def test_autotune_zfn_sharding_uses_backend_aware_fm_index_profile() -> None:
    """FM-index should use a more conservative worker/chunk profile than scan backends.

    For pyahocorasick the chunk_size_bp is a large-contig fallback threshold;
    the primary execution path (contig_first) does not use it for scheduling.
    """
    fm_index = _autotune_zfn_sharding(cores_budget=12, search_backend=ZFNSearchBackend.FM_INDEX)
    pyaho = _autotune_zfn_sharding(cores_budget=12, search_backend=ZFNSearchBackend.PYAHOCORASICK)

    assert fm_index.max_workers == 4
    assert fm_index.chunk_size_bp == 16_000_000
    assert pyaho.max_workers == 2
    # chunk_size_bp is retained as the large-contig chunked-fallback threshold.
    assert pyaho.chunk_size_bp == 8_000_000


def test_build_zfn_design_configuration_rejects_invalid_half_site_bases() -> None:
    """Invalid half-site bases should fail fast during typed parameter construction."""
    with pytest.raises(ValueError, match=r"Invalid base\(s\) in half-site"):
        _build_zfn_design_configuration(
            zfn_subfinger_mutation=[],
            zfn_max_mismatches_per_subfinger=None,
            zfn_max_substitutions_overall=None,
            zfn_left_half_site="BADSEQ",
            zfn_right_half_site="GCCCACGCG",
            zfn_search_space="ensembl_human_hg38_primary",
            zfn_search_space_index=None,
            zfn_search_backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON,
            zfn_algorithm=ZFNAlgorithm.ZFN_V2,
            zfn_dimer_mode=DimerMode.HETERODIMER_ONLY,
            zfn_spacer_lengths="5,6",
            zfn_max_mismatches=2,
            zfn_annotation=None,
        )


def test_build_zfn_design_configuration_rejects_invalid_spacer_lengths() -> None:
    """Spacer lengths outside valid biological range should be rejected."""
    with pytest.raises(ValueError, match="allowed_spacer_lengths"):
        _build_zfn_design_configuration(
            zfn_subfinger_mutation=[],
            zfn_max_mismatches_per_subfinger=None,
            zfn_max_substitutions_overall=None,
            zfn_left_half_site="GCGTGGGCG",
            zfn_right_half_site="GCCCACGCG",
            zfn_search_space="ensembl_human_hg38_primary",
            zfn_search_space_index=None,
            zfn_search_backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON,
            zfn_algorithm=ZFNAlgorithm.ZFN_V2,
            zfn_dimer_mode=DimerMode.HETERODIMER_ONLY,
            zfn_spacer_lengths="0",
            zfn_max_mismatches=2,
            zfn_annotation=None,
        )


def test_internal_zfn_build_search_index_delegates_and_prints_summary(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Internal index build command should pass typed args through and emit JSON summary."""
    genome_fasta = tmp_path / "genome.fa"
    genome_fasta.write_text(">chr1\nACGT\n", encoding="utf-8")

    seen: dict[str, object] = {}

    def _fake_build(*, backend: ZFNSearchBackend, genome_fasta: Path, output_dir: Path | None) -> dict[str, object]:
        seen["backend"] = backend
        seen["genome_fasta"] = genome_fasta
        seen["output_dir"] = output_dir
        return {
            "artifact": "/tmp/fm_bundle/multifm_index.pkl",
            "backend": backend.value,
            "bases": 4,
            "bundle_dir": "/tmp/fm_bundle",
            "contigs": 1,
            "genome_fasta": str(genome_fasta),
        }

    printed: dict[str, str] = {}

    def _capture_print(text: object) -> None:
        printed.setdefault("text", str(text))

    monkeypatch.setattr(cli_module, "build_zfn_search_index", _fake_build)
    monkeypatch.setattr(cli_module.console, "print", _capture_print)

    cli_module.internal_zfn_build_search_index(
        genome_fasta=genome_fasta,
        search_backend=ZFNSearchBackend.FM_INDEX,
        output_dir=None,
    )

    assert seen["backend"] == ZFNSearchBackend.FM_INDEX
    assert seen["genome_fasta"] == genome_fasta
    assert seen["output_dir"] is None
    payload = json.loads(printed["text"])
    assert payload["backend"] == "fm_index"
    assert payload["contigs"] == 1
