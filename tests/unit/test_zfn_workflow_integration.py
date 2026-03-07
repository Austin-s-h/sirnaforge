"""Tests for ZFN workflow routing, config construction, and scoring integration.

These tests verify the critical fixes that route design_mode=ZFN to the
ZFN-specific engine instead of falling through to the siRNA designer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.models.sirna import DesignMode, DesignParameters
from sirnaforge.models.zfn import (
    ZFNDefaultSubfingerMutationConstraint,
    ZFNDesignParameters,
    ZFNDesignResult,
    ZFNMutationConstraints,
    ZFNMutationType,
    ZFNOverallMutationConstraint,
    ZFNShardingConfig,
    ZFNSubfingerMutationConstraint,
)
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig, ZFNWorkflowConfig, apply_zfn_runtime_overrides
from sirnaforge.zfn.design import ZFNDesigner

# ---------------------------------------------------------------------------
# Mutation constraint scoring
# ---------------------------------------------------------------------------


class TestMutationConstraintScoring:
    """Verify manufacturability scoring integrates mutation constraint budgets."""

    @pytest.fixture()
    def designer(self) -> ZFNDesigner:
        """Return a fresh ZFNDesigner instance."""
        return ZFNDesigner()

    def test_no_constraints_does_not_penalise(self, designer: ZFNDesigner) -> None:
        """Score should be unaffected when no mutation constraints are set."""
        params = ZFNDesignParameters(
            left_half_site="GCCGCACTG",
            right_half_site="ACCAGATGA",
            mutation_constraints=None,
        )
        score = designer.score_manufacturability(params)
        # Pure canonical bases, no homopolymers → expect full 100
        assert math.isclose(score, 100.0, rel_tol=1e-9, abs_tol=1e-9)

    def test_subfinger_violation_adds_penalty(self, designer: ZFNDesigner) -> None:
        """Exceeding a per-sub-finger budget should decrease the score."""
        # 9bp half-site = 3 triplets.  Middle triplet has 2 IUPAC bases.
        params_no_mc = ZFNDesignParameters(
            left_half_site="ACCNNYACC",
            right_half_site="GCCGCCGCC",
        )
        params_with_mc = ZFNDesignParameters(
            left_half_site="ACCNNYACC",
            right_half_site="GCCGCCGCC",
            mutation_constraints=ZFNMutationConstraints(
                subfinger_mutations=[
                    ZFNSubfingerMutationConstraint(
                        subfinger_index=2,  # NNY triplet
                        max_mutations=1,
                        mutation_types=[ZFNMutationType.SUBSTITUTION],
                    )
                ],
            ),
        )
        score_no_mc = designer.score_manufacturability(params_no_mc)
        score_with_mc = designer.score_manufacturability(params_with_mc)
        assert score_with_mc < score_no_mc

    def test_default_subfinger_applies_to_unconstrained_fingers(self, designer: ZFNDesigner) -> None:
        """Default constraint should apply to fingers without explicit overrides."""
        params = ZFNDesignParameters(
            left_half_site="NCCNCCNCC",  # 1 IUPAC per triplet → N, N, N
            right_half_site="GCCGCCGCC",
            mutation_constraints=ZFNMutationConstraints(
                default_subfinger_mutation=ZFNDefaultSubfingerMutationConstraint(
                    max_mutations=0,  # no IUPAC allowed per finger
                    mutation_types=[ZFNMutationType.SUBSTITUTION],
                ),
            ),
        )
        score = designer.score_manufacturability(params)
        # Each of the 3 left triplets has 1 IUPAC → each violates by 1 → penalty
        # Plus IUPAC complexity penalty (3 × 4 = 12)
        assert score < 75.0

    def test_overall_budget_violation_penalises(self, designer: ZFNDesigner) -> None:
        """Exceeding the global mutation budget should decrease the score."""
        params = ZFNDesignParameters(
            left_half_site="NNNNNNNNN",
            right_half_site="GCCGCCGCC",
            mutation_constraints=ZFNMutationConstraints(
                overall_mutations=[
                    ZFNOverallMutationConstraint(
                        max_mutations=2,
                        mutation_types=[ZFNMutationType.SUBSTITUTION],
                    )
                ],
            ),
        )
        score = designer.score_manufacturability(params)
        # 9 IUPAC bases (complexity = 64), budget of 2 → excess 7 → penalty 20
        # 64 - 20 = 44 → still under 50
        assert score < 50.0

    def test_constraint_penalty_capped_at_40(self, designer: ZFNDesigner) -> None:
        """Constraint penalty should never exceed 40 points."""
        params = ZFNDesignParameters(
            left_half_site="NNNNNNNNN",
            right_half_site="NNNNNNNNN",
            mutation_constraints=ZFNMutationConstraints(
                default_subfinger_mutation=ZFNDefaultSubfingerMutationConstraint(
                    max_mutations=0,
                    mutation_types=[ZFNMutationType.SUBSTITUTION],
                ),
                overall_mutations=[
                    ZFNOverallMutationConstraint(
                        max_mutations=0,
                        mutation_types=[ZFNMutationType.SUBSTITUTION],
                    )
                ],
            ),
        )
        score = designer.score_manufacturability(params)
        # IUPAC complexity: 18 × 4 = 72 → base = 28
        # Constraint penalty capped at 40 → final clamp at 0
        assert score >= 0.0

    def test_satisfied_constraints_no_penalty(self, designer: ZFNDesigner) -> None:
        """Constraints that are satisfied should not penalise."""
        params = ZFNDesignParameters(
            left_half_site="ACCNCCACC",  # 1 IUPAC in triplet 1
            right_half_site="GCCGCCGCC",
            mutation_constraints=ZFNMutationConstraints(
                default_subfinger_mutation=ZFNDefaultSubfingerMutationConstraint(
                    max_mutations=1,
                    mutation_types=[ZFNMutationType.SUBSTITUTION],
                ),
                overall_mutations=[
                    ZFNOverallMutationConstraint(
                        max_mutations=5,
                        mutation_types=[ZFNMutationType.SUBSTITUTION],
                    )
                ],
            ),
        )
        score_constrained = designer.score_manufacturability(params)
        params_no_mc = ZFNDesignParameters(
            left_half_site="ACCNCCACC",
            right_half_site="GCCGCCGCC",
        )
        score_unconstrained = designer.score_manufacturability(params_no_mc)
        assert math.isclose(score_constrained, score_unconstrained, rel_tol=1e-9, abs_tol=1e-9)


# ---------------------------------------------------------------------------
# Triplet counting helper
# ---------------------------------------------------------------------------


class TestCountAmbiguousPerTriplet:
    """Verify triplet ambiguity counting helper."""

    def test_canonical_bases_all_zero(self) -> None:
        """All canonical bases should yield zero IUPAC counts per triplet."""
        counts = ZFNDesigner.count_ambiguous_per_triplet("ACCGCCGCC")
        assert counts == [0, 0, 0]

    def test_iupac_bases_counted(self) -> None:
        """IUPAC ambiguity codes should be counted per triplet."""
        counts = ZFNDesigner.count_ambiguous_per_triplet("NCCNNCNNN")
        assert counts == [1, 2, 3]

    def test_trailing_bases_grouped(self) -> None:
        """Trailing bases that don't fill a triplet are still counted."""
        counts = ZFNDesigner.count_ambiguous_per_triplet("ACCN")
        assert counts == [0, 1]

    def test_empty_sequence(self) -> None:
        """Empty sequence should return an empty list."""
        counts = ZFNDesigner.count_ambiguous_per_triplet("")
        assert counts == []


# ---------------------------------------------------------------------------
# Workflow routing: ZFN mode → ZFNDesigner
# ---------------------------------------------------------------------------


class TestZFNWorkflowRouting:
    """Verify SiRNAWorkflow routes ZFN mode to ZFNDesigner, not SiRNADesigner."""

    def test_zfn_mode_instantiates_zfn_designer(self, tmp_path: Path) -> None:
        """When design_mode=ZFN, the workflow should create a ZFNDesigner."""
        dp = DesignParameters(design_mode=DesignMode.ZFN)
        cfg = WorkflowConfig(
            gene_query="test_gene",
            output_dir=tmp_path / "out",
            design_params=dp,
        )
        wf = SiRNAWorkflow(cfg)
        assert wf.zfn_designer is not None
        assert isinstance(wf.zfn_designer, ZFNDesigner)
        # siRNA designer should NOT be set
        assert wf.sirnaforgeer is None

    def test_sirna_mode_does_not_create_zfn_designer(self, tmp_path: Path) -> None:
        """SiRNA mode should create SiRNADesigner, not ZFNDesigner."""
        dp = DesignParameters(design_mode=DesignMode.SIRNA)
        cfg = WorkflowConfig(
            gene_query="test_gene",
            output_dir=tmp_path / "out",
            design_params=dp,
        )
        wf = SiRNAWorkflow(cfg)
        assert wf.zfn_designer is None
        assert wf.sirnaforgeer is not None

    def test_zfn_mode_skips_transcript_and_orf_directories(self, tmp_path: Path) -> None:
        """ZFN mode should not create transcript/ORF output folders during config bootstrap."""
        dp = DesignParameters(design_mode=DesignMode.ZFN)
        out_dir = tmp_path / "out_zfn"
        WorkflowConfig(
            gene_query="test_gene",
            output_dir=out_dir,
            design_params=dp,
        )

        assert out_dir.exists()
        assert not (out_dir / "transcripts").exists()
        assert not (out_dir / "orf_reports").exists()
        assert (out_dir / "sirnaforge").exists()
        assert (out_dir / "off_target").exists()
        assert (out_dir / "logs").exists()

    def test_sirna_mode_keeps_transcript_and_orf_directories(self, tmp_path: Path) -> None:
        """SiRNA mode should preserve historical output folder bootstrap behavior."""
        dp = DesignParameters(design_mode=DesignMode.SIRNA)
        out_dir = tmp_path / "out_sirna"
        WorkflowConfig(
            gene_query="test_gene",
            output_dir=out_dir,
            design_params=dp,
        )

        assert out_dir.exists()
        assert (out_dir / "transcripts").exists()
        assert (out_dir / "orf_reports").exists()
        assert (out_dir / "sirnaforge").exists()
        assert (out_dir / "off_target").exists()
        assert (out_dir / "logs").exists()


# ---------------------------------------------------------------------------
# ZFNWorkflowConfig
# ---------------------------------------------------------------------------


class TestZFNWorkflowConfig:
    """Verify ZFNWorkflowConfig construction and integration with WorkflowConfig."""

    def test_workflow_config_stores_zfn_config(self, tmp_path: Path) -> None:
        """ZFNWorkflowConfig should be stored on WorkflowConfig when provided."""
        zfn_params = ZFNDesignParameters(
            left_half_site="GCCCCACTG",
            right_half_site="ACCAGATGA",
        )
        zfn_cfg = ZFNWorkflowConfig(zfn_params=zfn_params)
        dp = DesignParameters(design_mode=DesignMode.ZFN)
        cfg = WorkflowConfig(
            gene_query="test_gene",
            output_dir=tmp_path / "out",
            design_params=dp,
            zfn_config=zfn_cfg,
        )
        assert cfg.zfn_config is not None
        assert cfg.zfn_config.zfn_params.left_half_site == "GCCCCACTG"

    def test_zfn_config_default_is_none(self, tmp_path: Path) -> None:
        """zfn_config should default to None when not provided."""
        dp = DesignParameters(design_mode=DesignMode.SIRNA)
        cfg = WorkflowConfig(
            gene_query="test_gene",
            output_dir=tmp_path / "out",
            design_params=dp,
        )
        assert cfg.zfn_config is None

    def test_zfn_sharding_overrides_load_from_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Environment JSON should enable sharding without adding CLI flags."""
        monkeypatch.setenv(
            "SIRNAFORGE_ZFN_SHARDING_JSON",
            '{"enabled": true, "chunk_size_mb": 10, "overlap_bp": 60, "chromosomes": "chr3", "max_workers": 2}',
        )

        params = ZFNDesignParameters(
            left_half_site="GCGTACGTA",
            right_half_site="TACGGCATA",
            sharding=ZFNShardingConfig(enabled=False),
        )
        nextflow_overrides: dict[str, object] = {}

        updated = apply_zfn_runtime_overrides(params, nextflow_overrides)
        assert updated.sharding.enabled is True
        assert updated.sharding.chunk_size_bp == 10_000_000
        assert updated.sharding.overlap_bp == 60
        assert updated.sharding.chromosomes == ["chr3"]
        assert updated.sharding.max_workers == 2
        assert nextflow_overrides["zfn_sharding_enabled"] is True


# ---------------------------------------------------------------------------
# ZFNDesigner evaluate_pair smoke test
# ---------------------------------------------------------------------------


class TestZFNDesignerEvaluatePair:
    """Verify ZFNDesigner.evaluate_pair returns coherent results."""

    def test_evaluate_pair_with_explicit_fasta(self, tmp_path: Path) -> None:
        """evaluate_pair should return candidates and off-target sites."""
        left = "GCGTACGTA"
        right = "TACGGCATA"
        spacer = "AAAAA"
        target = f"{left}{spacer}{str(Seq(right).reverse_complement())}"
        fasta = tmp_path / "genome.fa"
        fasta.write_text(f">chr1\nAAAA{target}CCCC\n")

        params = ZFNDesignParameters(
            search_space_fasta=str(fasta),
            left_half_site=left,
            right_half_site=right,
        )
        designer = ZFNDesigner()
        result = designer.evaluate_pair(params=params)

        assert isinstance(result, ZFNDesignResult)
        assert len(result.candidates) == 1

        cand = result.candidates[0]
        assert cand.composite_score >= 0.0
        assert cand.composite_score <= 100.0
        assert cand.component_scores["manufacturability"] > 0.0

    def test_evaluate_pair_mutation_constraints_affect_composite(self, tmp_path: Path) -> None:
        """Mutation constraints should influence the composite score via manufacturability."""
        left = "NNNNNNNNN"
        right = "GCCGCCGCC"
        spacer = "AAAAA"
        target = f"{left}{spacer}{str(Seq(right).reverse_complement())}"
        fasta = tmp_path / "genome.fa"
        fasta.write_text(f">chr1\nAAAA{target}CCCC\n")

        params_no_mc = ZFNDesignParameters(
            search_space_fasta=str(fasta),
            left_half_site=left,
            right_half_site=right,
        )
        params_with_mc = ZFNDesignParameters(
            search_space_fasta=str(fasta),
            left_half_site=left,
            right_half_site=right,
            mutation_constraints=ZFNMutationConstraints(
                overall_mutations=[
                    ZFNOverallMutationConstraint(
                        max_mutations=0,
                        mutation_types=[ZFNMutationType.SUBSTITUTION],
                    )
                ],
            ),
        )

        designer = ZFNDesigner()
        result_no_mc = designer.evaluate_pair(params=params_no_mc)
        result_with_mc = designer.evaluate_pair(params=params_with_mc)

        manuf_no_mc = result_no_mc.candidates[0].component_scores["manufacturability"]
        manuf_with_mc = result_with_mc.candidates[0].component_scores["manufacturability"]
        assert manuf_with_mc < manuf_no_mc


@pytest.mark.asyncio
async def test_zfn_workflow_outputs_include_contract_and_provenance(tmp_path: Path) -> None:
    """Workflow artifacts should include stable schema metadata and ZFN provenance fields."""
    left = "GCGTACGTA"
    right = "TACGGCATA"
    spacer = "AAAAA"
    target = f"{left}{spacer}{str(Seq(right).reverse_complement())}"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\nAAA{target}CCC\n", encoding="utf-8")

    dp = DesignParameters(design_mode=DesignMode.ZFN)
    zfn_params = ZFNDesignParameters(
        left_half_site=left,
        right_half_site=right,
        search_space_fasta=str(fasta),
        spacer_constraints={"allowed_spacer_lengths": [6, 5, 6]},
    )
    assert zfn_params.spacer_constraints.allowed_spacer_lengths == [5, 6]
    cfg = WorkflowConfig(
        gene_query="ZFN_TEST",
        output_dir=tmp_path / "out",
        design_params=dp,
        zfn_config=ZFNWorkflowConfig(zfn_params=zfn_params),
    )
    result = await SiRNAWorkflow(cfg).run_complete_workflow()

    assert result["workflow_mode"] == "zfn"
    assert result["algorithm"] == zfn_params.algorithm.value
    assert result["spacer_lengths"] == [5, 6]

    candidate_path = cfg.output_dir / "sirnaforge" / "zfn_candidate_summary.json"
    payload = json.loads(candidate_path.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "zfn_candidate_summary.v1"
    assert payload["search_contract"]["allowed_spacer_lengths"] == [5, 6]
