"""The run manifest must be a complete record of the parameters applied.

Issue #78: `manifest.json` listed four design parameters, so thresholds that
decided the outcome -- `min_asymmetry_score`, the homopolymer gate -- appeared
nowhere in it.
"""

import json

import pytest

from sirnaforge import __version__
from sirnaforge.models.sirna import DesignParameters, FilterCriteria
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig


def _manifest(tmp_path, **filter_overrides) -> dict:
    """Build a manifest from a workflow that has not been run."""
    params = DesignParameters(filters=FilterCriteria(**filter_overrides))
    config = WorkflowConfig(output_dir=tmp_path / "out", gene_query="tp53", design_params=params)
    return SiRNAWorkflow(config)._build_fair_manifest(
        all_csv=tmp_path / "all.csv",
        pass_csv=tmp_path / "pass.csv",
        pass_fasta=tmp_path / "pass.fasta",
        orf_report=tmp_path / "orf.tsv",
    )


@pytest.mark.unit
def test_manifest_records_every_applied_threshold(tmp_path):
    """Filters that decide pass/fail must appear in the manifest."""
    manifest = _manifest(tmp_path, min_asymmetry_score=0.72, min_empirical_score=0.6, max_poly_runs=2)
    filters = manifest["design_parameters"]["filters"]

    assert filters["min_asymmetry_score"] == 0.72
    assert filters["min_empirical_score"] == 0.6
    assert filters["max_poly_runs"] == 2
    assert filters["max_paired_fraction"] == FilterCriteria().max_paired_fraction


@pytest.mark.unit
def test_manifest_records_mode_weights_and_version(tmp_path):
    """The manifest also carries the scoring weights, design mode and tool version."""
    manifest = _manifest(tmp_path)

    assert manifest["design_parameters"]["scoring"]["off_target"] == 0.30
    assert manifest["design_parameters"]["design_mode"] == "sirna"
    assert manifest["tool_version"] == __version__
    assert json.dumps(manifest), "manifest must stay JSON-serialisable"
