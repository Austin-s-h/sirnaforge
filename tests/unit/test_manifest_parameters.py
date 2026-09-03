"""The run manifest must be a complete record of the parameters applied.

Issue #78: `manifest.json` listed four design parameters, so thresholds that
decided the outcome -- `min_asymmetry_score`, the homopolymer gate -- appeared
nowhere in it.
"""

import json

import pytest
from pydantic import ValidationError

from sirnaforge import __version__
from sirnaforge.core.scoring import COMPOSITE_TERMS, SCORING_WEIGHT_SET_VERSION
from sirnaforge.models.sirna import DesignParameters, FilterCriteria, ScoringWeights
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

    assert manifest["design_parameters"]["scoring"]["off_target"] == 0.25
    assert manifest["design_parameters"]["design_mode"] == "sirna"
    assert manifest["tool_version"] == __version__
    assert json.dumps(manifest), "manifest must stay JSON-serialisable"


@pytest.mark.unit
def test_manifest_records_the_weight_set_version(tmp_path):
    """Issue #80: a score is meaningless without the weight definition that produced it.

    Two regimes of `composite_score` now exist (the pre-#80 five-term set and the seven-term
    set), so the manifest must name which one ran or results from the two get compared silently.
    """
    scoring = _manifest(tmp_path)["scoring"]

    assert scoring["weight_set_version"] == SCORING_WEIGHT_SET_VERSION
    assert scoring["active_terms"] == list(COMPOSITE_TERMS)
    # The recorded weights must be the ones actually applied, not a restatement of the defaults.
    assert scoring["weights"] == DesignParameters().scoring.model_dump(mode="json")


@pytest.mark.unit
def test_manifest_weights_track_a_custom_weight_set(tmp_path):
    """A run with reweighted scoring must record its own weights, not the defaults."""
    custom = ScoringWeights(
        asymmetry=0.10,
        gc_content=0.10,
        accessibility=0.10,
        empirical=0.10,
        off_target=0.40,
        isoform_coverage=0.10,
        conservation=0.10,
    )
    params = DesignParameters(scoring=custom)
    config = WorkflowConfig(output_dir=tmp_path / "custom", gene_query="tp53", design_params=params)
    manifest = SiRNAWorkflow(config)._build_fair_manifest(
        all_csv=tmp_path / "all.csv",
        pass_csv=tmp_path / "pass.csv",
        pass_fasta=tmp_path / "pass.fasta",
        orf_report=tmp_path / "orf.tsv",
    )

    assert manifest["scoring"]["weights"]["off_target"] == 0.40
    assert manifest["scoring"]["weights"]["asymmetry"] == 0.10


@pytest.mark.unit
def test_manifest_cannot_record_an_unnormalisable_weight_set(tmp_path):
    """Issue #80 story 20: a weight vector that cannot be normalised fails at configuration time.

    Catching it here rather than in the scorer is the point -- a misweighted run must never get
    far enough to produce a plausible-looking score and a manifest describing it.
    """
    with pytest.raises(ValidationError, match="sum to 1.0"):
        ScoringWeights(off_target=0.9)
