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
from sirnaforge.models.sirna import (
    DesignParameters,
    FilterCriteria,
    PostScreenSiRNAWeights,
    ScoringWeights,
)
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

    assert manifest["design_parameters"]["scoring"]["postscreen_sirna"]["off_target"] == 0.25
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
    assert scoring["scored_terms"] == list(COMPOSITE_TERMS)
    # Issue #96: each vector is recorded under its own NAME, because that name is stamped on every
    # candidate row -- a score without its vector cannot be traced to the weights that made it.
    assert set(scoring["vectors"]) == {"design_v4", "postscreen_sirna_v4", "postscreen_mirna_v4"}
    assert scoring["vectors"] == DesignParameters().scoring.as_manifest()
    assert scoring["vector_terms"]["design_v4"] == ["target_accessibility", "asymmetry", "gc_content"]
    # Terms that are computed and reported but score nothing must be named as such, or their
    # absence from the weights reads as an omission.
    assert "empirical" in scoring["reported_not_scored"]
    assert "conservation" in scoring["reported_not_scored"]


@pytest.mark.unit
def test_manifest_weights_track_a_custom_weight_set(tmp_path):
    """A run with reweighted scoring must record its own weights, not the defaults."""
    custom = ScoringWeights(
        postscreen_sirna=PostScreenSiRNAWeights(
            off_target=0.40,
            target_accessibility=0.30,
            asymmetry=0.20,
            gc_content=0.10,
        )
    )
    params = DesignParameters(scoring=custom)
    config = WorkflowConfig(output_dir=tmp_path / "custom", gene_query="tp53", design_params=params)
    manifest = SiRNAWorkflow(config)._build_fair_manifest(
        all_csv=tmp_path / "all.csv",
        pass_csv=tmp_path / "pass.csv",
        pass_fasta=tmp_path / "pass.fasta",
        orf_report=tmp_path / "orf.tsv",
    )

    assert manifest["scoring"]["vectors"]["postscreen_sirna_v4"]["off_target"] == 0.40
    assert manifest["scoring"]["vectors"]["postscreen_sirna_v4"]["gc_content"] == 0.10
    # The untouched vectors still record their own declared numbers.
    assert manifest["scoring"]["vectors"]["design_v4"]["target_accessibility"] == 0.35


@pytest.mark.unit
def test_manifest_cannot_record_an_unnormalisable_weight_set(tmp_path):
    """Issue #80 story 20: a weight vector that cannot be normalised fails at configuration time.

    Catching it here rather than in the scorer is the point -- a misweighted run must never get
    far enough to produce a plausible-looking score and a manifest describing it.
    """
    with pytest.raises(ValidationError, match="sum to exactly 1.0"):
        ScoringWeights(postscreen_sirna=PostScreenSiRNAWeights(off_target=0.9))
