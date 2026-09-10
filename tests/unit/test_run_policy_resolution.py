"""One resolution point for run mode, filter actions and thresholds (#99).

What is pinned here is the resolver's contract, not a set of numbers: the CLI, the Python API and
the off-target-only path must resolve the *same* configuration; precedence must be applied exactly
once; a preset must apply on provenance rather than on value equality; and a configuration error
must cost nothing.
"""

from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.config.run_policy import (
    BUILTIN_PROFILES,
    DEFAULT_PROFILE_NAME,
    LEGACY_DEFAULT_EXCEPTIONS,
    PASSTHROUGH_FIELDS,
    SELECTABLE_ACTIONS,
    SETTING_SPECS,
    EntryPoint,
    ResolvedRunPolicy,
    RunPolicyError,
    _model_default,
    declared_filter_ids,
    default_for,
    describe_parameters,
    resolve_run_policy,
    switchable_filter_ids,
)
from sirnaforge.models.policy import (
    FilterAction,
    FilterEvaluation,
    FilterStage,
    FilterVerdict,
    Requiredness,
    RunMode,
    ScreeningChannel,
    SettingSource,
    UnknownEvidenceAction,
)
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    FilterCriteria,
    MiRNADesignConfig,
    OffTargetFilterCriteria,
    PostScreenSiRNAWeights,
    ScoringWeights,
    SiRNACandidate,
    TargetAccessibilityConfig,
    build_candidate_row,
)

FASTA = ">t1\n" + "ATGCGCATGCATCGATCGATCGGCATCGATCGATCGACTAGCATCGACTGACTGCATCAGCATCAGCATCAGCTACGATCAG\n"


def _fasta(tmp_path: Path) -> Path:
    path = tmp_path / "toy.fa"
    path.write_text(FASTA)
    return path


def _toy_candidate() -> SiRNACandidate:
    """A minimal candidate, only so ``build_candidate_row`` can report its column set."""
    return SiRNACandidate(
        id="c1",
        transcript_id="t1",
        position=1,
        guide_sequence="A" * 21,
        passenger_sequence="T" * 21,
        gc_content=45.0,
        length=21,
        asymmetry_score=0.7,
    )


def _cli_policy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str) -> ResolvedRunPolicy:
    """Run the real `workflow` command and return the policy it handed to the workflow."""
    captured: dict[str, Any] = {}

    async def _fake_workflow(**kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return {"transcript_summary": {}, "design_summary": {}, "offtarget_summary": {}}

    monkeypatch.setattr("sirnaforge.cli.run_sirna_workflow", _fake_workflow)
    result = CliRunner().invoke(
        app,
        [
            "workflow",
            "TOY",
            "--input-fasta",
            str(_fasta(tmp_path)),
            "--output-dir",
            str(tmp_path / "out"),
            "--species",
            "human",
            *extra,
        ],
    )
    assert result.exit_code == 0, result.output
    policy = captured["resolved_policy"]
    assert isinstance(policy, ResolvedRunPolicy)
    return policy


def _design_cli(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str) -> DesignParameters:
    """Run the real `design` command and return the parameters the designer was constructed with."""
    captured: dict[str, Any] = {}

    class _StubDesigner:
        def __init__(self, parameters: DesignParameters) -> None:
            captured["parameters"] = parameters

        def design_from_file(self, _path: str) -> Any:
            raise SystemExit(0)

    monkeypatch.setattr("sirnaforge.cli.SiRNADesigner", _StubDesigner)
    monkeypatch.setattr("sirnaforge.core.design.MiRNADesigner", _StubDesigner)
    CliRunner().invoke(
        app,
        ["design", str(_fasta(tmp_path)), "--output", str(tmp_path / "out.csv"), *extra],
    )
    return captured["parameters"]


# --------------------------------------------------------------------------------------
# One resolved configuration across every entry point
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_cli_api_and_offtarget_only_resolve_the_same_configuration(tmp_path, monkeypatch):
    """The tranche criterion: one resolution point, three public surfaces, one answer."""
    cli_policy = _cli_policy(tmp_path, monkeypatch, "--gc-max", "58", "--max-off-targets", "9")
    api_policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        stated={"gc_max": 58.0, "max_off_target_count": 9},
        query_species="human",
        screen_species=["human"],
    )
    offtarget_policy = resolve_run_policy(
        entry_point=EntryPoint.OFFTARGET_ONLY,
        stated={"gc_max": 58.0, "max_off_target_count": 9},
        query_species="human",
        screen_species=["human"],
    )

    assert cli_policy.design_parameters == api_policy.design_parameters == offtarget_policy.design_parameters
    assert cli_policy.run_mode is api_policy.run_mode is offtarget_policy.run_mode is RunMode.QUALIFIED
    assert cli_policy.profile.content_hash == api_policy.profile.content_hash
    # Provenance agrees too: the entry point changes the default mode, never a threshold's authority.
    assert [record.model_dump() for record in cli_policy.resolved if record.key != "entry_point"] == [
        record.model_dump() for record in api_policy.resolved if record.key != "entry_point"
    ]
    assert cli_policy.source_of("gc_max") is SettingSource.EXPLICIT
    assert cli_policy.source_of("gc_min") is SettingSource.BUILTIN_PROFILE


@pytest.mark.unit
def test_the_offtarget_only_path_no_longer_gates_on_thresholds_from_nowhere(tmp_path):
    """It used to build a bare DesignParameters(), so its gates came from an unresolved default."""
    policy = resolve_run_policy(entry_point=EntryPoint.OFFTARGET_ONLY, stated={"max_off_target_count": 4})

    assert policy.descriptor("max_off_target_count").threshold == 4
    assert policy.source_of("max_off_target_count") is SettingSource.EXPLICIT
    assert policy.run_mode is RunMode.QUALIFIED


@pytest.mark.unit
def test_a_directly_built_parameter_object_is_described_without_changing_a_number(tmp_path):
    """The adapter for direct WorkflowConfig/API callers must not re-default anything.

    Driven with a *non-default* value on every model, because an all-default input satisfies
    ``policy.design_parameters == parameters`` even if the adapter dropped half of SETTING_SPECS.
    """
    parameters = DesignParameters(
        sirna_length=23,
        top_n=7,
        check_off_targets=False,
        apply_modifications=True,
        modification_pattern="none",
        default_overhang="UU",
        filters=FilterCriteria(gc_min=31.0, gc_max=61.0, max_poly_runs=4, min_isoform_coverage=0.5),
        offtarget_filters=OffTargetFilterCriteria(
            max_off_target_count=None, fail_on_high_risk_mirna=False, max_total_offtarget_hits=9
        ),
        target_accessibility=TargetAccessibilityConfig(window_size=120, max_bp_span=90, log_floor=-3.0),
    )
    policy = describe_parameters(parameters, entry_point=EntryPoint.SCREENING_WORKFLOW)

    assert policy.design_parameters == parameters
    assert all(record.source is SettingSource.EXPLICIT for record in policy.resolved if record.key in {"gc_min"})


@pytest.mark.unit
def test_every_model_field_the_adapter_must_carry_is_declared():
    """A field absent from SETTING_SPECS silently reverts to its own default in the adapter.

    ``describe_parameters`` rebuilds ``stated`` from SETTING_SPECS alone, so an undeclared field
    would be reported in the manifest as a value the run did not apply. #101 adds filter fields,
    which is exactly when this fires, so the coverage is pinned rather than assumed.
    """
    declared = {spec.field for spec in SETTING_SPECS}
    structural = {"design_mode", "filters", "offtarget_filters", "target_accessibility"}
    for model in (DesignParameters, FilterCriteria, OffTargetFilterCriteria, TargetAccessibilityConfig):
        missing = set(model.model_fields) - declared - set(PASSTHROUGH_FIELDS) - structural
        assert not missing, (
            f"{model.__name__} fields carried by neither SETTING_SPECS nor PASSTHROUGH_FIELDS: {missing}"
        )


# --------------------------------------------------------------------------------------
# Precedence, resolved exactly once
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_precedence_is_profile_then_preset_then_config_file_then_explicit(tmp_path):
    """Four layers, one answer each, and the winning authority is recorded per setting."""
    config = tmp_path / "policy.json"
    config.write_text('{"settings": {"gc_max": 55.0, "max_poly_runs": 4}}')

    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        design_mode="mirna",
        config_file=config,
        stated={"max_poly_runs": 5},
    )

    # profile only
    assert policy.value_of("max_paired_fraction") == default_for("max_paired_fraction")
    assert policy.source_of("max_paired_fraction") is SettingSource.BUILTIN_PROFILE
    # preset beats profile
    assert policy.value_of("default_overhang") == MiRNADesignConfig().overhang
    assert policy.source_of("default_overhang") is SettingSource.DESIGN_MODE_PRESET
    # config file beats the preset
    assert policy.value_of("gc_max") == 55.0
    assert policy.source_of("gc_max") is SettingSource.CONFIG_FILE
    # explicit beats the config file
    assert policy.value_of("max_poly_runs") == 5
    assert policy.source_of("max_poly_runs") is SettingSource.EXPLICIT


@pytest.mark.unit
def test_a_config_file_records_both_what_was_requested_and_what_won(tmp_path):
    """Requested and resolved are different records: an override that lost must still be visible."""
    config = tmp_path / "policy.toml"
    config.write_text("[settings]\ngc_max = 55.0\n")
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, config_file=config, stated={"gc_max": 62.0})

    requested = {(record.key, record.value, record.source) for record in policy.requested}
    assert ("gc_max", 55.0, SettingSource.CONFIG_FILE) in requested
    assert ("gc_max", 62.0, SettingSource.EXPLICIT) in requested
    assert policy.value_of("gc_max") == 62.0


@pytest.mark.unit
def test_a_run_mode_and_design_mode_from_a_config_file_are_attributed_to_the_file(tmp_path):
    """A file-supplied directive is a config-file override, not something the caller typed.

    Both directives were derived from the *arguments* rather than from where they came from, so a
    file-supplied run mode read as `explicit` and a file-supplied design mode read as
    `builtin_profile` with the detail "no design mode stated" -- false, and it had already moved
    gc_min/gc_max/default_overhang through the miRNA preset. Neither appeared in `requested` at all.
    """
    config = tmp_path / "policy.json"
    config.write_text('{"design_mode": "mirna", "run_mode": "exploratory", "settings": {"gc_max": 57.0}}')
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, config_file=config)

    assert policy.run_mode is RunMode.EXPLORATORY
    assert policy.design_mode is DesignMode.MIRNA
    assert policy.source_of("run_mode") is SettingSource.CONFIG_FILE
    assert policy.source_of("design_mode") is SettingSource.CONFIG_FILE
    design_mode_record = next(record for record in policy.resolved if record.key == "design_mode")
    assert design_mode_record.detail == "stated in the policy config file"

    requested = {(record.key, record.value, record.source) for record in policy.requested}
    assert ("run_mode", "exploratory", SettingSource.CONFIG_FILE) in requested
    assert ("design_mode", "mirna", SettingSource.CONFIG_FILE) in requested

    # A directive the caller passed directly still reads as explicit.
    typed = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, config_file=config, run_mode="qualified")
    assert typed.source_of("run_mode") is SettingSource.EXPLICIT


# --------------------------------------------------------------------------------------
# Field-set provenance, not value equality (#101 item 1)
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_an_explicit_gc_max_survives_mirna_mode(tmp_path, monkeypatch):
    """The regression: `--design-mode mirna --gc-max 60` was silently rewritten to 52.

    The old test was a value comparison against the siRNA default, which cannot tell an omitted
    option from one typed with that same value. Run through the real CLI so the fix is exercised
    where the defect lived.
    """
    policy = _cli_policy(tmp_path, monkeypatch, "--design-mode", "mirna", "--gc-max", "60")

    assert policy.design_parameters.filters.gc_max == 60.0
    assert policy.source_of("gc_max") is SettingSource.EXPLICIT


@pytest.mark.unit
def test_an_omitted_gc_max_still_takes_the_mirna_preset(tmp_path, monkeypatch):
    """The preset must keep working; only how it decides to apply has changed."""
    policy = _cli_policy(tmp_path, monkeypatch, "--design-mode", "mirna")

    assert policy.design_parameters.filters.gc_max == MiRNADesignConfig().gc_max
    assert policy.source_of("gc_max") is SettingSource.DESIGN_MODE_PRESET


@pytest.mark.unit
def test_an_explicit_overhang_and_modification_pattern_survive_mirna_mode(tmp_path, monkeypatch):
    """Same sentinel pattern, two more options: `--overhang dTdT` used to become UU."""
    policy = _cli_policy(
        tmp_path, monkeypatch, "--design-mode", "mirna", "--overhang", "dTdT", "--modifications", "none"
    )

    assert policy.design_parameters.default_overhang == "dTdT"
    assert policy.design_parameters.modification_pattern == "none"
    assert policy.design_parameters.apply_modifications is False


@pytest.mark.unit
def test_gc_max_65_is_a_supported_setting(tmp_path, monkeypatch):
    """Documented and supported in both modes: FilterCriteria.gc_max is already bounded 0-100."""
    for extra in (("--gc-max", "65"), ("--design-mode", "mirna", "--gc-max", "65")):
        policy = _cli_policy(tmp_path, monkeypatch, *extra)
        assert policy.design_parameters.filters.gc_max == 65.0


@pytest.mark.unit
def test_the_design_command_resolves_the_same_way(tmp_path, monkeypatch):
    """`design` shares the resolver, so the same provenance rule holds there."""
    assert _design_cli(tmp_path, monkeypatch, "--design-mode", "mirna", "--gc-max", "60").filters.gc_max == 60.0
    assert _design_cli(tmp_path, monkeypatch, "--design-mode", "mirna").filters.gc_max == MiRNADesignConfig().gc_max


# --------------------------------------------------------------------------------------
# Run modes
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_entry_points_carry_the_documented_default_run_mode():
    """Design defaults to design-only; the screening workflow and the off-target path to qualified."""
    assert resolve_run_policy(entry_point=EntryPoint.DESIGN_COMMAND).run_mode is RunMode.DESIGN_ONLY
    assert resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).run_mode is RunMode.QUALIFIED
    assert resolve_run_policy(entry_point=EntryPoint.OFFTARGET_ONLY).run_mode is RunMode.QUALIFIED


@pytest.mark.unit
def test_the_legacy_skip_flag_maps_visibly_to_design_only(tmp_path, monkeypatch):
    """Visibly, which is the requirement: the rule that produced the mode is named in the record."""
    policy = _cli_policy(tmp_path, monkeypatch, "--skip-off-targets")

    assert policy.run_mode is RunMode.DESIGN_ONLY
    assert policy.design_parameters.check_off_targets is False
    record = next(record for record in policy.resolved if record.key == "run_mode")
    assert record.source is SettingSource.RUN_MODE_RULE
    assert "skip" in (record.detail or "")


@pytest.mark.unit
def test_the_two_ways_to_ask_for_a_design_only_run_publish_the_same_reference_record(tmp_path, monkeypatch):
    """`--run-mode design_only` and `--skip-off-targets` are one run and must leave one record.

    Reference selection read the raw `skip_off_targets` flag, so `--run-mode design_only` resolved the
    four default Ensembl cDNA references and published them in `reference_summary.transcriptome` for a
    run that screened nothing. No download happened -- step 5 short-circuits earlier -- but the run
    record claimed references the run never touched.
    """
    captured: dict[str, Any] = {}

    async def _fake_workflow(**kwargs: Any) -> dict[str, Any]:
        captured[kwargs["gene_query"]] = kwargs["transcriptome_selection"]
        return {"transcript_summary": {}, "design_summary": {}, "offtarget_summary": {}}

    monkeypatch.setattr("sirnaforge.cli.run_sirna_workflow", _fake_workflow)
    for label, extra in (("SKIP", ["--skip-off-targets"]), ("MODE", ["--run-mode", "design_only"])):
        result = CliRunner().invoke(
            app, ["workflow", label, "--output-dir", str(tmp_path / label), "--species", "human", *extra]
        )
        assert result.exit_code == 0, result.output

    assert captured["SKIP"].to_metadata() == captured["MODE"].to_metadata()
    for selection in captured.values():
        assert selection.choices == ()
        assert selection.to_metadata()["enabled"] is False


@pytest.mark.unit
def test_qualified_with_screening_explicitly_off_is_rejected():
    """Both stated, and they contradict: qualified evidence cannot come from a run that screens nothing."""
    with pytest.raises(RunPolicyError, match="qualified requires off-target screening"):
        resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW,
            run_mode=RunMode.QUALIFIED,
            stated={"check_off_targets": False},
        )
    with pytest.raises(RunPolicyError, match="qualified requires off-target screening"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, run_mode="qualified", legacy_skip_screening=True)


@pytest.mark.unit
def test_screening_off_without_an_explicit_mode_resolves_to_design_only():
    """Not an error: the caller said nothing about the mode, so the mode follows the switch."""
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"check_off_targets": False})

    assert policy.run_mode is RunMode.DESIGN_ONLY
    assert policy.source_of("run_mode") is SettingSource.RUN_MODE_RULE


@pytest.mark.unit
def test_required_qualified_completeness_cannot_be_waived_while_staying_qualified():
    """Waiving it is exploratory by definition; asking for both spellings at once is an error."""
    with pytest.raises(RunPolicyError, match="cannot waive required screening completeness"):
        resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW,
            run_mode=RunMode.QUALIFIED,
            require_screening_completeness=False,
        )

    downgraded = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, require_screening_completeness=False)
    assert downgraded.run_mode is RunMode.EXPLORATORY


@pytest.mark.unit
def test_design_mode_stays_orthogonal_to_run_mode():
    """Choosing miRNA design says nothing about how complete the screening evidence must be."""
    for design_mode in (DesignMode.SIRNA, DesignMode.MIRNA):
        assert (
            resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, design_mode=design_mode).run_mode
            is RunMode.QUALIFIED
        )
        assert (
            resolve_run_policy(entry_point=EntryPoint.DESIGN_COMMAND, design_mode=design_mode).run_mode
            is RunMode.DESIGN_ONLY
        )


@pytest.mark.unit
def test_qualified_requires_the_query_species_transcriptome_and_nothing_else():
    """Requiredness is per channel and species; a secondary species is reportable, not disqualifying."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human", screen_species=["human", "mouse"]
    )
    requirements = policy.evidence_requirements

    assert requirements.required_pairs == frozenset({("transcriptome", "human")})
    assert requirements.requiredness_of(ScreeningChannel.TRANSCRIPTOME, "mouse") is Requiredness.EXPLORATORY
    assert requirements.requiredness_of(ScreeningChannel.MIRNA_SEED, "human") is Requiredness.EXPLORATORY
    assert requirements.unknown_evidence_action is UnknownEvidenceAction.FAIL


@pytest.mark.unit
def test_exploratory_requires_nothing_and_a_design_only_run_declares_no_channels():
    """A design-only run holds no screening evidence, so it cannot have a requirement about any."""
    exploratory = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, run_mode="exploratory")
    design_only = resolve_run_policy(entry_point=EntryPoint.DESIGN_COMMAND)

    assert exploratory.evidence_requirements.required_pairs == frozenset()
    assert exploratory.evidence_requirements.unknown_evidence_action is UnknownEvidenceAction.WARN
    assert design_only.evidence_requirements.channel_requirements == ()


# --------------------------------------------------------------------------------------
# Per-filter action, threshold, scope and missing-evidence policy
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_every_filter_with_a_clearable_threshold_can_be_disabled_independently():
    """Switching one gate off must actually clear its threshold, and disturb no other gate.

    Recording the action alone would not be enough: nothing in the 0.7.1 gate application reads a
    filter action, so a gate is switched off by removing the number it compares against -- which is
    the state ``_check_offtarget_filters`` already reads as "no gate".
    """
    baseline = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    switchable = switchable_filter_ids()
    assert len(switchable) == 10, switchable

    for filter_id in switchable:
        policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={filter_id: "off"})
        descriptor = policy.descriptor(filter_id)
        assert descriptor.action is FilterAction.OFF
        # The threshold is gone from the validated parameters, so the gate cannot fire.
        setting = next(resolved.setting_key for resolved in policy.filters if resolved.filter_id == filter_id)
        assert policy.value_of(setting) in (None, False)
        others = {
            resolved.filter_id: (resolved.descriptor.action, resolved.descriptor.threshold)
            for resolved in policy.filters
            if resolved.filter_id != filter_id
        }
        assert others == {
            resolved.filter_id: (resolved.descriptor.action, resolved.descriptor.threshold)
            for resolved in baseline.filters
            if resolved.filter_id != filter_id
        }


@pytest.mark.unit
def test_a_design_stage_gate_cannot_be_switched_off_in_0_7_1_and_says_so():
    """Reported rather than faked: the six design thresholds are floats with no absent value.

    Setting one to an inert extreme would record a threshold the user never chose, and making the
    action authoritative needs the gate application to read it, which belongs to the filter-verdict
    work rather than to the resolver.
    """
    unswitchable = set(declared_filter_ids()) - set(switchable_filter_ids())
    assert unswitchable == {
        "gc_content_min",
        "gc_content_max",
        "max_poly_runs",
        "max_paired_fraction",
        "min_asymmetry_score",
        "min_empirical_score",
    }

    for filter_id in sorted(unswitchable):
        # The message must name the filter the user typed, not the internal setting key it lands on.
        with pytest.raises(RunPolicyError, match=f"filter {filter_id!r} cannot be switched off"):
            resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={filter_id: "off"})


@pytest.mark.unit
def test_a_disabled_filter_is_not_evaluated_rather_than_passed():
    """The distinction the whole vocabulary exists for: off makes no claim about the candidate."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={"max_off_target_count": "off"}
    )
    descriptor = policy.descriptor("max_off_target_count")

    assert FilterVerdict(descriptor=descriptor, observed=40, evaluation=FilterEvaluation.NOT_EVALUATED)
    for evaluation in (FilterEvaluation.PASS, FilterEvaluation.FAIL, FilterEvaluation.UNKNOWN):
        with pytest.raises(ValidationError):
            FilterVerdict(descriptor=descriptor, observed=40, evaluation=evaluation)


@pytest.mark.unit
def test_a_post_screen_gate_is_off_in_a_design_only_run():
    """A design-only run holds no screening evidence, so its off-target gates are not evaluated."""
    policy = resolve_run_policy(entry_point=EntryPoint.DESIGN_COMMAND)

    post_screen = [resolved for resolved in policy.filters if resolved.descriptor.stage is FilterStage.POST_SCREEN]
    assert post_screen, "the registry must declare post-screen gates"
    assert all(resolved.descriptor.action is FilterAction.OFF for resolved in post_screen)
    assert all(resolved.descriptor.stage is FilterStage.DESIGN for resolved in policy.evaluated_filters)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("run_mode", "kwargs"),
    [
        ("exploratory", {"legacy_skip_screening": True}),
        ("exploratory", {"stated": {"check_off_targets": False}}),
        (None, {"stated": {"check_off_targets": False}}),
    ],
)
def test_a_post_screen_gate_is_off_whenever_the_run_screens_nothing(run_mode, kwargs):
    """The rule is the resolved screening switch, not the ``design_only`` label.

    ``--run-mode exploratory --skip-off-targets`` is accepted and resolves ``check_off_targets=False``,
    so step 5 short-circuits and no candidate gets an off-target verdict. Keying the post-screen OFF
    rule on ``DESIGN_ONLY`` alone reported six gates as ``action=fail, evaluated=true`` for that run.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, run_mode=run_mode, **kwargs)

    assert policy.design_parameters.check_off_targets is False
    post_screen = [resolved for resolved in policy.filters if resolved.descriptor.stage is FilterStage.POST_SCREEN]
    assert post_screen, "the registry must declare post-screen gates"
    still_enforced = [resolved.filter_id for resolved in post_screen if resolved.is_evaluated]
    assert not still_enforced, f"gates reported as enforced on a run that screens nothing: {still_enforced}"
    assert all(resolved.descriptor.action is FilterAction.OFF for resolved in post_screen)


@pytest.mark.unit
def test_a_post_screen_gate_cannot_be_switched_on_when_the_run_screens_nothing():
    """Nor can an explicit action put it back: there is no evidence for it to read."""
    with pytest.raises(RunPolicyError, match="holds no screening evidence"):
        resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW,
            run_mode="exploratory",
            legacy_skip_screening=True,
            filter_actions={"max_off_target_count": "fail"},
        )


@pytest.mark.unit
def test_an_undeclared_threshold_resolves_to_off_not_to_a_silent_pass():
    """max_transcriptome_seed_perfect and min_isoform_coverage ship as None, so they cannot decide."""
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)

    for filter_id in ("max_transcriptome_seed_perfect", "min_isoform_coverage"):
        descriptor = policy.descriptor(filter_id)
        assert descriptor.threshold is None
        assert descriptor.action is FilterAction.OFF

    opted_in = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"max_transcriptome_seed_perfect": 2}
    )
    assert opted_in.descriptor("max_transcriptome_seed_perfect").action is FilterAction.FAIL


@pytest.mark.unit
def test_a_filter_cannot_be_switched_on_without_a_threshold_to_compare_against():
    """An action of fail with no threshold is a gate that cannot be evaluated, reported as in force."""
    with pytest.raises(RunPolicyError, match="declares no threshold"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={"min_isoform_coverage": "fail"})


@pytest.mark.unit
def test_a_gate_no_code_reads_cannot_be_switched_on():
    """`fail` on a gate nothing applies would report an enforced limit that does not exist.

    `max_mirna_1mm_seed` carries a threshold of 10, so the "no threshold" guard did not catch it and
    `--filter-action max_mirna_1mm_seed=fail` resolved to action=fail, threshold=10, evaluated=true.
    """
    with pytest.raises(RunPolicyError, match="no code in 0.7.1 reads max_mirna_1mm_seed"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={"max_mirna_1mm_seed": "fail"})


@pytest.mark.unit
def test_a_boolean_gate_that_is_disabled_cannot_be_switched_on_by_its_action():
    """`fail_on_high_risk_mirna=False` means the gate does not run; an action cannot override that."""
    with pytest.raises(RunPolicyError, match="fail_on_high_risk_mirna is False"):
        resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW,
            stated={"fail_on_high_risk_mirna": False},
            filter_actions={"fail_on_high_risk_mirna": "fail"},
        )


@pytest.mark.unit
def test_warn_is_not_selectable_in_0_7_1_because_nothing_applies_it():
    """A `warn` in the manifest beside a rejected candidate is a claim the gate code does not hold.

    `FilterAction.WARN` stays in the vocabulary, but no 0.7.1 code path demotes a rejection to a
    label, so resolving it would emit `action: warn, evaluated: true` for a candidate the gate failed.
    """
    with pytest.raises(RunPolicyError, match="which 0.7.1 does not apply"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={"max_off_target_count": "warn"})

    assert FilterAction.WARN not in SELECTABLE_ACTIONS
    assert set(SELECTABLE_ACTIONS) == {FilterAction.OFF, FilterAction.FAIL}


@pytest.mark.unit
def test_the_gate_that_no_code_reads_resolves_to_off():
    """max_mirna_1mm_seed carries a threshold of 10 and is compared by nothing in 0.7.1."""
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    descriptor = policy.descriptor("max_mirna_1mm_seed")

    assert descriptor.threshold == 10
    assert descriptor.action is FilterAction.OFF


@pytest.mark.unit
def test_the_human_stratified_gates_declare_their_scope_and_their_missing_column():
    """Six gates read a human-stratified counter that no candidate column exports (#101's work)."""
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)

    stratified = [resolved for resolved in policy.filters if resolved.descriptor.scope.species == frozenset({"human"})]
    assert {resolved.filter_id for resolved in stratified} == {
        "max_transcriptome_hits_0mm",
        "max_transcriptome_hits_1mm",
        "max_transcriptome_hits_2mm",
        "max_mirna_perfect_seed",
        "fail_on_high_risk_mirna",
        "max_total_offtarget_hits",
    }
    # None of their counters is in the candidate CSV, so a client cannot re-apply them and get the
    # pipeline's answer -- said in the descriptor rather than discovered as a disagreement.
    assert all(resolved.evidence_exported is False for resolved in stratified)
    assert all("HUMAN-STRATIFIED" not in resolved.definition for resolved in stratified)
    assert all("all-species number" in resolved.definition for resolved in stratified)


@pytest.mark.unit
def test_evidence_exported_agrees_with_the_columns_the_candidate_row_actually_carries():
    """`evidence_exported` is the one flag #103 has to trust, and it lives away from the row builder.

    A column renamed in ``build_candidate_row`` would turn ``evidence_exported=True`` into a false
    claim that a client can re-apply the gate, with nothing to catch it.
    """
    published = set(build_candidate_row(_toy_candidate()))
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)

    for resolved in policy.filters:
        exported = resolved.descriptor.column in published
        assert resolved.evidence_exported == exported, (
            f"{resolved.filter_id} declares evidence_exported={resolved.evidence_exported} but its "
            f"column {resolved.descriptor.column!r} is {'in' if exported else 'not in'} the candidate row"
        )


@pytest.mark.unit
def test_each_human_stratified_gate_names_its_own_stratification_convention():
    """The six do not stratify the same way, and one shared sentence made three of them false.

    ``_process_nextflow_results`` counts a transcriptome hit as human when its species is human *or*
    the label is blank, but counts a miRNA hit only when ``is_human_species(label)`` is True -- and
    ``is_human_species(None)`` is False. ``max_total_offtarget_hits`` sums the two.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)

    for filter_id in ("max_transcriptome_hits_0mm", "max_transcriptome_hits_1mm", "max_transcriptome_hits_2mm"):
        definition = policy.filter(filter_id).definition
        assert "human or unlabelled" in definition, filter_id

    for filter_id in ("max_mirna_perfect_seed", "fail_on_high_risk_mirna"):
        definition = policy.filter(filter_id).definition
        assert "human or unlabelled" not in definition, f"{filter_id} counts human-labelled hits only"
        assert "labelled human only" in definition, filter_id

    combined = policy.filter("max_total_offtarget_hits").definition
    assert "different conventions" in combined


@pytest.mark.unit
def test_a_boolean_flag_is_expressed_as_a_ceiling_of_zero():
    """fail_on_high_risk_mirna has no threshold field, so as data it is "at most zero hits"."""
    on = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).descriptor("fail_on_high_risk_mirna")
    off = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"fail_on_high_risk_mirna": False}
    ).descriptor("fail_on_high_risk_mirna")

    assert on.threshold == 0
    assert on.action is FilterAction.FAIL
    assert on.comparator.passes(0, on.threshold) and not on.comparator.passes(1, on.threshold)
    assert off.action is FilterAction.OFF


@pytest.mark.unit
def test_a_resolved_descriptor_reproduces_the_verdict_it_reports():
    """#103 re-thresholds client-side, so the descriptor a run publishes must decide as the run did."""
    descriptor = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).descriptor("max_off_target_count")

    assert FilterVerdict(descriptor=descriptor, observed=15, evaluation=FilterEvaluation.PASS)
    assert FilterVerdict(descriptor=descriptor, observed=16, evaluation=FilterEvaluation.FAIL)
    with pytest.raises(ValidationError):
        FilterVerdict(descriptor=descriptor, observed=16, evaluation=FilterEvaluation.PASS)


# --------------------------------------------------------------------------------------
# Failing before expensive work
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_an_invalid_bound_fails_in_the_resolver():
    """Constructed, not copied, so every field bound and cross-field rule still applies."""
    with pytest.raises(RunPolicyError, match="gc_max"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"gc_min": 70.0, "gc_max": 40.0})
    with pytest.raises(RunPolicyError, match="max_bp_span"):
        resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"plfold_window": 40, "plfold_max_bp_span": 100}
        )
    with pytest.raises(RunPolicyError, match="min_empirical_score"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"min_empirical_score": 0.9})


@pytest.mark.unit
def test_an_undeclared_term_is_rejected_rather_than_ignored(tmp_path):
    """A silently ignored setting reads as an applied override, which is the worst of both."""
    with pytest.raises(RunPolicyError, match="does not have"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"gc_maximum": 60.0})

    config = tmp_path / "policy.json"
    config.write_text('{"settings": {"gc_maximum": 60.0}}')
    with pytest.raises(RunPolicyError, match="does not have"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, config_file=config)

    config.write_text('{"gc_max": 60.0}')
    with pytest.raises(RunPolicyError, match="unknown top-level keys"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, config_file=config)

    with pytest.raises(RunPolicyError, match="this build does not have"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={"max_gc": "off"})


@pytest.mark.unit
def test_a_non_finite_weight_fails_in_the_resolver():
    """A NaN weight must not survive to produce a NaN score at the end of a screen."""
    with pytest.raises(ValidationError):
        ScoringWeights(postscreen_sirna=PostScreenSiRNAWeights(off_target=float("nan")))

    broken = ScoringWeights()
    object.__setattr__(broken.postscreen_sirna, "off_target", float("inf"))
    parameters = DesignParameters(scoring=broken)
    with pytest.raises(RunPolicyError, match="non-finite weight"):
        describe_parameters(parameters, entry_point=EntryPoint.SCREENING_WORKFLOW)


@pytest.mark.unit
def test_the_cli_reports_a_cross_field_error_and_creates_no_output_directory(tmp_path):
    """#95 item 2: this escaped as an unhandled traceback, after the log directory had been made."""
    output = tmp_path / "never"
    result = CliRunner().invoke(
        app,
        [
            "workflow",
            "TOY",
            "--input-fasta",
            str(_fasta(tmp_path)),
            "--output-dir",
            str(output),
            "--plfold-window",
            "40",
            "--plfold-max-bp-span",
            "100",
        ],
    )

    assert result.exit_code == 1
    assert "max_bp_span" in result.output
    assert "Traceback" not in result.output
    assert not output.exists(), "validation must happen before the output tree is created"


@pytest.mark.unit
def test_the_design_command_reports_the_same_error_without_a_traceback(tmp_path):
    """The `design` path built its models outside any try block, so it raised the raw ValidationError."""
    result = CliRunner().invoke(
        app,
        [
            "design",
            str(_fasta(tmp_path)),
            "--output",
            str(tmp_path / "out.csv"),
            "--plfold-window",
            "40",
            "--plfold-max-bp-span",
            "100",
        ],
    )

    assert result.exit_code == 1
    assert "max_bp_span" in result.output
    assert "Traceback" not in result.output


@pytest.mark.unit
def test_a_bad_filter_action_is_reported_by_name(tmp_path):
    """A typo in --filter-action must name the filter and the vocabulary, not raise."""
    result = CliRunner().invoke(
        app,
        ["design", str(_fasta(tmp_path)), "--filter-action", "min_asymmetry_score=nope"],
    )

    assert result.exit_code == 1
    assert "min_asymmetry_score" in result.output


# --------------------------------------------------------------------------------------
# Profile identity, and defaults that cannot drift
# --------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_legacy_profile_is_derived_from_the_model_defaults_with_one_declared_exception():
    """This is the mechanism that stops documented and actual defaults drifting apart again."""
    profile = BUILTIN_PROFILES[DEFAULT_PROFILE_NAME]
    differing = {spec.key for spec in SETTING_SPECS if profile.baseline[spec.key] != _model_default(spec)}

    assert differing == set(LEGACY_DEFAULT_EXCEPTIONS)
    assert differing == {"gc_min"}
    assert profile.baseline["gc_min"] == 30.0
    assert profile.exceptions["gc_min"]
    assert profile.experimental is True


@pytest.mark.unit
def test_the_legacy_profile_carries_the_shipped_off_target_cap_and_design_weights():
    """The two documented-versus-actual drifts found at 8dce4ae, pinned on the code's side."""
    assert default_for("max_off_target_count") == 15
    weights = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).design_parameters.scoring.design
    assert (weights.asymmetry, weights.target_accessibility) == (0.40, 0.35)


@pytest.mark.unit
def test_the_profile_hash_changes_with_the_baseline_and_not_with_the_run():
    """Two runs on one profile compare; a profile edited in place must not claim the old identity."""
    first = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    second = resolve_run_policy(entry_point=EntryPoint.DESIGN_COMMAND, stated={"gc_max": 51.0})

    assert first.profile.content_hash == second.profile.content_hash
    assert first.profile.content_hash.startswith("sha256:")

    edited = BUILTIN_PROFILES[DEFAULT_PROFILE_NAME]
    mutated = type(edited)(
        name=edited.name,
        version=edited.version,
        description=edited.description,
        experimental=edited.experimental,
        baseline={**edited.baseline, "gc_min": 31.0},
        exceptions=edited.exceptions,
    )
    assert mutated.identity().content_hash != edited.identity().content_hash


@pytest.mark.unit
def test_an_unknown_profile_is_an_error():
    """Silently falling back to the default profile would misreport which numbers applied."""
    with pytest.raises(RunPolicyError, match="unknown profile"):
        resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, profile_name="calibrated")


@pytest.mark.unit
def test_the_resolved_policy_is_immutable():
    """Resolved once, read everywhere: a mutable policy could differ between two readers."""
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    with pytest.raises(ValidationError):
        policy.run_mode = RunMode.DESIGN_ONLY


@pytest.mark.unit
def test_resolving_twice_with_the_same_inputs_gives_the_same_answer():
    """Precedence is applied once per call and depends on nothing outside the arguments."""
    kwargs: dict[str, Any] = {
        "entry_point": EntryPoint.SCREENING_WORKFLOW,
        "design_mode": "mirna",
        "stated": {"gc_max": 58.0},
        "query_species": "human",
        "screen_species": ["human", "mouse"],
    }
    assert resolve_run_policy(**kwargs).as_manifest() == resolve_run_policy(**kwargs).as_manifest()
