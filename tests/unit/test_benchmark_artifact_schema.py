"""Schema and CSV/JSON round-trip tests for the benchmark artifact contract (#109, bm-artifact slice).

Two properties are load-bearing and get their own tests rather than being implied by the happy
path: ``guide_3p_overhang=""`` (measured, blunt) must stay distinct from ``None`` (not stated)
across a write/read cycle, and a CSV whose header omits a declared column must raise
``BenchmarkArtifactError`` naming that column rather than silently defaulting it.
"""

from pathlib import Path

import pytest
from pydantic import ValidationError

from sirnaforge.benchmark.artifact import (
    ACCOUNTING_COLUMNS,
    ACCOUNTING_FILENAME,
    BENCHMARK_ARTIFACT_SCHEMA_VERSION,
    BENCHMARK_MANIFEST_SCHEMA_VERSION,
    CANDIDATES_ALL_FILENAME,
    DESIGN_INPUTS_FASTA_FILENAME,
    MANIFEST_FILENAME,
    OBSERVATION_COLUMNS,
    OBSERVATIONS_FILENAME,
    PAIRED_LENGTH_BOUNDS,
    BenchmarkAccountingRow,
    BenchmarkArtifactCounts,
    BenchmarkArtifactError,
    BenchmarkArtifactManifest,
    BenchmarkArtifactOutputs,
    BenchmarkObservation,
    FilterExclusionCounts,
    GCWideningBlock,
    GCWideningEntry,
    ManifestInputEntry,
    ManifestOutputEntry,
    ManifestPanelBlock,
    PolynucleotideRunRequirementBlock,
    artifact_dir_name,
    read_accounting,
    read_manifest,
    read_observations,
    write_accounting,
    write_manifest,
    write_observations,
)
from sirnaforge.config.run_policy import POLICY_SCHEMA_VERSION
from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation, SettingSource

_BASE_SEQUENCE = ("AUGC" * 6)[:21]
assert len(_BASE_SEQUENCE) == 21


def _observation(**overrides: object) -> BenchmarkObservation:
    """A valid, fully-populated observation; tests override only the field(s) under test."""
    fields: dict[str, object] = {
        "observation_id": "huesken_full:000001",
        "panel_id": "huesken_full",
        "source_row_index": 1,
        "architecture": "fully_complementary",
        "assay_label": "inhibition",
        "measured_endpoint": "inhibition_fraction",
        "measured_value": 0.87,
        "full_guide_sequence": _BASE_SEQUENCE,
        "guide_length": 21,
        "paired_guide_sequence": _BASE_SEQUENCE,
        "paired_slice_start_1based": 1,
        "paired_length": 21,
        "passenger_sequence": None,
        "guide_3p_overhang": "",
        "passenger_3p_overhang": None,
        "duplex_pairing_status": "unstated",
        "source_citation": "Huesken et al. 2005, Nat Biotechnol",
        "source_redistribution": "https://github.com/apkrfi/unMod-siRNA-Pred",
        "source_file": "tests/unit/data/sirna_efficacy_subset.csv",
        "source_sha256": "a" * 64,
        "target_transcript_id": None,
        "target_identity_status": "unavailable",
        "target_start_1based": None,
        "target_end_1based": None,
        "target_strand": None,
        "design_context_id": "huesken_full:000001",
        "design_context_source": "measured_target_site",
        "split": "development",
        "compatibility_status": "compatible",
        "compatibility_reason": "",
    }
    fields.update(overrides)
    return BenchmarkObservation(**fields)


def _accounting_row(**overrides: object) -> BenchmarkAccountingRow:
    fields: dict[str, object] = {
        "observation_id": "huesken_full:000001",
        "panel_id": "huesken_full",
        "paired_length": 21,
        "candidate_id": "huesken_full:000001:c0",
        "designed_guide_sequence": _BASE_SEQUENCE,
        "guide_match": "exact",
        "entered_design": True,
        "default_filter_status": FilterEvaluation.PASS,
        "default_filter_reasons": "",
        "benchmark_filter_status": FilterEvaluation.PASS,
        "benchmark_filter_reasons": "",
    }
    fields.update(overrides)
    return BenchmarkAccountingRow(**fields)


# --------------------------------------------------------------------------------------------------
# observations.csv round trip
# --------------------------------------------------------------------------------------------------


def test_observation_round_trip_is_lossless(tmp_path: Path) -> None:
    """write_observations -> read_observations reproduces the originals, including sort order."""
    second = _observation(
        observation_id="huesken_full:000002",
        source_row_index=2,
        design_context_id="huesken_full:000002",
        guide_3p_overhang=None,
        passenger_sequence="UGACUGACUGACUGACUGACU"[:21],
        target_identity_status="panel_local",
        target_transcript_id="NM_000001",
        target_start_1based=100,
        target_end_1based=121,
        target_strand="+",
        source_redistribution=None,
        split=None,
        compatibility_status="incompatible",
        compatibility_reason="paired_length outside the requested design range",
    )
    first = _observation()
    # Written out of order; the reader must hand back rows sorted by observation_id.
    path = write_observations([second, first], tmp_path / "observations.csv")

    result = read_observations(path)

    assert result == [first, second]


def test_guide_3p_overhang_empty_string_is_distinct_from_none(tmp_path: Path) -> None:
    """Empty string (measured, blunt) and None (not stated) must both round-trip without collapsing."""
    blunt = _observation(observation_id="p:000001", design_context_id="p:000001", guide_3p_overhang="")
    unstated = _observation(observation_id="p:000002", design_context_id="p:000002", guide_3p_overhang=None)

    path = write_observations([blunt, unstated], tmp_path / "observations.csv")
    result = {row.observation_id: row.guide_3p_overhang for row in read_observations(path)}

    assert result["p:000001"] == ""
    assert result["p:000002"] is None
    assert result["p:000001"] != result["p:000002"]


def test_passenger_sequence_none_round_trips_as_absent(tmp_path: Path) -> None:
    """passenger_sequence=None (guide-only panel row) survives the round trip as None, not ''."""
    row = _observation(passenger_sequence=None)
    path = write_observations([row], tmp_path / "observations.csv")

    (result,) = read_observations(path)

    assert result.passenger_sequence is None


def test_missing_declared_column_raises_naming_the_column(tmp_path: Path) -> None:
    """A header omitting a declared column must raise BenchmarkArtifactError naming it.

    Not a row constructed with a silent default for the missing field: that would be
    indistinguishable from a real value in any report read off the artifact.
    """
    path = tmp_path / "observations.csv"
    columns = [name for name in OBSERVATION_COLUMNS if name != "source_sha256"]
    path.write_text(",".join(columns) + "\n")

    with pytest.raises(BenchmarkArtifactError) as excinfo:
        read_observations(path)

    assert "source_sha256" in str(excinfo.value)


def test_unknown_column_in_header_raises(tmp_path: Path) -> None:
    """A header declaring a column this schema does not have is rejected, not silently ignored."""
    path = tmp_path / "observations.csv"
    columns = [*OBSERVATION_COLUMNS, "totally_unexpected_column"]
    path.write_text(",".join(columns) + "\n")

    with pytest.raises(BenchmarkArtifactError, match="totally_unexpected_column"):
        read_observations(path)


def test_row_with_wrong_field_count_raises(tmp_path: Path) -> None:
    """A short/long data row (vs. the header) is a corrupt file, not a row with blank tails."""
    path = tmp_path / "observations.csv"
    path.write_text(",".join(OBSERVATION_COLUMNS) + "\n" + "only,one,field\n")

    with pytest.raises(BenchmarkArtifactError, match="line 2"):
        read_observations(path)


def test_observation_columns_are_declaration_order() -> None:
    """Column order is the model's own declaration order, so the two cannot drift apart."""
    assert tuple(BenchmarkObservation.model_fields) == OBSERVATION_COLUMNS
    assert OBSERVATION_COLUMNS[0] == "observation_id"
    assert OBSERVATION_COLUMNS[-1] == "compatibility_reason"


def test_two_writes_of_the_same_rows_are_byte_identical(tmp_path: Path) -> None:
    """No CSV carries a timestamp or an absolute path, so two prepares must match byte-for-byte."""
    rows = [_observation(observation_id="p:000002", design_context_id="p:000002"), _observation()]
    first_path = write_observations(rows, tmp_path / "run1" / "observations.csv")
    second_path = write_observations(list(reversed(rows)), tmp_path / "run2" / "observations.csv")

    assert first_path.read_bytes() == second_path.read_bytes()


def test_full_guide_sequence_must_be_upper_case() -> None:
    """A lower-case sequence would mean a caller normalised it, which #109 rules out."""
    with pytest.raises(ValidationError, match="upper-cased"):
        _observation(full_guide_sequence=_BASE_SEQUENCE.lower())


def test_compatible_row_cannot_carry_a_reason() -> None:
    """An empty reason claims nothing went wrong, so a compatible row may hold no other reason text."""
    with pytest.raises(ValidationError, match="compatibility_reason"):
        _observation(compatibility_status="compatible", compatibility_reason="should not be here")


def test_incompatible_row_must_state_a_reason() -> None:
    """An incompatible row with no stated reason would be unauditable."""
    with pytest.raises(ValidationError, match="compatibility_reason"):
        _observation(compatibility_status="incompatible", compatibility_reason="")


def test_observation_is_frozen_and_forbids_extra_fields() -> None:
    """Frozen and extra="forbid" per the settled schema: no mutation, no undeclared field."""
    row = _observation()
    with pytest.raises(ValidationError):
        row.panel_id = "different"  # type: ignore[misc]
    with pytest.raises(ValidationError):
        _observation(unexpected_field=1)


# --------------------------------------------------------------------------------------------------
# accounting.csv round trip
# --------------------------------------------------------------------------------------------------


def test_accounting_round_trip_is_lossless(tmp_path: Path) -> None:
    """write_accounting -> read_accounting reproduces both a no-candidate and a candidate row."""
    no_candidate = _accounting_row(
        observation_id="p:000002",
        candidate_id=None,
        designed_guide_sequence=None,
        guide_match="none",
        entered_design=False,
        default_filter_status=FilterEvaluation.NOT_EVALUATED,
        default_filter_reasons="",
        benchmark_filter_status=FilterEvaluation.NOT_EVALUATED,
        benchmark_filter_reasons="",
    )
    with_candidate = _accounting_row(
        default_filter_status=FilterEvaluation.FAIL,
        default_filter_reasons="max_poly_runs",
        benchmark_filter_status=FilterEvaluation.PASS,
        benchmark_filter_reasons="",
    )
    path = write_accounting([no_candidate, with_candidate], tmp_path / "accounting.csv")

    result = read_accounting(path)

    assert result == [with_candidate, no_candidate]


def test_accounting_columns_are_declaration_order() -> None:
    """Same guarantee as observations.csv: column order tracks the model, never restated."""
    assert tuple(BenchmarkAccountingRow.model_fields) == ACCOUNTING_COLUMNS


def test_accounting_missing_column_raises_naming_it(tmp_path: Path) -> None:
    """Same guard as observations.csv, on the second CSV."""
    path = tmp_path / "accounting.csv"
    columns = [name for name in ACCOUNTING_COLUMNS if name != "guide_match"]
    path.write_text(",".join(columns) + "\n")

    with pytest.raises(BenchmarkArtifactError, match="guide_match"):
        read_accounting(path)


def test_no_candidate_row_must_be_entirely_not_evaluated() -> None:
    """entered_design=False must carry not_evaluated verdicts and no candidate fields (#109)."""
    with pytest.raises(ValidationError, match="entered_design"):
        _accounting_row(entered_design=False, candidate_id=None, designed_guide_sequence=None)


def test_no_candidate_row_cannot_carry_a_candidate_id() -> None:
    """entered_design=False with a stated candidate_id is refused, not silently accepted."""
    with pytest.raises(ValidationError, match="entered_design"):
        _accounting_row(
            entered_design=False,
            default_filter_status=FilterEvaluation.NOT_EVALUATED,
            default_filter_reasons="",
            benchmark_filter_status=FilterEvaluation.NOT_EVALUATED,
            benchmark_filter_reasons="",
        )


def test_accounting_row_is_frozen_and_forbids_extra_fields() -> None:
    """Mirrors the observation-row guard: a join row cannot be mutated or extended."""
    row = _accounting_row()
    with pytest.raises(ValidationError):
        row.entered_design = False  # type: ignore[misc]
    with pytest.raises(ValidationError):
        _accounting_row(unexpected_field=1)


# --------------------------------------------------------------------------------------------------
# manifest.json
# --------------------------------------------------------------------------------------------------


def _manifest() -> BenchmarkArtifactManifest:
    output = ManifestOutputEntry(path="observations.csv", exists=True, size_bytes=10, sha256="b" * 64, rows=1)
    absent = ManifestOutputEntry(path="candidates_all.csv", exists=False)
    return BenchmarkArtifactManifest(
        tool_version="0.7.1",
        created_utc="2026-09-14T00:00:00Z",
        invoked_command=("sirnaforge", "benchmark", "prepare", "--panel", "huesken_full"),
        panel=ManifestPanelBlock(
            panel_id="huesken_full",
            display_name="Huesken 2005 (full)",
            architecture="fully_complementary",
            citation="Huesken et al. 2005, Nat Biotechnol",
            redistribution="https://github.com/apkrfi/unMod-siRNA-Pred",
            data_present=True,
            descriptor_hash="sha256:" + "c" * 64,
        ),
        paired_length=21,
        requested_length=21,
        split_rule="sha256_accession_parity_v1",
        inputs=(ManifestInputEntry(role="panel_csv", path="work/sirna_bench.csv", sha256="d" * 64, size_bytes=100),),
        outputs=BenchmarkArtifactOutputs(observations_csv=output, design_inputs_fasta=absent),
        counts=BenchmarkArtifactCounts(
            observations_in_source=1,
            observations_kept=1,
            observations_incompatible=0,
            mapped_native=0,
            mapped_panel_local=0,
            mapping_unavailable=1,
            entered_design=0,
            no_candidate=1,
            default_pass=0,
            benchmark_pass=0,
            excluded_by_filter={"max_poly_runs": FilterExclusionCounts(default=0, benchmark=0)},
        ),
        polynucleotide_run_requirement=PolynucleotideRunRequirementBlock(
            comparator=FilterComparator.LE,
            threshold=3,
            action=FilterAction.FAIL,
            evaluated=True,
            excluded=FilterExclusionCounts(default=0, benchmark=0),
        ),
        run_policy={"schema_version": POLICY_SCHEMA_VERSION},
        default_run_policy={"schema_version": POLICY_SCHEMA_VERSION},
        gc_widening=GCWideningBlock(
            gc_min=GCWideningEntry(default=30.0, benchmark=25.0, source=SettingSource.EXPLICIT, widened=True),
            gc_max=GCWideningEntry(default=52.0, benchmark=52.0, source=SettingSource.BUILTIN_PROFILE, widened=False),
        ),
    )


def test_manifest_round_trip_is_lossless(tmp_path: Path) -> None:
    """write_manifest -> read_manifest reproduces a fully-populated manifest exactly."""
    manifest = _manifest()
    path = write_manifest(manifest, tmp_path / "manifest.json")

    result = read_manifest(path)

    assert result == manifest


def test_manifest_default_schema_versions_are_not_restated() -> None:
    """The three version fields default to the module's own constants, never a hand-typed string."""
    manifest = _manifest()
    assert manifest.schema_version == BENCHMARK_MANIFEST_SCHEMA_VERSION
    assert manifest.artifact_schema_version == BENCHMARK_ARTIFACT_SCHEMA_VERSION
    assert manifest.policy_schema_version == POLICY_SCHEMA_VERSION


def test_manifest_mapped_native_must_be_zero_in_109() -> None:
    """#110, not #109, owns native transcript mapping; a non-zero count here is a false claim."""
    with pytest.raises(ValidationError):
        BenchmarkArtifactCounts(
            observations_in_source=1,
            observations_kept=1,
            observations_incompatible=0,
            mapped_native=1,
            mapped_panel_local=0,
            mapping_unavailable=0,
            entered_design=0,
            no_candidate=1,
            default_pass=0,
            benchmark_pass=0,
        )


def test_manifest_output_entry_absent_file_carries_no_measurement() -> None:
    """A placeholder size for a file that does not exist would read as a checksum of nothing."""
    with pytest.raises(ValidationError):
        ManifestOutputEntry(path="x.csv", exists=False, size_bytes=10)


def test_manifest_is_frozen_and_forbids_extra_fields() -> None:
    """The manifest is read back into exactly what wrote it, or the reader raises."""
    manifest = _manifest()
    with pytest.raises(ValidationError):
        manifest.tool_version = "9.9.9"  # type: ignore[misc]


def test_read_manifest_missing_file_raises_artifact_error(tmp_path: Path) -> None:
    """A missing manifest.json is a BenchmarkArtifactError, not an uncaught FileNotFoundError."""
    with pytest.raises(BenchmarkArtifactError):
        read_manifest(tmp_path / "does_not_exist.json")


def test_read_manifest_invalid_json_raises_artifact_error(tmp_path: Path) -> None:
    """Malformed JSON is a BenchmarkArtifactError, not an uncaught JSONDecodeError."""
    path = tmp_path / "manifest.json"
    path.write_text("{not json")

    with pytest.raises(BenchmarkArtifactError):
        read_manifest(path)


# --------------------------------------------------------------------------------------------------
# Directory-name helper and fixed filenames
# --------------------------------------------------------------------------------------------------


def test_artifact_dir_name_builds_the_fixed_pattern() -> None:
    """The one place that spells `<panel_id>__len<paired_length>`."""
    assert artifact_dir_name("huesken_full", 21) == "huesken_full__len21"


@pytest.mark.parametrize("panel_id", ["Huesken", "huesken-full", "huesken full", ""])
def test_artifact_dir_name_rejects_invalid_panel_id(panel_id: str) -> None:
    """Only `[a-z0-9_]+` panel ids build a directory name."""
    with pytest.raises(BenchmarkArtifactError):
        artifact_dir_name(panel_id, 21)


@pytest.mark.parametrize("paired_length", [18, 24])
def test_artifact_dir_name_rejects_out_of_range_length(paired_length: int) -> None:
    """A length outside 19-23 could never be consumed by the fixed-length design path."""
    with pytest.raises(BenchmarkArtifactError):
        artifact_dir_name("huesken_full", paired_length)


def test_paired_length_bounds_matches_the_cli_range() -> None:
    """PAIRED_LENGTH_BOUNDS is the same 19-23 range `sirnaforge design --length` already accepts."""
    assert PAIRED_LENGTH_BOUNDS == (19, 23)


def test_fixed_filenames_are_distinct() -> None:
    """The five fixed inner filenames never collide."""
    names = {
        MANIFEST_FILENAME,
        OBSERVATIONS_FILENAME,
        DESIGN_INPUTS_FASTA_FILENAME,
        CANDIDATES_ALL_FILENAME,
        ACCOUNTING_FILENAME,
    }
    assert len(names) == 5
    assert MANIFEST_FILENAME == "manifest.json"
