"""Tests for the SelectionState vocabulary and its schema enforcement (#100).

Before this module existed, workflow.py spelled three ``selection_state`` values as private
string literals while the model field's own default ("not_selected") was a fourth value its
docstring never mentioned, and nothing validated any of the four against ``candidates_all.csv``.
Config.strict=False plus add_missing_columns=True let an undeclared column -- or a misspelled
value in a declared one -- through validation silently.
"""

from typing import Any

import pandas as pd
import pandera
import pytest

from sirnaforge.models.policy import DECLARED_FILTER_IDS, FilterEvaluation
from sirnaforge.models.schemas import SiRNACandidateSchema
from sirnaforge.models.sirna import SelectionState, SiRNACandidate, build_candidate_row


def _minimal_row(**overrides: Any) -> dict[str, Any]:
    """One valid SiRNACandidateSchema row, with every non-nullable column filled in."""
    row: dict[str, Any] = {
        "id": "test_1",
        "transcript_id": "ENST00000000001",
        "position": 100,
        "guide_sequence": "TTTTTTTTTTTTTTTTTTTTT",
        "passenger_sequence": "AAAAAAAAAAAAAAAAAAAAA",
        "gc_content": 40.0,
        "asymmetry_score": 0.7,
        "paired_fraction": 0.3,
        "off_target_count": 0,
        "transcript_hit_count": 1,
        "transcript_hit_fraction": 1.0,
        "passes_filters": True,
        # Series[Any] columns hit a pandera add_missing_columns bug (dtype.try_coerce on a
        # None dtype) when left for auto-fill, so every row supplies them explicitly.
        "structure": None,
    }
    row.update(overrides)
    return row


def _validate(**overrides: Any) -> pd.DataFrame:
    df = pd.DataFrame([_minimal_row(**overrides)])
    return SiRNACandidateSchema.validate(df)


# A single-row DataFrame built straight from build_candidate_row's dict infers these
# all-None columns as `object`, which pandera's coercion does not resolve to the declared
# nullable dtype; test_schemas.py works around the same columns the same way. Pre-existing,
# unrelated to this slice's fix.
_NULLABLE_FLOAT_WORKAROUND_COLUMNS = (
    "mfe",
    "duplex_stability_dg",
    "duplex_stability_score",
    "dg_5p",
    "dg_3p",
    "delta_dg_end",
    "melting_temp_c",
)
_NULLABLE_INT_WORKAROUND_COLUMNS = ("seed_7mer_hits", "seed_8mer_hits")


def _coerce_known_nullable_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.astype(dict.fromkeys(_NULLABLE_FLOAT_WORKAROUND_COLUMNS, "float64"))
    return df.astype(dict.fromkeys(_NULLABLE_INT_WORKAROUND_COLUMNS, "Int64"))


class TestSelectionStateEnum:
    """The vocabulary itself."""

    def test_five_values(self):
        """The five states are exactly eligible/provisional/withheld/not_eligible/not_selected."""
        assert {s.value for s in SelectionState} == {
            "eligible",
            "provisional_incomplete_evidence",
            "withheld_incomplete_evidence",
            "not_eligible",
            "not_selected",
        }

    def test_model_default_is_not_selected(self):
        """The model field's default must be a real SelectionState member, not a bare literal."""
        assert SiRNACandidate.model_fields["selection_state"].default == SelectionState.NOT_SELECTED.value


class TestSchemaDeclaresSelectionState:
    """Acceptance: 'all reasons survive CSV/JSON export' -- the column must be declared, not merely tolerated as an undeclared extra."""

    def test_selection_state_is_a_declared_column(self):
        """Must go red without the fix: today the column is not declared at all."""
        assert "selection_state" in SiRNACandidateSchema.to_schema().columns

    def test_valid_selection_state_passes(self):
        """Every SelectionState member is an accepted value."""
        for state in SelectionState:
            result = _validate(selection_state=state.value)
            assert result.loc[0, "selection_state"] == state.value

    def test_misspelled_selection_state_rejected(self):
        """Must go red without the fix: today this silently validates."""
        with pytest.raises(pandera.errors.SchemaError):
            _validate(selection_state="elligible")

    def test_absent_selection_state_column_still_validates(self):
        """An older run's CSV without this column at all must not be broken by declaring it."""
        df = pd.DataFrame([_minimal_row()])
        result = SiRNACandidateSchema.validate(df)
        assert "selection_state" not in result.columns


class TestFilterVerdictAllowList:
    """Every present ``{filter_id}_verdict`` column is checked against FilterEvaluation, mirroring check_passes_filters_values."""

    def test_valid_verdict_values_pass(self):
        """A declared filter's verdict column accepts a genuine FilterEvaluation value."""
        row = _minimal_row(max_off_target_count_verdict=FilterEvaluation.PASS.value)
        result = _validate(**row)
        assert result.loc[0, "max_off_target_count_verdict"] == FilterEvaluation.PASS.value

    def test_misspelled_verdict_rejected(self):
        """Must go red without the fix: today this silently validates."""
        with pytest.raises(pandera.errors.SchemaError):
            _validate(max_off_target_count_verdict="passed")

    def test_every_declared_filter_id_checked(self):
        """No filter's verdict column is exempt from the allow-list."""
        for filter_id in DECLARED_FILTER_IDS:
            with pytest.raises(pandera.errors.SchemaError):
                _validate(**{f"{filter_id}_verdict": "bogus"})

    def test_unrelated_column_named_like_a_verdict_is_ignored(self):
        """Only columns keyed on a DECLARED_FILTER_IDS entry are checked."""
        result = _validate(some_other_column_verdict="bogus")
        assert result.loc[0, "some_other_column_verdict"] == "bogus"


class TestBuildCandidateRowRoundTrips:
    """The producer side: build_candidate_row must emit values the schema accepts."""

    def test_default_candidate_row_validates(self):
        """A freshly-built candidate's default selection_state and every NOT_EVALUATED verdict validate unchanged."""
        candidate = SiRNACandidate(
            id="c1",
            transcript_id="ENST00000000001",
            position=1,
            guide_sequence="TTTTTTTTTTTTTTTTTTTTT",
            passenger_sequence="AAAAAAAAAAAAAAAAAAAAA",
            gc_content=40.0,
            length=21,
            asymmetry_score=0.7,
        )
        row = build_candidate_row(candidate)
        assert row["selection_state"] == SelectionState.NOT_SELECTED.value
        for filter_id in DECLARED_FILTER_IDS:
            assert row[f"{filter_id}_verdict"] == FilterEvaluation.NOT_EVALUATED.value

        df = _coerce_known_nullable_dtypes(pd.DataFrame([row]))
        SiRNACandidateSchema.validate(df)

    def test_every_selection_state_round_trips_through_the_schema(self):
        """Each of the five states, set on a real candidate, survives build_candidate_row and validation unchanged."""
        for state in SelectionState:
            candidate = SiRNACandidate(
                id="c1",
                transcript_id="ENST00000000001",
                position=1,
                guide_sequence="TTTTTTTTTTTTTTTTTTTTT",
                passenger_sequence="AAAAAAAAAAAAAAAAAAAAA",
                gc_content=40.0,
                length=21,
                asymmetry_score=0.7,
                selection_state=state.value,
            )
            row = build_candidate_row(candidate)
            df = _coerce_known_nullable_dtypes(pd.DataFrame([row]))
            validated = SiRNACandidateSchema.validate(df)
            assert validated.loc[0, "selection_state"] == state.value


class TestObservedColumnsDeclaredAsNullableFloat:
    """'leave {filter_id}_observed as nullable float' -- declared, untyped vocabulary, just a type."""

    def test_every_declared_filter_id_has_an_observed_column(self):
        """No filter's observed column is missing from the declared schema."""
        columns = SiRNACandidateSchema.to_schema().columns
        for filter_id in DECLARED_FILTER_IDS:
            assert f"{filter_id}_observed" in columns

    def test_observed_column_accepts_a_float_and_null(self):
        """Observed columns carry no vocabulary check, only a float type."""
        result = _validate(max_off_target_count_observed=3.0)
        assert result.loc[0, "max_off_target_count_observed"] == 3.0

        result_null = _validate(max_off_target_count_observed=None)
        assert pd.isna(result_null.loc[0, "max_off_target_count_observed"])
