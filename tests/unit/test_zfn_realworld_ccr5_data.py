"""Validation tests for integrated CCR5 ZFN real-world study facts and templates."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def _data_dir() -> Path:
    """Return the unit-test data directory for ZFN benchmark fixtures."""
    return Path(__file__).parent / "data" / "zfn"


def test_ccr5_study_facts_include_expected_assertions() -> None:
    """Validate integrated CCR5 study facts include all expected assertion IDs."""
    payload = json.loads((_data_dir() / "ccr5_study_facts.json").read_text(encoding="utf-8"))

    assert payload["dataset_id"] == "prognos_ccr5_extractable_facts"
    assert payload["nuclease"] == "CCR5 ZFNs"

    facts = payload["facts"]
    ids = {fact["dataset_id"] for fact in facts}
    assert {
        "prognos_ccr5_validated_sites",
        "prognos_ccr5_prediction_recovery_3x",
        "prognos_ccr5_novel_site",
        "prognos_ranking_tiebreak",
        "prognos_site_sequence_orientation_rule",
    }.issubset(ids)

    validated = next(f for f in facts if f["dataset_id"] == "prognos_ccr5_validated_sites")
    assert "12" in validated["assertion"]

    tiebreak = next(f for f in facts if f["dataset_id"] == "prognos_ranking_tiebreak")
    assert "Exon > Promoter > Intron > Intergenic" in tiebreak["assertion"]


def test_ccr5_design_parameter_template_has_expected_columns() -> None:
    """Validate CCR5 design-parameter template schema for future table ingestion."""
    with (_data_dir() / "ccr5_zfn_design_parameters_template.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "dataset_id",
            "source",
            "left_half_site",
            "right_half_site",
            "max_mismatches",
            "seed_len_from_fokI",
            "seed_max_mismatches",
            "allowed_spacer_lengths",
            "dimer_mode",
            "algorithm",
            "notes",
        ]
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["dataset_id"] == "prognos_ccr5_zfn"


def test_ccr5_golden_sites_template_has_expected_columns() -> None:
    """Validate CCR5 golden-sites template schema for result alignment tests."""
    with (_data_dir() / "ccr5_zfn_golden_sites_template.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "dataset_id",
            "source",
            "site_id",
            "chrom",
            "start_1based",
            "end_1based",
            "orientation",
            "spacer_len",
            "sequence",
            "score",
            "region",
            "nearest_gene",
            "validation_smrt_percent",
            "p_value",
            "rank_homology",
            "rank_conserved_g",
            "rank_zfn_v2",
        ]
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["dataset_id"] == "prognos_ccr5_sites"


def test_ccr5_s10_visible_rows_fixture_schema_and_count() -> None:
    """Validate S10 visible-row fixture schema and row count from extracted text."""
    with (_data_dir() / "ccr5_s10_visible_rows.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "Closest gene",
            "Match type",
            "hg19 coordinate",
            "(+) half-site",
            "(−) half-site",
            "Empty indels",
            "Empty total",
            "Active indels",
            "Active total",
            "Active mutation freq",
            "p-value",
            "Notes",
        ]
        rows = list(reader)

    assert len(rows) == 17

    genes = {row["Closest gene"] for row in rows}
    assert {"CCR5", "LPCAT2", "CSNK1G3"}.issubset(genes)

    csnk1g3 = next(row for row in rows if row["Closest gene"] == "CSNK1G3")
    assert csnk1g3["Match type"] == "L-5-R"
    assert csnk1g3["hg19 coordinate"] == "chr5:123393701"
    assert csnk1g3["Active indels"] == "17"
    assert csnk1g3["Active mutation freq"] == "0.086%"
    assert csnk1g3["p-value"] == "0.000019"


def test_ccr5_s11_homology_visible_rows_fixture_schema_and_count() -> None:
    """Validate S11 homology fixture schema and row count from extracted text."""
    with (_data_dir() / "ccr5_s11_homology_visible_rows.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "Homology rank",
            "T mism",
            "+ mism",
            "− mism",
            "Interrogated by",
            "Closest gene",
            "Match type",
            "hg19 coordinate",
            "(+) half-site",
            "(−) half-site",
        ]
        rows = list(reader)

    assert len(rows) == 32
    first_row = rows[0]
    last_row = rows[-1]
    assert first_row["Homology rank"] == "2"
    assert last_row["Homology rank"] == "33"
    assert last_row["Closest gene"] == "CSNK1G3"
    assert last_row["Match type"] == "L-5-R"


def test_ccr5_annotation_model_contains_validated_benchmark_context() -> None:
    """Validate benchmark annotation metadata for real-world CCR5 context."""
    payload = json.loads((_data_dir() / "ccr5_benchmark_annotations.json").read_text(encoding="utf-8"))

    assert payload["dataset_id"] == "prognos_ccr5_benchmark_annotations"
    assert payload["genome_build"] == "hg19"
    assert payload["design"]["fingers_per_monomer"] == 4
    assert payload["performance"]["validated_offtargets_prior_methods"] == 12
    assert payload["performance"]["prognos_recovery_at_3x"] == "10/12"
    assert payload["performance"]["novel_site_near_gene"] == "CSNK1G3"
    assert payload["ranking_policy"]["tie_break_order"] == ["exon", "promoter", "intron", "intergenic"]
    assert payload["fixture_scope"]["s10_visible_rows"] == 17
    assert payload["fixture_scope"]["s11_homology_visible_rows"] == 32
