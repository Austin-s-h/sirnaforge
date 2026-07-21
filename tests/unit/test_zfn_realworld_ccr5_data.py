"""Validation tests for integrated CCR5 ZFN real-world study facts and templates."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from Bio.Seq import Seq

from sirnaforge.models.zfn import (
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNSearchBackend,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.design import ZFNDesigner


def _data_dir() -> Path:
    """Return the unit-test data directory for ZFN benchmark fixtures."""
    return Path(__file__).parent / "data" / "zfn"


def _match_type_spacer(match_type: str) -> int:
    """Extract spacer length from compact match-type labels like ``L-5-R``."""
    return int(match_type.split("-")[1])


def _genomic_site_from_visible_row(plus_half_site: str, minus_half_site: str, spacer_len: int) -> str:
    """Build genomic site string for search-space FASTA from S10 visible row halves.

    The ZFN searcher expects ``right_half_site`` in model orientation and compares
    against reverse-complemented genomic sequence on the right half.
    """
    plus = plus_half_site.upper()
    minus = minus_half_site.upper()
    right_genomic = str(Seq(minus).reverse_complement())
    return f"{plus}{'A' * spacer_len}{right_genomic}"


def test_00_ccr5_sirnaforge_run_sanity_against_visible_rows(tmp_path: Path) -> None:
    """Run sirnaforge ZFN design on synthetic CCR5 subset and compare to visible real-data rows."""
    s10_path = _data_dir() / "ccr5_s10_visible_rows.csv"
    with s10_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    ccr5_row = next(row for row in rows if row["Closest gene"] == "CCR5")
    csnk1g3_row = next(row for row in rows if row["Closest gene"] == "CSNK1G3")

    ccr5_spacer = _match_type_spacer(ccr5_row["Match type"])
    csnk1g3_spacer = _match_type_spacer(csnk1g3_row["Match type"])

    synthetic_genome = (
        "TTT"
        + _genomic_site_from_visible_row(ccr5_row["(+) half-site"], ccr5_row["(−) half-site"], ccr5_spacer)
        + "GGGG"
        + _genomic_site_from_visible_row(
            csnk1g3_row["(+) half-site"],
            csnk1g3_row["(−) half-site"],
            csnk1g3_spacer,
        )
        + "CCC"
    )

    fasta_path = tmp_path / "ccr5_visible_subset.fa"
    fasta_path.write_text(f">chr5\n{synthetic_genome}\n", encoding="utf-8")

    params = ZFNDesignParameters(
        left_half_site=ccr5_row["(+) half-site"],
        right_half_site=ccr5_row["(−) half-site"],
        search_space_fasta=str(fasta_path),
        search_backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=4),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5, 6]),
        top_n_sites=100,
    )

    result = ZFNDesigner().evaluate_pair(params)
    sites = result.off_target_sites
    assert sites, "Expected at least one predicted site from synthetic CCR5 subset"

    assert any(site.left_mismatches == 0 and site.right_mismatches == 0 for site in sites)

    assert any(site.spacer_len == 5 and site.left_mismatches >= 1 and site.right_mismatches >= 1 for site in sites)


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
