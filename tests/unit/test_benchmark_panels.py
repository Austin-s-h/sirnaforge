"""Benchmark panel registry: architectures, compatibility verdicts, honestly-labelled fixtures.

Issue #109. Fixture provenance and the checksum table live in
``tests/unit/data/benchmark/README.md``; the predeclared-split accession lists this file pins
``predeclared_split`` against live in ``tests/unit/data/README.md`` (issues #97, #102).
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from pydantic import ValidationError

from sirnaforge.benchmark.panels import (
    PANEL_REGISTRY,
    PREDECLARED_SPLIT_RULE_ID,
    AssayLabelSource,
    CompatibilityStatus,
    DuplexPairingStatus,
    PanelArchitecture,
    PanelColumnMapping,
    PanelDescriptor,
    derive_observation,
    describe_panel,
    predeclared_split,
)

DATA_DIR = Path(__file__).parent / "data"
BENCHMARK_DATA_DIR = DATA_DIR / "benchmark"

# Verbatim from tests/unit/data/README.md's "Predeclared splits" table (issues #97, #102). Copied,
# not re-derived, so this test catches a drift in *either* file against the other's authority.
_README_DEVELOPMENT_ACCESSIONS = [
    "NM_002559",
    "NM_003342",
    "NM_003344",
    "NM_004223",
    "NM_004359",
    "NM_005339",
    "NM_007019",
    "NM_012864",
    "NM_014501",
    "NM_016021",
    "NM_017346",
    "NM_020548",
    "NM_021988",
    "NM_022005",
    "NM_025237",
    "NM_031313",
    "NM_053656",
    "U47298",
    "U92436",
    "XM_214061",
    "XM_371822",
]
_README_HELD_OUT_ACCESSIONS = [
    "AH001498",
    "BD135193",
    "J03132",
    "M60857",
    "NM_001001481",
    "NM_001632",
    "NM_002046",
    "NM_003337",
    "NM_003340",
    "NM_003345",
    "NM_003347",
    "NM_003348",
    "NM_003969",
    "NM_005450",
    "NM_006357",
    "NM_014176",
    "NM_015213",
    "NM_016406",
    "U47296",
    "X75932",
]


def _read_rows(csv_path: Path) -> list[dict[str, str]]:
    """Raw rows from a fixture CSV, as ``derive_observation`` expects to receive them."""
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


# --------------------------------------------------------------------------------------------------
# predeclared_split
# --------------------------------------------------------------------------------------------------


def test_predeclared_split_rule_id_is_declared() -> None:
    """The rule id ``panels.py`` exports must match the one this test file pins against."""
    assert PREDECLARED_SPLIT_RULE_ID == "sha256_accession_parity_v1"


def test_predeclared_split_reproduces_development_accessions() -> None:
    """Every accession the README calls ``development`` must resolve the same way here."""
    for accession in _README_DEVELOPMENT_ACCESSIONS:
        assert predeclared_split(accession) == "development", accession


def test_predeclared_split_reproduces_held_out_accessions() -> None:
    """Every accession the README calls ``held_out`` must resolve the same way here."""
    for accession in _README_HELD_OUT_ACCESSIONS:
        assert predeclared_split(accession) == "held_out", accession


def test_predeclared_split_covers_every_readme_accession_exactly_once() -> None:
    """No accession moves sides, and no accession is claimed by both lists."""
    development = {a for a in _README_DEVELOPMENT_ACCESSIONS if predeclared_split(a) == "development"}
    held_out = {a for a in _README_HELD_OUT_ACCESSIONS if predeclared_split(a) == "held_out"}
    assert development == set(_README_DEVELOPMENT_ACCESSIONS)
    assert held_out == set(_README_HELD_OUT_ACCESSIONS)
    assert set(_README_DEVELOPMENT_ACCESSIONS).isdisjoint(_README_HELD_OUT_ACCESSIONS)


def test_predeclared_split_is_deterministic() -> None:
    """Computable from the accession alone -- calling it twice must not change the answer."""
    assert predeclared_split("NM_003969") == predeclared_split("NM_003969")


# --------------------------------------------------------------------------------------------------
# Registry contents and describe_panel
# --------------------------------------------------------------------------------------------------


def test_registry_has_exactly_the_six_named_panels() -> None:
    """Huesken (subset + full) plus the four issue-named panels; no extra, no missing."""
    assert set(PANEL_REGISTRY) == {
        "huesken_subset",
        "huesken_full",
        "ichihara",
        "martinelli",
        "shmushkovich",
        "oligogym",
    }


def test_only_huesken_subset_has_vendored_bytes() -> None:
    """The repository's second verified fact, pinned so it cannot silently drift.

    ``huesken_subset`` is the only descriptor with real bytes under ``tests/unit/data/``; every
    other panel the issue names is declared but not vendored.
    """
    present = {panel_id for panel_id, d in PANEL_REGISTRY.items() if d.data_present}
    assert present == {"huesken_subset"}


def test_describe_panel_shmushkovich_data_present_is_false() -> None:
    """Must-go-red bullet: named, declared, and honestly unvendored."""
    assert describe_panel("shmushkovich").data_present is False


def test_describe_panel_returns_the_registered_descriptor() -> None:
    """A lookup returns the same object the registry holds, not a copy."""
    assert describe_panel("huesken_subset") is PANEL_REGISTRY["huesken_subset"]


def test_describe_panel_unknown_id_names_the_valid_ids() -> None:
    """A typo gets a list of valid ids, not a bare ``KeyError``."""
    with pytest.raises(ValueError, match="huesken_subset"):
        describe_panel("not_a_real_panel")


@pytest.mark.parametrize("panel_id", sorted(PANEL_REGISTRY))
def test_every_descriptor_has_a_stable_content_hash(panel_id: str) -> None:
    """Same descriptor, called twice, must hash the same and use the declared prefix."""
    descriptor = describe_panel(panel_id)
    first = descriptor.content_hash()
    second = describe_panel(panel_id).content_hash()
    assert first == second
    assert first.startswith("sha256:")


def test_descriptor_content_hash_changes_when_the_descriptor_changes() -> None:
    """An edited descriptor must not report the same identity as its predecessor."""
    base = describe_panel("ichihara")
    edited = base.model_copy(update={"display_name": "a different display name"})
    assert base.content_hash() != edited.content_hash()


def test_asymmetric_descriptor_declares_no_paired_length() -> None:
    """An asymmetric panel has no fixed paired length to declare."""
    assert describe_panel("shmushkovich").declared_paired_length is None


@pytest.mark.parametrize("panel_id", ["huesken_subset", "huesken_full", "ichihara", "martinelli", "oligogym"])
def test_non_asymmetric_descriptors_declare_a_paired_length(panel_id: str) -> None:
    """Every non-asymmetric panel must declare the core its architecture implies."""
    assert describe_panel(panel_id).declared_paired_length is not None


def test_panel_descriptor_rejects_asymmetric_with_a_declared_paired_length() -> None:
    """An asymmetric panel declaring a paired length is a contradiction the model must catch."""
    with pytest.raises(ValidationError):
        PanelDescriptor(
            panel_id="bad_asymmetric",
            display_name="bad",
            citation="none",
            architecture=PanelArchitecture.ASYMMETRIC,
            declared_paired_length=19,
            columns=PanelColumnMapping(guide_column="guide_sequence"),
            assay_label=AssayLabelSource(constant="x"),
            measured_endpoint="x",
            data_present=False,
        )


def test_panel_descriptor_rejects_paired_core_without_a_declared_paired_length() -> None:
    """A non-asymmetric panel with no declared paired length is under-specified."""
    with pytest.raises(ValidationError):
        PanelDescriptor(
            panel_id="bad_paired_core",
            display_name="bad",
            citation="none",
            architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
            columns=PanelColumnMapping(guide_column="guide_sequence"),
            assay_label=AssayLabelSource(constant="x"),
            measured_endpoint="x",
            data_present=False,
        )


def test_assay_label_source_requires_exactly_one_of_column_or_constant() -> None:
    """Neither both nor neither: exactly one source must win."""
    with pytest.raises(ValidationError):
        AssayLabelSource()
    with pytest.raises(ValidationError):
        AssayLabelSource(column="a", constant="b")


# --------------------------------------------------------------------------------------------------
# derive_observation: paired_core_with_overhang (real vendored Huesken row + synthetic fixture)
# --------------------------------------------------------------------------------------------------


def test_huesken_subset_row_is_compatible_and_slices_a_19nt_core() -> None:
    """The real vendored Huesken row, not just a synthetic stand-in, must slice correctly too."""
    descriptor = describe_panel("huesken_subset")
    rows = _read_rows(DATA_DIR / "sirna_efficacy_subset.csv")
    row = rows[0]
    assert len(row["guide_sequence"]) == 21

    observed = derive_observation(descriptor, row, requested_paired_length=19)

    assert observed.compatibility_status is CompatibilityStatus.COMPATIBLE
    assert observed.full_guide_sequence == row["guide_sequence"].upper()
    assert len(observed.full_guide_sequence) == 21
    assert observed.paired_guide_sequence == observed.full_guide_sequence[:19]
    assert observed.guide_3p_overhang == observed.full_guide_sequence[19:21]
    assert len(observed.guide_3p_overhang) == 2
    assert observed.paired_length == 19
    assert observed.paired_slice_start_1based == 1
    assert observed.duplex_pairing_status is DuplexPairingStatus.PAIRED_CORE_WITH_OVERHANG
    assert observed.split in {"development", "held_out"}
    assert observed.split == predeclared_split(row["accession"])


def test_synthetic_paired_core_record_matches_the_must_pass_contract() -> None:
    """The slice's must-go-red assertion: full 21 nt, a 19 nt core, and the measured 2 nt tail."""
    descriptor = PanelDescriptor(
        panel_id="synthetic_paired_core",
        display_name="synthetic paired-core fixture",
        citation="synthetic; tests/unit/data/benchmark/README.md",
        architecture=PanelArchitecture.PAIRED_CORE_WITH_OVERHANG,
        declared_paired_length=19,
        columns=PanelColumnMapping(guide_column="guide_sequence", accession_column="accession"),
        assay_label=AssayLabelSource(constant="synthetic_paired_core"),
        measured_endpoint="inhibition_fraction",
        data_present=False,
    )
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_paired_core_with_overhang.csv")
    row = next(r for r in rows if r["accession"] == "SYN0001")

    observed = derive_observation(descriptor, row, requested_paired_length=19)

    assert observed.compatibility_status is CompatibilityStatus.COMPATIBLE
    assert len(observed.full_guide_sequence) == 21
    assert observed.paired_guide_sequence == observed.full_guide_sequence[:19]
    assert observed.guide_3p_overhang == observed.full_guide_sequence[19:21]
    assert observed.guide_3p_overhang == "TT"


def test_paired_core_guide_too_short_for_the_requested_core_is_incompatible() -> None:
    """A guide shorter than the requested core cannot supply one; never padded to fake it."""
    descriptor = describe_panel("ichihara")
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_paired_core_with_overhang.csv")
    row = next(r for r in rows if r["accession"] == "SYN0005")
    assert len(row["guide_sequence"]) == 18

    observed = derive_observation(descriptor, row, requested_paired_length=19)

    assert observed.compatibility_status is CompatibilityStatus.INCOMPATIBLE
    assert "18" in observed.compatibility_reason
    assert "19" in observed.compatibility_reason
    # Never coerced: the full sequence travels through unsliced.
    assert observed.paired_guide_sequence == observed.full_guide_sequence
    assert observed.paired_length == 18
    assert observed.guide_3p_overhang is None


# --------------------------------------------------------------------------------------------------
# derive_observation: fully_complementary
# --------------------------------------------------------------------------------------------------


def test_fully_complementary_record_is_compatible_and_blunt() -> None:
    """Equal-length guide and passenger: compatible, and the overhang is blunt, not absent."""
    descriptor = describe_panel("martinelli")
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_fully_complementary.csv")
    row = next(r for r in rows if r["accession"] == "FC0001")

    observed = derive_observation(descriptor, row, requested_paired_length=21)

    assert observed.compatibility_status is CompatibilityStatus.COMPATIBLE
    assert observed.paired_guide_sequence == observed.full_guide_sequence
    assert observed.guide_3p_overhang == ""
    assert observed.passenger_3p_overhang == ""
    assert observed.duplex_pairing_status is DuplexPairingStatus.FULLY_COMPLEMENTARY


def test_fully_complementary_length_mismatch_between_strands_is_incompatible() -> None:
    """A guide/passenger length mismatch contradicts "fully complementary"; both lengths are named."""
    descriptor = describe_panel("martinelli")
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_fully_complementary.csv")
    row = next(r for r in rows if r["accession"] == "FC0004")
    assert len(row["guide_sequence"]) == 21
    assert len(row["passenger_sequence"]) == 20

    observed = derive_observation(descriptor, row, requested_paired_length=21)

    assert observed.compatibility_status is CompatibilityStatus.INCOMPATIBLE
    assert "21" in observed.compatibility_reason
    assert "20" in observed.compatibility_reason
    assert observed.paired_guide_sequence == observed.full_guide_sequence


def test_fully_complementary_requested_length_shorter_than_measured_is_incompatible() -> None:
    """A fully complementary duplex cannot be sliced without inventing an overhang."""
    descriptor = PanelDescriptor(
        panel_id="synthetic_fully_complementary_guide_only",
        display_name="synthetic fully-complementary fixture, guide only",
        citation="synthetic; tests/unit/data/benchmark/README.md",
        architecture=PanelArchitecture.FULLY_COMPLEMENTARY,
        declared_paired_length=21,
        columns=PanelColumnMapping(guide_column="guide_sequence", accession_column="accession"),
        assay_label=AssayLabelSource(constant="synthetic_fully_complementary"),
        measured_endpoint="inhibition_fraction",
        data_present=False,
    )
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_fully_complementary.csv")
    row = next(r for r in rows if r["accession"] == "FC0001")

    observed = derive_observation(descriptor, row, requested_paired_length=19)

    assert observed.compatibility_status is CompatibilityStatus.INCOMPATIBLE
    assert observed.paired_guide_sequence == observed.full_guide_sequence
    assert len(observed.full_guide_sequence) == 21


# --------------------------------------------------------------------------------------------------
# derive_observation: asymmetric (the 15/20 must-go-red case)
# --------------------------------------------------------------------------------------------------


def test_asymmetric_15_20_record_is_incompatible_and_never_coerced() -> None:
    """Slice's core must-go-red assertion: incompatible, reason names the asymmetry, no coercion."""
    descriptor = describe_panel("shmushkovich")
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_asymmetric_15_20.csv")
    row = rows[0]
    assert len(row["guide_sequence"]) == 20
    assert len(row["passenger_sequence"]) == 15

    observed = derive_observation(descriptor, row, requested_paired_length=19)

    assert observed.compatibility_status is CompatibilityStatus.INCOMPATIBLE
    assert "asymmetric" in observed.compatibility_reason.lower()
    assert "20" in observed.compatibility_reason
    assert "15" in observed.compatibility_reason
    # Never coerced to a symmetric (19-23 nt) length: the guide travels through at its own,
    # unsliced, 20 nt length rather than being truncated or padded to the requested core.
    assert observed.paired_guide_sequence == observed.full_guide_sequence
    assert len(observed.paired_guide_sequence) == 20
    assert observed.paired_length == 20
    assert observed.guide_3p_overhang is None
    assert observed.passenger_3p_overhang is None
    assert observed.duplex_pairing_status is DuplexPairingStatus.ASYMMETRIC


@pytest.mark.parametrize("requested_paired_length", [15, 19, 20, 23])
def test_asymmetric_is_incompatible_regardless_of_requested_length(requested_paired_length: int) -> None:
    """#109 excludes asymmetric design entirely: no requested length can make this compatible."""
    descriptor = describe_panel("shmushkovich")
    rows = _read_rows(BENCHMARK_DATA_DIR / "synthetic_asymmetric_15_20.csv")
    row = rows[0]

    observed = derive_observation(descriptor, row, requested_paired_length=requested_paired_length)

    assert observed.compatibility_status is CompatibilityStatus.INCOMPATIBLE


# --------------------------------------------------------------------------------------------------
# Other derived fields: assay_label, measured_value, sequence validation
# --------------------------------------------------------------------------------------------------


def test_assay_label_falls_back_to_the_descriptor_constant() -> None:
    """Huesken's rows carry no per-row assay label, so the descriptor's constant applies."""
    descriptor = describe_panel("huesken_subset")
    rows = _read_rows(DATA_DIR / "sirna_efficacy_subset.csv")
    observed = derive_observation(descriptor, rows[0], requested_paired_length=19)
    assert observed.assay_label == "huesken_2005_knockdown_inhibition"


def test_measured_value_is_read_verbatim_never_rescaled() -> None:
    """The efficacy column round-trips exactly, no normalisation or sign flip."""
    descriptor = describe_panel("huesken_subset")
    rows = _read_rows(DATA_DIR / "sirna_efficacy_subset.csv")
    row = rows[0]
    observed = derive_observation(descriptor, row, requested_paired_length=19)
    assert observed.measured_value == float(row["efficacy"])


def test_measured_value_is_none_when_the_panel_has_no_value_column() -> None:
    """A panel with no numeric endpoint must report ``None``, not ``0.0`` or an empty string."""
    descriptor = describe_panel("huesken_subset").model_copy(
        update={"columns": describe_panel("huesken_subset").columns.model_copy(update={"measured_value_column": None})}
    )
    rows = _read_rows(DATA_DIR / "sirna_efficacy_subset.csv")
    observed = derive_observation(descriptor, rows[0], requested_paired_length=19)
    assert observed.measured_value is None


def test_full_guide_sequence_is_upper_cased_but_otherwise_verbatim() -> None:
    """Case-normalised for comparison, but no base is added, removed or reordered."""
    descriptor = describe_panel("huesken_subset")
    row = {"guide_sequence": "acgtacgtacgtacgtacgtt", "accession": "NM_003969", "efficacy": "0.5"}
    observed = derive_observation(descriptor, row, requested_paired_length=19)
    assert observed.full_guide_sequence == "ACGTACGTACGTACGTACGTT"


def test_non_nucleotide_characters_raise() -> None:
    """A malformed source row must fail loudly, not silently produce a bad observation."""
    descriptor = describe_panel("huesken_subset")
    row = {"guide_sequence": "ACGTACGTACGTACGTACGXX", "accession": "NM_003969", "efficacy": "0.5"}
    with pytest.raises(ValueError, match="non-nucleotide"):
        derive_observation(descriptor, row, requested_paired_length=19)


# --------------------------------------------------------------------------------------------------
# Fixtures: readable, well-formed, and honestly not real panel bytes
# --------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "filename",
    [
        "synthetic_paired_core_with_overhang.csv",
        "synthetic_fully_complementary.csv",
        "synthetic_asymmetric_15_20.csv",
    ],
)
def test_synthetic_fixtures_are_named_generically_not_after_a_real_panel(filename: str) -> None:
    """A synthetic fixture's filename must never be mistaken for a real panel's."""
    for real_panel_name in ("huesken", "ichihara", "martinelli", "shmushkovich", "oligogym"):
        assert real_panel_name not in filename


def test_synthetic_context_fasta_is_well_formed() -> None:
    """Two-plus records, each a valid nucleotide sequence."""
    text = (BENCHMARK_DATA_DIR / "synthetic_context.fa").read_text()
    records = [line for line in text.splitlines() if line]
    headers = [line for line in records if line.startswith(">")]
    assert len(headers) >= 2
    for line in records:
        if not line.startswith(">"):
            assert set(line.upper()) <= set("ACGTN")


def test_benchmark_fixtures_readme_exists_and_disclaims_the_missing_panels() -> None:
    """The fixtures README must name, and disclaim, every panel with no vendored bytes."""
    text = (BENCHMARK_DATA_DIR / "README.md").read_text()
    for missing_panel in ("Ichihara", "Martinelli", "Shmushkovich", "OligoGym"):
        assert missing_panel in text
    assert "not vendored" in text.lower() or "no bytes" in text.lower()
