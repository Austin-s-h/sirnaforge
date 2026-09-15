"""Benchmark panel registry: architectures, compatibility verdicts, honestly-labelled fixtures.

Issue #109. Two authorities this file pins ``panels.py`` against, and never restates from memory:

* ``tests/data/benchmarks/README.md`` and ``tests/data/benchmarks/oligogym/manifest.json`` for the
  vendored OligoGym-derived bytes (row counts, the guide convention, the synthetic-flank rule). Every
  count and geometry the registry's comments claim is re-derived here from ``records.csv`` itself, so
  a regenerated table cannot leave a stale number standing in a docstring.
* ``tests/unit/data/README.md`` for the predeclared-split accession lists (issues #97, #102), and
  ``tests/unit/data/benchmark/README.md`` for the synthetic fixtures' provenance and checksums.

The synthetic fixtures are exercised through descriptors declared *in this file* rather than through
registry entries: since ``f4beab7`` the registry's Ichihara/Martinelli/Shmushkovich descriptors are
mapped onto real vendored bytes, so borrowing one to read a four-row synthetic CSV would test a
contract against a panel whose real rows are right there to test it against instead.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pytest
from Bio.Seq import Seq
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
    PanelRowSelector,
    TargetIdentityStatus,
    derive_observation,
    describe_panel,
    predeclared_split,
)

DATA_DIR = Path(__file__).parent / "data"
BENCHMARK_DATA_DIR = DATA_DIR / "benchmark"
#: ``tests/unit/`` -> ``tests/`` -> repo root, the same walk ``prepare.py`` does from ``src/``.
REPO_ROOT = Path(__file__).resolve().parents[2]
OLIGOGYM_DIR = REPO_ROOT / "tests" / "data" / "benchmarks" / "oligogym"

#: Panels whose bytes are the vendored OligoGym-derived table, and the id ``oligogym`` aggregates.
_OLIGOGYM_MEMBER_PANELS = ("ichihara", "martinelli", "shmushkovich")


@lru_cache(maxsize=1)
def _records_rows() -> tuple[dict[str, str], ...]:
    """Every row of the vendored ``oligogym/records.csv``, read once for the whole module."""
    with (OLIGOGYM_DIR / "records.csv").open(newline="") as handle:
        return tuple(dict(row) for row in csv.DictReader(handle))


def _rows_for(panel_id: str) -> list[dict[str, str]]:
    """The vendored rows one registry panel selects, via its own selector."""
    descriptor = describe_panel(panel_id)
    return [row for row in _records_rows() if descriptor.selects_row(row)]


def _synthetic_paired_core_descriptor() -> PanelDescriptor:
    """A paired-core-with-overhang descriptor for the synthetic fixtures, declared locally.

    Deliberately not ``describe_panel("ichihara")``: that descriptor now maps records.csv's
    ``context_type``/``site_start``/``site_end`` columns, which a synthetic fixture does not have.
    """
    return PanelDescriptor(
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


def _synthetic_asymmetric_descriptor() -> PanelDescriptor:
    """An asymmetric descriptor for the two-row synthetic fixture, declared locally.

    The real 356-row Shmushkovich panel is exercised against its vendored bytes in
    ``test_shmushkovich_...``; this keeps the fixture's boundary cases readable without borrowing a
    descriptor that now expects records.csv's columns.
    """
    return PanelDescriptor(
        panel_id="synthetic_asymmetric",
        display_name="synthetic asymmetric 15/20 fixture",
        citation="synthetic; tests/unit/data/benchmark/README.md",
        architecture=PanelArchitecture.ASYMMETRIC,
        declared_paired_length=None,
        columns=PanelColumnMapping(
            guide_column="guide_sequence",
            passenger_column="passenger_sequence",
            accession_column="accession",
        ),
        assay_label=AssayLabelSource(constant="synthetic_asymmetric"),
        measured_endpoint="inhibition_fraction",
        data_present=False,
    )


def _synthetic_fully_complementary_descriptor(*, with_passenger: bool = True) -> PanelDescriptor:
    """The fully-complementary architecture's exemplar, now that no vendored panel measures one.

    All 907 vendored Martinelli rows are 21 nt strands whose 5' 19 nt pair (see
    ``test_martinelli_declared_geometry_matches_the_vendored_bytes``), so
    :attr:`PanelArchitecture.FULLY_COMPLEMENTARY` has no registry panel left to exercise it. It stays
    part of the vocabulary -- a blunt duplex is a real architecture a future panel will declare -- so
    its slicing and mismatch branches are pinned here against the synthetic fixture instead.
    """
    return PanelDescriptor(
        panel_id="synthetic_fully_complementary",
        display_name="synthetic fully-complementary fixture",
        citation="synthetic; tests/unit/data/benchmark/README.md",
        architecture=PanelArchitecture.FULLY_COMPLEMENTARY,
        declared_paired_length=21,
        columns=PanelColumnMapping(
            guide_column="guide_sequence",
            passenger_column="passenger_sequence" if with_passenger else None,
            accession_column="accession",
        ),
        assay_label=AssayLabelSource(constant="synthetic_fully_complementary"),
        measured_endpoint="inhibition_fraction",
        data_present=False,
    )


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


def test_registry_has_the_six_named_panels_plus_one_id_per_architecture() -> None:
    """Huesken (subset + full), the four issue-named panels, and the three ``user_supplied_*`` ids.

    The six named panels are #109's own list. The three added ids are architecture-level rather than
    panel-level -- they name a duplex geometry, ship no bytes, and read the table a caller passes with
    ``--panel-csv`` -- so there must be exactly one per :class:`PanelArchitecture` member: a caller
    holding an asymmetric or a blunt fully-complementary table had no reachable descriptor at all once
    the OligoGym bytes landed and every published panel began refusing ``--panel-csv``.
    """
    assert set(PANEL_REGISTRY) == {
        "huesken_subset",
        "huesken_full",
        "ichihara",
        "martinelli",
        "shmushkovich",
        "oligogym",
        "user_supplied_paired_core_with_overhang",
        "user_supplied_fully_complementary",
        "user_supplied_asymmetric",
    }
    user_supplied = {panel_id: d for panel_id, d in PANEL_REGISTRY.items() if panel_id.startswith("user_supplied_")}
    assert {d.architecture for d in user_supplied.values()} == set(PanelArchitecture)
    assert len(user_supplied) == len(PanelArchitecture)
    # The id names the geometry it reads a table under, so a caller cannot pick one by accident.
    assert all(panel_id == f"user_supplied_{d.architecture.value}" for panel_id, d in user_supplied.items())


def test_vendored_panels_are_huesken_subset_plus_the_four_oligogym_ids() -> None:
    """``data_present`` must name exactly the panels this repository ships bytes for.

    ``f4beab7`` vendored ``tests/data/benchmarks/oligogym/records.csv``, so the four OligoGym-derived
    ids joined ``huesken_subset``; ``huesken_full``'s table (``work/sirna_bench.csv``) is untracked
    here and stays ``False``, and the three ``user_supplied_*`` ids vendor nothing by construction --
    that is what makes ``--panel-csv`` reachable at all, since a vendored panel refuses one. Both
    halves are pinned, so this goes red the day a panel's bytes are added or removed without its
    descriptor being told, in either direction.
    """
    present = {panel_id for panel_id, d in PANEL_REGISTRY.items() if d.data_present}
    absent = set(PANEL_REGISTRY) - present
    assert present == {"huesken_subset", "ichihara", "martinelli", "shmushkovich", "oligogym"}
    assert absent == {
        "huesken_full",
        "user_supplied_paired_core_with_overhang",
        "user_supplied_fully_complementary",
        "user_supplied_asymmetric",
    }


@pytest.mark.parametrize("panel_id", sorted(PANEL_REGISTRY))
def test_data_present_panels_name_bytes_that_exist(panel_id: str) -> None:
    """A ``data_present`` claim must resolve to a real file, and its absence must name none."""
    descriptor = describe_panel(panel_id)
    if descriptor.data_present:
        assert descriptor.vendored_csv is not None
        assert (REPO_ROOT / descriptor.vendored_csv).is_file(), descriptor.vendored_csv
    else:
        assert descriptor.vendored_csv is None


def test_no_descriptor_still_carries_a_placeholder_citation() -> None:
    """The four vendored panels' citations must state what the repository can show, not a placeholder.

    #109 shipped ``martinelli``/``shmushkovich``/``oligogym`` with citations reading "no primary
    citation has been independently verified in this repository" and endpoints named
    ``efficacy_unspecified_placeholder``, which was true before the bytes landed and is false now.
    Every panel that reads a DOI out of its own rows must record it (see
    ``test_every_vendored_panel_records_the_doi_its_own_rows_carry``).
    """
    for panel_id, descriptor in PANEL_REGISTRY.items():
        combined = f"{descriptor.citation} {descriptor.measured_endpoint} {descriptor.assay_label.constant or ''}"
        assert "placeholder" not in combined.lower(), panel_id


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


@pytest.mark.parametrize("panel_id", ["shmushkovich", "user_supplied_asymmetric"])
def test_asymmetric_descriptor_declares_no_paired_length(panel_id: str) -> None:
    """An asymmetric panel has no fixed paired length to declare."""
    assert describe_panel(panel_id).declared_paired_length is None


@pytest.mark.parametrize(
    "panel_id",
    [
        "huesken_subset",
        "huesken_full",
        "ichihara",
        "martinelli",
        "oligogym",
        "user_supplied_paired_core_with_overhang",
        "user_supplied_fully_complementary",
    ],
)
def test_non_asymmetric_descriptors_declare_a_paired_length(panel_id: str) -> None:
    """Every non-asymmetric panel must declare the core its architecture implies."""
    assert describe_panel(panel_id).declared_paired_length is not None


# --------------------------------------------------------------------------------------------------
# The user-supplied architecture-level ids: one column set, no claim about anyone's rows
# --------------------------------------------------------------------------------------------------

#: The ``user_supplied_*`` ids' column set as ``docs/benchmark_artifacts.md`` publishes it: the one
#: mandatory column, then the three optional ones. Spelled here independently of
#: ``panels.USER_TABLE_COLUMNS`` so the two must agree.
_DOCUMENTED_USER_TABLE_COLUMNS = ("guide_sequence", "passenger_sequence", "accession", "efficacy")

_USER_SUPPLIED_IDS = (
    "user_supplied_paired_core_with_overhang",
    "user_supplied_fully_complementary",
    "user_supplied_asymmetric",
)


def test_the_user_supplied_ids_share_one_column_set_with_one_mandatory_column() -> None:
    """All three read the same column names, and only ``guide_sequence`` is fatal if absent.

    One fixed set is the whole mechanism -- there is no per-run column-mapping option, because a
    mapping supplied at run time would not be covered by the ``descriptor_hash`` the manifest records,
    so two artifacts could share a panel id and a hash while having read different columns. The caller
    renames their columns instead, and this pins the names ``docs/benchmark_artifacts.md`` tells them
    to rename to.
    """
    mappings = {describe_panel(panel_id).columns for panel_id in _USER_SUPPLIED_IDS}
    assert len(mappings) == 1  # frozen models, so equal mappings collapse to one entry
    columns = mappings.pop()
    assert (
        columns.guide_column,
        columns.passenger_column,
        columns.accession_column,
        columns.measured_value_column,
    ) == _DOCUMENTED_USER_TABLE_COLUMNS
    for panel_id in _USER_SUPPLIED_IDS:
        assert describe_panel(panel_id).required_source_columns() == ("guide_sequence",)


@pytest.mark.parametrize(
    "filename",
    [
        "synthetic_paired_core_with_overhang.csv",
        "synthetic_fully_complementary.csv",
        "synthetic_asymmetric_15_20.csv",
    ],
)
def test_every_synthetic_fixture_uses_the_documented_user_table_column_names(filename: str) -> None:
    """The fixtures are prepared under these ids, so their headers must be a subset of the declared set.

    A fixture carrying a column the descriptors do not map would be silently ignored, and one missing
    ``guide_sequence`` would be refused by ``prepare``'s header check -- either way the fixture, not
    the required set, is what has to move.
    """
    header = tuple((BENCHMARK_DATA_DIR / filename).read_text().splitlines()[0].split(","))
    assert set(header) <= set(_DOCUMENTED_USER_TABLE_COLUMNS), header
    assert "guide_sequence" in header


def test_the_user_supplied_ids_vouch_for_nothing_and_locate_nothing() -> None:
    """No bytes, no split, and no target coordinate: the honest content of "the caller's own table".

    ``target_identity_status`` is ``UNAVAILABLE`` because this repository knows nothing about a
    caller's target context -- not even whether one exists -- and ``split_rule_id`` is ``None`` because
    ``predeclared_split``'s id is pinned to the accession lists in ``tests/unit/data/README.md``, so
    stamping it on a caller's accessions would report an un-audited partition under an audited rule's
    name.
    """
    for panel_id in _USER_SUPPLIED_IDS:
        descriptor = describe_panel(panel_id)
        assert descriptor.data_present is False
        assert descriptor.vendored_csv is None
        assert descriptor.row_selector is None
        assert descriptor.split_rule_id is None
        assert descriptor.target_identity_status is TargetIdentityStatus.UNAVAILABLE
        assert descriptor.aggregate_of is None
        assert "vouches for nothing" in descriptor.citation
        assert "docs/benchmark_artifacts.md" in (descriptor.redistribution or "")


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
    descriptor = _synthetic_paired_core_descriptor()
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
    descriptor = _synthetic_paired_core_descriptor()
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
    descriptor = _synthetic_fully_complementary_descriptor()
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
    descriptor = _synthetic_fully_complementary_descriptor()
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
    descriptor = _synthetic_fully_complementary_descriptor(with_passenger=False)
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
    descriptor = _synthetic_asymmetric_descriptor()
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
    descriptor = _synthetic_asymmetric_descriptor()
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


def test_benchmark_fixtures_readme_names_the_panels_and_labels_its_own_fixtures_synthetic() -> None:
    """The fixtures README must name every panel it stands in for, and admit that it stands in.

    Its "Ichihara, Martinelli, Shmushkovich and OligoGym are not vendored anywhere in this
    repository" claim was true when #109 was settled and was falsified by ``f4beab7``; rewriting that
    prose belongs to the docs pass, so this assertion no longer depends on the sentence. What it still
    pins is the part that stays true: the files in that directory are synthetic stand-ins, and the
    authority for the vendored bytes is ``tests/data/benchmarks/README.md``.
    """
    text = (BENCHMARK_DATA_DIR / "README.md").read_text()
    for panel_name in ("Ichihara", "Martinelli", "Shmushkovich", "OligoGym"):
        assert panel_name in text
    assert "synthetic" in text.lower()
    assert (REPO_ROOT / "tests" / "data" / "benchmarks" / "README.md").is_file()


# --------------------------------------------------------------------------------------------------
# The vendored OligoGym-derived bytes (f4beab7). Every number the registry's comments claim is
# re-derived here from tests/data/benchmarks/oligogym/ rather than restated, so a regenerated
# records.csv cannot leave a stale count or geometry standing in panels.py.
# --------------------------------------------------------------------------------------------------


def test_vendored_records_csv_matches_its_own_manifest() -> None:
    """``records.csv`` and ``manifest.json`` must agree on how many rows each dataset has.

    ``manifest.json`` is what ``tests/data/benchmarks/README.md`` points a reader at, so it is the
    number a descriptor's comment is allowed to quote; this is the check that the number is real.
    """
    manifest = json.loads((OLIGOGYM_DIR / "manifest.json").read_text())
    rows = _records_rows()

    counted = Counter(row["dataset"] for row in rows)
    assert dict(counted) == manifest["datasets"]
    assert len(rows) == manifest["total_records"] == 4113
    assert manifest["site_start_1based"] == 71
    assert manifest["flank_nt"] == 70


def test_row_selectors_partition_the_vendored_table_with_no_row_claimed_twice() -> None:
    """The three member panels must tile ``records.csv`` exactly: 2,850 + 907 + 356 = 4,113.

    One file holds four datasets, so a selector that overreached would silently ingest a sibling
    panel's rows under this panel's architecture -- and one that under-reached would drop measured
    observations #109 promises never to drop. Both failures are visible only as counts.
    """
    selected = {panel_id: _rows_for(panel_id) for panel_id in _OLIGOGYM_MEMBER_PANELS}
    assert [len(selected[panel_id]) for panel_id in _OLIGOGYM_MEMBER_PANELS] == [2850, 907, 356]

    ids_by_panel = {panel_id: {row["benchmark_id"] for row in rows} for panel_id, rows in selected.items()}
    assert sum(len(ids) for ids in ids_by_panel.values()) == len(set().union(*ids_by_panel.values())) == 4113
    # The aggregate selects everything, which is exactly why it is not ingestible (see below).
    assert len(_rows_for("oligogym")) == 4113


@pytest.mark.parametrize("panel_id", ["ichihara", "martinelli", "shmushkovich", "oligogym"])
def test_every_vendored_panel_records_the_doi_its_own_rows_carry(panel_id: str) -> None:
    """A citation must be readable out of the bytes it describes, not asserted alongside them.

    ``records.csv``'s ``source_url`` is the only citation this repository holds for Martinelli and
    Shmushkovich, so the descriptor has to carry that exact DOI -- this is the assertion that would
    have gone red against #109's "no primary citation has been independently verified" placeholders.
    """
    descriptor = describe_panel(panel_id)
    dois = {row["source_url"] for row in _rows_for(panel_id)}
    assert dois, panel_id
    for doi in dois:
        assert doi in descriptor.citation, f"{panel_id} does not record {doi}"


def test_vendored_panels_disclose_the_synthetic_flanks_and_the_absent_upstream_extract() -> None:
    """Redistribution text must say the context is fabricated and the upstream extract is missing.

    Both are things a reader of ``observations.csv`` cannot see for themselves: the rows look like
    ordinary target sites, and ``source_file`` names ``tests/data/external/oligogym/*.csv.gz``, which
    is not in this tree at all.
    """
    for panel_id in (*_OLIGOGYM_MEMBER_PANELS, "oligogym"):
        redistribution = describe_panel(panel_id).redistribution or ""
        assert "synthetic" in redistribution, panel_id
        assert "not vendored" in redistribution.lower(), panel_id
        assert "oligogym_native_design" in redistribution, panel_id
    assert not (REPO_ROOT / "tests" / "data" / "external").exists()


def test_ichihara_declared_geometry_matches_all_2850_vendored_rows() -> None:
    """The registry claims a 19 nt paired core plus a measured 2 nt 3' guide overhang. Check every row.

    ``passenger == revcomp(guide[:19])`` is what makes ``PAIRED_CORE_WITH_OVERHANG`` at 19 the
    measured geometry rather than a convention borrowed from Huesken: the passenger the panel recorded
    pairs the first 19 nt of the guide and nothing beyond it.
    """
    descriptor = describe_panel("ichihara")
    assert descriptor.architecture is PanelArchitecture.PAIRED_CORE_WITH_OVERHANG
    assert descriptor.declared_paired_length == 19

    rows = _rows_for("ichihara")
    for row in rows:
        guide = row["guide_sequence"]
        assert len(guide) == 21
        assert row["passenger_sequence"] == str(Seq(guide[:19]).reverse_complement())

    observed = [derive_observation(descriptor, row, requested_paired_length=19) for row in rows]
    assert all(o.compatibility_status is CompatibilityStatus.COMPATIBLE for o in observed)
    assert all(o.paired_guide_sequence == o.full_guide_sequence[:19] for o in observed)
    assert all(o.guide_3p_overhang == o.full_guide_sequence[19:21] for o in observed)
    # "" means measured and blunt (artifact.py's convention), which is what a 19 nt passenger against
    # a 19 nt core is; None would claim the panel stated no overhang at all.
    assert all(o.passenger_3p_overhang == "" for o in observed)


def test_martinelli_declared_geometry_matches_the_vendored_bytes_not_the_issue_text() -> None:
    """#109 declared Martinelli fully complementary at 21; its own bytes say 19-nt core plus 2/2 tails.

    The counts pinned here are the evidence the re-declaration rests on. Under the old
    ``FULLY_COMPLEMENTARY``/21 declaration, 889 of the 907 rows would have been written with
    ``guide_3p_overhang=""`` -- "measured and blunt" -- which their bytes contradict.
    """
    descriptor = describe_panel("martinelli")
    assert descriptor.architecture is PanelArchitecture.PAIRED_CORE_WITH_OVERHANG
    assert descriptor.declared_paired_length == 19

    rows = _rows_for("martinelli")
    core_19 = sum(
        1 for row in rows if row["passenger_sequence"][:19] == str(Seq(row["guide_sequence"][:19]).reverse_complement())
    )
    blunt_21 = sum(
        1 for row in rows if row["passenger_sequence"] == str(Seq(row["guide_sequence"]).reverse_complement())
    )
    assert (core_19, blunt_21, len(rows)) == (858, 18, 907)
    assert all(len(row["guide_sequence"]) == len(row["passenger_sequence"]) == 21 for row in rows)

    observed = [derive_observation(descriptor, row, requested_paired_length=19) for row in rows]
    assert all(o.paired_guide_sequence == o.full_guide_sequence[:19] for o in observed)
    # A measured 2 nt tail on each strand, never the empty string a blunt duplex would report.
    assert all(o.guide_3p_overhang == o.full_guide_sequence[19:21] for o in observed)
    assert all(o.passenger_3p_overhang and len(o.passenger_3p_overhang) == 2 for o in observed)


def test_shmushkovich_keeps_all_356_vendored_rows_incompatible_and_unsliced() -> None:
    """#109 keeps asymmetric out of scope, so every one of the 356 hsiRNA rows is refused, not adapted.

    The count is the point: a panel that quietly sliced its 20 nt guide to 19 would still report
    "incompatible" nowhere and would look like a working benchmark. All 356 rows are a 20 nt guide
    against a 15 nt passenger equal to ``revcomp(guide[:15])`` -- a 15 nt core with a 5 nt tail.
    """
    descriptor = describe_panel("shmushkovich")
    assert descriptor.architecture is PanelArchitecture.ASYMMETRIC
    assert descriptor.declared_paired_length is None

    rows = _rows_for("shmushkovich")
    assert len(rows) == 356
    for row in rows:
        assert (len(row["guide_sequence"]), len(row["passenger_sequence"])) == (20, 15)
        assert row["passenger_sequence"] == str(Seq(row["guide_sequence"][:15]).reverse_complement())

    observed = [derive_observation(descriptor, row, requested_paired_length=19) for row in rows]
    assert [o.compatibility_status for o in observed] == [CompatibilityStatus.INCOMPATIBLE] * 356
    assert all("asymmetric" in o.compatibility_reason for o in observed)
    assert all("#110" in o.compatibility_reason for o in observed)
    assert all(o.paired_guide_sequence == o.full_guide_sequence for o in observed)
    assert all(o.paired_length == 20 for o in observed)
    assert all(o.duplex_pairing_status is DuplexPairingStatus.ASYMMETRIC for o in observed)


def test_shmushkovich_measured_value_is_the_measured_direction_not_the_derived_flip() -> None:
    """``label_processed`` (percent remaining), never records.csv's derived ``100 - label`` column.

    The adapter wrote ``efficacy_higher_is_better`` so ranking comparisons could pool datasets; for
    this panel alone that column is a transformation rather than a measurement, and recording it as
    ``measured_value`` under an endpoint named "remaining" would invert the panel's own direction.
    """
    descriptor = describe_panel("shmushkovich")
    assert descriptor.columns.measured_value_column == "label_processed"
    assert "remaining" in descriptor.measured_endpoint
    assert "lower_is_better" in descriptor.measured_endpoint

    for row in _rows_for("shmushkovich"):
        observed = derive_observation(descriptor, row, requested_paired_length=19)
        assert observed.measured_value == float(row["label_processed"])
        assert observed.measured_value == pytest.approx(100.0 - float(row["efficacy_higher_is_better"]))
        assert row["label_direction"] == "lower_is_better"


def test_ichihara_and_martinelli_measured_value_is_the_one_value_all_three_label_columns_hold() -> None:
    """For these two panels ``efficacy_higher_is_better`` is not a transform: it equals the raw label.

    That equality is what makes reading the pooled column honest here and dishonest for Shmushkovich,
    so it is checked rather than assumed -- on all 3,757 rows.
    """
    for panel_id in ("ichihara", "martinelli"):
        descriptor = describe_panel(panel_id)
        assert descriptor.columns.measured_value_column == "efficacy_higher_is_better"
        rows = _rows_for(panel_id)
        for row in rows:
            assert row["label_raw"] == row["label_processed"] == row["efficacy_higher_is_better"]
            assert row["label_direction"] == "higher_is_better"
        observed = derive_observation(descriptor, rows[0], requested_paired_length=19)
        assert observed.measured_value == float(rows[0]["label_raw"])
    # Percent, not Huesken's 0-1 fraction: named differently so no reader pools the two panels.
    assert "percent" in describe_panel("ichihara").measured_endpoint
    assert describe_panel("huesken_subset").measured_endpoint == "inhibition_fraction"


# --------------------------------------------------------------------------------------------------
# Target identity: a real coordinate in a fabricated context
# --------------------------------------------------------------------------------------------------


def test_target_identity_vocabulary_cannot_express_confirmed() -> None:
    """#110 owns native mapping, so #109's vocabulary has no ``confirmed`` member to emit by accident."""
    assert {status.value for status in TargetIdentityStatus} == {
        "unavailable",
        "panel_local",
        "synthetic_context_local",
    }


@pytest.mark.parametrize(
    ("panel_id", "expected_span"),
    [("ichihara", (71, 91)), ("martinelli", (71, 91)), ("shmushkovich", (71, 90))],
)
def test_vendored_sites_are_recorded_as_synthetic_context_local(panel_id: str, expected_span: tuple[int, int]) -> None:
    """Every vendored row's site is exact, reproducible, and in a context this repository fabricated.

    ``site_start`` is 71 on all 4,113 rows because the adapter wrote ``70 x A + revcomp(guide) +
    70 x A``. Recording that as ``panel_local`` would read as the panel's own target coordinate, and
    ``unavailable`` would throw away a coordinate the artifact can honestly reproduce -- so #109 says
    which context it is, and never that the identity is confirmed.
    """
    descriptor = describe_panel(panel_id)
    assert descriptor.target_identity_status is TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL

    for row in _rows_for(panel_id):
        observed = derive_observation(descriptor, row, requested_paired_length=19)
        assert observed.target_identity_status is TargetIdentityStatus.SYNTHETIC_CONTEXT_LOCAL
        assert observed.target_context_type == "synthetic_neutral_flanks"
        assert (observed.target_start_1based, observed.target_end_1based) == expected_span


def test_a_row_claiming_a_native_context_is_refused_rather_than_relabelled() -> None:
    """The descriptor's context claim is checked against the row, not stamped onto it.

    This is the guard that keeps #110's work and #109's apart: if a native-context table were ever
    fed to one of these panels, its real coordinates would be written under a synthetic label (or a
    synthetic 71 under a native one), and nothing downstream could tell the difference afterwards.
    """
    descriptor = describe_panel("ichihara")
    row = dict(_rows_for("ichihara")[0])
    row["context_type"] = "native_transcript"

    with pytest.raises(ValueError, match="synthetic_neutral_flanks"):
        derive_observation(descriptor, row, requested_paired_length=19)


def test_no_oligogym_panel_declares_a_split_rule() -> None:
    """No accession in records.csv means no predeclared split -- not a new rule wearing its id.

    ``predeclared_split`` is pinned to the accession lists in ``tests/unit/data/README.md``; splitting
    these panels on ``target`` (a gene symbol, empty on all 1,263 Martinelli/Shmushkovich rows) would
    be a different rule reported under ``sha256_accession_parity_v1``.
    """
    for panel_id in (*_OLIGOGYM_MEMBER_PANELS, "oligogym"):
        descriptor = describe_panel(panel_id)
        assert descriptor.split_rule_id is None, panel_id
        assert descriptor.columns.accession_column is None, panel_id

    for panel_id in _OLIGOGYM_MEMBER_PANELS:
        descriptor = describe_panel(panel_id)
        observed = derive_observation(descriptor, _rows_for(panel_id)[0], requested_paired_length=19)
        assert observed.split is None, panel_id

    unnamed_targets = sum(1 for row in _records_rows() if not row["target"])
    assert unnamed_targets >= 1263


# --------------------------------------------------------------------------------------------------
# The aggregate, and the descriptor invariants the new fields carry
# --------------------------------------------------------------------------------------------------


def test_oligogym_is_an_aggregate_and_refuses_to_derive_a_row() -> None:
    """One table, three geometries: the aggregate names its members instead of picking one.

    405 of the 4,113 rows (Martinelli's 49 non-core-19 rows plus all 356 asymmetric Shmushkovich rows)
    are not the paired-core-19 geometry the majority shows, so ingesting the file as one panel would
    stamp the wrong architecture on them.
    """
    descriptor = describe_panel("oligogym")
    assert descriptor.aggregate_of == ("ichihara", "martinelli", "shmushkovich")
    assert descriptor.data_present is True

    with pytest.raises(ValueError, match="aggregate") as excinfo:
        derive_observation(descriptor, dict(_records_rows()[0]), requested_paired_length=19)
    for member in _OLIGOGYM_MEMBER_PANELS:
        assert member in str(excinfo.value)


def test_data_present_and_vendored_csv_must_agree() -> None:
    """A claim about bytes must name them, and a named path must not be disowned."""
    payload = describe_panel("huesken_subset").model_dump()

    with pytest.raises(ValidationError, match="names no vendored_csv"):
        PanelDescriptor.model_validate({**payload, "vendored_csv": None})
    with pytest.raises(ValidationError, match="claims no data"):
        PanelDescriptor.model_validate({**payload, "data_present": False})


def test_a_declared_site_status_requires_the_columns_that_locate_it() -> None:
    """``synthetic_context_local`` without coordinates -- or coordinates without a status -- is refused."""
    payload = describe_panel("ichihara").model_dump()

    without_columns = {
        **payload,
        "columns": {**payload["columns"], "target_site_start_column": None, "target_site_end_column": None},
    }
    with pytest.raises(ValidationError, match="maps no site start/end columns"):
        PanelDescriptor.model_validate(without_columns)

    unlabelled = {**payload, "target_identity_status": TargetIdentityStatus.UNAVAILABLE}
    with pytest.raises(ValidationError, match="must label"):
        PanelDescriptor.model_validate(unlabelled)

    without_context = {**payload, "columns": {**payload["columns"], "context_type_column": None}}
    with pytest.raises(ValidationError, match="context_type_column"):
        PanelDescriptor.model_validate(without_context)


def test_row_selector_is_a_closed_set_and_fails_loudly_on_a_missing_column() -> None:
    """A table that cannot say which panel a row belongs to must raise, not contribute silently."""
    selector = PanelRowSelector(column="dataset", values=("ichihara_2007_1",))
    assert selector.matches({"dataset": "ichihara_2007_1"}) is True
    # Closed set, not a prefix: a regenerated dataset id must be declared, never swept in by name.
    assert selector.matches({"dataset": "ichihara_2007_1_v2"}) is False
    with pytest.raises(KeyError):
        selector.matches({"benchmark_id": "x"})


def test_a_panel_without_a_selector_claims_every_row() -> None:
    """``huesken_subset``'s table is its own, so it declares no selector and selects everything."""
    descriptor = describe_panel("huesken_subset")
    assert descriptor.row_selector is None
    assert descriptor.selects_row({"anything": "at all"}) is True
