"""Native OligoGym mapping imports retain every exact transcript hit and source provenance."""

import csv
import json
from pathlib import Path

import pytest

from sirnaforge.benchmark.artifact import read_manifest, read_native_mappings, read_observations
from sirnaforge.benchmark.design import design_artifact
from sirnaforge.benchmark.prepare import BenchmarkPrepareError, prepare_artifact

NATIVE = Path(__file__).parents[1] / "data" / "benchmarks" / "oligogym_native_design"


@pytest.mark.unit
def test_prepare_imports_all_native_hits_and_keeps_measured_overhangs(tmp_path: Path) -> None:
    """Imported evidence keeps exact coordinates and both measured OligoGym tails."""
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"species": "human", "assembly": "GRCh38", "release": "Ensembl 110", "sha256": "a" * 64}))

    manifest = prepare_artifact(
        panel_id="ichihara",
        paired_length=19,
        out_dir=tmp_path / "out",
        native_observations_csv=NATIVE / "observations.csv",
        native_mappings_csv=NATIVE / "mappings.csv",
        native_reference_manifest=reference,
    )
    artifact = tmp_path / "out" / "ichihara__len19"
    observations = read_observations(artifact / "observations.csv")
    mappings = read_native_mappings(artifact / "native_mappings.csv")
    by_id = {row.observation_id: row for row in observations}

    assert manifest.artifact_schema_version == "1.1"
    assert manifest.native_reference.species == "human"
    assert manifest.counts.mapped_native == sum(row.target_identity_status == "confirmed" for row in observations)
    assert manifest.counts.mapped_native_ambiguous == sum(row.target_identity_status == "ambiguous" for row in observations)
    assert len(mappings) == sum(row.native_mapping_count for row in observations)
    assert any(row.guide_3p_overhang == "AT" for row in observations)
    for hit in mappings:
        obs = by_id[hit.observation_id]
        assert hit.end_1based - hit.start_1based + 1 == obs.paired_length
        assert hit.assembly == manifest.native_reference.assembly
        assert obs.target_identity_status in {"confirmed", "ambiguous"}


@pytest.mark.unit
def test_native_import_rejects_mapping_outside_measured_paired_core(tmp_path: Path) -> None:
    """A mapping whose reported span exceeds its sequence cannot become confirmed evidence."""
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"species": "human", "assembly": "GRCh38", "release": "Ensembl 110", "sha256": "a" * 64}))
    altered = tmp_path / "mappings.csv"
    with (NATIVE / "mappings.csv").open(newline="") as source, altered.open("w", newline="") as target:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames)
        writer.writeheader()
        changed = False
        for row in reader:
            if not changed and row["mapping_status"] == "mapped":
                row["site_end_1based"] = str(int(row["site_end_1based"]) + 1)
                changed = True
            writer.writerow(row)

    with pytest.raises(BenchmarkPrepareError, match="span site_sequence exactly"):
        prepare_artifact(
            panel_id="ichihara",
            paired_length=19,
            out_dir=tmp_path / "out",
            native_observations_csv=NATIVE / "observations.csv",
            native_mappings_csv=altered,
            native_reference_manifest=reference,
        )


@pytest.mark.unit
def test_martinelli_mixed_measured_lengths_partition_into_fixed_length_artifacts(tmp_path: Path) -> None:
    """Measured 19-nt and 21-nt rows enter only their matching fixed-length artifacts."""
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"species": "human", "assembly": "GRCh38", "release": "Ensembl 110", "sha256": "a" * 64}))
    results = {}
    for requested_length in (19, 21):
        manifest = prepare_artifact(
            panel_id="martinelli",
            paired_length=requested_length,
            out_dir=tmp_path / f"out{requested_length}",
            native_observations_csv=NATIVE / "observations.csv",
            native_mappings_csv=NATIVE / "mappings.csv",
            native_reference_manifest=reference,
        )
        artifact = tmp_path / f"out{requested_length}" / f"martinelli__len{requested_length}"
        rows = read_observations(artifact / "observations.csv")
        results[requested_length] = (manifest, rows)

    manifest19, rows19 = results[19]
    manifest21, rows21 = results[21]
    assert manifest19.counts.observations_kept == 858
    assert manifest19.counts.observations_incompatible == 49
    assert manifest21.counts.observations_kept == 18
    assert manifest21.counts.observations_incompatible == 889
    assert sum(row.compatibility_status == "compatible" and row.paired_length == 19 for row in rows19) == 858
    assert sum(row.compatibility_status == "compatible" and row.paired_length == 21 for row in rows21) == 18
    assert sum(row.duplex_pairing_status == "non_complementary" for row in rows19) == 30
    assert sum(row.duplex_pairing_status == "partial_unresolved" for row in rows19) == 1
    assert all(row.target_identity_status in {"confirmed", "ambiguous"} for row in rows21 if row.native_mapping_count)
    # A 21-nt native site remains mappable in the 19-nt design artifact while review-only for design.
    long_rows = [row for row in rows19 if row.paired_length == 21]
    assert len(long_rows) == 18
    assert all(row.compatibility_status == "incompatible" and row.native_mapping_count > 0 for row in long_rows)
    designed = design_artifact(artifact_dir=tmp_path / "out21" / "martinelli__len21")
    designed_manifest = read_manifest(designed.manifest_path)
    assert designed_manifest.outputs.native_mappings_csv is not None
    assert designed_manifest.counts.observations_kept == 18


@pytest.mark.unit
def test_shmushkovich_native_rows_are_preserved_as_review_only(tmp_path: Path) -> None:
    """Unsupported 15/20 measured architectures stay auditable and never enter design FASTA."""
    reference = tmp_path / "reference.json"
    reference.write_text(json.dumps({"species": "human", "assembly": "GRCh38", "release": "test", "sha256": "a" * 64}))
    manifest = prepare_artifact(
        panel_id="shmushkovich",
        paired_length=19,
        out_dir=tmp_path / "out",
        native_observations_csv=NATIVE / "observations.csv",
        native_mappings_csv=NATIVE / "mappings.csv",
        native_reference_manifest=reference,
    )
    artifact = tmp_path / "out" / "shmushkovich__len19"
    rows = read_observations(artifact / "observations.csv")

    assert manifest.counts.observations_in_source == 356
    assert manifest.counts.observations_incompatible == 356
    assert all(row.architecture == "unsupported" for row in rows)
    assert all(row.compatibility_status == "incompatible" for row in rows)
    assert all(row.compatibility_reason.startswith("review_only_") for row in rows)
    assert (artifact / "design_inputs.fasta").read_text() == ""
