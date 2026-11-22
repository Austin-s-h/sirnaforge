"""Deterministic four-species exact mapping against tiny, explicitly pinned Ensembl-like FASTAs."""

import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest

from scripts.map_oligo_benchmarks_to_ensembl import main as mapper_main
from scripts.native_mapping_generator import (
    SUPPORTED_SPECIES,
    NativeMappingInputError,
    generate_species_mappings,
    load_reference_manifests,
    write_species_sidecars,
)

QUERY = "ACGTTGCAACGATTCGATC"
GUIDE = "GATCGAATCGTTGCAACGT"


def _reference(tmp_path: Path, species: str, assembly: str, fasta_text: str) -> Path:
    fasta = tmp_path / f"{species}.fa"
    fasta.write_text(fasta_text)
    manifest = tmp_path / f"{species}.json"
    manifest.write_text(
        json.dumps(
            {
                "species": species,
                "assembly": assembly,
                "release": "Ensembl test fixture 1",
                "fasta": fasta.name,
                "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
            }
        )
    )
    return manifest


def _observations() -> list[dict[str, object]]:
    return [
        {
            "benchmark_id": "dataset__row_000001",
            "dataset": "dataset",
            "source_row": 1,
            "full_guide_sequence": GUIDE,
            "paired_guide_sequence": GUIDE,
            "passenger_sequence": QUERY,
            "paired_target_sequence": QUERY,
            "paired_length": 19,
            "duplex_pairing_status": "fully_complementary",
            "guide_overhang_sequence": "",
            "passenger_overhang_sequence": "",
            "efficacy_higher_is_better": 91.0,
        },
        {
            "benchmark_id": "dataset__row_000002",
            "dataset": "dataset",
            "source_row": 2,
            "full_guide_sequence": GUIDE[:15],
            "paired_guide_sequence": GUIDE[:15],
            "passenger_sequence": QUERY[:15],
            "paired_target_sequence": QUERY[:15],
            "paired_length": 15,
            "duplex_pairing_status": "paired_core_with_overhang",
            "guide_overhang_sequence": "",
            "passenger_overhang_sequence": "",
            "efficacy_higher_is_better": 9.0,
        },
    ]


@pytest.mark.unit
def test_four_species_sidecars_keep_projection_status_and_all_exact_hits(tmp_path: Path) -> None:
    """All-hit locations and species-specific status stay explicit in generated sidecars."""
    manifests = [
        _reference(
            tmp_path,
            "human",
            "GRCh38",
            f">ENSTH cdna chromosome:GRCh38:1:100:500:1 gene:ENSGH gene_symbol:HG\nGG{QUERY}CCC\n",
        ),
        _reference(
            tmp_path,
            "mouse",
            "GRCm39",
            f">ENSTM1 cdna chromosome:GRCm39:2:100:500:-1 gene:ENSGM1 gene_symbol:MG1\nA{QUERY}TT\n"
            f">ENSTM2 cdna chromosome:GRCm39:3:100:500:1 gene:ENSGM2 gene_symbol:MG2\n{QUERY}\n",
        ),
        _reference(
            tmp_path,
            "rat",
            "GRCr8",
            ">ENSRT cdna primary_assembly:GRCr8:4:100:500:1 gene:ENSRG gene_symbol:RG\nTTTTTTTTTT\n",
        ),
        _reference(
            tmp_path,
            "macaque",
            "Mmul_10",
            f">ENSMM cdna primary_assembly:Mmul_10:5:100:500:-1 gene:ENSGMM gene_symbol:MM\n{QUERY}\n",
        ),
    ]
    references = load_reference_manifests(manifests)
    summaries = write_species_sidecars(_observations(), references, tmp_path / "out")

    assert tuple(references) == SUPPORTED_SPECIES
    assert summaries["human"]["mapping_status_counts"] == {"mapped": 1, "unavailable_unsupported_geometry": 1}
    assert (
        summaries["mouse"]["mappings_sha256"]
        == hashlib.sha256((tmp_path / "out" / "mouse" / "mappings.csv").read_bytes()).hexdigest()
    )
    assert summaries["mouse"]["mapping_status_counts"] == {"mapped": 1, "unavailable_unsupported_geometry": 1}
    assert summaries["rat"]["mapping_status_counts"] == {"no_exact_site": 1, "unavailable_unsupported_geometry": 1}
    assert summaries["macaque"]["mapping_status_counts"] == {"mapped": 1, "unavailable_unsupported_geometry": 1}

    with (tmp_path / "out" / "mouse" / "mappings.csv").open(newline="") as handle:
        mouse_hits = list(csv.DictReader(handle))
    assert len(mouse_hits) == 2
    assert {hit["mapping_transcript_id"] for hit in mouse_hits} == {"ENSTM1", "ENSTM2"}
    assert {hit["mapping_gene_id"] for hit in mouse_hits} == {"ENSGM1", "ENSGM2"}
    assert {hit["mapping_evidence_type"] for hit in mouse_hits} == {"cross_species_sequence_projection"}
    assert {hit["transcript_genomic_strand"] for hit in mouse_hits} == {"+", "-"}
    assert all(int(hit["site_end_1based"]) - int(hit["site_start_1based"]) + 1 == len(QUERY) for hit in mouse_hits)

    with (tmp_path / "out" / "mouse" / "observations.csv").open(newline="") as handle:
        projected_rows = list(csv.DictReader(handle))
    assert all("efficacy_higher_is_better" not in row for row in projected_rows)
    assert all(row["assay_species"] == "human" for row in projected_rows)
    assert projected_rows[0]["native_mapping_status"] == "mapped"
    assert projected_rows[1]["native_mapping_status"] == "unavailable_unsupported_geometry"


@pytest.mark.unit
def test_four_species_mapping_requires_all_pinned_manifests_and_verified_sha(tmp_path: Path) -> None:
    """An incomplete reference set or changed FASTA cannot be reported as unmapped biology."""
    human = _reference(
        tmp_path,
        "human",
        "GRCh38",
        f">ENSTH cdna chromosome:GRCh38:1:100:500:1 gene:ENSGH gene_symbol:HG\n{QUERY}\n",
    )
    with pytest.raises(NativeMappingInputError, match="requires reference manifests for: macaque, mouse, rat"):
        load_reference_manifests([human])

    fasta_path = tmp_path / "human.fa"
    fasta_path.write_text(fasta_path.read_text() + "A\n")
    others = [
        _reference(tmp_path, species, assembly, ">tx\nACGT\n")
        for species, assembly in (("mouse", "GRCm39"), ("rat", "GRCr8"), ("macaque", "Mmul_10"))
    ]
    with pytest.raises(NativeMappingInputError, match="checksum mismatch for human"):
        load_reference_manifests([human, *others])


@pytest.mark.unit
def test_native_search_rejects_a_paired_target_that_disagrees_with_the_strands(tmp_path: Path) -> None:
    """A malformed query cannot turn an unrelated exact site into a native mapping."""
    manifest = _reference(
        tmp_path,
        "human",
        "GRCh38",
        f">ENSTH cdna chromosome:GRCh38:1:100:500:1 gene:ENSGH\n{QUERY}\n",
    )
    reference = load_reference_manifests(
        [
            manifest,
            *[
                _reference(tmp_path, species, assembly, ">tx\nACGT\n")
                for species, assembly in (("mouse", "GRCm39"), ("rat", "GRCr8"), ("macaque", "Mmul_10"))
            ],
        ]
    )["human"]
    row = _observations()[0]
    row["paired_target_sequence"] = "T" + QUERY[1:]
    with pytest.raises(NativeMappingInputError, match="geometry inconsistent"):
        generate_species_mappings([row], reference)


@pytest.mark.unit
def test_four_species_mapper_cli_writes_reproducible_species_directories(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """The mapper CLI emits isolated all-hit sidecars for all four declared references."""
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    with gzip.open(source_dir / "ichihara_2007_1.csv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("fasta", "y", "y_raw", "targets"))
        writer.writeheader()
        writer.writerow({"fasta": f"{QUERY}.{GUIDE}", "y": "91", "y_raw": "91", "targets": "GENE1"})

    manifests = [
        _reference(tmp_path, "human", "GRCh38", f">ENSTH cdna chromosome:GRCh38:1:100:500:1 gene:ENSGH\n{QUERY}\n"),
        _reference(tmp_path, "mouse", "GRCm39", f">ENSTM cdna chromosome:GRCm39:2:100:500:-1 gene:ENSGM\n{QUERY}\n"),
        _reference(tmp_path, "rat", "GRCr8", f">ENSRT cdna primary_assembly:GRCr8:4:100:500:1 gene:ENSRG\n{QUERY}\n"),
        _reference(
            tmp_path, "macaque", "Mmul_10", f">ENSMM cdna primary_assembly:Mmul_10:5:100:500:-1 gene:ENSGMM\n{QUERY}\n"
        ),
    ]
    arguments = [
        "map_oligo_benchmarks_to_ensembl",
        "--source-dir",
        str(source_dir),
        "--output-dir",
        str(tmp_path / "mapped"),
        "--dataset",
        "ichihara_2007_1",
    ]
    for reference_manifest in manifests:
        arguments.extend(("--species-reference-manifest", str(reference_manifest)))
    monkeypatch.setattr(sys, "argv", arguments)
    mapper_main()
    capsys.readouterr()
    manifest = json.loads((tmp_path / "mapped" / "manifest.json").read_text())

    assert set(manifest["species"]) == set(SUPPORTED_SPECIES)
    assay_fields = next(csv.DictReader((tmp_path / "mapped" / "assay_observations.csv").open(newline="")))
    assert assay_fields["efficacy_higher_is_better"] == "91.0"
    for species in SUPPORTED_SPECIES:
        species_dir = tmp_path / "mapped" / "species" / species
        with (species_dir / "observations.csv").open(newline="") as handle:
            species_rows = list(csv.DictReader(handle))
        with (species_dir / "mappings.csv").open(newline="") as handle:
            species_hits = list(csv.DictReader(handle))
        assert len(species_rows) == 1
        assert len(species_hits) == 1
        assert species_rows[0]["native_reference_species"] == species
        assert "efficacy_higher_is_better" not in species_rows[0]
