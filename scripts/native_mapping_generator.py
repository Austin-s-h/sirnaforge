"""Generate per-species OligoGym mapping sidecars from explicitly pinned cDNA references."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping
from pathlib import Path
from typing import Any, TextIO

import ahocorasick

SUPPORTED_SPECIES = ("human", "mouse", "rat", "macaque")
SPECIES_EVIDENCE = {
    "human": "measured_species_native",
    "mouse": "cross_species_sequence_projection",
    "rat": "cross_species_sequence_projection",
    "macaque": "cross_species_sequence_projection",
}
OBSERVATION_FIELDS = (
    "benchmark_id",
    "dataset",
    "source_row",
    "full_guide_sequence",
    "paired_guide_sequence",
    "passenger_sequence",
    "paired_length",
    "duplex_pairing_status",
    "guide_overhang_sequence",
    "passenger_overhang_sequence",
    "native_reference_species",
    "native_reference_assembly",
    "native_reference_release",
    "native_reference_sha256",
    "native_mapping_count",
    "native_mapping_status",
    "target_identity_status",
    "mapping_evidence_type",
    "assay_species",
)
MAPPING_FIELDS = (
    "mapping_id",
    "benchmark_id",
    "dataset",
    "query_sequence",
    "mapping_status",
    "mapping_transcript_id",
    "mapping_transcript_header",
    "mapping_gene_id",
    "mapping_gene_symbol",
    "mapping_assembly",
    "mapping_release",
    "site_start_1based",
    "site_end_1based",
    "site_sequence",
    "mapping_strand",
    "transcript_genomic_strand",
    "genomic_contig",
    "mapping_rank",
    "mapping_evidence_type",
)


class NativeMappingInputError(ValueError):
    """Invalid pinned reference or observation input for multi-species mapping."""


def file_sha256(path: Path) -> str:
    """Return the SHA-256 of a reference file without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_reference_manifests(paths: Iterable[Path | str]) -> dict[str, dict[str, str]]:
    """Read all four required, checksum-verified species reference declarations."""
    references: dict[str, dict[str, str]] = {}
    for raw_path in paths:
        manifest_path = Path(raw_path)
        try:
            payload = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise NativeMappingInputError(f"cannot read reference manifest {manifest_path}: {exc}") from exc
        if not isinstance(payload, dict):
            raise NativeMappingInputError(f"reference manifest {manifest_path} must be a JSON object")
        required = ("species", "assembly", "release", "fasta", "sha256")
        if any(not str(payload.get(field, "")).strip() for field in required):
            raise NativeMappingInputError(
                f"reference manifest {manifest_path} must declare species, assembly, release, fasta, and sha256"
            )
        species = str(payload["species"]).lower()
        if species not in SUPPORTED_SPECIES:
            raise NativeMappingInputError(f"unsupported reference species {species!r}")
        if species in references:
            raise NativeMappingInputError(f"duplicate reference manifest for species {species!r}")
        fasta = Path(str(payload["fasta"]))
        if not fasta.is_absolute():
            fasta = manifest_path.parent / fasta
        if not fasta.is_file():
            raise NativeMappingInputError(f"reference FASTA for {species} does not exist: {fasta}")
        expected_digest = str(payload["sha256"]).lower().removeprefix("sha256:")
        if not re.fullmatch(r"[0-9a-f]{64}", expected_digest):
            raise NativeMappingInputError(f"reference manifest for {species} has an invalid SHA-256")
        observed_digest = file_sha256(fasta)
        if expected_digest != observed_digest:
            raise NativeMappingInputError(
                f"reference FASTA checksum mismatch for {species}: declared {expected_digest}, found {observed_digest}"
            )
        references[species] = {
            "species": species,
            "assembly": str(payload["assembly"]).strip(),
            "release": str(payload["release"]).strip(),
            "fasta": str(fasta),
            "sha256": observed_digest,
            "manifest_path": str(manifest_path),
            "manifest_sha256": file_sha256(manifest_path),
        }
    missing = sorted(set(SUPPORTED_SPECIES) - references.keys())
    if missing:
        raise NativeMappingInputError(f"four-species mapping requires reference manifests for: {', '.join(missing)}")
    return {species: references[species] for species in SUPPORTED_SPECIES}


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield transcript id, full FASTA description, and sequence from plain or gzipped FASTA."""
    header: str | None = None
    chunks: list[str] = []
    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(chunks).upper().replace("U", "T")
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header.split()[0], header, "".join(chunks).upper().replace("U", "T")


def header_metadata(description: str) -> dict[str, str]:
    """Extract Ensembl gene, symbol, genome assembly/contig, and genomic strand metadata."""
    gene = re.search(r"\bgene:([^ ]+)", description)
    symbol = re.search(r"\bgene_symbol:([^ ]+)", description)
    location = re.search(
        r"(?:chromosome|scaffold|primary_assembly):([^: ]+):([^: ]+):(-?\d+):(-?\d+):(-?1|1)(?:\s|$)",
        description,
    )
    return {
        "gene_id": gene.group(1) if gene else "",
        "gene_symbol": symbol.group(1) if symbol else "",
        "header_assembly": location.group(1) if location else "",
        "genomic_contig": location.group(2) if location else "",
        "transcript_genomic_strand": ("-" if location and location.group(5) == "-1" else "+") if location else "",
    }


def _automaton(patterns: Iterable[str]) -> ahocorasick.Automaton:
    automaton = ahocorasick.Automaton()
    for pattern in sorted(set(patterns)):
        automaton.add_word(pattern, pattern)
    automaton.make_automaton()
    return automaton


def _dna(sequence: Any) -> str:
    return "".join(str(sequence or "").upper().replace("U", "T").split())


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def generate_species_mappings(
    observation_rows: Iterable[Mapping[str, Any]],
    reference: Mapping[str, str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Map exact paired passenger targets to every transcript hit in one pinned species reference.

    Cross-species outputs carry only source identity/sequences and mapping evidence. Measured efficacy
    labels remain in the source `observations.csv` and are never copied into these projected rows.
    """
    species = str(reference.get("species", "")).lower()
    if species not in SUPPORTED_SPECIES:
        raise NativeMappingInputError(f"unsupported mapping species {species!r}")
    observations = [dict(row) for row in observation_rows]
    evidence_type = SPECIES_EVIDENCE[species]
    queries: dict[str, list[dict[str, Any]]] = defaultdict(list)
    identities: set[str] = set()
    for row in observations:
        benchmark_id = str(row.get("benchmark_id", ""))
        if not benchmark_id or benchmark_id in identities:
            raise NativeMappingInputError(f"empty or duplicate benchmark observation identity {benchmark_id!r}")
        identities.add(benchmark_id)
        query = _dna(row.get("paired_target_sequence", ""))
        try:
            paired_length = int(row.get("paired_length", 0))
        except (TypeError, ValueError):
            paired_length = 0
        pairing = str(row.get("duplex_pairing_status", ""))
        if query and len(query) == paired_length and paired_length >= 19 and pairing in {
            "fully_complementary",
            "paired_core_with_overhang",
        }:
            paired_guide = _dna(row.get("paired_guide_sequence", ""))
            full_guide = _dna(row.get("full_guide_sequence", ""))
            passenger = _dna(row.get("passenger_sequence", ""))
            if (
                paired_guide != full_guide[:paired_length]
                or query != passenger[:paired_length]
                or _reverse_complement(paired_guide) != query
            ):
                raise NativeMappingInputError(
                    f"observation {benchmark_id!r} has paired-target geometry inconsistent with its measured strands"
                )
            queries[query].append(row)

    auto = _automaton(queries)
    hits_by_observation: dict[str, list[dict[str, Any]]] = defaultdict(list)
    serial = 0
    for transcript_id, header, raw_sequence in iter_fasta(Path(reference["fasta"])):
        sequence = raw_sequence.upper().replace("U", "T")
        metadata = header_metadata(header)
        declared_assembly = str(reference["assembly"])
        if metadata["header_assembly"] and metadata["header_assembly"] != declared_assembly:
            raise NativeMappingInputError(
                f"transcript {transcript_id} header assembly {metadata['header_assembly']} conflicts with "
                f"reference manifest {declared_assembly}"
            )
        for end_index, query in auto.iter(sequence):
            if (
                not metadata["gene_id"]
                or not metadata["header_assembly"]
                or not metadata["genomic_contig"]
                or not metadata["transcript_genomic_strand"]
            ):
                raise NativeMappingInputError(
                    f"mapped transcript {transcript_id} lacks Ensembl gene, assembly, contig, or strand metadata"
                )
            start = end_index - len(query) + 2
            serial += 1
            row_ids = queries[query]
            for obs in row_ids:
                benchmark_id = str(obs.get("benchmark_id", ""))
                hit = {
                    "mapping_id": f"{species}_mapping_{serial:012d}_{benchmark_id}",
                    "benchmark_id": benchmark_id,
                    "dataset": str(obs.get("dataset", "")),
                    "query_sequence": query,
                    "mapping_status": "mapped",
                    "mapping_transcript_id": transcript_id,
                    "mapping_transcript_header": header,
                    "mapping_gene_id": metadata["gene_id"],
                    "mapping_gene_symbol": metadata["gene_symbol"],
                    "mapping_assembly": declared_assembly,
                    "mapping_release": str(reference["release"]),
                    "site_start_1based": start,
                    "site_end_1based": end_index + 1,
                    "site_sequence": query,
                    # Coordinates are 1-based within the transcript sequence; this strand states
                    # the exact query orientation within the reference transcript.
                    "mapping_strand": "+",
                    "transcript_genomic_strand": metadata["transcript_genomic_strand"],
                    "genomic_contig": metadata["genomic_contig"],
                    "mapping_rank": 1 + len(hits_by_observation[benchmark_id]),
                    "mapping_evidence_type": evidence_type,
                }
                hits_by_observation[benchmark_id].append(hit)

    observation_results: list[dict[str, Any]] = []
    mapping_results: list[dict[str, Any]] = []
    for obs in observations:
        benchmark_id = str(obs.get("benchmark_id", ""))
        hits = hits_by_observation.get(benchmark_id, [])
        paired_length = int(obs.get("paired_length", 0) or 0)
        is_searchable = any(obs is query_row for query_row in queries.get(
            _dna(obs.get("paired_target_sequence", "")), []
        ))
        status = "mapped" if hits else "no_exact_site" if is_searchable else "unavailable_unsupported_geometry"
        observation_results.append(
            {
                "benchmark_id": benchmark_id,
                "dataset": str(obs.get("dataset", "")),
                "source_row": obs.get("source_row", ""),
                "full_guide_sequence": obs.get("full_guide_sequence", ""),
                "paired_guide_sequence": obs.get("paired_guide_sequence", ""),
                "passenger_sequence": obs.get("passenger_sequence", ""),
                "paired_length": paired_length,
                "duplex_pairing_status": obs.get("duplex_pairing_status", ""),
                "guide_overhang_sequence": obs.get("guide_overhang_sequence", ""),
                "passenger_overhang_sequence": obs.get("passenger_overhang_sequence", ""),
                "native_reference_species": species,
                "native_reference_assembly": reference["assembly"],
                "native_reference_release": reference["release"],
                "native_reference_sha256": reference["sha256"],
                "native_mapping_count": len(hits),
                "native_mapping_status": status,
                "target_identity_status": "confirmed" if len(hits) == 1 else "ambiguous" if hits else "unavailable",
                "mapping_evidence_type": evidence_type,
                "assay_species": "human",
            }
        )
        mapping_results.extend(hits)
    return observation_results, mapping_results


def write_species_sidecars(
    observation_rows: Iterable[Mapping[str, Any]],
    references: Mapping[str, Mapping[str, str]],
    output_dir: Path,
) -> dict[str, dict[str, Any]]:
    """Write per-species observations and all-hit mapping CSVs, plus a provenance manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    source_rows = [dict(row) for row in observation_rows]
    missing = sorted(set(SUPPORTED_SPECIES) - references.keys())
    if missing:
        raise NativeMappingInputError(f"four-species output requires references for: {', '.join(missing)}")
    summaries: dict[str, dict[str, Any]] = {}
    for species in SUPPORTED_SPECIES:
        reference = references[species]
        obs_rows, hit_rows = generate_species_mappings(source_rows, reference)
        species_dir = output_dir / species
        species_dir.mkdir(parents=True, exist_ok=True)
        with (species_dir / "observations.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=OBSERVATION_FIELDS)
            writer.writeheader()
            writer.writerows(obs_rows)
        with (species_dir / "mappings.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MAPPING_FIELDS)
            writer.writeheader()
            writer.writerows(hit_rows)
        statuses: dict[str, int] = defaultdict(int)
        for row in obs_rows:
            statuses[str(row["native_mapping_status"])] += 1
        summaries[species] = {
            "reference": dict(reference),
            "mapping_evidence_type": SPECIES_EVIDENCE[species],
            "observations": len(obs_rows),
            "mapping_rows": len(hit_rows),
            "mapping_status_counts": dict(sorted(statuses.items())),
            "observations_csv": f"{species}/observations.csv",
            "mappings_csv": f"{species}/mappings.csv",
            "observations_sha256": file_sha256(species_dir / "observations.csv"),
            "mappings_sha256": file_sha256(species_dir / "mappings.csv"),
        }
    manifest = {
        "schema_version": "1.0",
        "description": "Four-species exact paired-target projections against caller-pinned Ensembl cDNA references",
        "assay_species": "human",
        "label_transfer": "none; efficacy labels remain only on the human assay source observations",
        "species": summaries,
        "mapping_coordinates": "1-based inclusive offsets within the mapped transcript sequence",
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return summaries
