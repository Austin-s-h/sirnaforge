#!/usr/bin/env python
"""Build measured review records without inventing strand or endpoint evidence.

Input extracts preserve source fields, chemistry and measurement strings. Core-only
records are suitable for sequence regression, not full-duplex or native-accessibility
validation. No model, fold or alignment is run by this builder.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any


def canonical_sequence(value: str) -> str:
    """Normalize the comparison alphabet, preserving literal notation on the source row."""
    sequence = "".join(value.split()).upper().replace("U", "T")
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("canonical sequence must contain only A/C/G/T/U")
    return sequence


def reverse_complement(sequence: str) -> str:
    """Reverse-complement a validated DNA-alphabet comparison sequence."""
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _geometry(source: dict[str, Any]) -> dict[str, Any]:
    """Check only the explicitly supplied core, retaining absent full strands as null."""
    length = int(source["paired_length"])
    representation = source["representation"]
    if representation == "paired_core_only":
        sequence = canonical_sequence(str(source["sequence"]))
        orientation = source["orientation"]
        if orientation not in {"guide", "sense"}:
            raise ValueError("source orientation must be verified as guide or sense")
        core = sequence if orientation == "guide" else reverse_complement(sequence)
        if len(core) != length:
            raise ValueError("reported core length disagrees with its sequence")
        return {
            "paired_guide_sequence": core,
            "full_guide_sequence": None,
            "full_passenger_sequence": None,
            "guide_overhang": None,
            "passenger_overhang": None,
            "compatibility_status": "compatible_core_only" if 19 <= length <= 23 else "incompatible_length",
            "geometry_evidence": "declared_paired_core_only",
        }
    if representation not in {"full_duplex", "guide_with_passenger_core", "guide_with_declared_core"}:
        raise ValueError(f"unrecognized sequence representation: {representation}")
    guide = canonical_sequence(str(source["guide_sequence"]))
    start = int(source.get("guide_core_start", 0))
    core = guide[start : start + length]
    passenger = source.get("passenger_sequence")
    passenger = canonical_sequence(str(passenger)) if passenger else None
    p_start = int(source.get("passenger_core_start", 0))
    if start < 0 or p_start < 0 or length < 1 or len(core) != length:
        raise ValueError("reported paired span is outside its strand")
    matched = passenger is None or reverse_complement(passenger[p_start : p_start + length]) == core
    status = "compatible_full_duplex" if representation == "full_duplex" else "compatible_declared_core"
    if not 19 <= length <= 23:
        status = "incompatible_length"
    elif not matched or (representation == "full_duplex" and passenger is None):
        status = "unresolved_pairing"
    return {
        "paired_guide_sequence": core,
        "full_guide_sequence": guide,
        "full_passenger_sequence": passenger if representation == "full_duplex" else None,
        "guide_overhang": guide[start + length :],
        "passenger_overhang": passenger[p_start + length :] if representation == "full_duplex" and passenger else None,
        "compatibility_status": status,
        "geometry_evidence": "verified_complementary_core" if passenger and matched else "source_declared_core",
    }


def _measurement(raw: str, kind: str) -> dict[str, Any]:
    """Preserve label direction, units and censoring; never convert between endpoints."""
    if kind != "measured":
        return {"measurement_value": None, "measurement_relation": None, "label_exclusion": "excluded_not_measured"}
    if raw in {"Effective", "Ineffective"}:
        return {"measurement_value": None, "measurement_relation": "category", "label_exclusion": None}
    match = re.fullmatch(r"\s*(>=|<=|>|<|=)?\s*(.+?)\s*", raw)
    if match is None:
        return {
            "measurement_value": None,
            "measurement_relation": None,
            "label_exclusion": "excluded_invalid_measurement",
        }
    try:
        number = Decimal(match[2])
    except InvalidOperation:
        number = Decimal("NaN")
    if not number.is_finite():
        return {
            "measurement_value": None,
            "measurement_relation": None,
            "label_exclusion": "excluded_invalid_measurement",
        }
    return {"measurement_value": match[2], "measurement_relation": match[1] or "=", "label_exclusion": None}


def normalize_record(source: dict[str, Any]) -> dict[str, Any]:
    """Normalize measured evidence without upgrading its geometry or scientific endpoint."""
    row = dict(source)
    row.update(_geometry(source))
    label = _measurement(str(source["measurement_raw"]), str(source["label_kind"]))
    exclusion = label.pop("label_exclusion")
    row.update(label)
    if exclusion:
        row["compatibility_status"] = exclusion
    row["target_identity_status"] = "not_established"
    row["sequence_cluster"] = hashlib.sha256(row["paired_guide_sequence"].encode()).hexdigest()
    return row


_CSV_COLUMNS = (
    "record_id",
    "collection_id",
    "source_record_id",
    "paired_guide_sequence",
    "paired_length",
    "full_guide_sequence",
    "full_passenger_sequence",
    "representation",
    "geometry_evidence",
    "measurement_raw",
    "measurement_value",
    "measurement_relation",
    "endpoint",
    "units",
    "direction",
    "chemistry",
    "concentration",
    "target_gene",
    "target_accession",
    "compatibility_status",
    "sequence_cluster",
    "existing_sequence_overlap",
    "target_identity_status",
)


def _existing_cores(path: Path | None) -> set[str]:
    """Index demonstrably paired sites in the historical table; never assume asymmetry away."""
    if path is None:
        return set()
    cores: set[str] = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            guide = canonical_sequence(row["guide_sequence"])
            passenger = canonical_sequence(row["passenger_sequence"])
            if guide == reverse_complement(passenger) and 19 <= len(guide) <= 23:
                cores.add(guide)
            elif len(guide) >= 19 and len(passenger) >= 19 and guide[:19] == reverse_complement(passenger[:19]):
                cores.add(guide[:19])
    return cores


def _load_sources(source_dir: Path, existing_cores: set[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Validate source counts/identities before creating any output directory."""
    catalog = json.loads((source_dir / "collections.json").read_text())
    rows: list[dict[str, Any]] = []
    inputs: list[dict[str, Any]] = []
    for collection in catalog["collections"]:
        path = source_dir / collection["extract"]
        if not path.resolve().is_relative_to(source_dir.resolve()):
            raise ValueError("extract path must remain inside its source directory")
        raw_rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        if len(raw_rows) != collection["expected_rows"]:
            raise ValueError(f"source count disagrees for {collection['id']}")
        inputs.append({"path": collection["extract"], "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        for source in raw_rows:
            if source["collection_id"] != collection["id"]:
                raise ValueError("extract row disagrees with its collection")
            row = normalize_record(source)
            row["record_id"] = f"{collection['id']}:{source['source_record_id']}"
            row["existing_sequence_overlap"] = row["paired_guide_sequence"] in existing_cores
            rows.append(row)
    ids = [row["record_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate observation identity; retain dose/replicate in source_record_id")
    return sorted(rows, key=lambda row: row["record_id"]), inputs


def build_dataset(source_dir: Path, out_dir: Path, existing_records_csv: Path | None = None) -> dict[str, Any]:
    """Write deterministic evidence, compatible records and complete exclusion counts."""
    rows, inputs = _load_sources(source_dir, _existing_cores(existing_records_csv))
    accepted = [row for row in rows if row["compatibility_status"].startswith("compatible_")]
    excluded = [row for row in rows if not row["compatibility_status"].startswith("compatible_")]
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, records in (("observations.jsonl", rows), ("excluded.jsonl", excluded)):
        (out_dir / name).write_text(
            "".join(json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n" for row in records)
        )
    with (out_dir / "records.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(accepted)
    manifest = {
        "schema_version": "1.0.0",
        "counts": {
            "source_rows": len(rows),
            "compatible_rows": len(accepted),
            "excluded_rows": len(excluded),
            "sequence_clusters": len({row["sequence_cluster"] for row in accepted}),
        },
        "collections": {
            name: dict(Counter(row["compatibility_status"] for row in rows if row["collection_id"] == name))
            for name in sorted({row["collection_id"] for row in rows})
        },
        "inputs": inputs,
        "outputs": {
            name: hashlib.sha256((out_dir / name).read_bytes()).hexdigest()
            for name in ("observations.jsonl", "excluded.jsonl", "records.csv")
        },
        "interpretation": "Source rows are overlapping assay views, not independent new measurements. No validation split assigned.",
        "existing_sequence_overlap_rows": sum(row["existing_sequence_overlap"] for row in accepted),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    """Build local extracts only; never download, fit a model or run screening."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--existing-records-csv", type=Path)
    args = parser.parse_args()
    manifest = build_dataset(args.source_dir, args.out_dir, args.existing_records_csv)
    print(json.dumps(manifest["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
