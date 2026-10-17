#!/usr/bin/env python3
"""Build memory-safe representative references from native Ensembl cDNA FASTA files.

The complete Ensembl cDNA files are valid native references but the human file is too large for
BWA-MEM2 indexing on a 26.8 GiB host.  This profile retains every protein-coding transcript for
genes listed in the benchmark mapping table and the longest protein-coding transcript for every
other gene.  It is therefore suitable for a broad gene-level screen, while its reduced isoform
scope is explicit and auditable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def _iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    header = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.upper().replace("U", "T"))
    if header is not None:
        yield header.split()[0], header, "".join(chunks)


def _metadata(header: str) -> dict[str, str]:
    def find(name: str) -> str:
        match = re.search(rf"\b{name}:([^\s]+)", header)
        return match.group(1) if match else ""

    return {
        "gene_id": find("gene"),
        "gene_symbol": find("gene_symbol"),
        "gene_biotype": find("gene_biotype"),
        "transcript_biotype": find("transcript_biotype"),
    }


def _protein_coding(metadata: dict[str, str]) -> bool:
    return "protein_coding" in {metadata["gene_biotype"], metadata["transcript_biotype"]}


def _preserved_gene_ids(mapping_path: Path | None) -> set[str]:
    if mapping_path is None:
        return set()
    ids: set[str] = set()
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("mapping_gene_id"):
                ids.add(row["mapping_gene_id"])
    return ids


def _write_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def prepare(input_fasta: Path, output_fasta: Path, mapping_path: Path | None = None) -> dict[str, object]:
    """Select protein-coding gene representatives while preserving mapped transcripts."""
    preserved_genes = _preserved_gene_ids(mapping_path)
    best_by_gene: dict[str, tuple[int, str, str, str, dict[str, str]]] = {}
    total = 0
    protein_coding = 0
    for record_id, header, sequence in _iter_fasta(input_fasta):
        total += 1
        metadata = _metadata(header)
        if not metadata["gene_id"] or not _protein_coding(metadata):
            continue
        protein_coding += 1
        current = (len(sequence), record_id, header, sequence, metadata)
        previous = best_by_gene.get(metadata["gene_id"])
        if previous is None or (len(sequence), record_id) > (previous[0], previous[1]):
            best_by_gene[metadata["gene_id"]] = current

    selected = {record[1]: record for record in best_by_gene.values()}
    preserved_transcripts = 0
    if preserved_genes:
        for record_id, header, sequence in _iter_fasta(input_fasta):
            metadata = _metadata(header)
            if metadata["gene_id"] in preserved_genes and _protein_coding(metadata):
                selected[record_id] = (len(sequence), record_id, header, sequence, metadata)
                preserved_transcripts += 1

    records = sorted(selected.values(), key=lambda item: item[1])
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as handle:
        for _, _, header, sequence, _ in records:
            _write_record(handle, header, sequence)

    digest = hashlib.sha256(output_fasta.read_bytes()).hexdigest()
    return {
        "input_fasta": str(input_fasta),
        "output_fasta": str(output_fasta),
        "input_records": total,
        "input_protein_coding_records": protein_coding,
        "selected_records": len(records),
        "selected_genes": len(best_by_gene),
        "preserved_mapping_genes": len(preserved_genes),
        "preserved_mapping_transcripts_seen": preserved_transcripts,
        "profile": "all protein-coding transcripts for mapped genes; longest protein-coding transcript per other gene",
        "sha256": digest,
    }


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--mapping-path", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    result = prepare(args.input_fasta, args.output_fasta, args.mapping_path)
    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
