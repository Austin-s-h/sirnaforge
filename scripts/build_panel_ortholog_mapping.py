#!/usr/bin/env python3
"""Resolve and save offline Ensembl Compara mappings for a multi-gene benchmark panel."""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
from collections import defaultdict
from pathlib import Path

from sirnaforge.data.orthology import resolve_orthologues


async def build(mapping_path: Path, dataset: str, query_key: str, target_species: set[str]) -> tuple[dict, dict]:
    """Resolve panel gene orthologues and return merged mappings with provenance."""
    gene_ids = set()
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("dataset") == dataset and row.get("mapping_gene_id"):
                gene_ids.add(row["mapping_gene_id"])
    merged: dict[str, set[str]] = defaultdict(set)
    per_gene: dict[str, dict[str, list[str]]] = {}
    for gene_id in sorted(gene_ids):
        result = await resolve_orthologues(
            {gene_id},
            "human",
            target_species,
            budget=60.0,
        )
        per_gene[gene_id] = {species: sorted(ids) for species, ids in sorted(result.gene_ids_by_species.items())}
        for species, ids in result.gene_ids_by_species.items():
            merged[species].update(ids)
        print(f"{dataset} {gene_id}: {result.summary()}", flush=True)
    panel_mapping = {query_key: {species: sorted(ids) for species, ids in sorted(merged.items())}}
    metadata = {
        "dataset": dataset,
        "query_key": query_key,
        "query_species": "human",
        "target_species": sorted(target_species),
        "source_gene_ids": sorted(gene_ids),
        "per_gene": per_gene,
        "merged_mapping": panel_mapping[query_key],
        "provenance": "Ensembl Compara via sirnaforge.data.orthology.resolve_orthologues",
    }
    return panel_mapping, metadata


def main() -> None:
    """Parse command-line arguments and write the requested benchmark artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping-path", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--query-key", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    mapping, metadata = asyncio.run(build(args.mapping_path, args.dataset, args.query_key, {"mouse", "rat", "macaque"}))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(mapping, indent=2) + "\n", encoding="utf-8")
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
