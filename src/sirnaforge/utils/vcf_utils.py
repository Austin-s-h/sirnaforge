"""Lightweight VCF reader helpers for tests."""

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass
class VcfRecord:
    chrom: str
    pos: int
    ref: str
    alts: list[str]


def read_vcf_records(path: Path) -> Iterator[VcfRecord]:
    """Yield VCF records with minimal fields."""
    vcf = pysam.VariantFile(str(path))
    for rec in vcf:
        yield VcfRecord(rec.chrom if rec.chrom.startswith("chr") else f"chr{rec.chrom}", rec.pos, rec.ref, list(rec.alts))
    vcf.close()
