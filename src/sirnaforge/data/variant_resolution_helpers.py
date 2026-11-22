"""Strict identity checks shared by public variant adapters.

GRCh38 chromosome accessions are fixed by NCBI assembly GCF_000001405.40:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
"""

import re

from sirnaforge.models.variant import VariantQuery, VariantQueryType, VariantRecord

_CHROMOSOME_VERSIONS = [11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11, 11, 10]
GRCH38_ACCESSIONS = {
    f"NC_{number:06d}.{version}": f"chr{number if number < 23 else ('X' if number == 23 else 'Y')}"
    for number, version in enumerate(_CHROMOSOME_VERSIONS, start=1)
}
GRCH38_ACCESSIONS["NC_012920.1"] = "chrMT"


class VariantResolutionAmbiguity(ValueError):
    """A query identifies multiple distinct genomic alleles and needs clarification."""


def normalize_chromosome(value: str) -> str | None:
    """Normalize human primary chromosome names, rejecting contigs and invalid IDs."""
    chrom = re.sub(r"^chr", "", value, flags=re.IGNORECASE).upper()
    if chrom in {"M", "MT"}:
        return "chrMT"
    if chrom in {"X", "Y"}:
        return f"chr{chrom}"
    if chrom.isdigit() and 1 <= int(chrom) <= 22:
        return f"chr{int(chrom)}"
    return None


def exact_query_match(record: VariantRecord, query: VariantQuery) -> bool:
    """Validate assembly and, when present, every explicit coordinate/allele field."""
    if record.assembly != "GRCh38" or normalize_chromosome(record.chr) is None:
        return False
    if query.query_type != VariantQueryType.COORDINATE:
        return True
    return (
        normalize_chromosome(record.chr) == normalize_chromosome(query.chr or "")
        and record.pos == query.pos
        and record.ref.upper() == (query.ref or "").upper()
        and record.alt.upper() == (query.alt or "").upper()
    )


def unique_variant(records: list[VariantRecord], query: VariantQuery) -> VariantRecord | None:
    """Select a proven exact allele, never an arbitrary first mapping or alternate."""
    matches = {record.to_vcf_style(): record for record in records if exact_query_match(record, query)}
    if len(matches) > 1:
        raise VariantResolutionAmbiguity(
            f"{query.raw_input} resolves to multiple GRCh38 alleles: {', '.join(sorted(matches))}; "
            "supply chromosome:position:reference:alternate or an unambiguous nucleotide HGVS identifier"
        )
    return next(iter(matches.values()), None)
