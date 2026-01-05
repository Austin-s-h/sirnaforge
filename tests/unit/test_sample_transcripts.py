"""Validates synthetic sample transcripts and toy variants."""

from pathlib import Path

import pytest
from Bio import SeqIO

from sirnaforge.data.variant_resolver import VariantResolver


@pytest.fixture(scope="module")
def sample_sequences():
    """Load sample transcripts from FASTA file."""
    fasta_path = Path("examples/sample_transcripts.fasta")
    return {rec.id.split()[0]: rec.seq for rec in SeqIO.parse(fasta_path, "fasta")}


def test_variants_match_reference(sample_sequences):
    """Ensure toy variant refs match bases in the synthetic sample transcripts."""
    vcf_path = Path("examples/variant_demo.vcf")
    resolver = VariantResolver(min_af=0.0)  # keep low-AF demo variants
    records = resolver.read_vcf(vcf_path)
    assert records, "variant_demo.vcf should contain records"

    seq_map = {k.split("_")[0]: v for k, v in sample_sequences.items()}

    for rec in records:
        chrom = rec.chr
        assert chrom in seq_map, f"No transcript for {chrom}"
        seq = seq_map[chrom]
        assert 1 <= rec.pos <= len(seq)
        ref_base = seq[rec.pos - 1].upper()
        assert ref_base == rec.ref.upper(), f"{chrom}:{rec.pos} expected {rec.ref}, found {ref_base}"
