"""Pandera schemas for siRNAforge data validation.

This module defines pandera schemas for validating the structure and content
of various table-like outputs from the siRNAforge pipeline.

Use schemas directly: MySchema.validate(df) - let validation errors bubble up.
"""

import pandas as pd
import pandera.pandas as pa
from pandera.pandas import Column, DataFrameSchema

# Schema for siRNA candidate results (CSV output)
SiRNACandidateSchema = DataFrameSchema(
    {
        "id": Column(str, description="Unique identifier for siRNA candidate"),
        "transcript_id": Column(str, description="Source transcript identifier"),
        "position": Column(int, checks=pa.Check.ge(1), description="1-based position in transcript"),
        "guide_sequence": Column(
            str,
            checks=[pa.Check.str_length(19, 23), pa.Check.str_matches(r"^[ATCGU]+$")],
            description="Guide (antisense) sequence",
        ),
        "passenger_sequence": Column(
            str,
            checks=[pa.Check.str_length(19, 23), pa.Check.str_matches(r"^[ATCGU]+$")],
            description="Passenger (sense) sequence",
        ),
        "gc_content": Column(float, checks=pa.Check.in_range(0.0, 100.0), description="GC content percentage"),
        "asymmetry_score": Column(
            float, checks=pa.Check.in_range(0.0, 1.0), description="Thermodynamic asymmetry score"
        ),
        "paired_fraction": Column(
            float, checks=pa.Check.in_range(0.0, 1.0), description="Fraction of paired bases in secondary structure"
        ),
        "off_target_count": Column(int, checks=pa.Check.ge(0), description="Number of potential off-targets"),
        "transcript_hit_count": Column(
            int, checks=pa.Check.ge(0), description="Number of input transcripts containing this guide"
        ),
        "transcript_hit_fraction": Column(
            float, checks=pa.Check.in_range(0.0, 1.0), description="Fraction of input transcripts hit by this guide"
        ),
        "composite_score": Column(float, checks=pa.Check.in_range(0.0, 100.0), description="Final composite score"),
        "passes_filters": Column(bool, description="Whether candidate passes quality filters"),
    },
    strict=True,
    coerce=True,
)


# Schema for ORF validation report (tab-delimited output)
ORFValidationSchema = DataFrameSchema(
    {
        "transcript_id": Column(str, description="Transcript identifier"),
        "sequence_length": Column(int, checks=pa.Check.ge(1), description="Sequence length in nucleotides"),
        "gc_content": Column(float, checks=pa.Check.in_range(0.0, 100.0), description="Overall GC content percentage"),
        "orfs_found": Column(int, checks=pa.Check.ge(0), description="Number of ORFs detected"),
        "has_valid_orf": Column(bool, description="Whether transcript has valid ORF"),
        "longest_orf_start": Column(
            pd.Int64Dtype(), checks=pa.Check.ge(1), nullable=True, description="Start position of longest ORF (1-based)"
        ),
        "longest_orf_end": Column(
            pd.Int64Dtype(), checks=pa.Check.ge(1), nullable=True, description="End position of longest ORF (1-based)"
        ),
        "longest_orf_length": Column(
            pd.Int64Dtype(), checks=pa.Check.ge(0), nullable=True, description="Length of longest ORF in nucleotides"
        ),
        "longest_orf_frame": Column(
            pd.Int64Dtype(),
            checks=pa.Check.in_range(0, 2),
            nullable=True,
            description="Reading frame of longest ORF (0, 1, or 2)",
        ),
        "start_codon": Column(
            str, checks=pa.Check.str_matches(r"^(ATG|)$"), nullable=True, description="Start codon of longest ORF"
        ),
        "stop_codon": Column(
            str,
            checks=pa.Check.str_matches(r"^(TAA|TAG|TGA|)$"),
            nullable=True,
            description="Stop codon of longest ORF",
        ),
        "orf_gc_content": Column(
            float, checks=pa.Check.in_range(0.0, 100.0), nullable=True, description="GC content of longest ORF"
        ),
    },
    strict=True,
    coerce=True,
)


# Schema for off-target analysis results (TSV output)
OffTargetHitsSchema = DataFrameSchema(
    {
        "qname": Column(str, description="Query sequence name/ID"),
        "target_id": Column(str, nullable=True, description="Target sequence identifier"),
        "species": Column(str, nullable=True, description="Target species"),
        "chromosome": Column(str, nullable=True, description="Chromosome/contig"),
        "position": Column(pd.Int64Dtype(), checks=pa.Check.ge(1), nullable=True, description="Genomic position"),
        "strand": Column(str, checks=pa.Check.str_matches(r"^[+-]$"), nullable=True, description="Strand orientation"),
        "mismatches": Column(pd.Int64Dtype(), checks=pa.Check.ge(0), nullable=True, description="Number of mismatches"),
        "alignment_score": Column(float, nullable=True, description="Alignment score"),
        "offtarget_score": Column(
            float, checks=pa.Check.ge(0.0), nullable=True, description="Off-target penalty score"
        ),
        "target_sequence": Column(
            str,
            checks=pa.Check.str_matches(r"^[ATCGUN-]*$"),
            nullable=True,
            description="Target sequence with alignment",
        ),
    },
    strict=False,  # More lenient for external tool outputs
    coerce=True,
)
