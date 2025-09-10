"""Test script for ORF analysis functionality."""

import asyncio
import json
from pathlib import Path

import pytest

from sirnaforge.data import TranscriptInfo
from sirnaforge.data.base import DatabaseType
from sirnaforge.data.orf_analysis import analyze_multiple_transcript_orfs
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


async def _test_orf_analysis():
    """Async implementation of ORF analysis on TP53 transcripts using saved test data."""

    logger.info("=== Testing ORF Analysis on TP53 Transcripts (Test Data) ===")

    # Load test data from saved file
    test_data_path = Path(__file__).parent / "data" / "tp53_test_data.json"

    with test_data_path.open() as f:
        test_data = json.load(f)

    logger.info("Loaded test data from file")

    if not test_data["success"]:
        logger.error("Test data indicates failure")
        return

    # Convert test data to TranscriptInfo objects
    transcripts = []

    for transcript_data in test_data["transcripts"]:
        transcript = TranscriptInfo(
            transcript_id=transcript_data["transcript_id"],
            transcript_name=transcript_data["transcript_name"],
            transcript_type=transcript_data["transcript_type"],
            gene_id=transcript_data["gene_id"],
            gene_name=transcript_data["gene_name"],
            sequence=transcript_data["sequence"],
            database=DatabaseType.ENSEMBL,
            is_canonical=transcript_data.get("is_canonical", False),
            length=transcript_data.get("length"),
        )
        transcripts.append(transcript)

    logger.info(f"Found {len(transcripts)} transcripts for TP53")

    # Filter to protein-coding transcripts with sequences
    protein_coding = [t for t in transcripts if t.transcript_type == "protein_coding" and t.sequence]

    logger.info(f"Analyzing {len(protein_coding)} protein-coding transcripts...")

    # Perform ORF analysis on all transcripts
    analyses = await analyze_multiple_transcript_orfs(protein_coding)

    # Summary report
    logger.info("\n=== ORF Analysis Summary ===")
    for transcript_id, analysis in analyses.items():
        logger.info(f"\nTranscript: {transcript_id}")
        logger.info(f"  Sequence Type: {analysis.sequence_type}")
        logger.info(f"  Sequence Length: {analysis.sequence_length} bp")
        logger.info(f"  Valid ORF: {analysis.has_valid_orf}")
        if analysis.longest_orf:
            logger.info(f"  Longest ORF: {analysis.longest_orf.length} bp")
        if analysis.cds_sequence:
            logger.info(f"  Known CDS Length: {len(analysis.cds_sequence)} bp")
        if analysis.protein_sequence:
            logger.info(f"  Protein Length: {len(analysis.protein_sequence)} aa")


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
def test_orf_analysis():
    """Synchronous pytest-compatible wrapper that runs the async ORF analysis test."""
    asyncio.run(_test_orf_analysis())


if __name__ == "__main__":
    asyncio.run(_test_orf_analysis())
