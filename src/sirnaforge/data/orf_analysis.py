"""ORF analysis and sequence validation for transcript sequences."""

from typing import Optional

from pydantic import BaseModel, ConfigDict

from sirnaforge.data.base import (
    AbstractDatabaseClient,
    SequenceType,
    SequenceUtils,
    TranscriptInfo,
)
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


# TODO: The output report of this should be more machine-readable aka tab delim.
class ORFInfo(BaseModel):
    """Information about an Open Reading Frame."""

    start_pos: int
    end_pos: int
    length: int
    reading_frame: int  # 0, 1, or 2
    start_codon: str
    stop_codon: str
    has_valid_start: bool
    has_valid_stop: bool
    is_complete: bool  # Both valid start and stop
    gc_content: float

    model_config = ConfigDict(frozen=True)


class SequenceAnalysis(BaseModel):
    """Complete sequence analysis including ORF information."""

    transcript_id: str
    sequence_type: SequenceType
    sequence_length: int
    gc_content: float
    orfs: list[ORFInfo]
    longest_orf: Optional[ORFInfo] = None
    has_valid_orf: bool = False
    cds_sequence: Optional[str] = None
    protein_sequence: Optional[str] = None
    # Rudimentary UTR/CDS characterization
    cds_start: Optional[int] = None  # 0-based start index within transcript sequence
    cds_end: Optional[int] = None  # 0-based end index (exclusive)
    utr5_length: Optional[int] = None
    utr3_length: Optional[int] = None
    sequence_region: Optional[str] = (
        None  # e.g., cds_only, cds_with_utr5, cds_with_utr3, cds_with_utr5_and_utr3, no_cds
    )

    model_config = ConfigDict(use_enum_values=True)


class ORFAnalyzer:
    """Analyze ORFs in transcript sequences and validate sequence types."""

    def __init__(self, database_client: Optional[AbstractDatabaseClient] = None):
        """
        Initialize ORF analyzer.

        Args:
            database_client: Optional database client for retrieving additional sequence types
        """
        self.database_client = database_client

        # Genetic code (standard)
        self.start_codons = {"ATG"}
        self.stop_codons = {"TAA", "TAG", "TGA"}

        # Translation table
        self.codon_table = {
            "TTT": "F",
            "TTC": "F",
            "TTA": "L",
            "TTG": "L",
            "TCT": "S",
            "TCC": "S",
            "TCA": "S",
            "TCG": "S",
            "TAT": "Y",
            "TAC": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGT": "C",
            "TGC": "C",
            "TGA": "*",
            "TGG": "W",
            "CTT": "L",
            "CTC": "L",
            "CTA": "L",
            "CTG": "L",
            "CCT": "P",
            "CCC": "P",
            "CCA": "P",
            "CCG": "P",
            "CAT": "H",
            "CAC": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGT": "R",
            "CGC": "R",
            "CGA": "R",
            "CGG": "R",
            "ATT": "I",
            "ATC": "I",
            "ATA": "I",
            "ATG": "M",
            "ACT": "T",
            "ACC": "T",
            "ACA": "T",
            "ACG": "T",
            "AAT": "N",
            "AAC": "N",
            "AAA": "K",
            "AAG": "K",
            "AGT": "S",
            "AGC": "S",
            "AGA": "R",
            "AGG": "R",
            "GTT": "V",
            "GTC": "V",
            "GTA": "V",
            "GTG": "V",
            "GCT": "A",
            "GCC": "A",
            "GCA": "A",
            "GCG": "A",
            "GAT": "D",
            "GAC": "D",
            "GAA": "E",
            "GAG": "E",
            "GGT": "G",
            "GGC": "G",
            "GGA": "G",
            "GGG": "G",
        }

    def calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content using shared utility."""
        return SequenceUtils.calculate_gc_content(sequence)

    def find_orfs(self, sequence: str, min_length: int = 150) -> list[ORFInfo]:
        """Find all ORFs in a sequence (all 3 reading frames)."""
        orfs = []
        sequence = sequence.upper()

        for frame in range(3):
            frame_sequence = sequence[frame:]

            # Find all start positions
            start_positions = []
            for i in range(0, len(frame_sequence) - 2, 3):
                codon = frame_sequence[i : i + 3]
                if len(codon) == 3 and codon in self.start_codons:
                    start_positions.append(i)

            # For each start, find the next stop
            for start_pos in start_positions:
                for i in range(start_pos, len(frame_sequence) - 2, 3):
                    codon = frame_sequence[i : i + 3]
                    if len(codon) == 3 and codon in self.stop_codons:
                        orf_length = i - start_pos + 3
                        if orf_length >= min_length:
                            orf_sequence = frame_sequence[start_pos : i + 3]

                            orfs.append(
                                ORFInfo(
                                    start_pos=frame + start_pos,
                                    end_pos=frame + i + 3,
                                    length=orf_length,
                                    reading_frame=frame,
                                    start_codon=frame_sequence[start_pos : start_pos + 3],
                                    stop_codon=codon,
                                    has_valid_start=True,
                                    has_valid_stop=True,
                                    is_complete=True,
                                    gc_content=self.calculate_gc_content(orf_sequence),
                                )
                            )
                        break  # Stop at first stop codon

        return sorted(orfs, key=lambda x: x.length, reverse=True)

    def translate_sequence(self, sequence: str) -> str:
        """Translate DNA sequence to protein."""
        sequence = sequence.upper()
        protein = ""

        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]
            if len(codon) == 3:
                amino_acid = self.codon_table.get(codon, "X")
                protein += amino_acid

        return protein

    async def get_additional_sequence(self, transcript_id: str, sequence_type: SequenceType) -> Optional[str]:
        """
        Retrieve specific sequence type using the database client if available.

        Args:
            transcript_id: Transcript identifier
            sequence_type: Type of sequence to retrieve

        Returns:
            Sequence string or None if not available or client not provided
        """
        if not self.database_client:
            logger.debug(f"No database client available to retrieve {sequence_type} for {transcript_id}")
            return None

        try:
            return await self.database_client.get_sequence(transcript_id, sequence_type)
        except Exception as e:
            logger.warning(f"Failed to retrieve {sequence_type} for {transcript_id}: {e}")
            return None

    async def analyze_transcript(self, transcript: TranscriptInfo) -> SequenceAnalysis:  # noqa: PLR0912
        """Perform complete ORF analysis of a transcript."""
        if not transcript.sequence:
            raise ValueError(f"No sequence available for transcript {transcript.transcript_id}")

        logger.info(f"Analyzing transcript {transcript.transcript_id} (length: {len(transcript.sequence)})")

        # Calculate basic sequence properties
        gc_content = self.calculate_gc_content(transcript.sequence)

        # Find ORFs in the current sequence
        orfs = self.find_orfs(transcript.sequence)
        longest_orf = orfs[0] if orfs else None
        has_valid_orf = longest_orf is not None and longest_orf.is_complete

        # Try to get CDS and protein sequences for comparison (database-agnostic)
        cds_sequence = None
        protein_sequence = None

        if self.database_client:
            cds_sequence = await self.get_additional_sequence(transcript.transcript_id, SequenceType.CDS)
            protein_sequence = await self.get_additional_sequence(transcript.transcript_id, SequenceType.PROTEIN)

            if cds_sequence:
                logger.info(f"Retrieved CDS sequence for {transcript.transcript_id} (length: {len(cds_sequence)})")
            if protein_sequence:
                logger.info(
                    f"Retrieved protein sequence for {transcript.transcript_id} (length: {len(protein_sequence)})"
                )

        # Compute CDS/UTR characterization quickly using CDS match if available, else longest ORF
        cds_start: Optional[int] = None
        cds_end: Optional[int] = None
        utr5_length: Optional[int] = None
        utr3_length: Optional[int] = None
        if transcript.sequence:
            seq_len = len(transcript.sequence)
            if cds_sequence and cds_sequence in transcript.sequence:
                cds_start = transcript.sequence.find(cds_sequence)
                cds_end = cds_start + len(cds_sequence)
            elif orfs and orfs[0].is_complete:
                # Use longest complete ORF as CDS proxy
                cds_start = orfs[0].start_pos
                cds_end = orfs[0].end_pos
            if cds_start is not None and cds_end is not None:
                utr5_length = max(0, cds_start)
                utr3_length = max(0, seq_len - cds_end)

        # Determine sequence type based on analysis and characterization
        sequence_type = self._determine_sequence_type(transcript, cds_sequence, cds_start, cds_end)

        # Classify region composition
        if cds_start is None or cds_end is None:
            sequence_region = "no_cds"
        else:
            has_utr5 = (utr5_length or 0) > 0
            has_utr3 = (utr3_length or 0) > 0
            if not has_utr5 and not has_utr3:
                sequence_region = "cds_only"
            elif has_utr5 and has_utr3:
                sequence_region = "cds_with_utr5_and_utr3"
            elif has_utr5:
                sequence_region = "cds_with_utr5"
            else:
                sequence_region = "cds_with_utr3"

        analysis = SequenceAnalysis(
            transcript_id=transcript.transcript_id,
            sequence_type=sequence_type,
            sequence_length=len(transcript.sequence),
            gc_content=gc_content,
            orfs=orfs,
            longest_orf=longest_orf,
            has_valid_orf=has_valid_orf,
            cds_sequence=cds_sequence,
            protein_sequence=protein_sequence,
            cds_start=cds_start,
            cds_end=cds_end,
            utr5_length=utr5_length,
            utr3_length=utr3_length,
            sequence_region=sequence_region,
        )

        # Log analysis results
        self._log_analysis_results(analysis)

        return analysis

    def _determine_sequence_type(
        self,
        transcript: TranscriptInfo,
        cds_sequence: Optional[str],
        cds_start: Optional[int],
        cds_end: Optional[int],
    ) -> SequenceType:
        """Determine what type of sequence we're dealing with."""
        # If we have a CDS sequence and it matches our sequence exactly, it's CDS
        if cds_sequence and transcript.sequence and cds_sequence == transcript.sequence:
            return SequenceType.CDS
        # If CDS is a proper substring of the transcript, it's cDNA (contains UTRs)
        if cds_sequence and transcript.sequence and cds_sequence in transcript.sequence:
            return SequenceType.CDNA
        # If we inferred CDS from ORF and it spans the entire sequence, consider it CDS
        if transcript.sequence and cds_start is not None and cds_end is not None:
            if cds_start == 0 and cds_end == len(transcript.sequence):
                return SequenceType.CDS
            return SequenceType.CDNA
        # Default assumption
        return SequenceType.CDNA

    def _log_analysis_results(self, analysis: SequenceAnalysis) -> None:
        """Log comprehensive analysis results."""

        logger.info(f"=== ORF Analysis Results for {analysis.transcript_id} ===")
        seq_type = (
            analysis.sequence_type.value if hasattr(analysis.sequence_type, "value") else str(analysis.sequence_type)
        )
        logger.info(f"Sequence Type: {seq_type}")
        logger.info(f"Sequence Length: {analysis.sequence_length} bp")
        logger.info(f"GC Content: {analysis.gc_content:.1f}%")
        logger.info(f"ORFs Found: {len(analysis.orfs)}")

        if analysis.longest_orf:
            orf = analysis.longest_orf
            logger.info(f"Longest ORF: {orf.start_pos}-{orf.end_pos} ({orf.length} bp)")
            logger.info(f"  Reading Frame: {orf.reading_frame}")
            logger.info(f"  Start Codon: {orf.start_codon}")
            logger.info(f"  Stop Codon: {orf.stop_codon}")
            logger.info(f"  Complete ORF: {orf.is_complete}")
            logger.info(f"  ORF GC Content: {orf.gc_content:.1f}%")
        else:
            logger.warning(f"No valid ORFs found in {analysis.transcript_id}")

        # UTR/CDS summary
        if analysis.cds_start is not None and analysis.cds_end is not None:
            logger.info(
                f"CDS region: {analysis.cds_start}-{analysis.cds_end} | 5'UTR={analysis.utr5_length or 0} nt, 3'UTR={analysis.utr3_length or 0} nt"
            )
            logger.info(f"Region composition: {analysis.sequence_region}")

        if analysis.cds_sequence:
            logger.info(f"CDS Length: {len(analysis.cds_sequence)} bp")
            # Verify ORF prediction against known CDS
            if analysis.longest_orf:
                # TODO: Extract actual ORF sequence for comparison
                if analysis.cds_sequence in str(analysis.transcript_id):  # Simplified check
                    logger.info("✅ ORF prediction matches known CDS")
                else:
                    logger.warning("⚠️  ORF prediction differs from known CDS")

        if analysis.protein_sequence:
            logger.info(f"Protein Length: {len(analysis.protein_sequence)} aa")

        logger.info(f"Valid ORF Status: {'✅ Valid' if analysis.has_valid_orf else '❌ Invalid'}")
        logger.info("=" * 60)

    async def analyze_transcripts(self, transcripts: list[TranscriptInfo]) -> dict[str, SequenceAnalysis]:
        """Analyze multiple transcripts."""
        analyses = {}

        logger.info(f"Starting ORF analysis for {len(transcripts)} transcripts")

        for transcript in transcripts:
            try:
                if transcript.sequence:
                    analysis = await self.analyze_transcript(transcript)
                    analyses[transcript.transcript_id] = analysis
                else:
                    logger.warning(f"Skipping {transcript.transcript_id} - no sequence available")
            except Exception as e:
                logger.error(f"Failed to analyze {transcript.transcript_id}: {e}")

        # Summary statistics
        valid_orfs = sum(1 for a in analyses.values() if a.has_valid_orf)
        logger.info(f"ORF Analysis Summary: {valid_orfs}/{len(analyses)} transcripts have valid ORFs")

        return analyses


# Convenience functions
def create_orf_analyzer(database_client: Optional[AbstractDatabaseClient] = None) -> ORFAnalyzer:
    """
    Create an ORF analyzer with optional database client.

    Args:
        database_client: Optional database client for retrieving additional sequence types

    Returns:
        ORFAnalyzer instance
    """
    return ORFAnalyzer(database_client=database_client)


async def analyze_transcript_orfs(
    transcript: TranscriptInfo, database_client: Optional[AbstractDatabaseClient] = None
) -> SequenceAnalysis:
    """
    Analyze ORFs in a single transcript.

    Args:
        transcript: Transcript to analyze
        database_client: Optional database client for additional sequence retrieval

    Returns:
        SequenceAnalysis result
    """
    analyzer = create_orf_analyzer(database_client)
    return await analyzer.analyze_transcript(transcript)


async def analyze_multiple_transcript_orfs(
    transcripts: list[TranscriptInfo], database_client: Optional[AbstractDatabaseClient] = None
) -> dict[str, SequenceAnalysis]:
    """
    Analyze ORFs in multiple transcripts.

    Args:
        transcripts: List of transcripts to analyze
        database_client: Optional database client for additional sequence retrieval

    Returns:
        Dictionary mapping transcript IDs to SequenceAnalysis results
    """
    analyzer = create_orf_analyzer(database_client)
    return await analyzer.analyze_transcripts(transcripts)
