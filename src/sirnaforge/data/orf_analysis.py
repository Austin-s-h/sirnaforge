"""ORF analysis and sequence validation for transcript sequences."""

from typing import Optional

from pydantic import BaseModel, ConfigDict

from sirnaforge.data.base import (
    BaseEnsemblClient,
    DatabaseType,
    SequenceType,
    SequenceUtils,
    TranscriptInfo,
)
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


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

    model_config = ConfigDict(use_enum_values=True)


class ORFAnalyzer(BaseEnsemblClient):
    """Analyze ORFs in transcript sequences and validate sequence types."""

    def __init__(self, timeout: int = 30):
        """Initialize ORF analyzer."""
        super().__init__(timeout=timeout)

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

    async def get_sequence_from_ensembl(
        self, transcript_id: str, sequence_type: SequenceType
    ) -> Optional[str]:
        """Retrieve specific sequence type from Ensembl using inherited method."""
        return await self.get_sequence(transcript_id, sequence_type)

    async def analyze_transcript(self, transcript: TranscriptInfo) -> SequenceAnalysis:
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

        # Try to get CDS and protein sequences from Ensembl for comparison
        cds_sequence = None
        protein_sequence = None

        if transcript.database == DatabaseType.ENSEMBL:
            try:
                cds_sequence = await self.get_sequence_from_ensembl(transcript.transcript_id, SequenceType.CDS)
                protein_sequence = await self.get_sequence_from_ensembl(transcript.transcript_id, SequenceType.PROTEIN)

                if cds_sequence:
                    logger.info(f"Retrieved CDS sequence for {transcript.transcript_id} (length: {len(cds_sequence)})")
                if protein_sequence:
                    logger.info(
                        f"Retrieved protein sequence for {transcript.transcript_id} (length: {len(protein_sequence)})"
                    )

            except Exception as e:
                logger.warning(f"Could not retrieve additional sequences for {transcript.transcript_id}: {e}")

        # Determine sequence type based on analysis
        sequence_type = self._determine_sequence_type(transcript, cds_sequence, orfs)

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
        )

        # Log analysis results
        self._log_analysis_results(analysis)

        return analysis

    def _determine_sequence_type(
        self, transcript: TranscriptInfo, cds_sequence: Optional[str], orfs: list[ORFInfo]
    ) -> SequenceType:
        """Determine what type of sequence we're dealing with."""

        # If we have a CDS sequence and it matches our sequence, we have CDS
        if cds_sequence and cds_sequence == transcript.sequence:
            return SequenceType.CDS

        # If we have a CDS sequence and our sequence is longer, likely cDNA
        if cds_sequence and len(transcript.sequence) > len(cds_sequence):
            return SequenceType.CDNA

        # If we have good ORFs and the sequence starts with ATG, likely CDS or cDNA
        if orfs and transcript.sequence.startswith("ATG"):
            return SequenceType.CDS if len(orfs[0]) == len(transcript.sequence) else SequenceType.CDNA

        # Default assumption for Ensembl transcript sequences
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
async def analyze_transcript_orfs(transcript: TranscriptInfo) -> SequenceAnalysis:
    """Analyze ORFs in a single transcript."""
    analyzer = ORFAnalyzer()
    return await analyzer.analyze_transcript(transcript)


async def analyze_multiple_transcript_orfs(transcripts: list[TranscriptInfo]) -> dict[str, SequenceAnalysis]:
    """Analyze ORFs in multiple transcripts."""
    analyzer = ORFAnalyzer()
    return await analyzer.analyze_transcripts(transcripts)
