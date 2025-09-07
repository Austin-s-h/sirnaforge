"""Shared base classes and utilities for genomic data analysis."""

import re
from enum import Enum
from pathlib import Path
from typing import Optional, Union

import aiohttp
from pydantic import BaseModel, ConfigDict, field_validator

from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class DatabaseType(str, Enum):
    """Supported genomic databases."""

    ENSEMBL = "ensembl"
    REFSEQ = "refseq"
    GENCODE = "gencode"


class SequenceType(str, Enum):
    """Types of sequence data that can be retrieved."""

    CDNA = "cdna"  # Complete cDNA sequence (includes UTRs)
    CDS = "cds"  # Coding sequence only (ORF)
    PROTEIN = "protein"  # Translated protein sequence
    GENOMIC = "genomic"  # Genomic sequence with introns


class GeneInfo(BaseModel):
    """Gene information model."""

    gene_id: str
    gene_name: Optional[str] = None
    gene_type: Optional[str] = None
    chromosome: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[int] = None
    description: Optional[str] = None
    database: DatabaseType

    model_config = ConfigDict(use_enum_values=True)


class TranscriptInfo(BaseModel):
    """Transcript information model."""

    transcript_id: str
    transcript_name: Optional[str] = None
    transcript_type: Optional[str] = None
    gene_id: str
    gene_name: Optional[str] = None
    sequence: Optional[str] = None
    length: Optional[int] = None
    database: DatabaseType
    is_canonical: bool = False

    model_config = ConfigDict(use_enum_values=True)

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v):
        """Validate RNA sequence."""
        if v is not None:
            # Convert to uppercase and check for valid RNA bases
            v = v.upper()
            if not re.match(r"^[ACGTU]*$", v):
                raise ValueError("Sequence contains invalid RNA bases")
        return v


class BaseEnsemblClient:
    """Base class for Ensembl API interactions."""

    def __init__(self, timeout: int = 30, base_url: str = "https://rest.ensembl.org"):
        """Initialize Ensembl client."""
        self.timeout = timeout
        self.base_url = base_url
        self.species = "homo_sapiens"

    async def get_sequence(
        self,
        identifier: str,
        sequence_type: SequenceType = SequenceType.CDNA,
        headers: Optional[dict] = None
    ) -> Optional[str]:
        """
        Get sequence from Ensembl REST API.

        Args:
            identifier: Gene ID, transcript ID, etc.
            sequence_type: Type of sequence to retrieve
            headers: Optional HTTP headers

        Returns:
            Sequence string or None if not found
        """
        # Map sequence type to Ensembl API parameter
        type_mapping = {
            SequenceType.CDNA: "cdna",
            SequenceType.CDS: "cds",
            SequenceType.PROTEIN: "protein",
            SequenceType.GENOMIC: "genomic"
        }

        seq_type = type_mapping.get(sequence_type, "cdna")
        url = f"{self.base_url}/sequence/id/{identifier}?species={self.species}&type={seq_type}"

        if headers is None:
            headers = {"Content-Type": "text/plain"}

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, headers=headers) as response
            ):
                if response.status == 200:
                    sequence = await response.text()
                    # Remove FASTA header if present
                    if sequence.startswith(">"):
                        sequence = "\n".join(sequence.split("\n")[1:])
                    return sequence.replace("\n", "").upper()
                logger.debug(f"Failed to get {seq_type} for {identifier}: HTTP {response.status}")
                return None
        except Exception as e:
            logger.debug(f"Error fetching {seq_type} sequence for {identifier}: {e}")
            return None

    async def lookup_gene(
        self,
        query: str,
        expand: bool = False
    ) -> Optional[dict]:
        """
        Look up gene information from Ensembl.

        Args:
            query: Gene ID, gene name, or transcript ID
            expand: Whether to expand transcript information

        Returns:
            Gene data dictionary or None if not found
        """
        headers = {"Content-Type": "application/json"}

        # Try different lookup endpoints
        lookup_urls = [
            f"{self.base_url}/lookup/id/{query}?species={self.species}",
            f"{self.base_url}/lookup/symbol/{self.species}/{query}",
        ]

        if expand:
            lookup_urls = [url + "&expand=1" for url in lookup_urls]

        async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session:
            for url in lookup_urls:
                try:
                    async with session.get(url, headers=headers) as response:
                        if response.status == 200:
                            return await response.json()
                except Exception as e:
                    logger.debug(f"Failed lookup at {url}: {e}")
                    continue

        return None


class SequenceUtils:
    """Utility functions for sequence analysis."""

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculate GC content of a sequence."""
        if not sequence:
            return 0.0

        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100.0

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence.upper()))

    @staticmethod
    def transcribe_dna_to_rna(sequence: str) -> str:
        """Convert DNA sequence to RNA (T -> U)."""
        return sequence.upper().replace('T', 'U')

    @staticmethod
    def reverse_transcribe_rna_to_dna(sequence: str) -> str:
        """Convert RNA sequence to DNA (U -> T)."""
        return sequence.upper().replace('U', 'T')


class FastaUtils:
    """Utility functions for FASTA file operations."""

    @staticmethod
    def save_sequences_fasta(
        sequences: list[tuple[str, str]],
        output_path: Union[str, Path],
        line_length: int = 80
    ) -> None:
        """
        Save sequences to FASTA format.

        Args:
            sequences: List of (header, sequence) tuples
            output_path: Output file path
            line_length: Maximum line length for sequence
        """
        output_path = Path(output_path)

        with output_path.open("w") as f:
            for header, sequence in sequences:
                # Ensure header starts with >
                output_header = header if header.startswith(">") else ">" + header

                f.write(output_header + "\n")

                # Write sequence with line wrapping
                for i in range(0, len(sequence), line_length):
                    f.write(sequence[i:i + line_length] + "\n")

        logger.info(f"Saved {len(sequences)} sequences to {output_path}")

    @staticmethod
    def read_fasta(file_path: Union[str, Path]) -> list[tuple[str, str]]:
        """
        Read sequences from FASTA file.

        Args:
            file_path: Path to FASTA file

        Returns:
            List of (header, sequence) tuples
        """
        sequences = []
        current_header = None
        current_sequence = []

        with Path(file_path).open() as f:
            for file_line in f:
                line = file_line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_sequence)))
                    current_header = line[1:]  # Remove >
                    current_sequence = []
                elif current_header is not None:
                    current_sequence.append(line.upper())

            # Add last sequence
            if current_header is not None:
                sequences.append((current_header, "".join(current_sequence)))

        return sequences


def get_database_display_name(database: DatabaseType) -> str:
    """Get display name for database, handling both enum and string values."""
    if hasattr(database, 'value'):
        return database.value
    return str(database)


def get_sequence_type_display_name(seq_type: SequenceType) -> str:
    """Get display name for sequence type, handling both enum and string values."""
    if hasattr(seq_type, 'value'):
        return seq_type.value
    return str(seq_type)
