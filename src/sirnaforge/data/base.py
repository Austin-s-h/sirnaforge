"""Shared base classes and utilities for genomic data analysis."""

import asyncio
import re
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path
from typing import Callable, Optional, TypeVar, Union, cast

import aiohttp
from pydantic import BaseModel, ConfigDict, field_validator

from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class DatabaseError(Exception):
    """Base exception for database-related errors."""

    def __init__(self, message: str, database: Optional[str] = None):
        super().__init__(message)
        self.database = database


class DatabaseAccessError(DatabaseError):
    """Exception for network/access issues (firewall, timeout, server down)."""

    pass


class GeneNotFoundError(DatabaseError):
    """Exception for when a gene is not found in the database."""

    def __init__(self, query: str, database: Optional[str] = None):
        super().__init__(f"Gene '{query}' not found", database)
        self.query = query


# mypy-friendly typed alias for pydantic field_validator
F = TypeVar("F", bound=Callable[..., object])
FieldValidatorFactory = Callable[..., Callable[[F], F]]
field_validator_typed: FieldValidatorFactory = field_validator


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

    @field_validator_typed("sequence")
    @classmethod
    def validate_sequence(cls, v: Optional[str]) -> Optional[str]:
        """Validate RNA sequence."""
        if v is not None:
            # Convert to uppercase and check for valid RNA bases
            v = v.upper()
            if not re.match(r"^[ACGTU]*$", v):
                raise ValueError("Sequence contains invalid RNA bases")
        return v


class AbstractDatabaseClient(ABC):
    """Abstract base class for database clients."""

    def __init__(self, timeout: int = 30):
        """Initialize database client."""
        self.timeout = timeout

    @abstractmethod
    async def search_gene(
        self, query: str, include_sequence: bool = True
    ) -> tuple[Optional[GeneInfo], list[TranscriptInfo]]:
        """
        Search for a gene and return gene info and transcripts.

        Args:
            query: Gene ID, gene name, or transcript ID
            include_sequence: Whether to fetch transcript sequences

        Returns:
            Tuple of (gene_info, transcripts)

        Raises:
            DatabaseAccessError: For network/server access issues
            GeneNotFoundError: When gene is not found in database
        """
        pass

    @abstractmethod
    async def get_sequence(self, identifier: str, sequence_type: SequenceType = SequenceType.CDNA) -> str:
        """
        Get sequence for a specific identifier.

        Args:
            identifier: Gene ID, transcript ID, etc.
            sequence_type: Type of sequence to retrieve

        Returns:
            Sequence string

        Raises:
            DatabaseAccessError: For network/server access issues
            GeneNotFoundError: When identifier is not found in database
        """
        pass

    @property
    @abstractmethod
    def database_type(self) -> DatabaseType:
        """Return the database type this client handles."""
        pass


class EnsemblClient(AbstractDatabaseClient):
    """Client for Ensembl REST API interactions."""

    def __init__(self, timeout: int = 30, base_url: str = "https://rest.ensembl.org"):
        """Initialize Ensembl client."""
        super().__init__(timeout)
        self.base_url = base_url
        self.species = "homo_sapiens"

    @property
    def database_type(self) -> DatabaseType:
        """Return the database type this client handles."""
        return DatabaseType.ENSEMBL

    async def search_gene(
        self, query: str, include_sequence: bool = True
    ) -> tuple[Optional[GeneInfo], list[TranscriptInfo]]:
        """Search for a gene and return gene info and transcripts."""
        # First, try to resolve the query to a gene
        gene_info = await self._lookup_gene(query)

        # Get all transcripts for the gene
        transcripts = await self._get_transcripts(gene_info.gene_id, include_sequence)

        return gene_info, transcripts

    async def get_sequence(
        self, identifier: str, sequence_type: SequenceType = SequenceType.CDNA, headers: Optional[dict] = None
    ) -> str:
        """
        Get sequence from Ensembl REST API.

        Args:
            identifier: Gene ID, transcript ID, etc.
            sequence_type: Type of sequence to retrieve
            headers: Optional HTTP headers

        Returns:
            Sequence string

        Raises:
            DatabaseAccessError: For network/server access issues
            GeneNotFoundError: When identifier is not found in database
        """
        # Map sequence type to Ensembl API parameter
        type_mapping = {
            SequenceType.CDNA: "cdna",
            SequenceType.CDS: "cds",
            SequenceType.PROTEIN: "protein",
            SequenceType.GENOMIC: "genomic",
        }

        seq_type = type_mapping.get(sequence_type, "cdna")
        url = f"{self.base_url}/sequence/id/{identifier}?species={self.species}&type={seq_type}"

        if headers is None:
            headers = {"Content-Type": "text/plain"}

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, headers=headers) as response,
            ):
                if response.status == 200:
                    sequence_text: str = str(await response.text())
                    # Remove FASTA header if present
                    if sequence_text.startswith(">"):
                        sequence_text = "\n".join(sequence_text.split("\n")[1:])
                    return sequence_text.replace("\n", "").upper()
                if response.status == 404:
                    raise GeneNotFoundError(identifier, "Ensembl")
                if response.status in (403, 502, 503, 504):
                    # Server errors or access denied - likely firewall/access issue
                    raise DatabaseAccessError(f"HTTP {response.status}: Access denied or server unavailable", "Ensembl")
                logger.debug(f"Failed to get {seq_type} for {identifier}: HTTP {response.status}")
                raise DatabaseAccessError(f"HTTP {response.status}", "Ensembl")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "Ensembl") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "Ensembl") from e
        except (DatabaseAccessError, GeneNotFoundError):
            # Re-raise our custom exceptions
            raise
        except Exception as e:
            logger.debug(f"Error fetching {seq_type} sequence for {identifier}: {e}")
            raise DatabaseAccessError(f"Unexpected error: {e}", "Ensembl") from e

    async def _lookup_gene(self, query: str) -> GeneInfo:
        """Look up gene information from Ensembl."""
        data = await self._lookup_gene_data(query)

        return GeneInfo(
            gene_id=data.get("id", query),
            gene_name=data.get("display_name"),
            gene_type=data.get("biotype"),
            chromosome=data.get("seq_region_name"),
            start=data.get("start"),
            end=data.get("end"),
            strand=data.get("strand"),
            description=data.get("description"),
            database=DatabaseType.ENSEMBL,
        )

    async def _lookup_gene_data(self, query: str, expand: bool = False) -> dict:
        """
        Look up gene information from Ensembl.

        Args:
            query: Gene ID, gene name, or transcript ID
            expand: Whether to expand transcript information

        Returns:
            Gene data dictionary

        Raises:
            DatabaseAccessError: For network/server access issues
            GeneNotFoundError: When gene is not found in database
        """
        headers = {"Content-Type": "application/json"}

        # Try different lookup endpoints
        lookup_urls = [
            f"{self.base_url}/lookup/id/{query}?species={self.species}",
            f"{self.base_url}/lookup/symbol/{self.species}/{query}",
        ]

        if expand:
            lookup_urls = [url + "&expand=1" for url in lookup_urls]

        last_error = None

        async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session:
            for url in lookup_urls:
                try:
                    async with session.get(url, headers=headers) as response:
                        if response.status == 200:
                            return cast(dict, await response.json())
                        if response.status == 404:
                            # Continue to next URL, gene might be found there
                            continue
                        if response.status in (403, 502, 503, 504):
                            # Access denied or server error - raise immediately
                            raise DatabaseAccessError(
                                f"HTTP {response.status}: Access denied or server unavailable", "Ensembl"
                            )
                except aiohttp.ClientConnectorError as e:
                    last_error = DatabaseAccessError(f"Connection failed: {e}", "Ensembl")
                except asyncio.TimeoutError as e:
                    last_error = DatabaseAccessError(f"Request timeout: {e}", "Ensembl")
                except DatabaseAccessError:
                    # Re-raise access errors immediately
                    raise
                except Exception as e:
                    logger.debug(f"Failed lookup at {url}: {e}")
                    last_error = DatabaseAccessError(f"Unexpected error: {e}", "Ensembl")

        # If we had access errors, raise them
        if last_error:
            raise last_error

        # If no results from any URL, gene not found
        raise GeneNotFoundError(query, "Ensembl")

    async def _get_transcripts(self, gene_id: str, include_sequence: bool) -> list[TranscriptInfo]:
        """Get all transcripts for a gene from Ensembl."""
        transcripts: list[TranscriptInfo] = []

        try:
            # Get transcript list with expansion
            data = await self._lookup_gene_data(gene_id, expand=True)

            transcript_data = data.get("Transcript", [])

            for transcript in transcript_data:
                transcript_id = transcript.get("id")
                if not transcript_id:
                    continue

                sequence = None
                if include_sequence:
                    try:
                        sequence = await self.get_sequence(transcript_id)
                    except (DatabaseAccessError, GeneNotFoundError):
                        # If we can't get the sequence, log and continue without it
                        logger.warning(f"Could not retrieve sequence for transcript {transcript_id}")
                        sequence = None

                transcripts.append(
                    TranscriptInfo(
                        transcript_id=transcript_id,
                        transcript_name=transcript.get("display_name"),
                        transcript_type=transcript.get("biotype"),
                        gene_id=gene_id,
                        gene_name=data.get("display_name"),
                        sequence=sequence,
                        length=len(sequence) if sequence else None,
                        database=DatabaseType.ENSEMBL,
                        is_canonical=transcript.get("is_canonical", False),
                    )
                )

        except (DatabaseAccessError, GeneNotFoundError):
            # Propagate access and not-found errors
            raise
        except Exception as e:
            logger.error(f"Failed to get transcripts for {gene_id}: {e}")
            raise DatabaseAccessError(f"Failed to get transcripts: {e}", "Ensembl") from e

        return transcripts


class RefSeqClient(AbstractDatabaseClient):
    """Client for RefSeq database via NCBI E-utilities API."""

    def __init__(self, timeout: int = 30, base_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"):
        """Initialize RefSeq client."""
        super().__init__(timeout)
        self.base_url = base_url
        self.email = "sirnaforge@example.com"  # Required by NCBI
        self.tool = "sirnaforge"

    @property
    def database_type(self) -> DatabaseType:
        """Return the database type this client handles."""
        return DatabaseType.REFSEQ

    async def search_gene(
        self, query: str, include_sequence: bool = True
    ) -> tuple[Optional[GeneInfo], list[TranscriptInfo]]:
        """Search for a gene and return gene info and transcripts."""
        # Search for gene in NCBI Gene database
        gene_id = await self._search_gene_id(query)

        # Get gene information
        gene_info = await self._get_gene_info(gene_id, query)

        # Get transcripts for this gene
        transcripts = await self._get_transcripts(gene_id, gene_info, include_sequence)

        return gene_info, transcripts

    async def get_sequence(self, identifier: str, _sequence_type: SequenceType = SequenceType.CDNA) -> str:
        """Get sequence for a specific identifier from NCBI."""
        url = f"{self.base_url}/efetch.fcgi"
        params = {
            "db": "nucleotide",
            "id": identifier,
            "rettype": "fasta",
            "retmode": "text",
            "email": self.email,
            "tool": self.tool,
        }

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, params=params) as response,
            ):
                if response.status == 200:
                    fasta_text = await response.text()
                    # Remove FASTA header and extract sequence
                    lines = fasta_text.strip().split("\n")
                    if lines[0].startswith(">"):
                        sequence = "".join(lines[1:])
                        return sequence.replace("\n", "").upper()
                    raise GeneNotFoundError(identifier, "RefSeq")
                if response.status == 404:
                    raise GeneNotFoundError(identifier, "RefSeq")
                raise DatabaseAccessError(f"HTTP {response.status}", "RefSeq")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "RefSeq") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "RefSeq") from e
        except (DatabaseAccessError, GeneNotFoundError):
            raise
        except Exception as e:
            raise DatabaseAccessError(f"Unexpected error: {e}", "RefSeq") from e

    async def _search_gene_id(self, query: str) -> str:
        """Search for gene ID using NCBI esearch."""
        url = f"{self.base_url}/esearch.fcgi"
        params = {
            "db": "gene",
            "term": f"{query}[Gene Name] AND Homo sapiens[Organism]",
            "retmode": "json",
            "email": self.email,
            "tool": self.tool,
        }

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, params=params) as response,
            ):
                if response.status == 200:
                    data = await response.json()
                    id_list = data.get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        return str(id_list[0])
                    raise GeneNotFoundError(query, "RefSeq")
                raise DatabaseAccessError(f"HTTP {response.status}", "RefSeq")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "RefSeq") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "RefSeq") from e
        except (DatabaseAccessError, GeneNotFoundError):
            raise
        except Exception as e:
            raise DatabaseAccessError(f"Unexpected error: {e}", "RefSeq") from e

    async def _get_gene_info(self, gene_id: str, original_query: str) -> GeneInfo:
        """Get detailed gene information from NCBI."""
        url = f"{self.base_url}/esummary.fcgi"
        params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "json",
            "email": self.email,
            "tool": self.tool,
        }

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, params=params) as response,
            ):
                if response.status == 200:
                    data = await response.json()
                    gene_data = data.get("result", {}).get(gene_id, {})

                    return GeneInfo(
                        gene_id=gene_id,
                        gene_name=gene_data.get("name", original_query),
                        gene_type=gene_data.get("genetype", "unknown"),
                        chromosome=gene_data.get("chromosome", None),
                        start=(
                            gene_data.get("genomicinfo", [{}])[0].get("chrstart")
                            if gene_data.get("genomicinfo")
                            else None
                        ),
                        end=(
                            gene_data.get("genomicinfo", [{}])[0].get("chrstop")
                            if gene_data.get("genomicinfo")
                            else None
                        ),
                        strand=None,  # Not readily available in summary
                        description=gene_data.get("summary", ""),
                        database=DatabaseType.REFSEQ,
                    )
                raise DatabaseAccessError(f"HTTP {response.status}", "RefSeq")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "RefSeq") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "RefSeq") from e
        except (DatabaseAccessError, GeneNotFoundError):
            raise
        except Exception as e:
            raise DatabaseAccessError(f"Unexpected error: {e}", "RefSeq") from e

    async def _get_transcripts(self, gene_id: str, gene_info: GeneInfo, include_sequence: bool) -> list[TranscriptInfo]:
        """Get transcripts for a gene from RefSeq using NCBI E-utilities."""
        transcripts: list[TranscriptInfo] = []

        try:
            # Step 1: Use elink to find associated nucleotide records (mRNAs/transcripts)
            transcript_ids = await self._find_linked_transcripts(gene_id)

            if not transcript_ids:
                logger.info(f"No linked transcripts found for gene {gene_id}")
                return transcripts

            logger.info(f"Found {len(transcript_ids)} linked transcript(s) for gene {gene_id}")

            # Step 2: Get transcript metadata using esummary
            transcript_metadata = await self._get_transcript_metadata(transcript_ids)

            # Step 3: Build TranscriptInfo objects
            for transcript_id, metadata in transcript_metadata.items():
                sequence = None
                if include_sequence:
                    try:
                        sequence = await self.get_sequence(transcript_id)
                    except (DatabaseAccessError, GeneNotFoundError):
                        logger.warning(f"Could not retrieve sequence for transcript {transcript_id}")
                        sequence = None

                # Extract information from metadata
                title = metadata.get("title", "")
                accession = metadata.get("accessionversion", transcript_id)

                # Parse transcript type from title (RefSeq convention)
                transcript_type = self._parse_transcript_type(title, accession)

                # Determine if this is a canonical transcript (NM_ prefixes are typically canonical)
                is_canonical = accession.startswith("NM_")

                transcripts.append(
                    TranscriptInfo(
                        transcript_id=accession,
                        transcript_name=title.split(",")[0] if title else None,  # First part of title
                        transcript_type=transcript_type,
                        gene_id=gene_id,
                        gene_name=gene_info.gene_name,
                        sequence=sequence,
                        length=len(sequence) if sequence else metadata.get("slen"),
                        database=DatabaseType.REFSEQ,
                        is_canonical=is_canonical,
                    )
                )

            logger.info(f"Successfully processed {len(transcripts)} transcript(s) for gene {gene_id}")

        except (DatabaseAccessError, GeneNotFoundError):
            # Propagate access and not-found errors
            raise
        except Exception as e:
            logger.error(f"Failed to get transcripts for gene {gene_id}: {e}")
            raise DatabaseAccessError(f"Failed to get transcripts: {e}", "RefSeq") from e

        return transcripts

    async def _find_linked_transcripts(self, gene_id: str) -> list[str]:
        """Use elink to find nucleotide records linked to a gene."""
        url = f"{self.base_url}/elink.fcgi"
        params = {
            "dbfrom": "gene",
            "db": "nucleotide",
            "id": gene_id,
            "retmode": "json",
            "email": self.email,
            "tool": self.tool,
        }

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, params=params) as response,
            ):
                if response.status == 200:
                    data = await response.json()
                    linksets = data.get("linksets", [])

                    transcript_ids = []
                    for linkset in linksets:
                        if linkset.get("dbto") == "nucleotide":
                            for link in linkset.get("linksetdbs", []):
                                if link.get("dbto") == "nucleotide":
                                    transcript_ids.extend(link.get("links", []))

                    return transcript_ids
                raise DatabaseAccessError(f"HTTP {response.status}", "RefSeq")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "RefSeq") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "RefSeq") from e
        except (DatabaseAccessError, GeneNotFoundError):
            raise
        except Exception as e:
            raise DatabaseAccessError(f"Unexpected error: {e}", "RefSeq") from e

    async def _get_transcript_metadata(self, transcript_ids: list[str]) -> dict[str, dict]:
        """Get metadata for multiple transcripts using esummary."""
        if not transcript_ids:
            return {}

        # NCBI recommends batching requests, but limit to reasonable size
        batch_size = 200
        all_metadata = {}

        for i in range(0, len(transcript_ids), batch_size):
            batch_ids = transcript_ids[i : i + batch_size]
            batch_metadata = await self._get_transcript_metadata_batch(batch_ids)
            all_metadata.update(batch_metadata)

        return all_metadata

    async def _get_transcript_metadata_batch(self, transcript_ids: list[str]) -> dict[str, dict]:
        """Get metadata for a batch of transcripts."""
        url = f"{self.base_url}/esummary.fcgi"
        params = {
            "db": "nucleotide",
            "id": ",".join(transcript_ids),
            "retmode": "json",
            "email": self.email,
            "tool": self.tool,
        }

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, params=params) as response,
            ):
                if response.status == 200:
                    data = await response.json()
                    result = data.get("result", {})

                    # Filter out the 'uids' key which is metadata about the result
                    return {k: v for k, v in result.items() if k != "uids" and isinstance(v, dict)}
                raise DatabaseAccessError(f"HTTP {response.status}", "RefSeq")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "RefSeq") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "RefSeq") from e
        except (DatabaseAccessError, GeneNotFoundError):
            raise
        except Exception as e:
            raise DatabaseAccessError(f"Unexpected error: {e}", "RefSeq") from e

    def _parse_transcript_type(self, title: str, accession: str) -> str:
        """Parse transcript type from RefSeq title and accession."""
        # RefSeq accession prefixes indicate type
        accession_types = {
            "NM_": "protein_coding",  # mRNA
            "NR_": "non_coding",  # non-coding RNA
            "XM_": "protein_coding",  # predicted mRNA
            "XR_": "non_coding",  # predicted non-coding RNA
        }

        # Check accession prefix first
        for prefix, transcript_type in accession_types.items():
            if accession.startswith(prefix):
                return transcript_type

        # Fall back to parsing title
        title_lower = title.lower()
        if "mrna" in title_lower or "protein" in title_lower:
            return "protein_coding"
        if any(term in title_lower for term in ["ncrna", "lncrna", "lincrna", "mirna", "snrna", "snorna"]):
            return "non_coding"

        # Default fallback
        return "unknown"


class GencodeClient(AbstractDatabaseClient):
    """Client for GENCODE database."""

    def __init__(self, timeout: int = 30):
        """Initialize GENCODE client."""
        super().__init__(timeout)
        self.base_url = "https://www.gencodegenes.org"
        self.version = "44"  # GENCODE version

    @property
    def database_type(self) -> DatabaseType:
        """Return the database type this client handles."""
        return DatabaseType.GENCODE

    async def search_gene(
        self,
        query: str,
        include_sequence: bool = True,  # noqa: ARG002
    ) -> tuple[Optional[GeneInfo], list[TranscriptInfo]]:
        """Search for a gene and return gene info and transcripts."""
        # GENCODE doesn't have a simple REST API like Ensembl
        # This would typically require parsing GTF/GFF files or using their FTP download
        # For now, this is a placeholder implementation
        logger.info(f"GENCODE search for '{query}:{include_sequence}' not implemented")
        raise DatabaseAccessError("GENCODE search not yet implemented - requires GTF file parsing", "GENCODE")

    async def get_sequence(self, _identifier: str, _sequence_type: SequenceType = SequenceType.CDNA) -> str:
        """Get sequence for a specific identifier from GENCODE."""
        # GENCODE sequences are typically accessed via FASTA files
        # This would require downloading and indexing GENCODE FASTA files
        raise DatabaseAccessError("GENCODE sequence retrieval not yet implemented", "GENCODE")


class SequenceUtils:
    """Utility functions for sequence analysis."""

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculate GC content of a sequence."""
        if not sequence:
            return 0.0

        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100.0

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(complement.get(base, "N") for base in reversed(sequence.upper()))

    @staticmethod
    def transcribe_dna_to_rna(sequence: str) -> str:
        """Convert DNA sequence to RNA (T -> U)."""
        return sequence.upper().replace("T", "U")

    @staticmethod
    def reverse_transcribe_rna_to_dna(sequence: str) -> str:
        """Convert RNA sequence to DNA (U -> T)."""
        return sequence.upper().replace("U", "T")


class FastaUtils:
    """Utility functions for FASTA file operations."""

    @staticmethod
    def save_sequences_fasta(
        sequences: list[tuple[str, str]], output_path: Union[str, Path], line_length: int = 80
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
                    f.write(sequence[i : i + line_length] + "\n")

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
        sequences: list[tuple[str, str]] = []
        current_header: Optional[str] = None
        current_sequence: list[str] = []

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

    @staticmethod
    def parse_fasta_to_dict(file_path: Union[str, Path]) -> dict[str, str]:
        """
        Parse FASTA file into a dictionary.

        Args:
            file_path: Path to FASTA file

        Returns:
            Dictionary mapping sequence names to sequences
        """
        sequences_list = FastaUtils.read_fasta(file_path)

        # Convert to dictionary
        sequences_dict = {}
        for header, sequence in sequences_list:
            # Clean header (remove > if present)
            clean_header = header.lstrip(">")
            sequences_dict[clean_header] = sequence.upper().replace("U", "T")

        logger.info(f"Parsed {len(sequences_dict)} sequences from {file_path}")
        return sequences_dict

    @staticmethod
    def write_dict_to_fasta(sequences: dict[str, str], output_path: Union[str, Path]) -> None:
        """
        Write sequences dictionary to FASTA format.

        Args:
            sequences: Dictionary of sequence name -> sequence
            output_path: Output file path
        """
        # Convert to list format
        fasta_sequences = list(sequences.items())

        # Use existing save method
        FastaUtils.save_sequences_fasta(fasta_sequences, output_path)

        logger.info(f"Wrote {len(sequences)} sequences to {output_path}")

    @staticmethod
    def validate_sirna_sequences(sequences: dict[str, str], expected_length: int = 21) -> dict[str, str]:
        """
        Validate siRNA sequences for correct length and nucleotide content.

        Args:
            sequences: Dictionary of sequence name -> sequence
            expected_length: Expected siRNA length

        Returns:
            Dictionary of valid sequences
        """
        valid_sequences = {}
        invalid_count = 0

        for name, seq in sequences.items():
            # Clean sequence
            clean_seq = seq.upper().replace("U", "T")

            # Validate length and nucleotide content
            if len(clean_seq) == expected_length and all(base in "ATCG" for base in clean_seq):
                valid_sequences[name] = clean_seq
            else:
                invalid_count += 1
                logger.debug(f"Invalid sequence {name}: length={len(clean_seq)}, sequence={clean_seq}")

        logger.info(f"Validation complete: {len(valid_sequences)} valid, {invalid_count} invalid sequences")

        if len(valid_sequences) == 0:
            raise ValueError("No valid sequences found after validation")

        return valid_sequences


def get_database_display_name(database: DatabaseType) -> str:
    """Get display name for database, handling both enum and string values."""
    if hasattr(database, "value"):
        return database.value
    return str(database)


def get_sequence_type_display_name(seq_type: SequenceType) -> str:
    """Get display name for sequence type, handling both enum and string values."""
    if hasattr(seq_type, "value"):
        return seq_type.value
    return str(seq_type)
