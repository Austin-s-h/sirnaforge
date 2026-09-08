"""Shared base classes and utilities for genomic data analysis."""

from __future__ import annotations

import asyncio
import re
from abc import ABC, abstractmethod
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import aiohttp
from pydantic import BaseModel, ConfigDict

from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.typed_decorators import field_validator_typed

if TYPE_CHECKING:
    from sirnaforge.config.reference_policy import ReferenceChoice
    from sirnaforge.models.transcript_annotation import TranscriptAnnotationBundle

logger = get_logger(__name__)


class DatabaseError(Exception):
    """Base exception for database-related errors."""

    def __init__(self, message: str, database: str | None = None):
        """Initialize database error."""
        super().__init__(message)
        self.database = database


class DatabaseAccessError(DatabaseError):
    """Exception for network/access issues (firewall, timeout, server down)."""

    pass


class GeneNotFoundError(DatabaseError):
    """Exception for when a gene is not found in the database."""

    def __init__(self, query: str, database: str | None = None):
        """Initialize gene not found error."""
        super().__init__(f"Gene '{query}' not found", database)
        self.query = query


#: Ensembl REST accepts at most 50 identifiers per POST body on its id endpoints.
ENSEMBL_POST_CHUNK_SIZE = 50
#: 429 is Ensembl's per-IP rate limit and the 5xx family is a saturated or wedged backend.
#: All of them are transient -- a plain sequence fetch draws 503s under load -- so retry them.
ENSEMBL_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
ENSEMBL_MAX_ATTEMPTS = 3


@asynccontextmanager
async def ensembl_session(session: aiohttp.ClientSession | None, timeout: int) -> AsyncIterator[aiohttp.ClientSession]:
    """Yield ``session`` when the caller owns one, else open one for this operation.

    Retrieving a gene means several requests to one host, so they share a connection
    instead of paying a TLS handshake each.
    """
    if session is not None:
        yield session
        return

    async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(timeout), trust_env=True) as owned:
        yield owned


def _ensembl_retry_delay(response: aiohttp.ClientResponse, attempt: int, cap: int) -> float:
    """Seconds to wait before retrying, honouring Retry-After when Ensembl sends it."""
    retry_after = response.headers.get("Retry-After")
    if retry_after:
        try:
            return min(float(retry_after), float(cap))
        except ValueError:
            pass
    return float(2**attempt)


async def ensembl_request_json(
    session: aiohttp.ClientSession,
    method: str,
    url: str,
    *,
    headers: dict[str, str],
    json_body: dict[str, Any] | None = None,
    retry_cap: int = 30,
) -> Any:
    """Issue a JSON request to Ensembl, retrying its transient statuses.

    Args:
        session: Session to issue the request on
        method: HTTP method ("GET" or "POST")
        url: Fully-formed request URL
        headers: Request headers
        json_body: JSON body for POST requests
        retry_cap: Upper bound in seconds on a Retry-After the server asks for

    Returns:
        Parsed JSON payload, or None when Ensembl answers 404.

    Raises:
        DatabaseAccessError: For network errors, and for 429/5xx that outlive the retries
    """
    last_status: int | None = None

    for attempt in range(ENSEMBL_MAX_ATTEMPTS):
        try:
            async with session.request(method, url, headers=headers, json=json_body) as response:
                if response.status == 200:
                    return await response.json()
                if response.status == 404:
                    return None

                last_status = response.status
                if response.status in ENSEMBL_RETRY_STATUSES and attempt < ENSEMBL_MAX_ATTEMPTS - 1:
                    delay = _ensembl_retry_delay(response, attempt, retry_cap)
                    logger.debug(f"Ensembl returned HTTP {response.status} for {url}; retrying in {delay:.0f}s")
                    await asyncio.sleep(delay)
                    continue

                raise DatabaseAccessError(f"HTTP {response.status}", "Ensembl")
        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "Ensembl") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "Ensembl") from e

    raise DatabaseAccessError(f"HTTP {last_status} after {ENSEMBL_MAX_ATTEMPTS} attempts", "Ensembl")


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
    gene_name: str | None = None
    gene_type: str | None = None
    chromosome: str | None = None
    start: int | None = None
    end: int | None = None
    strand: int | None = None
    description: str | None = None
    database: DatabaseType

    model_config = ConfigDict(use_enum_values=True)


class TranscriptInfo(BaseModel):
    """Transcript information model."""

    transcript_id: str
    transcript_name: str | None = None
    transcript_type: str | None = None
    gene_id: str
    gene_name: str | None = None
    sequence: str | None = None
    length: int | None = None
    database: DatabaseType
    is_canonical: bool = False

    model_config = ConfigDict(use_enum_values=True)

    @field_validator_typed("sequence")
    @classmethod
    def validate_sequence(cls, v: str | None) -> str | None:
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
    ) -> tuple[GeneInfo | None, list[TranscriptInfo]]:
        """Search for a gene and return gene info and transcripts.

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
        """Get sequence for a specific identifier.

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

    @property
    def query_species(self) -> str:
        """Canonical organism whose transcripts a gene query against this database returns.

        Every concrete client here is wired to exactly one organism -- ``EnsemblClient`` pins
        ``species=homo_sapiens`` on every lookup, ``RefSeqClient`` appends
        ``AND Homo sapiens[Organism]`` to every esearch term, and ``GencodeClient`` targets the
        human GENCODE release -- so the organism of a retrieved transcript set is a property of
        *where the transcripts came from*. Callers that need the "query species" must read it from
        here rather than infer it from the off-target species list, which is an unordered set of
        genomes to screen against and says nothing about the target.
        """
        return "human"


class AbstractTranscriptAnnotationClient(ABC):
    """Abstract base class for transcript annotation clients.

    **Purpose and Scope:**
    Provides genomic annotation metadata (exon/CDS structure, coordinates, biotype)
    WITHOUT fetching full transcript sequences. This is complementary to, not overlapping
    with, AbstractDatabaseClient which focuses on sequence retrieval.

    **Key Differences from GeneSearcher/AbstractDatabaseClient:**

    1. **Focus**: Structural annotations (exons, CDS intervals, genomic coordinates)
       vs. sequence data (cDNA, CDS, protein sequences)

    2. **Use Case**: Enriching existing transcript metadata with genomic context
       vs. discovering and retrieving transcripts with sequences

    3. **Query Patterns**:
       - By stable IDs: fetch_by_ids(['ENST00000269305'])
       - By genomic regions: fetch_by_regions(['17:7661779-7687550'])
       vs. GeneSearcher which queries by gene name/symbol

    4. **Caching Strategy**: In-memory LRU cache with TTL for transient annotation data
       vs. ReferenceManager's persistent file cache for large sequence datasets

    **When to Use:**
    - Need exon/CDS boundaries for visualization or analysis
    - Need genomic coordinates for variant mapping
    - Need biotype information without full sequence download
    - Need to query multiple transcripts in a genomic region

    **When to Use GeneSearcher Instead:**
    - Need transcript sequences for siRNA design
    - Need to discover transcripts by gene name/symbol
    - Need protein sequences or translations
    """

    def __init__(self, timeout: int = 30):
        """Initialize transcript annotation client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout

    @abstractmethod
    async def fetch_by_ids(
        self, ids: list[str], *, species: str, reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        """Fetch transcript annotations by stable IDs.

        Args:
            ids: List of transcript or gene IDs (e.g., ENST00000269305, TP53)
            species: Species name (e.g., 'homo_sapiens', 'human')
            reference: Reference assembly/release choice

        Returns:
            TranscriptAnnotationBundle containing resolved annotations

        Raises:
            DatabaseAccessError: For network/server access issues
        """
        pass

    @abstractmethod
    async def fetch_by_regions(
        self, regions: list[str], *, species: str, reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        """Fetch transcript annotations by genomic regions.

        Args:
            regions: List of genomic regions in format 'chr:start-end' (e.g., '17:7661779-7687550')
            species: Species name (e.g., 'homo_sapiens', 'human')
            reference: Reference assembly/release choice

        Returns:
            TranscriptAnnotationBundle containing all transcripts overlapping regions

        Raises:
            DatabaseAccessError: For network/server access issues
        """
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

    @property
    def query_species(self) -> str:
        """Canonical form of the single Ensembl species this client queries."""
        return normalize_species_name(self.species)

    async def search_gene(
        self, query: str, include_sequence: bool = True
    ) -> tuple[GeneInfo | None, list[TranscriptInfo]]:
        """Search for a gene and return gene info and transcripts."""
        # First, try to resolve the query to a gene
        gene_info = await self._lookup_gene(query)

        # Get all transcripts for the gene
        transcripts = await self._get_transcripts(gene_info.gene_id, include_sequence)

        return gene_info, transcripts

    async def get_sequences(
        self,
        identifiers: list[str],
        sequence_type: SequenceType = SequenceType.CDNA,
        session: aiohttp.ClientSession | None = None,
    ) -> dict[str, str]:
        """Get sequences for many identifiers, one POST per chunk of identifiers.

        ``POST /sequence/id`` serves exactly what the per-id GET serves, but in one round
        trip: all 40 TP53 transcripts come back in ~13s, where fetching them one at a time
        takes ~350s and loses several to transient 503s.

        Args:
            identifiers: Identifiers to fetch
            sequence_type: Type of sequence to retrieve
            session: Session to reuse; one is opened for this call when omitted

        Returns:
            Mapping of identifier to sequence. Identifiers Ensembl does not know are absent
            from the mapping rather than raising -- the POST endpoint simply omits them.

        Raises:
            DatabaseAccessError: For network/server access issues
        """
        if not identifiers:
            return {}

        seq_type = self._sequence_type_param(sequence_type)
        url = f"{self.base_url}/sequence/id?species={self.species}&type={seq_type}"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        sequences: dict[str, str] = {}

        async with ensembl_session(session, self.timeout) as active:
            for start in range(0, len(identifiers), ENSEMBL_POST_CHUNK_SIZE):
                chunk = identifiers[start : start + ENSEMBL_POST_CHUNK_SIZE]
                payload = await ensembl_request_json(
                    active, "POST", url, headers=headers, json_body={"ids": chunk}, retry_cap=self.timeout
                )

                for record in payload if isinstance(payload, list) else []:
                    identifier = record.get("query") or record.get("id")
                    sequence = record.get("seq")
                    if identifier and sequence:
                        sequences[identifier] = str(sequence).replace("\n", "").upper()

        return sequences

    @staticmethod
    def _sequence_type_param(sequence_type: SequenceType) -> str:
        """Map a sequence type onto Ensembl's ``type`` query parameter."""
        type_mapping = {
            SequenceType.CDNA: "cdna",
            SequenceType.CDS: "cds",
            SequenceType.PROTEIN: "protein",
            SequenceType.GENOMIC: "genomic",
        }
        return type_mapping.get(sequence_type, "cdna")

    async def get_sequence(
        self,
        identifier: str,
        sequence_type: SequenceType = SequenceType.CDNA,
        headers: dict | None = None,
        session: aiohttp.ClientSession | None = None,
    ) -> str:
        """Get sequence from Ensembl REST API.

        Args:
            identifier: Gene ID, transcript ID, etc.
            sequence_type: Type of sequence to retrieve
            headers: Optional HTTP headers
            session: Session to reuse; one is opened for this call when omitted

        Returns:
            Sequence string

        Raises:
            DatabaseAccessError: For network/server access issues
            GeneNotFoundError: When identifier is not found in database
        """
        seq_type = self._sequence_type_param(sequence_type)
        url = f"{self.base_url}/sequence/id/{identifier}?species={self.species}&type={seq_type}"

        if headers is None:
            headers = {"Content-Type": "text/plain"}

        try:
            async with (
                ensembl_session(session, self.timeout) as active,
                active.get(url, headers=headers) as response,
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
        """Look up gene information from Ensembl.

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

        async with ensembl_session(None, self.timeout) as session:
            for url in lookup_urls:
                try:
                    # None means 404 here: continue to the next URL, the gene may be found there.
                    payload = await ensembl_request_json(session, "GET", url, headers=headers, retry_cap=self.timeout)
                    if payload is not None:
                        return cast(dict, payload)
                except DatabaseAccessError as e:
                    # Keep the error but let the remaining endpoint try; a transient failure on
                    # one lookup form should not hide a working answer from the other.
                    last_error = e
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

            sequences: dict[str, str] = {}
            if include_sequence:
                transcript_ids = [t["id"] for t in transcript_data if t.get("id")]
                sequences = await self._get_sequences_with_fallback(transcript_ids)

            for transcript in transcript_data:
                transcript_id = transcript.get("id")
                if not transcript_id:
                    continue

                sequence = sequences.get(transcript_id)

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

    async def _get_sequences_with_fallback(self, identifiers: list[str]) -> dict[str, str]:
        """Batch-fetch sequences, then retry per id for whatever the batch did not return.

        The batch is the fast path; the per-id retry covers ids the POST endpoint omits and
        the case where the batch request itself fails. An id that fails both ways is dropped
        with a warning, as before -- the transcript is still reported, without a sequence.
        """
        if not identifiers:
            return {}

        async with ensembl_session(None, self.timeout) as session:
            try:
                sequences = await self.get_sequences(identifiers, session=session)
            except DatabaseAccessError as e:
                logger.warning(f"Batch sequence request failed ({e}); falling back to per-transcript requests")
                sequences = {}

            for identifier in [i for i in identifiers if i not in sequences]:
                try:
                    sequences[identifier] = await self.get_sequence(identifier, session=session)
                except (DatabaseAccessError, GeneNotFoundError):
                    logger.warning(f"Could not retrieve sequence for transcript {identifier}")

        return sequences


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
    ) -> tuple[GeneInfo | None, list[TranscriptInfo]]:
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
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout), trust_env=True) as session,
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
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout), trust_env=True) as session,
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
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout), trust_env=True) as session,
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
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout), trust_env=True) as session,
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
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout), trust_env=True) as session,
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
    ) -> tuple[GeneInfo | None, list[TranscriptInfo]]:
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
    def save_sequences_fasta(sequences: list[tuple[str, str]], output_path: str | Path, line_length: int = 80) -> None:
        """Save sequences to FASTA format.

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
    def read_fasta(file_path: str | Path) -> list[tuple[str, str]]:
        """Read sequences from FASTA file.

        Args:
            file_path: Path to FASTA file

        Returns:
            List of (header, sequence) tuples
        """
        sequences: list[tuple[str, str]] = []
        current_header: str | None = None
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
    def parse_fasta_to_dict(file_path: str | Path) -> dict[str, str]:
        """Parse FASTA file into a dictionary.

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
    def write_dict_to_fasta(sequences: dict[str, str], output_path: str | Path) -> None:
        """Write sequences dictionary to FASTA format.

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
        """Validate siRNA sequences for correct length and nucleotide content.

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
