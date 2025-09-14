"""Gene search and sequence retrieval from multiple databases."""

import asyncio
from pathlib import Path
from typing import Optional, Union

from pydantic import BaseModel, ConfigDict, Field

from sirnaforge.data.base import (
    AbstractDatabaseClient,
    DatabaseAccessError,
    DatabaseType,
    EnsemblClient,
    FastaUtils,
    GencodeClient,
    GeneInfo,
    GeneNotFoundError,
    RefSeqClient,
    TranscriptInfo,
    get_database_display_name,
)
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class GeneSearchResult(BaseModel):
    """Complete gene search result."""

    query: str
    database: DatabaseType
    gene_info: Optional[GeneInfo] = None
    transcripts: list[TranscriptInfo] = Field(default_factory=list)
    error: Optional[str] = None
    is_access_error: bool = False  # True if error was due to access/network issues

    model_config = ConfigDict(use_enum_values=True)

    @property
    def success(self) -> bool:
        """Check if search was successful."""
        return self.gene_info is not None and len(self.transcripts) > 0


class GeneSearcher:
    """Search genes and retrieve sequences from genomic databases using multiple clients."""

    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """
        Initialize gene searcher with database clients.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum retry attempts
        """
        self.timeout = timeout
        self.max_retries = max_retries

        # Initialize database clients
        self.clients: dict[DatabaseType, AbstractDatabaseClient] = {
            DatabaseType.ENSEMBL: EnsemblClient(timeout=timeout),
            DatabaseType.REFSEQ: RefSeqClient(timeout=timeout),
            DatabaseType.GENCODE: GencodeClient(timeout=timeout),
        }

    def get_client(self, database: DatabaseType) -> AbstractDatabaseClient:
        """Get the client for a specific database."""
        return self.clients[database]

    async def search_gene_with_fallback(self, query: str, include_sequence: bool = True) -> GeneSearchResult:
        """
        Search for a gene with automatic fallback to other databases.

        Tries databases in order: Ensembl -> RefSeq -> GENCODE
        Falls back to next database only if access is blocked (not if gene is not found).

        Args:
            query: Gene ID, gene name, or transcript ID
            include_sequence: Whether to fetch transcript sequences

        Returns:
            GeneSearchResult from the first accessible database
        """
        databases_to_try = [DatabaseType.ENSEMBL, DatabaseType.REFSEQ, DatabaseType.GENCODE]

        for db in databases_to_try:
            logger.info(f"Searching for '{query}' in {get_database_display_name(db)}")

            try:
                result = await self.search_gene(query, db, include_sequence)
                if result.success:
                    logger.info(f"Successfully found '{query}' in {get_database_display_name(db)}")
                    return result
                if result.is_access_error:
                    # Access error - try next database
                    logger.warning(f"Access blocked to {get_database_display_name(db)}: {result.error}")
                    continue
                if result.error and "not yet implemented" in result.error:
                    # Implementation missing - try next database
                    logger.info(f"{get_database_display_name(db)} search not implemented, trying next database")
                    continue
                # Gene not found in this database, but database is accessible
                # Continue to next database to see if gene exists there
                logger.info(f"Gene '{query}' not found in {get_database_display_name(db)}, trying next database")
                continue

            except DatabaseAccessError as e:
                logger.warning(f"Access blocked to {get_database_display_name(db)}: {e}")
                # Try next database
                continue
            except Exception as e:
                logger.error(f"Unexpected error searching {get_database_display_name(db)}: {e}")
                # Try next database
                continue

        # If we get here, all databases failed or gene not found anywhere
        return GeneSearchResult(
            query=query,
            database=DatabaseType.ENSEMBL,  # Default to first attempted
            error=f"Gene '{query}' not found in any accessible database",
        )

    async def search_gene(
        self, query: str, database: Optional[DatabaseType] = None, include_sequence: bool = True
    ) -> GeneSearchResult:
        """
        Search for a gene and retrieve its isoforms.

        Args:
            query: Gene ID, gene name, or transcript ID
            database: Database to search (defaults to Ensembl)
            include_sequence: Whether to fetch transcript sequences

        Returns:
            GeneSearchResult with gene info and transcripts
        """
        db = database or DatabaseType.ENSEMBL

        logger.info(f"Searching for '{query}' in {get_database_display_name(db)}")

        try:
            client = self.get_client(db)
            gene_info, transcripts = await client.search_gene(query, include_sequence)

            return GeneSearchResult(query=query, database=db, gene_info=gene_info, transcripts=transcripts)

        except DatabaseAccessError as e:
            logger.warning(f"Access error for '{query}' in {get_database_display_name(db)}: {e}")
            return GeneSearchResult(query=query, database=db, error=str(e), is_access_error=True)
        except GeneNotFoundError as e:
            logger.info(f"Gene '{query}' not found in {get_database_display_name(db)}")
            return GeneSearchResult(query=query, database=db, error=str(e), is_access_error=False)
        except Exception as e:
            logger.error(f"Search failed for '{query}' in {get_database_display_name(db)}: {e}")
            return GeneSearchResult(query=query, database=db, error=str(e), is_access_error=False)

    async def search_multiple_databases(
        self, query: str, databases: Optional[list[DatabaseType]] = None, include_sequence: bool = True
    ) -> list[GeneSearchResult]:
        """
        Search across multiple databases.

        Args:
            query: Gene ID, gene name, or transcript ID
            databases: List of databases to search
            include_sequence: Whether to fetch sequences

        Returns:
            List of search results from each database
        """
        if databases is None:
            databases = list(DatabaseType)

        tasks = [self.search_gene(query, db, include_sequence) for db in databases]

        results = await asyncio.gather(*tasks, return_exceptions=True)

        # Convert exceptions to error results
        processed_results = []
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                processed_results.append(GeneSearchResult(query=query, database=databases[i], error=str(result)))
            elif isinstance(result, GeneSearchResult):
                processed_results.append(result)
            else:
                # Handle unexpected result type
                processed_results.append(
                    GeneSearchResult(query=query, database=databases[i], error="Unexpected result type")
                )

        return processed_results

    def save_transcripts_fasta(
        self, transcripts: list[TranscriptInfo], output_path: Union[str, Path], include_metadata: bool = True
    ) -> None:
        """
        Save transcripts to FASTA format using shared utility.

        Args:
            transcripts: List of transcript information
            output_path: Output file path
            include_metadata: Include metadata in FASTA headers
        """
        # Convert transcripts to (header, sequence) tuples
        fasta_sequences = []

        for transcript in transcripts:
            if not transcript.sequence:
                logger.warning(f"No sequence for transcript {transcript.transcript_id}")
                continue

            # Create FASTA header
            header = transcript.transcript_id
            if include_metadata:
                metadata = []
                if transcript.gene_name:
                    metadata.append(f"gene_name:{transcript.gene_name}")
                if transcript.transcript_type:
                    metadata.append(f"type:{transcript.transcript_type}")
                if transcript.is_canonical:
                    metadata.append("canonical:true")
                if transcript.length:
                    metadata.append(f"length:{transcript.length}")

                if metadata:
                    header += " " + " ".join(metadata)

            fasta_sequences.append((header, transcript.sequence))

        # Use shared FastaUtils to save
        FastaUtils.save_sequences_fasta(fasta_sequences, output_path)


# Convenience functions for synchronous usage
def search_gene_sync(
    query: str, database: DatabaseType = DatabaseType.ENSEMBL, include_sequence: bool = True
) -> GeneSearchResult:
    """Synchronous wrapper for gene search."""
    searcher = GeneSearcher()
    return asyncio.run(searcher.search_gene(query, database, include_sequence))


def search_gene_with_fallback_sync(query: str, include_sequence: bool = True) -> GeneSearchResult:
    """Synchronous wrapper for gene search with fallback."""
    searcher = GeneSearcher()
    return asyncio.run(searcher.search_gene_with_fallback(query, include_sequence))


def search_multiple_databases_sync(
    query: str, databases: Optional[list[DatabaseType]] = None, include_sequence: bool = True
) -> list[GeneSearchResult]:
    """Synchronous wrapper for multi-database search."""
    searcher = GeneSearcher()
    return asyncio.run(searcher.search_multiple_databases(query, databases, include_sequence))
