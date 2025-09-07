"""Gene search and sequence retrieval from multiple databases."""

import asyncio
from pathlib import Path
from typing import Optional, Union

from pydantic import BaseModel, ConfigDict, Field

from sirnaforge.data.base import (
    BaseEnsemblClient,
    DatabaseType,
    FastaUtils,
    GeneInfo,
    SequenceType,
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

    model_config = ConfigDict(use_enum_values=True)

    @property
    def success(self) -> bool:
        """Check if search was successful."""
        return self.gene_info is not None and len(self.transcripts) > 0


class GeneSearcher(BaseEnsemblClient):
    """Search genes and retrieve sequences from genomic databases."""

    def __init__(self, preferred_db: DatabaseType = DatabaseType.ENSEMBL, timeout: int = 30, max_retries: int = 3):
        """
        Initialize gene searcher.

        Args:
            preferred_db: Preferred database for searches
            timeout: Request timeout in seconds
            max_retries: Maximum retry attempts
        """
        super().__init__(timeout=timeout)
        self.preferred_db = preferred_db
        self.max_retries = max_retries

        # Database-specific configurations
        self.db_configs = {
            DatabaseType.ENSEMBL: {"base_url": "https://rest.ensembl.org", "species": "homo_sapiens"},
            DatabaseType.REFSEQ: {
                "base_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
                "database": "nucleotide",
            },
            DatabaseType.GENCODE: {"base_url": "https://www.gencodegenes.org", "version": "44"},
        }

    async def search_gene(
        self, query: str, database: Optional[DatabaseType] = None, include_sequence: bool = True
    ) -> GeneSearchResult:
        """
        Search for a gene and retrieve its isoforms.

        Args:
            query: Gene ID, gene name, or transcript ID
            database: Database to search (defaults to preferred_db)
            include_sequence: Whether to fetch transcript sequences

        Returns:
            GeneSearchResult with gene info and transcripts
        """
        db = database or self.preferred_db

        logger.info(f"Searching for '{query}' in {get_database_display_name(db)}")

        try:
            if db == DatabaseType.ENSEMBL:
                return await self._search_ensembl(query, include_sequence)
            if db == DatabaseType.REFSEQ:
                return await self._search_refseq(query, include_sequence)
            if db == DatabaseType.GENCODE:
                return await self._search_gencode(query, include_sequence)
            raise ValueError(f"Unsupported database: {db}")

        except Exception as e:
            logger.error(f"Search failed for '{query}' in {get_database_display_name(db)}: {e}")
            return GeneSearchResult(query=query, database=db, error=str(e))

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

    async def _search_ensembl(self, query: str, include_sequence: bool) -> GeneSearchResult:
        """Search Ensembl database."""
        # First, try to resolve the query to a gene
        gene_info = await self._ensembl_lookup_gene(query)

        if not gene_info:
            return GeneSearchResult(query=query, database=DatabaseType.ENSEMBL, error=f"Gene not found: {query}")

        # Get all transcripts for the gene
        transcripts = await self._ensembl_get_transcripts(gene_info.gene_id, include_sequence)

        return GeneSearchResult(
            query=query, database=DatabaseType.ENSEMBL, gene_info=gene_info, transcripts=transcripts
        )

    async def _ensembl_lookup_gene(self, query: str) -> Optional[GeneInfo]:
        """Look up gene information in Ensembl using inherited method."""
        data = await self.lookup_gene(query)

        if not data:
            return None

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

    async def _ensembl_get_transcripts(self, gene_id: str, include_sequence: bool) -> list[TranscriptInfo]:
        """Get all transcripts for a gene from Ensembl using inherited lookup method."""
        transcripts: list[TranscriptInfo] = []

        try:
            # Get transcript list with expansion
            data = await self.lookup_gene(gene_id, expand=True)

            if not data:
                logger.error(f"Failed to get transcript data for {gene_id}")
                return transcripts

            transcript_data = data.get("Transcript", [])

            for transcript in transcript_data:
                transcript_id = transcript.get("id")
                if not transcript_id:
                    continue

                sequence = None
                if include_sequence:
                    sequence = await self._ensembl_get_sequence(transcript_id)

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

        except Exception as e:
            logger.error(f"Failed to get transcripts for {gene_id}: {e}")

        return transcripts

    async def _ensembl_get_sequence(self, transcript_id: str) -> Optional[str]:
        """Get transcript sequence from Ensembl using inherited method."""
        return await self.get_sequence(transcript_id, SequenceType.CDNA)

    async def _search_refseq(self, query: str, include_sequence: bool) -> GeneSearchResult:  # noqa: ARG002
        """Search RefSeq database via NCBI E-utilities."""
        # This is a placeholder - would implement NCBI E-utilities API calls
        # include_sequence parameter kept for interface consistency
        return GeneSearchResult(query=query, database=DatabaseType.REFSEQ, error="RefSeq search not yet implemented")

    async def _search_gencode(self, query: str, include_sequence: bool) -> GeneSearchResult:  # noqa: ARG002
        """Search GENCODE database."""
        # This is a placeholder - would implement GENCODE API calls
        # include_sequence parameter kept for interface consistency
        return GeneSearchResult(query=query, database=DatabaseType.GENCODE, error="GENCODE search not yet implemented")

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


def search_multiple_databases_sync(
    query: str, databases: Optional[list[DatabaseType]] = None, include_sequence: bool = True
) -> list[GeneSearchResult]:
    """Synchronous wrapper for multi-database search."""
    searcher = GeneSearcher()
    return asyncio.run(searcher.search_multiple_databases(query, databases, include_sequence))
