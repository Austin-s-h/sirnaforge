"""Per-species transcript→gene index for multi-species hit classification.

This module provides species-keyed transcript-to-gene lookups, preventing the
regression where multiple species' annotations shared a single dict and collided.
Each species' index is built from its Ensembl cDNA FASTA headers and stored
separately.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcriptome_filter import TranscriptFilter
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


def _strip_version(identifier: str) -> str:
    """Strip Ensembl version suffixes (e.g., .9) for comparison.

    Args:
        identifier: Ensembl ID with optional version suffix

    Returns:
        ID without version suffix

    Note:
        This is a local implementation to keep the module dependency-light.
        An equivalent exists in workflow.py as _normalize_transcript_id.
    """
    return re.sub(r"\.\d+$", "", identifier.strip())


@dataclass(frozen=True)
class TranscriptRecord:
    """Parsed transcript metadata from an Ensembl cDNA header.

    Attributes:
        transcript_id: Version-stripped transcript identifier
        gene_id: Version-stripped gene identifier, or None if missing
        gene_symbol: Uppercase gene symbol, or None if missing
        biotype: Transcript biotype, preferring transcript_biotype over gene_biotype
    """

    transcript_id: str
    gene_id: str | None
    gene_symbol: str | None
    biotype: str | None


@dataclass
class SpeciesTranscriptIndex:
    """Transcript→gene index for a single species.

    Attributes:
        species: Canonical species name (via normalize_species_name)
        records: Mapping from version-stripped transcript ID to TranscriptRecord
        transcript_count: Total transcripts indexed
        missing_symbol_count: Count of records parsed without a gene_symbol
    """

    species: str
    records: dict[str, TranscriptRecord] = field(default_factory=dict)
    missing_symbol_count: int = 0
    _by_symbol: dict[str, set[str]] = field(default_factory=dict, repr=False)

    @property
    def transcript_count(self) -> int:
        """Number of distinct transcripts indexed."""
        return len(self.records)

    def add(self, record: TranscriptRecord) -> None:
        """Index one record, keeping the symbol reverse-map in step.

        Re-indexing the same transcript ID replaces the earlier record rather than
        double-counting it, so counts stay consistent with a FASTA carrying duplicates.
        """
        previous = self.records.get(record.transcript_id)
        if previous is not None:
            if previous.gene_symbol is None:
                self.missing_symbol_count -= 1
            else:
                self._by_symbol.get(previous.gene_symbol, set()).discard(previous.transcript_id)

        self.records[record.transcript_id] = record
        if record.gene_symbol is None:
            self.missing_symbol_count += 1
        else:
            self._by_symbol.setdefault(record.gene_symbol, set()).add(record.transcript_id)

    def gene_id_for(self, transcript_id: str) -> str | None:
        """Look up the gene ID for a transcript.

        Args:
            transcript_id: Transcript identifier (version suffix optional)

        Returns:
            Version-stripped gene ID, or None if not found
        """
        stripped_tid = _strip_version(transcript_id)
        record = self.records.get(stripped_tid)
        return record.gene_id if record else None

    def symbol_for(self, transcript_id: str) -> str | None:
        """Look up the gene symbol for a transcript.

        Args:
            transcript_id: Transcript identifier (version suffix optional)

        Returns:
            Uppercase gene symbol, or None if not found
        """
        stripped_tid = _strip_version(transcript_id)
        record = self.records.get(stripped_tid)
        return record.gene_symbol if record else None

    def biotype_for(self, transcript_id: str) -> str | None:
        """Look up the biotype for a transcript.

        Args:
            transcript_id: Transcript identifier (version suffix optional)

        Returns:
            Biotype string, or None if not found
        """
        stripped_tid = _strip_version(transcript_id)
        record = self.records.get(stripped_tid)
        return record.biotype if record else None

    def transcripts_for_symbol(self, symbol: str) -> frozenset[str]:
        """Find all transcripts for a given gene symbol.

        Args:
            symbol: Gene symbol (case-insensitive)

        Returns:
            Frozen set of version-stripped transcript IDs
        """
        return frozenset(self._by_symbol.get(symbol.upper(), ()))


class TranscriptGeneIndex:
    """Multi-species transcript→gene index manager.

    Maintains separate per-species indices to prevent collision when the same
    transcript ID string exists in multiple species' annotations.
    """

    def __init__(self) -> None:
        """Initialize empty index."""
        self._indices: dict[str, SpeciesTranscriptIndex] = {}

    def build(self, species: str, fasta_path: Path) -> SpeciesTranscriptIndex:
        """Build or retrieve the index for a species from its cDNA FASTA.

        This method is idempotent: calling it multiple times for the same species
        returns the already-built index without re-parsing.

        Args:
            species: Species name (will be normalized to canonical form)
            fasta_path: Path to Ensembl cDNA FASTA file

        Returns:
            SpeciesTranscriptIndex for the requested species

        Note:
            If the FASTA is missing or unreadable, logs a warning and returns
            an empty index. Never raises.
        """
        canonical_species = normalize_species_name(species)

        # Return existing index if already built
        if canonical_species in self._indices:
            return self._indices[canonical_species]

        index = SpeciesTranscriptIndex(species=canonical_species)

        # Handle missing/unreadable FASTA gracefully
        if not fasta_path.exists():
            logger.warning(f"FASTA not found for species '{canonical_species}' at {fasta_path}; returning empty index")
            self._indices[canonical_species] = index
            return index

        try:
            # Stream headers only, never load entire file
            with fasta_path.open("r") as fh:
                for raw_line in fh:
                    line = raw_line.strip()
                    if not line.startswith(">"):
                        continue

                    # Strip leading '>' and parse header
                    description = line[1:]
                    metadata = TranscriptFilter.parse_ensembl_header(description)

                    # Extract and strip versions from IDs
                    raw_transcript_id = description.split()[0]
                    transcript_id = _strip_version(raw_transcript_id)

                    gene_id_raw = metadata.get("gene")
                    gene_id = _strip_version(gene_id_raw) if gene_id_raw else None

                    # Uppercase symbol for case-insensitive comparison
                    gene_symbol_raw = metadata.get("gene_symbol")
                    gene_symbol = gene_symbol_raw.upper() if gene_symbol_raw else None

                    # Prefer transcript_biotype, fall back to gene_biotype
                    biotype = metadata.get("transcript_biotype") or metadata.get("gene_biotype")

                    index.add(
                        TranscriptRecord(
                            transcript_id=transcript_id,
                            gene_id=gene_id,
                            gene_symbol=gene_symbol,
                            biotype=biotype,
                        )
                    )

            logger.info(
                f"Built transcript index for {canonical_species}: "
                f"{index.transcript_count} transcripts, "
                f"{index.missing_symbol_count} without gene_symbol"
            )

        except Exception as e:
            # Discard whatever was parsed: a silently partial index misclassifies hits,
            # which is worse than having no index at all.
            logger.warning(
                f"Failed to read FASTA for species '{canonical_species}' at {fasta_path}: {e}; returning empty index"
            )
            index = SpeciesTranscriptIndex(species=canonical_species)

        self._indices[canonical_species] = index
        return index

    def for_species(self, species: str) -> SpeciesTranscriptIndex | None:
        """Retrieve the index for a species.

        Args:
            species: Species name (will be normalized to canonical form)

        Returns:
            SpeciesTranscriptIndex if previously built, None otherwise
        """
        canonical_species = normalize_species_name(species)
        return self._indices.get(canonical_species)

    def species(self) -> tuple[str, ...]:
        """Get all species with built indices.

        Returns:
            Tuple of canonical species names
        """
        return tuple(self._indices.keys())
