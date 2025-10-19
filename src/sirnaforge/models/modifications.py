"""Data models for siRNA chemical modifications and metadata.

This module provides structured representations for chemical modifications,
overhangs, and provenance metadata associated with siRNA strands.
"""

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field, HttpUrl


class ConfirmationStatus(str, Enum):
    """Confirmation status for siRNA sequence data."""

    PENDING = "pending"
    CONFIRMED = "confirmed"


class SourceType(str, Enum):
    """Source type for siRNA provenance."""

    PATENT = "patent"
    PUBLICATION = "publication"
    CLINICAL_TRIAL = "clinical_trial"
    DATABASE = "database"
    DESIGNED = "designed"
    OTHER = "other"


class Provenance(BaseModel):
    """Provenance information for siRNA sequences.

    Tracks the origin and validation status of siRNA sequences.
    """

    source_type: SourceType = Field(description="Type of source for this siRNA")
    identifier: str = Field(description="Source identifier (e.g., patent number, DOI, PubMed ID)")
    url: Optional[str] = Field(default=None, description="URL to the source document")

    def to_header_string(self) -> str:
        """Convert provenance to FASTA header format.

        Returns:
            Formatted string like "Patent:US10060921B2"
        """
        return f"{self.source_type.value.replace('_', ' ').title()}:{self.identifier}"


class ChemicalModification(BaseModel):
    """Chemical modification annotation for siRNA strands.

    Represents a specific type of chemical modification and the positions
    where it occurs in the sequence.
    """

    type: str = Field(
        description="Modification type (e.g., '2OMe', '2F', 'PS')",
        examples=["2OMe", "2F", "PS", "LNA"],
    )
    positions: list[int] = Field(
        default_factory=list,
        description="1-based positions in the sequence where this modification occurs",
    )

    def to_header_string(self) -> str:
        """Convert modification to FASTA header format.

        Returns:
            Formatted string like "2OMe(1,4,6,11,13,16,19)" or "2F()" for no positions
        """
        if self.positions:
            pos_str = ",".join(str(p) for p in sorted(self.positions))
            return f"{self.type}({pos_str})"
        return f"{self.type}()"


class StrandRole(str, Enum):
    """Role of the siRNA strand in the duplex."""

    GUIDE = "guide"
    SENSE = "sense"
    ANTISENSE = "antisense"
    PASSENGER = "passenger"


class StrandMetadata(BaseModel):
    """Complete metadata for a single siRNA strand.

    This model captures all relevant information about a siRNA strand
    including sequence, modifications, overhangs, and provenance.
    """

    id: str = Field(description="Unique identifier for this strand")
    sequence: str = Field(description="RNA or DNA sequence")
    overhang: Optional[str] = Field(
        default=None,
        description="Overhang sequence (e.g., 'dTdT' for DNA, 'UU' for RNA)",
    )
    chem_mods: list[ChemicalModification] = Field(
        default_factory=list,
        description="List of chemical modifications applied to this strand",
    )
    notes: Optional[str] = Field(default=None, description="Additional notes or comments")
    provenance: Optional[Provenance] = Field(default=None, description="Source and validation information")
    confirmation_status: ConfirmationStatus = Field(
        default=ConfirmationStatus.PENDING,
        description="Experimental confirmation status",
    )

    def to_fasta_header(self, target_gene: Optional[str] = None, strand_role: Optional[StrandRole] = None) -> str:
        """Generate FASTA header with embedded metadata.

        Args:
            target_gene: Target gene name
            strand_role: Role of this strand in the duplex

        Returns:
            FASTA header string with key-value pairs
        """
        parts = [f">{self.id}"]
        
        if target_gene:
            parts.append(f"Target={target_gene}")
        
        if strand_role:
            parts.append(f"Role={strand_role.value}")
        
        parts.append(f"Confirmed={self.confirmation_status.value}")
        
        if self.overhang:
            parts.append(f"Overhang={self.overhang}")
        
        if self.chem_mods:
            mods_str = "|".join(mod.to_header_string() for mod in self.chem_mods)
            parts.append(f"ChemMods={mods_str}")
        
        if self.provenance:
            parts.append(f"Provenance={self.provenance.to_header_string()}")
            if self.provenance.url:
                parts.append(f"URL={self.provenance.url}")
        
        # Join with semicolons and space for readability
        return " ".join([parts[0]] + ["; ".join(parts[1:])])


class SequenceRecord(BaseModel):
    """Complete sequence record with strand metadata.

    Associates a strand with its target and role information.
    """

    target_gene: str = Field(description="Target gene symbol")
    strand_role: StrandRole = Field(description="Role of this strand (guide, sense, antisense, passenger)")
    metadata: StrandMetadata = Field(description="Complete strand metadata including modifications")

    def to_fasta(self) -> str:
        """Generate complete FASTA record.

        Returns:
            Multi-line FASTA string with header and sequence
        """
        header = self.metadata.to_fasta_header(
            target_gene=self.target_gene,
            strand_role=self.strand_role
        )
        return f"{header}\n{self.metadata.sequence}\n"
