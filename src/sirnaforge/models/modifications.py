"""Data models for siRNA chemical modifications and metadata.

This module provides structured representations for chemical modifications,
overhangs, and provenance metadata associated with siRNA strands.
"""

from collections.abc import Iterable
from enum import Enum
from typing import Any, Callable, Optional, TypeVar

from pydantic import BaseModel, Field, field_validator, model_validator

F = TypeVar("F", bound=Callable[..., Any])
FieldValidatorFactory = Callable[..., Callable[[F], F]]
ModelValidatorFactory = Callable[..., Callable[[F], F]]

field_validator_typed: FieldValidatorFactory = field_validator
model_validator_typed: ModelValidatorFactory = model_validator


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

    @field_validator_typed("type")
    @classmethod
    def validate_type(cls, v: str) -> str:
        """Validate modification type is not empty."""
        if not v or not v.strip():
            raise ValueError("Modification type cannot be empty")
        return v.strip()

    @field_validator_typed("positions")
    @classmethod
    def validate_positions(cls, v: list[int]) -> list[int]:
        """Validate positions are positive integers."""
        if any(p < 1 for p in v):
            raise ValueError("All positions must be >= 1 (1-based indexing)")
        # Remove duplicates and sort
        return sorted(set(v))

    def to_header_string(self) -> str:
        """Convert modification to FASTA header format.

        Returns:
            Formatted string like "2OMe(1,4,6,11,13,16,19)" or "2F()" for no positions
        """
        if self.positions:
            pos_str = ",".join(str(p) for p in self.positions)
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

    @field_validator_typed("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        """Validate sequence contains only valid nucleotides."""
        if not v:
            raise ValueError("Sequence cannot be empty")
        valid_bases = set("ATCGUatcgu")
        if not all(c in valid_bases for c in v):
            raise ValueError("Sequence contains invalid characters. Valid: A, T, C, G, U")
        return v.upper()

    @model_validator_typed(mode="after")
    def validate_modification_positions(self) -> "StrandMetadata":
        """Validate that modification positions don't exceed sequence length."""
        seq_len = len(self.sequence)
        for mod in self.chem_mods:
            if mod.positions and max(mod.positions) > seq_len:
                raise ValueError(
                    f"Modification {mod.type} has position {max(mod.positions)} but sequence length is only {seq_len}"
                )
        return self

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
            mods_str = "+".join(mod.to_header_string() for mod in self.chem_mods)
            parts.append(f"ChemMods={mods_str}")

        if self.provenance:
            parts.append(f"Provenance={self.provenance.to_header_string()}")
            if self.provenance.url:
                parts.append(f"URL={self.provenance.url}")

        if len(parts) == 1:
            return parts[0]

        metadata_segment = "; ".join(parts[1:])
        return f"{parts[0]} {metadata_segment}"

    # Convenience mapping-style access for backward compatibility with dict usage
    def __getitem__(self, item: str) -> Any:
        fields = type(self).model_fields
        if item in fields:
            return getattr(self, item)
        raise KeyError(item)

    def get(self, item: str, default: Any = None) -> Any:
        try:
            return self[item]
        except KeyError:
            return default

    def __contains__(self, item: object) -> bool:
        if not isinstance(item, str):
            return False
        return item in type(self).model_fields

    def keys(self) -> Iterable[str]:
        return tuple(type(self).model_fields.keys())

    def items(self) -> Iterable[tuple[str, Any]]:
        field_names = tuple(type(self).model_fields.keys())
        return ((key, getattr(self, key)) for key in field_names)


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
        header = self.metadata.to_fasta_header(target_gene=self.target_gene, strand_role=self.strand_role)
        return f"{header}\n{self.metadata.sequence}\n"
