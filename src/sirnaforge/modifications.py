"""Helper functions for working with siRNA chemical modifications metadata.

This module provides utilities for:
- Parsing FASTA headers to extract modification metadata
- Loading metadata from JSON sidecar files
- Encoding/decoding modification annotations
"""

import contextlib
import json
import re
from pathlib import Path
from typing import Any, Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sirnaforge.models.modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SourceType,
    StrandMetadata,
    StrandRole,
)


def parse_chem_mods(chem_mods_str: str) -> list[ChemicalModification]:
    """Parse ChemMods field from FASTA header.

    Args:
        chem_mods_str: String like "2OMe(1,4,6,11)+2F()"

    Returns:
        List of ChemicalModification objects
    """
    if not chem_mods_str:
        return []

    modifications = []
    # Support both historical '|' separator and the preferred '+' delimiter
    for raw_part in re.split(r"[|+]", chem_mods_str):
        mod_part = raw_part.strip()
        if not mod_part:
            continue
        # Parse pattern: TYPE(pos1,pos2,...)
        match = re.match(r"([^\(]+)\(([\d,\s]*)\)", mod_part)
        if match:
            mod_type = match.group(1).strip()
            pos_str = match.group(2)
            positions = [int(p.strip()) for p in pos_str.split(",") if p.strip()]
            modifications.append(ChemicalModification(type=mod_type, positions=positions))

    return modifications


def parse_provenance(prov_str: str, url: Optional[str] = None) -> Optional[Provenance]:
    """Parse Provenance field from FASTA header.

    Args:
        prov_str: String like "Patent:US10060921B2"
        url: Optional URL string

    Returns:
        Provenance object or None
    """
    if not prov_str:
        return None

    # Parse pattern: SourceType:Identifier
    match = re.match(r"([^:]+):(.+)", prov_str)
    if not match:
        return None

    source_type_str = match.group(1).strip().lower().replace(" ", "_")
    identifier = match.group(2).strip()

    # Map source type string to enum
    try:
        source_type = SourceType(source_type_str)
    except ValueError:
        source_type = SourceType.OTHER

    return Provenance(source_type=source_type, identifier=identifier, url=url)


def parse_header(record: SeqRecord) -> dict[str, Any]:
    """Parse FASTA header to extract metadata.

    Args:
        record: BioPython SeqRecord from FASTA file

    Returns:
        Dictionary with parsed metadata fields
    """
    # Parse header description for key-value pairs
    # Format: >id Target=TTR;Role=guide;Confirmed=pending;Overhang=dTdT;...
    record_id = record.id or ""
    description = record.description or ""

    metadata: dict[str, Any] = {
        "id": record_id,
        "sequence": str(record.seq),
    }

    if not description or description == record_id:
        return metadata

    # Remove the ID from description if it's at the start
    desc = description
    if record_id and desc.startswith(record_id):
        desc = desc[len(record_id) :].strip()

    # Split by semicolon for key=value pairs
    pairs = [pair.strip() for pair in desc.split(";")]

    parsed_fields = {}
    for pair in pairs:
        if "=" in pair:
            key, value = pair.split("=", 1)
            parsed_fields[key.strip()] = value.strip()

    # Extract known fields
    metadata["target_gene"] = parsed_fields.get("Target")

    role_str = parsed_fields.get("Role")
    if role_str:
        with contextlib.suppress(ValueError):
            metadata["strand_role"] = StrandRole(role_str.lower())

    metadata["overhang"] = parsed_fields.get("Overhang")

    # Parse chemical modifications
    chem_mods_str = parsed_fields.get("ChemMods")
    if chem_mods_str:
        metadata["chem_mods"] = parse_chem_mods(chem_mods_str)

    # Parse provenance
    prov_str = parsed_fields.get("Provenance")
    url_str = parsed_fields.get("URL")
    if prov_str:
        metadata["provenance"] = parse_provenance(prov_str, url_str)

    # Parse confirmation status
    confirmed_str = parsed_fields.get("Confirmed", "pending")
    try:
        metadata["confirmation_status"] = ConfirmationStatus(confirmed_str.lower())
    except ValueError:
        metadata["confirmation_status"] = ConfirmationStatus.PENDING

    metadata["notes"] = parsed_fields.get("Notes")

    return metadata


def _load_raw_metadata(path: Path) -> dict[str, Any]:
    """Load raw metadata mapping from JSON file without validation."""

    if not path.exists():
        return {}

    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("Metadata JSON root must be an object")

    raw_metadata = data.get("modifications", data)
    if not isinstance(raw_metadata, dict):
        raise ValueError("Metadata JSON 'modifications' entry must be an object")

    return raw_metadata


def _validate_metadata_entries(raw_metadata: dict[str, Any]) -> dict[str, StrandMetadata]:
    """Validate raw metadata payload into StrandMetadata instances."""

    metadata_dict: dict[str, StrandMetadata] = {}
    for strand_id, meta_data in raw_metadata.items():
        metadata_dict[strand_id] = StrandMetadata.model_validate(meta_data)
    return metadata_dict


def load_metadata(json_path: Union[str, Path]) -> dict[str, StrandMetadata]:
    """Load and validate metadata from JSON sidecar file using Pydantic.

    Args:
        json_path: Path to JSON file containing metadata

    Returns:
        Dictionary mapping strand IDs to StrandMetadata objects

    Raises:
        ValidationError: If JSON data doesn't match StrandMetadata schema

    Example JSON structure:
        {
            "modifications": {
                "patisiran_ttr_guide": {
                    "id": "patisiran_ttr_guide",
                    "sequence": "AUGGAAUACUCUUGGUUAC",
                    "overhang": "dTdT",
                    "chem_mods": [
                        {"type": "2OMe", "positions": [1, 4, 6, 11, 13, 16, 19]}
                    ],
                    "provenance": {
                        "source_type": "patent",
                        "identifier": "US10060921B2"
                    }
                }
            }
        }
    """
    path = Path(json_path)
    raw_metadata = _load_raw_metadata(path)
    return _validate_metadata_entries(raw_metadata)


def merge_metadata_into_fasta(
    fasta_path: Union[str, Path],
    metadata_path: Union[str, Path],
    output_path: Union[str, Path],
) -> int:
    """Merge metadata from JSON into FASTA headers.

    Uses Pydantic for automatic validation of metadata.

    Args:
        fasta_path: Input FASTA file
        metadata_path: JSON file with metadata
        output_path: Output FASTA file with updated headers

    Returns:
        Number of sequences with metadata applied

    Raises:
        ValidationError: If metadata doesn't match StrandMetadata schema
    """
    metadata_path = Path(metadata_path)
    raw_metadata = _load_raw_metadata(metadata_path)
    metadata_dict = _validate_metadata_entries(raw_metadata)

    # Read FASTA
    records = list(SeqIO.parse(fasta_path, "fasta"))

    updated_count = 0
    output_records = []

    for record in records:
        seq_id = record.id

        if seq_id in metadata_dict:
            # Metadata is already a validated StrandMetadata object
            strand_meta = metadata_dict[seq_id]

            # Extract extra fields preserved in the raw metadata payload
            raw_entry = raw_metadata.get(seq_id, {}) if isinstance(raw_metadata, dict) else {}
            target_gene = raw_entry.get("target_gene") if isinstance(raw_entry, dict) else None
            strand_role_value = None
            if isinstance(raw_entry, dict):
                strand_role_value = raw_entry.get("strand_role")

            strand_role = None
            if isinstance(strand_role_value, StrandRole):
                strand_role = strand_role_value
            elif isinstance(strand_role_value, str):
                try:
                    strand_role = StrandRole(strand_role_value.lower())
                except ValueError:
                    strand_role = None

            # Generate new header
            new_header = strand_meta.to_fasta_header(target_gene=target_gene, strand_role=strand_role)

            # Create new record with updated header
            # Remove the '>' from header for SeqRecord
            new_desc = new_header[1:] if new_header.startswith(">") else new_header
            new_record = SeqRecord(
                record.seq,
                id=seq_id,
                description=new_desc,
            )
            output_records.append(new_record)
            updated_count += 1
        else:
            # Keep original record
            output_records.append(record)

    # Write output
    SeqIO.write(output_records, output_path, "fasta")

    return updated_count


def save_metadata_json(
    metadata_dict: dict[str, StrandMetadata],
    output_path: Union[str, Path],
) -> None:
    """Save strand metadata to JSON file using Pydantic serialization.

    Args:
        metadata_dict: Dictionary mapping strand IDs to StrandMetadata objects
        output_path: Path to output JSON file
    """
    # Use Pydantic's model_dump with json mode for proper serialization
    output_data = {
        "modifications": {
            strand_id: meta.model_dump(mode="json", exclude_none=True) for strand_id, meta in metadata_dict.items()
        }
    }
    # Write using Pydantic's JSON serialization
    path = Path(output_path)

    # Use json module for pretty printing
    path.write_text(json.dumps(output_data, indent=2))
