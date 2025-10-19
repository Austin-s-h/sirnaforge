"""Helper functions for working with siRNA chemical modifications metadata.

This module provides utilities for:
- Parsing FASTA headers to extract modification metadata
- Loading metadata from JSON sidecar files
- Encoding/decoding modification annotations
"""

import re
from pathlib import Path
from typing import Any, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pydantic import TypeAdapter, ValidationError

from sirnaforge.models.modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SequenceRecord,
    SourceType,
    StrandMetadata,
    StrandRole,
)


def parse_chem_mods(chem_mods_str: str) -> list[ChemicalModification]:
    """Parse ChemMods field from FASTA header.

    Args:
        chem_mods_str: String like "2OMe(1,4,6,11)|2F()"

    Returns:
        List of ChemicalModification objects
    """
    if not chem_mods_str:
        return []
    
    modifications = []
    # Split by pipe for multiple modification types
    for mod_part in chem_mods_str.split("|"):
        # Parse pattern: TYPE(pos1,pos2,...)
        match = re.match(r"(\w+)\(([\d,]*)\)", mod_part)
        if match:
            mod_type = match.group(1)
            pos_str = match.group(2)
            positions = [int(p) for p in pos_str.split(",") if p]
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
    metadata: dict[str, Any] = {
        "id": record.id,
        "sequence": str(record.seq),
    }
    
    if not record.description or record.description == record.id:
        return metadata
    
    # Remove the ID from description if it's at the start
    desc = record.description
    if desc.startswith(record.id):
        desc = desc[len(record.id):].strip()
    
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
        try:
            metadata["strand_role"] = StrandRole(role_str.lower())
        except ValueError:
            pass
    
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


def load_metadata(json_path: str | Path) -> dict[str, StrandMetadata]:
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
    if not path.exists():
        return {}
    
    # Read JSON file
    json_text = path.read_text()
    
    # Parse JSON to get the raw data
    import json
    data = json.loads(json_text)
    
    # Support both direct dict and nested under "modifications" key
    if "modifications" in data:
        raw_metadata = data["modifications"]
    else:
        raw_metadata = data
    
    # Use Pydantic TypeAdapter to validate each metadata entry
    metadata_dict: dict[str, StrandMetadata] = {}
    for strand_id, meta_data in raw_metadata.items():
        # Validate and parse using Pydantic
        metadata_dict[strand_id] = StrandMetadata.model_validate(meta_data)
    
    return metadata_dict


def merge_metadata_into_fasta(
    fasta_path: str | Path,
    metadata_path: str | Path,
    output_path: str | Path,
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
    # Load and validate metadata using Pydantic
    metadata_dict = load_metadata(metadata_path)
    
    # Read FASTA
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    updated_count = 0
    output_records = []
    
    for record in records:
        seq_id = record.id
        
        if seq_id in metadata_dict:
            # Metadata is already a validated StrandMetadata object
            strand_meta = metadata_dict[seq_id]
            
            # Extract target_gene and strand_role if present in the original metadata
            # These might be in the JSON but not part of StrandMetadata model
            target_gene = None
            strand_role = None
            
            # Try to parse from the metadata dict if it was loaded from JSON
            # For now, we'll use None defaults since StrandMetadata doesn't store these
            
            # Generate new header
            new_header = strand_meta.to_fasta_header(
                target_gene=target_gene,
                strand_role=strand_role
            )
            
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
    output_path: str | Path,
) -> None:
    """Save strand metadata to JSON file using Pydantic serialization.

    Args:
        metadata_dict: Dictionary mapping strand IDs to StrandMetadata objects
        output_path: Path to output JSON file
    """
    # Use Pydantic's model_dump with json mode for proper serialization
    output_data = {
        "modifications": {
            strand_id: meta.model_dump(mode="json", exclude_none=True)
            for strand_id, meta in metadata_dict.items()
        }
    }
    
    # Write using Pydantic's JSON serialization
    path = Path(output_path)
    
    # Use json module for pretty printing
    import json
    path.write_text(json.dumps(output_data, indent=2))
