"""Helper functions for working with siRNA chemical modifications metadata.

This module provides utilities for:
- Parsing FASTA headers to extract modification metadata
- Loading metadata from JSON sidecar files
- Encoding/decoding modification annotations
"""

import json
import re
from pathlib import Path
from typing import Any, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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


def load_metadata(json_path: str | Path) -> dict[str, dict[str, Any]]:
    """Load metadata from JSON sidecar file.

    Args:
        json_path: Path to JSON file containing metadata

    Returns:
        Dictionary mapping strand IDs to metadata dictionaries

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
    
    with open(path, "r") as f:
        data = json.load(f)
    
    # Support both direct dict and nested under "modifications" key
    if "modifications" in data:
        return data["modifications"]
    return data


def merge_metadata_into_fasta(
    fasta_path: str | Path,
    metadata_path: str | Path,
    output_path: str | Path,
) -> int:
    """Merge metadata from JSON into FASTA headers.

    Args:
        fasta_path: Input FASTA file
        metadata_path: JSON file with metadata
        output_path: Output FASTA file with updated headers

    Returns:
        Number of sequences with metadata applied
    """
    # Load metadata
    metadata_dict = load_metadata(metadata_path)
    
    # Read FASTA
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    updated_count = 0
    output_records = []
    
    for record in records:
        seq_id = record.id
        
        if seq_id in metadata_dict:
            # Build StrandMetadata from JSON
            meta_data = metadata_dict[seq_id]
            
            # Handle chem_mods if present
            chem_mods = []
            if "chem_mods" in meta_data:
                for mod_dict in meta_data["chem_mods"]:
                    chem_mods.append(ChemicalModification(**mod_dict))
            
            # Handle provenance if present
            provenance = None
            if "provenance" in meta_data:
                provenance = Provenance(**meta_data["provenance"])
            
            # Build StrandMetadata
            strand_meta = StrandMetadata(
                id=seq_id,
                sequence=str(record.seq),
                overhang=meta_data.get("overhang"),
                chem_mods=chem_mods,
                notes=meta_data.get("notes"),
                provenance=provenance,
                confirmation_status=ConfirmationStatus(
                    meta_data.get("confirmation_status", "pending")
                ),
            )
            
            # Generate new header
            target_gene = meta_data.get("target_gene")
            strand_role = None
            if "strand_role" in meta_data:
                strand_role = StrandRole(meta_data["strand_role"])
            
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
    """Save strand metadata to JSON file.

    Args:
        metadata_dict: Dictionary mapping strand IDs to StrandMetadata objects
        output_path: Path to output JSON file
    """
    # Convert to dict format
    output_data = {
        "modifications": {
            strand_id: meta.model_dump(mode="json", exclude_none=True)
            for strand_id, meta in metadata_dict.items()
        }
    }
    
    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)
