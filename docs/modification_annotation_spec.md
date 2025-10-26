# Chemical Modification Annotation Specification

**Version:** 1.0
**Date:** 2025-10-19
**Status:** Implemented

## Overview

This specification defines how chemical modifications, overhangs, and provenance metadata are represented and managed in siRNAforge. The system provides a first-class data model for annotating modifications alongside siRNA sequences, ensuring annotations survive all pipeline stages without manual parsing.

## Data Models

### ChemicalModification

Represents a specific type of chemical modification and the positions where it occurs in the sequence.

**Fields:**
- `type` (str): Modification type (e.g., "2OMe", "2F", "PS", "LNA")
- `positions` (list[int]): 1-based positions in the sequence where this modification occurs

**Example:**
```python
from sirnaforge.models.modifications import ChemicalModification

mod = ChemicalModification(type="2OMe", positions=[1, 4, 6, 11, 13, 16, 19])
```

### Provenance

Tracks the origin and validation status of siRNA sequences.

**Fields:**
- `source_type` (SourceType): Type of source (patent, publication, clinical_trial, database, designed, other)
- `identifier` (str): Source identifier (e.g., patent number, DOI, PubMed ID)
- `url` (str, optional): URL to the source document

**Example:**
```python
from sirnaforge.models.modifications import Provenance, SourceType

prov = Provenance(
    source_type=SourceType.PATENT,
    identifier="US10060921B2",
    url="https://patents.google.com/patent/US10060921B2"
)
```

### StrandMetadata

Complete metadata for a single siRNA strand including sequence, modifications, overhangs, and provenance.

**Fields:**
- `id` (str): Unique identifier for this strand
- `sequence` (str): RNA or DNA sequence
- `overhang` (str, optional): Overhang sequence (e.g., "dTdT" for DNA, "UU" for RNA)
- `chem_mods` (list[ChemicalModification]): List of chemical modifications
- `notes` (str, optional): Additional notes or comments
- `provenance` (Provenance, optional): Source and validation information
- `confirmation_status` (ConfirmationStatus): Experimental confirmation status (pending, confirmed)

**Example:**
```python
from sirnaforge.models.modifications import (
    StrandMetadata,
    ChemicalModification,
    Provenance,
    SourceType,
    ConfirmationStatus
)

metadata = StrandMetadata(
    id="patisiran_ttr_guide",
    sequence="AUGGAAUACUCUUGGUUAC",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 4, 6, 11, 13, 16, 19])
    ],
    provenance=Provenance(
        source_type=SourceType.PATENT,
        identifier="US10060921B2"
    ),
    confirmation_status=ConfirmationStatus.CONFIRMED
)
```

### SequenceRecord

Associates a strand with its target and role information.

**Fields:**
- `target_gene` (str): Target gene symbol
- `strand_role` (StrandRole): Role of this strand (guide, sense, antisense, passenger)
- `metadata` (StrandMetadata): Complete strand metadata

## FASTA Header Encoding

The system uses standardized key-value pairs in FASTA headers to encode metadata:

```
>patisiran_ttr_guide Target=TTR; Role=guide; Confirmed=confirmed; Overhang=dTdT; ChemMods=2OMe(1,4,6,11,13,16,19); Provenance=Patent:US10060921B2; URL=https://patents.google.com/patent/US10060921B2
AUGGAAUACUCUUGGUUAC
```

**Format Rules:**
- Fields are separated by `;` (semicolon with space)
- Key-value pairs use `=` (equals sign)
- `ChemMods` syntax: `TYPE(pos1,pos2,...)` with multiple types separated by `|`
- Multiple modifications: `ChemMods=2OMe(1,4,6)|2F(2,5,8)|PS()`
- Empty positions allowed: `2F()` means modification type is annotated but no specific positions

## JSON Sidecar Format

Metadata can be stored in separate JSON files for easier curation:

```json
{
  "patisiran_ttr_guide": {
    "id": "patisiran_ttr_guide",
    "sequence": "AUGGAAUACUCUUGGUUAC",
    "target_gene": "TTR",
    "strand_role": "guide",
    "overhang": "dTdT",
    "chem_mods": [
      {
        "type": "2OMe",
        "positions": [1, 4, 6, 11, 13, 16, 19]
      }
    ],
    "provenance": {
      "source_type": "patent",
      "identifier": "US10060921B2",
      "url": "https://patents.google.com/patent/US10060921B2"
    },
    "confirmation_status": "confirmed",
    "notes": "Alnylam's patisiran (Onpattro) - FDA approved"
  }
}
```

## Python API

### Loading Metadata

```python
from sirnaforge.modifications import load_metadata

# Load from JSON file
metadata_dict = load_metadata("path/to/metadata.json")

# Access specific strand
strand_meta = metadata_dict["patisiran_ttr_guide"]
print(strand_meta["overhang"])  # "dTdT"
print(strand_meta["chem_mods"])  # List of modifications
```

### Parsing FASTA Headers

```python
from Bio import SeqIO
from sirnaforge.modifications import parse_header

# Parse FASTA file
records = SeqIO.parse("sequences.fasta", "fasta")
for record in records:
    metadata = parse_header(record)
    print(f"ID: {metadata['id']}")
    print(f"Target: {metadata.get('target_gene')}")
    print(f"Modifications: {metadata.get('chem_mods')}")
```

### Merging Metadata into FASTA

```python
from sirnaforge.modifications import merge_metadata_into_fasta

# Merge JSON metadata into FASTA headers
count = merge_metadata_into_fasta(
    fasta_path="sequences.fasta",
    metadata_path="metadata.json",
    output_path="sequences_annotated.fasta"
)
print(f"Updated {count} sequences")
```

### Creating Metadata Programmatically

```python
from sirnaforge.models.modifications import (
    StrandMetadata,
    ChemicalModification,
    SequenceRecord,
    StrandRole
)

# Create strand metadata
metadata = StrandMetadata(
    id="my_sirna_guide",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 3, 5])
    ]
)

# Create full sequence record
record = SequenceRecord(
    target_gene="KRAS",
    strand_role=StrandRole.GUIDE,
    metadata=metadata
)

# Generate FASTA
fasta = record.to_fasta()
print(fasta)
```

## CLI Usage

### Show Sequences with Metadata

Display sequences in table format:
```bash
sirnaforge sequences show sequences.fasta
```

Display specific sequence:
```bash
sirnaforge sequences show sequences.fasta --id patisiran_ttr_guide
```

Output as JSON:
```bash
sirnaforge sequences show sequences.fasta --format json
```

Output as FASTA:
```bash
sirnaforge sequences show sequences.fasta --format fasta
```

### Annotate FASTA with Metadata

Merge metadata from JSON into FASTA headers:
```bash
sirnaforge sequences annotate sequences.fasta metadata.json -o annotated.fasta
```

The annotate command:
1. Reads the input FASTA file
2. Loads metadata from the JSON file
3. Matches sequences by ID
4. Generates updated FASTA headers with metadata
5. Writes annotated sequences to output file

## Integration with SiRNACandidate

The `SiRNACandidate` model includes optional metadata fields:

```python
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.models.modifications import StrandMetadata

# Create candidate with metadata
candidate = SiRNACandidate(
    id="sirna_001",
    transcript_id="ENST00000123456",
    position=100,
    guide_sequence="AUCGAUCGAUCGAUCGAUCGA",
    passenger_sequence="UCGAUCGAUCGAUCGAUCGAU",
    gc_content=52.4,
    length=21,
    asymmetry_score=0.75,
    composite_score=85.2,
    guide_metadata=StrandMetadata(
        id="sirna_001_guide",
        sequence="AUCGAUCGAUCGAUCGAUCGA",
        overhang="dTdT"
    )
)

# Generate FASTA with metadata
fasta = candidate.to_fasta(include_metadata=True)
```

## Modification Type Reference

Common chemical modifications:

| Type | Full Name | Description |
|------|-----------|-------------|
| `2OMe` | 2'-O-methyl | Ribose 2' position methylation (nuclease resistance) |
| `2F` | 2'-fluoro | Ribose 2' position fluorination (enhanced stability) |
| `PS` | Phosphorothioate | Phosphate backbone sulfur substitution (nuclease resistance) |
| `LNA` | Locked Nucleic Acid | Bicyclic ribose analog (enhanced binding affinity) |
| `PNA` | Peptide Nucleic Acid | Peptide backbone (high stability) |
| `MOE` | 2'-O-methoxyethyl | Ribose modification (improved pharmacokinetics) |

**Note:** The system accepts free-text modification types, allowing for proprietary or novel modifications to be annotated.

## Position Numbering

- All positions are **1-based** (first nucleotide is position 1)
- Positions are relative to the **5' end** of the strand
- Example: For sequence "AUCGAUCG", position 1 is "A", position 4 is "G"

## Backward Compatibility

- FASTA files without metadata annotations work seamlessly
- Missing metadata fields default to `None` or empty values
- Legacy workflows continue to function without modification
- The system gracefully handles partial metadata

## Best Practices

1. **Consistent IDs**: Use consistent sequence IDs between FASTA and JSON files
2. **Validate Sources**: Include provenance information for all curated sequences
3. **Version Control**: Store metadata JSON files in version control alongside FASTA
4. **Documentation**: Use the `notes` field to document important details
5. **Confirmation Status**: Mark experimentally validated sequences as "confirmed"

## Example Workflow

1. **Create Initial FASTA**:
   ```bash
   # Start with basic sequences
   cat > sequences.fasta << EOF
   >sirna_001
   AUCGAUCGAUCGAUCGAUCGA
   EOF
   ```

2. **Curate Metadata**:
   ```bash
   # Create JSON with modifications
   cat > metadata.json << EOF
   {
     "sirna_001": {
       "id": "sirna_001",
       "target_gene": "BRCA1",
       "strand_role": "guide",
       "overhang": "dTdT",
       "chem_mods": [{"type": "2OMe", "positions": [1, 4, 6]}]
     }
   }
   EOF
   ```

3. **Merge Metadata**:
   ```bash
   sirnaforge sequences annotate sequences.fasta metadata.json -o annotated.fasta
   ```

4. **View Results**:
   ```bash
   sirnaforge sequences show annotated.fasta
   ```

5. **Use in Workflows**:
   ```python
   from sirnaforge.modifications import parse_header
   from Bio import SeqIO

   for record in SeqIO.parse("annotated.fasta", "fasta"):
       metadata = parse_header(record)
       # Process with full metadata available
   ```

## Future Extensions

Potential future enhancements:

- Validation rules for position-specific chemistry compatibility
- Support for duplex-level annotations (e.g., lipid conjugates)
- Delivery vehicle and tropism annotations
- Integration with chemical synthesis planning tools
- Visualization of modification patterns

## References

- siRNAforge GitHub: https://github.com/austin-s-h/sirnaforge
- Issue #[number]: Implement best practice system for metadata of chemical modifications
