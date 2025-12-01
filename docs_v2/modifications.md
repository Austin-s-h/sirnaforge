# Chemical Modifications

siRNAforge supports annotation and tracking of chemical modifications for therapeutic siRNA development.

## Overview

Chemical modifications improve siRNA stability, reduce off-target effects, and enhance delivery. siRNAforge provides:

- **Modification metadata model** - Track 2'-O-methyl, phosphorothioate, etc.
- **Pattern library** - Common modification patterns
- **FASTA annotation** - Embed modifications in sequence headers
- **JSON sidecar files** - Detailed modification records

## Quick Start

```python
from sirnaforge.models.modifications import StrandMetadata, ChemicalModification

# Create modification metadata
metadata = StrandMetadata(
    id="my_sirna_001",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9])
    ]
)
```

## Common Modification Types

| Type | Name | Purpose |
|------|------|---------|
| `2OMe` | 2'-O-methyl | Nuclease resistance, reduced immunogenicity |
| `2F` | 2'-fluoro | Enhanced binding affinity, stability |
| `PS` | Phosphorothioate | Nuclease resistance, protein binding |
| `LNA` | Locked Nucleic Acid | Very high affinity (use sparingly) |

## Modification Patterns

### Standard (Balanced)
Alternating 2'-O-methyl for stability without over-modification:
```
Positions: 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21
Type: 2OMe
```

### Minimal (Cost-Optimized)
Terminal modifications only:
```
Positions: 19, 20, 21 (3' end)
Type: 2OMe
```

### Therapeutic (FDA-approved pattern)
Based on Onpattro (patisiran):
```
Guide: 2OMe at positions 1, 4, 6, 11, 13, 16, 19
Passenger: Full modification
```

## CLI Usage

### View Modifications
```bash
sirnaforge sequences show candidates.fasta
sirnaforge sequences show candidates.fasta --format json
```

### Annotate FASTA
```bash
sirnaforge sequences annotate candidates.fasta modifications.json -o annotated.fasta
```

## JSON Metadata Format

```json
{
  "sirna_001": {
    "id": "sirna_001",
    "sequence": "AUCGAUCGAUCGAUCGAUCGA",
    "overhang": "dTdT",
    "chem_mods": [
      {"type": "2OMe", "positions": [1, 3, 5, 7, 9]}
    ],
    "provenance": {
      "source_type": "designed",
      "identifier": "sirnaforge_v1.0_TP53"
    }
  }
}
```

## FASTA Header Format

Modifications can be encoded in FASTA headers:

```
>sirna_001 Overhang=dTdT; ChemMods=2OMe(1,3,5,7,9)
AUCGAUCGAUCGAUCGAUCGA
```

## Integration with Workflow

Modifications are typically added post-design:

```bash
# 1. Design candidates
sirnaforge workflow TP53 --output-dir results/

# 2. Add modifications
sirnaforge sequences annotate \
  results/sirnaforge/TP53_pass.fasta \
  my_modifications.json \
  -o results/sirnaforge/TP53_modified.fasta
```

## Position Numbering

- All positions are **1-based** (first nucleotide = position 1)
- Positions are relative to **5' end**
- Example: `AUCG...` → position 1 = A, position 2 = U

## See Also

- [Python API](api.md) - Programmatic access to modification models
