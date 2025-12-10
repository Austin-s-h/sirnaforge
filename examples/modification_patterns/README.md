# Chemical Modification Pattern Library

This directory contains example chemical modification patterns for siRNA design and synthesis planning.

> **📖 Full Documentation:** See [`docs/tutorials/modification_patterns_reference.md`](../../docs/tutorials/modification_patterns_reference.md) for comprehensive coverage of all patterns, modification types, synthesis guidance, and best practices.

## Pattern Files

### Standard Patterns

| File | Pattern | Use Case | Cost | Stability |
|------|---------|----------|------|-----------|
| `minimal_terminal.json` | Terminal modifications only | In vitro screening, cost-sensitive | Low (1.1x) | Moderate |
| `standard_2ome.json` | Alternating 2'-O-methyl | General use, balanced | Medium (1.5x) | High |
| `maximal_stability.json` | Full modification + PS linkages | In vivo, therapeutics | High (3x) | Very High |

### FDA-Approved Examples

| File | Drug | Target | Approval | Indication |
|------|------|--------|----------|------------|
| `fda_approved_onpattro.json` | Patisiran (Onpattro) | TTR | 2018 | hATTR amyloidosis |

> **For detailed information** on each pattern, modification rationale, synthesis costs, and application-specific recommendations, see [`docs/tutorials/modification_patterns_reference.md`](../../docs/tutorials/modification_patterns_reference.md).

## Quick Start

## Quick Start

### 1. Apply Pattern to siRNA Design

```python
from sirnaforge.modifications import load_metadata, save_metadata_json
from sirnaforge.models.modifications import StrandMetadata, ChemicalModification

# Load pattern
pattern = load_metadata("examples/modification_patterns/standard_2ome.json")

# Apply to your candidate
candidate_metadata = {
    "my_sirna_001": StrandMetadata(
        id="my_sirna_001",
        sequence="AUCGAUCGAUCGAUCGAUCGA",
        overhang="dTdT",
        chem_mods=[
            ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
        ]
    )
}

# Save for use
save_metadata_json(candidate_metadata, "my_modifications.json")
```

### 2. Annotate FASTA with Modifications

```bash
# Merge modification metadata into FASTA headers
sirnaforge sequences annotate \
  my_candidates.fasta \
  my_modifications.json \
  -o my_candidates_annotated.fasta
```

### 3. View Annotated Sequences

```bash
# Display sequences with modification info
sirnaforge sequences show my_candidates_annotated.fasta
```

---

## Learn More

**Comprehensive Documentation:** [`docs/tutorials/modification_patterns_reference.md`](../../docs/tutorials/modification_patterns_reference.md)

Topics covered:
- Detailed pattern descriptions and rationale
- Modification chemistry reference (2OMe, PS, LNA, etc.)
- Application-specific recommendations
- Synthesis vendor information and cost estimates
- Custom pattern design principles
- FDA-approved therapeutic examples
- Troubleshooting and best practices

**Integration Guide:** [`docs/modification_integration_guide.md`](../../docs/modification_integration_guide.md)
- Python API usage
- Workflow integration
- Custom pattern creation
- FASTA annotation methods

---

**Last Updated:** December 2025
