# Chemical Modification Pattern Library

This directory contains example chemical modification patterns for siRNA design and synthesis planning.

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

## Usage

### 1. Apply Pattern to Designed siRNAs

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

## Pattern Details

### Minimal Terminal (`minimal_terminal.json`)

**Strategy:**
- Guide: 2'-O-methyl at positions 19, 20, 21 (3' terminal)
- Passenger: 2'-O-methyl at positions 1, 2 (5' terminal)

**Best For:**
- Initial in vitro screening
- Cost-sensitive experiments
- Transient knockdown studies

**Limitations:**
- Limited serum stability (~6 hours)
- Not suitable for in vivo work
- Requires regular media changes in cell culture

### Standard 2'-O-Methyl (`standard_2ome.json`)

**Strategy:**
- Guide: Alternating 2'-O-methyl (positions 1,3,5,7,9,11,13,15,17,19)
- Passenger: Offset alternating pattern

**Best For:**
- General research use
- In vitro efficacy studies
- Initial in vivo feasibility

**Properties:**
- Good nuclease resistance (~24 hour serum half-life)
- Maintains RISC loading efficiency
- Reasonable synthesis cost (1.5x unmodified)
- Industry-standard pattern

### Maximal Stability (`maximal_stability.json`)

**Strategy:**
- Complete 2'-O-methyl modification on all positions
- Phosphorothioate linkages at terminal dinucleotides
- Similar to FDA-approved therapeutics

**Best For:**
- In vivo efficacy studies
- Therapeutic development
- Long-term knockdown experiments

**Properties:**
- Excellent serum stability (~72 hours)
- Very high nuclease resistance
- Higher synthesis cost (3x unmodified)
- May require specialized delivery (LNP, GalNAc)

### FDA-Approved Onpattro (`fda_approved_onpattro.json`)

**Description:**
Complete metadata for Patisiran (Onpattro), the first FDA-approved RNAi therapeutic targeting TTR for hereditary transthyretin amyloidosis.

**Key Features:**
- Validated in clinical trials
- Proven efficacy and safety
- Strategic modification pattern
- Lipid nanoparticle delivery

**Use Cases:**
- Reference for therapeutic development
- Benchmark for stability studies
- Template for liver-targeted siRNAs
- Educational example

## Customizing Patterns

You can create custom patterns by modifying these templates:

```json
{
  "my_custom_pattern": {
    "id": "my_custom_pattern",
    "sequence": "AUCGAUCGAUCGAUCGAUCGA",
    "overhang": "dTdT",
    "chem_mods": [
      {
        "type": "2OMe",
        "positions": [1, 5, 10, 15, 20]
      },
      {
        "type": "PS",
        "positions": []
      }
    ],
    "notes": "Custom pattern optimized for my specific application"
  }
}
```

## Modification Types Reference

| Type | Full Name | Typical Positions | Stability Gain | Cost Factor |
|------|-----------|-------------------|----------------|-------------|
| 2OMe | 2'-O-methyl | Alternating, all | ++ | 1.2-1.5x |
| 2F | 2'-fluoro | Pyrimidines | +++ | 1.5-2x |
| PS | Phosphorothioate | Terminal linkages | +++ | 1.3-1.8x |
| LNA | Locked Nucleic Acid | Sparse (every 3-4) | ++++ | 2-3x |
| MOE | 2'-O-methoxyethyl | Alternating | ++ | 1.5-2x |

## Best Practices

1. **Start Simple:** Begin with minimal or standard patterns
2. **Validate In Vitro:** Test efficacy before adding modifications
3. **Consider Cost:** More modifications ≠ better results
4. **Match Application:** Choose pattern based on experimental needs
5. **Document Decisions:** Use `notes` field to explain choices
6. **Include Provenance:** Track source of sequences/patterns
7. **Version Control:** Store modification files with sequences

## Synthesis Notes

### Ordering Modified siRNAs

When ordering from synthesis vendors:

1. **Provide Modification Map:** Use FASTA format with ChemMods annotation
2. **Specify Purity:** HPLC purification recommended for in vivo
3. **Check Compatibility:** Verify vendor supports your modification types
4. **Request QC:** Ask for mass spec confirmation
5. **Plan Ahead:** Modified synthesis takes 2-4 weeks

### Typical Vendors

- Integrated DNA Technologies (IDT)
- Thermo Fisher (Dharmacon)
- Sigma-Aldrich
- GenePharma
- Biomers.net

### Cost Expectations

Unmodified siRNA (21bp duplex): ~$200-400
+ 2'-O-methyl modifications: +20-50% per modification
+ Phosphorothioate: +30-80% per linkage
+ LNA: +100-200% per residue
+ HPLC purification: +$100-300
+ Bulk discount: 30-50% for multiple sequences

**Example:**
- Standard 2OMe pattern (10 mods): ~$400-600
- Maximal stability (21 mods + PS): ~$800-1200

## References

1. **Alnylam Pharmaceuticals** - Industry leader in RNAi therapeutics
2. **Tuschl et al.** - Original siRNA design rules
3. **Watts et al. (2018)** - Patisiran clinical trial results
4. **FDA Drug Approval Package** - Onpattro regulatory documents

## Contributing

To add new patterns:

1. Follow JSON format in existing examples
2. Include comprehensive notes and rationale
3. Provide references if based on publications
4. Test with sirnaforge annotation tools
5. Submit via pull request

## License

Examples based on published patents and scientific literature are provided for educational and research purposes. Commercial use may require licensing. Consult legal counsel for therapeutic development.

---

**Last Updated:** 2025-10-24  
**Maintainer:** siRNAforge Team
