# Chemical Modification Patterns Reference

Comprehensive guide to chemical modification patterns for siRNA stability, delivery, and therapeutic applications.

## Overview

Chemical modifications enhance siRNA stability and reduce off-target effects while maintaining on-target potency. This guide covers built-in patterns, FDA-approved examples, and best practices for custom pattern design.

## Built-In Modification Patterns

### Standard 2'-O-Methyl (Recommended)

**File:** `examples/modification_patterns/standard_2ome.json`

**Strategy:**
- **Guide strand:** Alternating 2'-O-methyl at odd positions (1, 3, 5, 7, 9, 11, 13, 15, 17, 19)
- **Passenger strand:** Offset alternating at even positions (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
- **Overhang:** dTdT (DNA) on both 3' ends

**Properties:**
- Serum stability: ~24 hours
- Synthesis cost: 1.5× unmodified
- Nuclease resistance: Moderate-high
- RISC loading: High efficiency

**Best For:**
- General research applications
- In vitro efficacy studies
- Initial in vivo feasibility
- Industry-standard baseline

**Example:**
```bash
sirnaforge workflow TP53 \
  --modification-file examples/modification_patterns/standard_2ome.json \
  --output-dir tp53_standard
```

---

### Minimal Terminal (Cost-Optimized)

**File:** `examples/modification_patterns/minimal_terminal.json`

**Strategy:**
- **Guide strand:** 3' terminal only (19, 20, 21)
- **Passenger strand:** 5' terminal only (1, 2)
- **Overhang:** dTdT (DNA)

**Properties:**
- Serum stability: ~6 hours
- Synthesis cost: 1.1× unmodified
- Nuclease resistance: Low-moderate
- RISC loading: High efficiency

**Best For:**
- High-throughput screening
- Cost-sensitive experiments
- In vitro only (cell culture with serum-free media)
- Transient knockdown studies

**Limitations:**
- Not suitable for in vivo work
- Requires frequent media changes
- Limited serum stability

**Example:**
```bash
sirnaforge design input.fasta \
  --modification-file examples/modification_patterns/minimal_terminal.json \
  --output minimal_cost.csv
```

---

### Maximal Stability (Therapeutic-Grade)

**File:** `examples/modification_patterns/maximal_stability.json`

**Strategy:**
- **Guide strand:** Complete 2'-O-methyl at all positions (1-21) + PS linkages at terminal dinucleotides
- **Passenger strand:** Complete 2'-O-methyl (1-21) + 5' PS linkages
- **Overhang:** dTdT with optional PS modifications

**Properties:**
- Serum stability: ~72 hours
- Synthesis cost: 3.0× unmodified
- Nuclease resistance: Very high
- RISC loading: High efficiency

**Best For:**
- In vivo efficacy studies
- Therapeutic development programs
- Preclinical toxicology
- Long-term knockdown experiments

**Requirements:**
- Specialized delivery (LNP or GalNAc conjugates)
- HPLC purification recommended
- Extended synthesis time (2-4 weeks)

**Example:**
```bash
sirnaforge workflow TTR \
  --modification-file examples/modification_patterns/maximal_stability.json \
  --genome-species human \
  --output-dir ttr_therapeutic
```

**FDA-Approved References:**
- **Patisiran (Onpattro)** - TTR, 2018, hATTR amyloidosis
- **Givosiran (Givlaari)** - ALAS1, 2019, Acute hepatic porphyria

---

## FDA-Approved Example: Patisiran (Onpattro)

**File:** `examples/modification_patterns/fda_approved_onpattro.json`

First FDA-approved RNAi therapeutic (2018) targeting transthyretin (TTR) for hereditary transthyretin-mediated amyloidosis.

### Guide Strand (Patisiran)

**Sequence:** `AUGGAAUACUCUUGGUUAC`

**Modifications:**
- 2'-O-methyl at positions: 1, 4, 6, 11, 13, 16, 19
- Overhang: dTdT
- Strategic pattern optimized for stability + RISC loading

**Rationale:**
- Balances nuclease resistance with guide strand activity
- Maintains A/U-rich 5' seed region accessibility
- 7 modifications provide excellent stability without over-modification

### Passenger Strand (Patisiran)

**Sequence:** `GUAACCAAGAGUAUUCCAU`

**Modifications:**
- 2'-O-methyl at positions: 3, 8, 10, 15
- Overhang: dTdT
- Limited modifications to promote degradation

**Rationale:**
- Fewer modifications than guide (intentional)
- Promotes preferential RISC loading of guide strand
- Reduces passenger strand activity/off-targets

### Clinical Success

**Efficacy:**
- Significant TTR reduction in Phase III trials (80% knockdown)
- FDA approved August 2018
- First-in-class RNAi therapeutic

**Delivery:**
- Lipid nanoparticle (LNP) formulation
- IV infusion: 0.3 mg/kg every 3 weeks
- Hepatocyte-targeted delivery

**Use This Pattern For:**
- Liver-targeted therapeutic development
- Benchmarking modification strategies
- Regulatory submission templates
- Educational/training purposes

**Example:**
```bash
# Design TTR-targeting siRNA using Patisiran pattern
sirnaforge workflow TTR \
  --modification-file examples/modification_patterns/fda_approved_onpattro.json \
  --genome-species human \
  --output-dir ttr_patisiran_template
```

**References:**
TODO: Review
- Adams et al. (2018) NEJM 379:11-21. PMID: 30145929
- US Patent US10060921B2
- FDA Drug Approval Package (2018)

---

## Modification Types Reference

### 2'-O-Methyl (2OMe)

**Most Common Choice**

- **Chemistry:** Methyl group at 2' position of ribose
- **Stability gain:** ++ (good)
- **Cost factor:** 1.2-1.5× per modification
- **RISC compatibility:** Excellent
- **Typical positions:** Alternating or custom patterns

**Pros:**
- Industry standard
- Well-characterized
- Compatible with all synthesis platforms
- Maintains Watson-Crick pairing

**Cons:**
- Moderate cost increase
- May require 10+ modifications for therapeutic stability

---

### Phosphorothioate (PS)

**Internucleotide Linkages**

- **Chemistry:** Sulfur replaces non-bridging oxygen in phosphate backbone
- **Stability gain:** +++ (very good)
- **Cost factor:** 1.3-1.8× per linkage
- **Typical positions:** Terminal dinucleotides (5' and 3')

**Pros:**
- Excellent nuclease resistance
- Enhances protein binding (albumin)
- Can improve pharmacokinetics

**Cons:**
- Potential for non-specific binding
- Synthesis complexity increases
- Can affect duplex stability

**Best Practice:**
- Limit to 2-4 linkages per strand
- Focus on terminal positions
- Often combined with 2'-O-methyl

---

### 2'-Fluoro (2F)

**Pyrimidine-Specific**

- **Chemistry:** Fluorine at 2' position (C and U only)
- **Stability gain:** +++ (very good)
- **Cost factor:** 1.5-2.0×
- **Typical positions:** All pyrimidines

**Pros:**
- Superior nuclease resistance vs 2OMe
- Maintains duplex stability
- Good RISC compatibility

**Cons:**
- Higher synthesis cost
- Pyrimidine-restricted (can't modify A/G)
- Less commonly used than 2OMe

---

### Locked Nucleic Acid (LNA)

**High-Affinity Modifications**

- **Chemistry:** Methylene bridge locks ribose in C3'-endo conformation
- **Stability gain:** ++++ (excellent)
- **Cost factor:** 2-3× per residue
- **Typical positions:** Sparse (every 3-4 nt)

**Pros:**
- Extremely high binding affinity
- Superior nuclease resistance
- Very effective at low frequency

**Cons:**
- Expensive
- Can inhibit RISC loading if over-used
- Requires careful positioning
- Not all synthesis providers offer

**Best Practice:**
- Use sparingly (2-3 per strand maximum)
- Avoid seed region (positions 2-8)
- Often combined with 2'-O-methyl

---

### 2'-O-Methoxyethyl (MOE)

**Alternative to 2OMe**

- **Chemistry:** Methoxyethyl group at 2' position
- **Stability gain:** ++ (good)
- **Cost factor:** 1.5-2.0×

**Pros:**
- Good nuclease resistance
- Reduced immunogenicity vs 2OMe
- Used in some antisense applications

**Cons:**
- Less common than 2OMe for siRNA
- Higher cost
- Limited track record in approved therapeutics

---

## Custom Pattern Design

### Design Principles

**1. Start Conservative**
- Begin with standard patterns
- Add modifications incrementally
- Test efficacy at each step

**2. Balance Competing Goals**
- Stability ↔ Cost
- Nuclease resistance ↔ RISC loading
- On-target potency ↔ Off-target reduction

**3. Consider Application**
- In vitro: Minimal modifications sufficient
- In vivo (research): Standard patterns
- Therapeutic: Maximal stability required

**4. Strand Asymmetry**
- More modifications on guide strand
- Fewer modifications on passenger strand
- Promotes guide strand RISC loading

### Pattern Template

Create custom JSON files following this structure:

```json
{
  "pattern_name": "my_custom_pattern",
  "description": "Brief description of strategy",
  "reference": "Literature or internal study reference",

  "guide_modifications": {
    "2OMe": {
      "positions": [1, 3, 5, 7, 9],
      "strategy": "alternating",
      "rationale": "Why these positions?"
    },
    "PS": {
      "positions": [],
      "internucleotide_linkages": [[1, 2], [20, 21]],
      "rationale": "Terminal linkages for stability"
    }
  },

  "passenger_modifications": {
    "2OMe": {
      "positions": [2, 4, 6],
      "strategy": "sparse",
      "rationale": "Minimal to promote degradation"
    }
  },

  "overhang": {
    "guide_3prime": "dTdT",
    "passenger_3prime": "dTdT",
    "rationale": "Standard DNA overhangs"
  },

  "estimated_properties": {
    "serum_stability_half_life_hours": 24,
    "relative_synthesis_cost": 1.5,
    "nuclease_resistance": "moderate-high",
    "risc_loading_efficiency": "high",
    "recommended_for": "in_vitro, initial_in_vivo"
  },

  "notes": [
    "Additional context",
    "Testing recommendations",
    "Known limitations"
  ]
}
```

### Position Selection Strategies

**Alternating (Balanced)**
```
Guide:     [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
Passenger: [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
```
- Standard industry approach
- Good stability/cost balance
- Maintains duplex structure

**Terminal-Heavy (Cost-Optimized)**
```
Guide:     [1, 2, 19, 20, 21]
Passenger: [1, 2]
```
- Minimal synthesis cost
- Protects most vulnerable positions
- Suitable for in vitro only

**Seed-Region Preservation**
```
Guide:     [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
Passenger: [full complement]
```
- Avoids positions 2-8 (seed region)
- May improve target recognition
- Based on some miRNA-mimetic designs

**Complete (Therapeutic)**
```
Guide:     [1-21] all positions
Passenger: [1-21] all positions
+ PS linkages at terminals
```
- Maximum stability
- Therapeutic-grade
- Highest cost

---

## Application-Specific Recommendations

### In Vitro Screening (96-well, 384-well)

**Recommended Pattern:** Minimal Terminal

**Rationale:**
- Cost is primary concern
- Short exposure time (<72 hours)
- Serum-free or low-serum media
- High-throughput compatible

**Example:**
```bash
sirnaforge design library.fasta \
  --modification-file examples/modification_patterns/minimal_terminal.json \
  --top-n 500 \
  --output hts_library.csv
```

---

### In Vitro Validation Studies

**Recommended Pattern:** Standard 2'-O-Methyl

**Rationale:**
- Industry standard for publications
- Good stability in culture media
- Reproducible results
- Comparable to literature

**Example:**
```bash
sirnaforge workflow GENE \
  --modification-file examples/modification_patterns/standard_2ome.json \
  --top-n 20
```

---

### In Vivo Proof-of-Concept

**Recommended Pattern:** Standard or Maximal Stability

**Rationale:**
- Need extended half-life
- Systemic delivery challenges
- Nuclease-rich environment
- Justifies higher cost

**Example:**
```bash
sirnaforge workflow TARGET \
  --modification-file examples/modification_patterns/maximal_stability.json \
  --genome-species mouse \
  --output-dir invivo_poc
```

---

### Therapeutic Development

**Recommended Pattern:** Maximal Stability (Custom Optimized)

**Rationale:**
- Regulatory requirements
- Long-term efficacy needed
- Safety/tox studies
- Cost justified by value

**Requirements:**
- GMP-grade synthesis
- HPLC purification
- Mass spec confirmation
- Batch-to-batch QC

**Consider:**
- GalNAc conjugation for hepatocytes
- Lipid nanoparticle formulation
- Antibody conjugates for targeting
- Patent landscape review

---

## Synthesis and Ordering

### Synthesis Vendors

**Major Providers:**
- Integrated DNA Technologies (IDT)
- Thermo Fisher (Dharmacon)
- Sigma-Aldrich
- GenePharma
- Biomers.net
- TriLink BioTechnologies

### Cost Estimates (2025)

**Base siRNA (21bp duplex, unmodified):** ~$200-400

**Modifications (per strand):**
- 2'-O-methyl: +$20-50 per modification
- Phosphorothioate: +$30-80 per linkage
- 2'-Fluoro: +$40-100 per modification
- LNA: +$100-200 per residue

**Additional Costs:**
- HPLC purification: +$100-300
- Mass spec QC: +$50-150
- Bulk synthesis (5+ sequences): -30-50% discount
- GMP-grade: 2-5× standard pricing

**Pattern Cost Examples:**
- Minimal terminal: ~$250-450
- Standard 2OMe: ~$400-600
- Maximal stability: ~$800-1,200
- Therapeutic + conjugate: $2,000-5,000+

### Ordering Checklist

**Required Information:**
1. ✅ Sequences (guide + passenger)
2. ✅ Modification map (positions + types)
3. ✅ Overhang specification
4. ✅ Purification method (HPLC recommended)
5. ✅ Scale (nmol) - typically 100-200 nmol research scale
6. ✅ QC requirements (mass spec, HPLC traces)

**Recommended Documentation:**
- Provide modification JSON file
- Include synthesis notes from vendor
- Archive lot numbers and QC data
- Store at -20°C or -80°C

### Timeline Expectations

- Unmodified: 5-10 business days
- Standard modifications: 2-3 weeks
- Extensive modifications: 3-4 weeks
- GMP/therapeutic grade: 6-12 weeks
- Custom conjugates: 4-8 weeks

---

## Best Practices

### 1. Design Before You Modify

✅ **Do:** Design unmodified siRNAs first, validate efficacy, then apply modifications

❌ **Don't:** Start with heavily modified sequences without unmodified baseline

### 2. Test Incrementally

✅ **Do:** Compare minimal → standard → maximal patterns

❌ **Don't:** Jump to maximal modifications without intermediate validation

### 3. Document Everything

✅ **Do:** Record modification patterns, lot numbers, QC data, storage conditions

❌ **Don't:** Rely on memory or informal notes

### 4. Match Pattern to Application

✅ **Do:** Use minimal for screening, standard for validation, maximal for in vivo

❌ **Don't:** Over-spend on stability you don't need

### 5. Consider Delivery Method

✅ **Do:** Design modifications compatible with your delivery system (LNP, electroporation, etc.)

❌ **Don't:** Ignore how modifications affect delivery vehicle interactions

### 6. Archive Modification Files

✅ **Do:** Version control JSON files with sequences

❌ **Don't:** Recreate modification schemes from memory

---

## Integration with siRNAforge

### Command-Line Usage

```bash
# Apply pattern during design
sirnaforge design input.fasta \
  --modification-file pattern.json \
  --output designed.csv

# Apply pattern during workflow
sirnaforge workflow GENE \
  --modification-file pattern.json \
  --output-dir results
```

### Python API

```python
from sirnaforge.modifications import load_metadata, save_metadata_json
from sirnaforge.models.modifications import StrandMetadata, ChemicalModification

# Load built-in pattern
pattern = load_metadata("examples/modification_patterns/standard_2ome.json")

# Create custom metadata
metadata = StrandMetadata(
    id="custom_sirna_001",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9])
    ]
)

# Save for later use
save_metadata_json({"custom_sirna_001": metadata}, "my_mods.json")
```

### Annotate Existing FASTA

```bash
# Merge modification metadata into FASTA headers
sirnaforge sequences annotate \
  candidates.fasta \
  modifications.json \
  -o candidates_annotated.fasta
```

---

## Troubleshooting

### Low Efficacy After Modification

**Possible Causes:**
- Over-modification in seed region (positions 2-8)
- Excessive passenger modifications reducing RISC loading
- Delivery method incompatibility

**Solutions:**
- Reduce modifications in seed region
- Ensure guide strand preference (asymmetric modification)
- Test unmodified sequence as control

### High Cost Synthesis

**Possible Causes:**
- Too many modifications
- Exotic modification types (LNA, custom)
- Small scale orders

**Solutions:**
- Use minimal or standard patterns
- Order multiple sequences together (bulk discount)
- Request quotes from multiple vendors

### Inconsistent Results

**Possible Causes:**
- Vendor variability
- Storage degradation
- Batch differences

**Solutions:**
- Specify HPLC purification
- Request lot-to-lot QC data
- Store properly (-20°C or -80°C)
- Use consistent vendor/synthesis method

---

## References

### Key Publications

1. **Alnylam Platform Papers**
   - Foster et al. (2018) "Advanced siRNA designs further improve in vivo performance" - Mol Ther 26:708-717

2. **FDA-Approved Therapeutics**
   - Adams et al. (2018) "Patisiran, an RNAi therapeutic..." - NEJM 379:11-21
   - Balwani et al. (2020) "Givosiran for acute hepatic porphyria" - NEJM 382:2289-2301

3. **Modification Chemistry**
   - Deleavey & Damha (2012) "Designing chemically modified oligonucleotides" - Chem Biol 19:937-954

4. **Design Principles**
   - Khvorova & Watts (2017) "The chemical evolution of oligonucleotide therapies" - Nat Biotechnol 35:238-248

### Regulatory Guidance

- FDA Guidance for Industry: siRNA-Based Therapeutics (2020)
- EMA Guideline on Quality of Oligonucleotide-Based Medicinal Products (2021)

### Patents

- US10060921B2: Patisiran (Onpattro) - Alnylam
- US9605264B2: RNAi modifications - Alnylam
- Multiple patents covering specific modification patterns

---

## Contributing Patterns

To contribute new patterns to siRNAforge:

1. **Create JSON File:** Follow the template structure
2. **Validate Format:** Test with siRNAforge tools
3. **Document Rationale:** Include references and use cases
4. **Add Examples:** Provide working command-line examples
5. **Submit PR:** Via GitHub with detailed description

**Example Contribution:**
```bash
# Test your pattern
sirnaforge design test.fasta \
  --modification-file my_new_pattern.json \
  --output test_output.csv

# Verify it works
cat test_output.csv | head -10
```

---

**Document Version:** 1.0
**Last Updated:** December 2025
**Maintainer:** siRNAforge Team

For questions or contributions, see [GitHub repository](https://github.com/Austin-s-h/sirnaforge).
