# Chemical Modification Integration Guide for Developers

This guide is for developers who want to integrate chemical modification metadata into their siRNAforge workflows, scripts, or applications.

## Quick Start

### 1. Basic Workflow

```python
from sirnaforge.models.modifications import StrandMetadata, ChemicalModification
from sirnaforge.modifications import save_metadata_json

# Create metadata for your siRNA
metadata = StrandMetadata(
    id="my_sirna_001",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    overhang="dTdT",
    chem_mods=[
        ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9])
    ]
)

# Save to JSON file
save_metadata_json({"my_sirna_001": metadata}, "modifications.json")
```

### 2. Load Existing Patterns

```python
from pathlib import Path
import json

# Load example pattern
pattern_file = Path("examples/modification_patterns/standard_2ome.json")
with pattern_file.open() as f:
    pattern = json.load(f)

# Extract modification positions
guide_positions = pattern["guide_modifications"]["2OMe"]["positions"]
```

### 3. Annotate FASTA Files

```python
from sirnaforge.modifications import merge_metadata_into_fasta

# Merge metadata into FASTA headers
merge_metadata_into_fasta(
    fasta_path="candidates.fasta",
    metadata_path="modifications.json",
    output_path="candidates_annotated.fasta"
)
```

## Integration Patterns

### Pattern A: Post-Processing (Recommended for Most Cases)

Add modifications AFTER the main design workflow completes.

```python
from sirnaforge.workflow import run_sirna_workflow
from sirnaforge.modifications import save_metadata_json
from sirnaforge.models.modifications import StrandMetadata, ChemicalModification

async def design_with_modifications(gene, output_dir):
    # Step 1: Run standard design
    results = await run_sirna_workflow(
        gene_query=gene,
        output_dir=output_dir,
        top_n_candidates=20
    )

    # Step 2: Apply modifications to top candidates
    metadata = {}
    for candidate in results["design_summary"]["top_candidates"]:
        guide_metadata = StrandMetadata(
            id=f"{candidate.id}_guide",
            sequence=candidate.guide_sequence,
            overhang="dTdT",
            chem_mods=apply_modification_pattern(
                candidate.guide_sequence,
                pattern="standard_2ome"
            )
        )
        metadata[f"{candidate.id}_guide"] = guide_metadata

    # Step 3: Save metadata
    json_path = Path(output_dir) / "sirnaforge" / "modifications.json"
    save_metadata_json(metadata, json_path)

    return results, metadata

def apply_modification_pattern(sequence, pattern="standard_2ome"):
    """Apply a modification pattern to a sequence."""
    if pattern == "standard_2ome":
        # Alternating 2'-O-methyl
        positions = [i for i in range(1, len(sequence)+1) if i % 2 == 1]
        return [ChemicalModification(type="2OMe", positions=positions)]
    elif pattern == "minimal":
        # Terminal positions only
        positions = [len(sequence)-2, len(sequence)-1, len(sequence)]
        return [ChemicalModification(type="2OMe", positions=positions)]
    else:
        return []
```

**Pros:**
- Clean separation of design and modification logic
- No workflow slowdown
- Easy to experiment with different patterns

**Cons:**
- Requires separate step
- Modifications not considered during scoring

### Pattern B: Inline Integration

Add modifications during candidate creation.

```python
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.models.modifications import StrandMetadata, ChemicalModification

def create_candidate_with_modifications(
    guide_seq, passenger_seq, position, transcript_id, score
):
    """Create candidate with embedded modification metadata."""

    # Create modification metadata
    guide_metadata = StrandMetadata(
        id=f"candidate_{position}_guide",
        sequence=guide_seq,
        overhang="dTdT",
        chem_mods=[
            ChemicalModification(
                type="2OMe",
                positions=[i for i in range(1, len(guide_seq)+1) if i % 2 == 1]
            )
        ]
    )

    # Create candidate with metadata
    return SiRNACandidate(
        id=f"candidate_{position}",
        transcript_id=transcript_id,
        position=position,
        guide_sequence=guide_seq,
        passenger_sequence=passenger_seq,
        gc_content=calculate_gc_content(guide_seq),
        length=len(guide_seq),
        asymmetry_score=0.75,
        composite_score=score,
        guide_metadata=guide_metadata  # Embedded metadata
    )
```

**Pros:**
- Metadata travels with candidate
- Available in to_fasta() output
- Single data structure

**Cons:**
- More memory usage
- Harder to experiment with patterns later

### Pattern C: Hybrid (Best for Flexibility)

Store minimal info during design, enrich during export.

```python
class ModificationAnnotator:
    """Helper class for adding modifications to candidates."""

    def __init__(self, pattern_library_path):
        self.patterns = self._load_patterns(pattern_library_path)

    def _load_patterns(self, path):
        """Load modification patterns from JSON."""
        patterns = {}
        for pattern_file in Path(path).glob("*.json"):
            with pattern_file.open() as f:
                data = json.load(f)
                if "pattern_name" in data:
                    patterns[data["pattern_name"]] = data
        return patterns

    def annotate_candidates(
        self,
        candidates: list[SiRNACandidate],
        pattern_name: str = "standard_2ome"
    ) -> dict[str, StrandMetadata]:
        """Apply modification pattern to list of candidates."""

        pattern = self.patterns.get(pattern_name)
        if not pattern:
            raise ValueError(f"Pattern {pattern_name} not found")

        metadata = {}
        for candidate in candidates:
            # Create guide metadata
            guide_metadata = self._create_metadata_from_pattern(
                candidate.guide_sequence,
                candidate.id,
                pattern["guide_modifications"]
            )
            metadata[f"{candidate.id}_guide"] = guide_metadata

        return metadata

    def _create_metadata_from_pattern(self, sequence, sid, pattern_spec):
        """Create StrandMetadata from pattern specification."""
        modifications = []
        for mod_type, spec in pattern_spec.items():
            positions = self._apply_position_strategy(
                sequence,
                spec.get("strategy", "custom"),
                spec.get("positions", [])
            )
            modifications.append(
                ChemicalModification(type=mod_type, positions=positions)
            )

        return StrandMetadata(
            id=f"{sid}_guide",
            sequence=sequence,
            overhang="dTdT",
            chem_mods=modifications
        )

    def _apply_position_strategy(self, sequence, strategy, positions):
        """Calculate positions based on strategy."""
        if strategy == "alternating":
            return [i for i in range(1, len(sequence)+1) if i % 2 == 1]
        elif strategy == "terminal_3prime":
            return [len(sequence)-2, len(sequence)-1, len(sequence)]
        elif strategy == "custom":
            return positions
        else:
            return []

# Usage
annotator = ModificationAnnotator("examples/modification_patterns")
metadata = annotator.annotate_candidates(
    candidates=design_results.top_candidates,
    pattern_name="standard_2ome"
)
```

**Pros:**
- Flexible pattern management
- Easy to add new patterns
- Clean API

**Cons:**
- More code to maintain
- Requires pattern library

## Best Practices

### 1. Validation

Always validate modification positions against sequence length:

```python
from pydantic import ValidationError

try:
    metadata = StrandMetadata(
        id="test",
        sequence="AUCGAUCG",  # 8 nucleotides
        chem_mods=[
            ChemicalModification(type="2OMe", positions=[1, 2, 10])  # Position 10 > 8!
        ]
    )
except ValidationError as e:
    print(f"Invalid modification: {e}")
```

### 2. Pattern Library Organization

Organize patterns by use case:

```
modification_patterns/
├── cost_optimized/
│   ├── minimal.json
│   └── basic_terminal.json
├── standard/
│   ├── alternating_2ome.json
│   └── balanced.json
├── therapeutic/
│   ├── maximal_stability.json
│   └── fda_approved_onpattro.json
└── custom/
    └── my_pattern.json
```

### 3. Provenance Tracking

Always include provenance for traceability:

```python
from sirnaforge.models.modifications import Provenance, SourceType

metadata = StrandMetadata(
    id="my_candidate",
    sequence="AUCGAUCGAUCGAUCGAUCGA",
    provenance=Provenance(
        source_type=SourceType.DESIGNED,
        identifier=f"sirnaforge_v{version}_{gene}_{timestamp}",
        url="https://your-lab.org/experiments/exp-123"
    )
)
```

### 4. Version Control

Store modification metadata in version control alongside code:

```bash
git add modifications.json
git commit -m "Add modification pattern for TP53 candidates"
git tag -a "v1.0-modifications" -m "First modification pattern set"
```

### 5. Testing

Test your modification integration:

```python
def test_modification_application():
    """Test that modifications are correctly applied."""
    sequence = "AUCGAUCGAUCGAUCGAUCGA"

    # Apply pattern
    metadata = apply_pattern(sequence, "standard_2ome")

    # Validate
    assert len(metadata.chem_mods) > 0
    assert all(pos <= len(sequence) for mod in metadata.chem_mods for pos in mod.positions)
    assert metadata.sequence == sequence
```

## Common Use Cases

### Use Case 1: Batch Annotation

Annotate multiple FASTA files with the same pattern:

```python
from pathlib import Path

def batch_annotate(input_dir, pattern_name="standard_2ome"):
    """Annotate all FASTA files in directory."""
    for fasta_file in Path(input_dir).glob("*.fasta"):
        # Create metadata for sequences
        metadata = create_metadata_from_fasta(fasta_file, pattern_name)

        # Save JSON
        json_path = fasta_file.with_suffix(".json")
        save_metadata_json(metadata, json_path)

        # Create annotated FASTA
        output_path = fasta_file.parent / f"{fasta_file.stem}_annotated.fasta"
        merge_metadata_into_fasta(fasta_file, json_path, output_path)
```

### Use Case 2: Synthesis Planning

Generate synthesis order from annotated candidates:

```python
def generate_synthesis_order(metadata_dict, vendor="IDT"):
    """Generate synthesis order format for vendor."""
    order = []

    for strand_id, metadata in metadata_dict.items():
        entry = {
            "name": strand_id,
            "sequence": metadata.sequence,
            "scale": "25nm",
            "purification": "HPLC",
            "modifications": []
        }

        # Format modifications for vendor
        for mod in metadata.chem_mods:
            for pos in mod.positions:
                entry["modifications"].append({
                    "position": pos,
                    "type": mod.type,
                    "base": metadata.sequence[pos-1]
                })

        order.append(entry)

    return order
```

### Use Case 3: Cost Estimation

Estimate synthesis costs based on modifications:

```python
def estimate_cost(metadata: StrandMetadata) -> float:
    """Estimate synthesis cost based on modifications."""
    base_cost = 300.0  # Base cost for unmodified 21-mer

    # Cost per modification type
    costs = {
        "2OMe": 20.0,
        "2F": 30.0,
        "PS": 25.0,
        "LNA": 50.0
    }

    total_cost = base_cost
    for mod in metadata.chem_mods:
        mod_cost = costs.get(mod.type, 30.0)
        total_cost += len(mod.positions) * mod_cost

    # HPLC purification
    total_cost += 150.0

    return total_cost
```

## Troubleshooting

### Issue: Metadata not appearing in FASTA headers

**Solution:** Ensure you're using `merge_metadata_into_fasta()`:

```python
# Wrong - this won't add metadata
shutil.copy(input_fasta, output_fasta)

# Right - this merges metadata
merge_metadata_into_fasta(input_fasta, metadata_json, output_fasta)
```

### Issue: Validation errors for positions

**Solution:** Check sequence length and 1-based indexing:

```python
sequence = "AUCGAUCG"  # 8 nucleotides
positions = [1, 2, 3, 8]  # Valid (1-based, all <= 8)
positions = [0, 1, 2]     # Invalid (0 is not allowed)
positions = [1, 2, 9]     # Invalid (9 > 8)
```

### Issue: Pattern not found

**Solution:** Verify pattern file exists and has correct structure:

```python
from pathlib import Path
import json

pattern_file = Path("examples/modification_patterns/my_pattern.json")
if not pattern_file.exists():
    print(f"Pattern file not found: {pattern_file}")
else:
    with pattern_file.open() as f:
        data = json.load(f)
        print(f"Pattern name: {data.get('pattern_name')}")
```

## Advanced Topics

### Custom Modification Types

You can use any modification type string:

```python
# Standard modifications
ChemicalModification(type="2OMe", positions=[1, 3, 5])

# Proprietary modifications
ChemicalModification(type="CustomMod_v2", positions=[7, 14, 21])

# Delivery conjugates
ChemicalModification(type="GalNAc", positions=[21])  # Conjugated at 3' end
```

### Multi-Strand Complexes

For siRNAs with guide and passenger:

```python
duplex_metadata = {
    "sirna_001_guide": StrandMetadata(
        id="sirna_001_guide",
        sequence=guide_seq,
        chem_mods=[ChemicalModification(type="2OMe", positions=guide_positions)]
    ),
    "sirna_001_passenger": StrandMetadata(
        id="sirna_001_passenger",
        sequence=passenger_seq,
        chem_mods=[ChemicalModification(type="2OMe", positions=passenger_positions)]
    )
}
```

### Integration with Workflow Outputs

Add modifications to workflow outputs:

```python
# After workflow completes
output_dir = Path("results/TP53")
candidates_csv = output_dir / "sirnaforge" / "TP53_pass.csv"
candidates_fasta = output_dir / "sirnaforge" / "TP53_pass.fasta"

# Load candidates from CSV
df = pd.read_csv(candidates_csv)

# Create metadata for top candidates
metadata = {}
for _, row in df.head(10).iterrows():
    metadata[row['id']] = create_metadata_for_candidate(row)

# Save and annotate
json_path = output_dir / "sirnaforge" / "modifications.json"
save_metadata_json(metadata, json_path)

annotated_fasta = output_dir / "sirnaforge" / "TP53_pass_annotated.fasta"
merge_metadata_into_fasta(candidates_fasta, json_path, annotated_fasta)
```

## References

- **Specification:** `docs/modification_annotation_spec.md`
- **API Reference:** `docs/api_reference.rst`
- **Examples:** `examples/chemical_modifications_example.py`
- **Tests:** `tests/unit/test_modifications.py`

---

**Last Updated:** 2025-10-24
**Maintainer:** siRNAforge Development Team
