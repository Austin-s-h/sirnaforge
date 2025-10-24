#!/usr/bin/env python3
"""Demonstration of chemical modification annotation workflow.

This script shows how to:
1. Load designed siRNA candidates
2. Apply modification patterns
3. Export annotated sequences
4. Validate the results

Run this script to see a complete example workflow.
"""

import json
from pathlib import Path

from sirnaforge.models.modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SequenceRecord,
    SourceType,
    StrandMetadata,
    StrandRole,
)
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.modifications import merge_metadata_into_fasta, save_metadata_json


def create_example_candidates() -> list[SiRNACandidate]:
    """Create example siRNA candidates (simulating design output)."""
    return [
        SiRNACandidate(
            id="TP53_candidate_001",
            transcript_id="ENST00000269305",
            position=542,
            guide_sequence="GUAAUCUACUGGGACGGAACU",
            passenger_sequence="UUCCGUCCCAGUAGAUUACUU",
            gc_content=47.6,
            length=21,
            asymmetry_score=0.72,
            composite_score=87.3,
        ),
        SiRNACandidate(
            id="TP53_candidate_002",
            transcript_id="ENST00000269305",
            position=789,
            guide_sequence="CAGCACAUGACGGAGGCUGCC",
            passenger_sequence="GCAGCCUCCGUCAUGUGCUGC",
            gc_content=61.9,
            length=21,
            asymmetry_score=0.68,
            composite_score=84.1,
        ),
        SiRNACandidate(
            id="TP53_candidate_003",
            transcript_id="ENST00000269305",
            position=1203,
            guide_sequence="ACGUGUGUGAGGCUCUCUGGG",
            passenger_sequence="CAGAGAGCCUCACACACGUAG",
            gc_content=57.1,
            length=21,
            asymmetry_score=0.75,
            composite_score=89.2,
        ),
    ]


def apply_standard_2ome_pattern(sequence: str) -> list[ChemicalModification]:
    """Apply standard alternating 2'-O-methyl pattern."""
    # Alternating positions (1, 3, 5, 7, ...)
    positions = [i for i in range(1, len(sequence) + 1) if i % 2 == 1]
    return [ChemicalModification(type="2OMe", positions=positions)]


def apply_minimal_pattern(sequence: str) -> list[ChemicalModification]:
    """Apply minimal terminal modifications."""
    seq_len = len(sequence)
    # Last 3 positions for 3' terminal protection
    return [ChemicalModification(type="2OMe", positions=[seq_len - 2, seq_len - 1, seq_len])]


def demonstrate_workflow():
    """Run complete demonstration workflow."""
    print("=" * 80)
    print("Chemical Modification Annotation Workflow Demo")
    print("=" * 80)
    print()

    # Step 1: Create or load designed candidates
    print("Step 1: Loading designed siRNA candidates...")
    candidates = create_example_candidates()
    print(f"  Loaded {len(candidates)} candidates")
    for c in candidates:
        print(f"    - {c.id}: score={c.composite_score:.1f}, GC={c.gc_content:.1f}%")
    print()

    # Step 2: Apply modification patterns
    print("Step 2: Applying modification patterns...")
    print("  Using standard 2'-O-methyl pattern for top candidates")

    metadata_dict = {}
    for i, candidate in enumerate(candidates):
        # Apply standard pattern to top 2 candidates, minimal to others
        if i < 2:
            pattern_name = "standard_2ome"
            guide_mods = apply_standard_2ome_pattern(candidate.guide_sequence)
            passenger_mods = apply_standard_2ome_pattern(candidate.passenger_sequence)
        else:
            pattern_name = "minimal_terminal"
            guide_mods = apply_minimal_pattern(candidate.guide_sequence)
            passenger_mods = apply_minimal_pattern(candidate.passenger_sequence)

        # Create metadata for guide strand
        guide_metadata = StrandMetadata(
            id=f"{candidate.id}_guide",
            sequence=candidate.guide_sequence,
            overhang="dTdT",
            chem_mods=guide_mods,
            provenance=Provenance(
                source_type=SourceType.DESIGNED,
                identifier=f"sirnaforge_demo_{candidate.id}",
                url="https://github.com/Austin-s-h/sirnaforge",
            ),
            confirmation_status=ConfirmationStatus.PENDING,
            notes=f"Top candidate targeting TP53, pattern: {pattern_name}",
        )

        # Create metadata for passenger strand
        passenger_metadata = StrandMetadata(
            id=f"{candidate.id}_passenger",
            sequence=candidate.passenger_sequence,
            overhang="dTdT",
            chem_mods=passenger_mods,
            provenance=Provenance(
                source_type=SourceType.DESIGNED,
                identifier=f"sirnaforge_demo_{candidate.id}",
                url="https://github.com/Austin-s-h/sirnaforge",
            ),
            confirmation_status=ConfirmationStatus.PENDING,
            notes=f"Passenger strand, pattern: {pattern_name}",
        )

        # Store both strands
        metadata_dict[f"{candidate.id}_guide"] = guide_metadata
        metadata_dict[f"{candidate.id}_passenger"] = passenger_metadata

        print(f"    - {candidate.id}: applied {pattern_name} pattern")
        print(f"      Guide: {len(guide_mods[0].positions)} modifications")
        print(f"      Passenger: {len(passenger_mods[0].positions)} modifications")

    print()

    # Step 3: Export to JSON
    print("Step 3: Exporting metadata to JSON...")
    output_dir = Path("/tmp/sirnaforge_demo")
    output_dir.mkdir(exist_ok=True)

    json_path = output_dir / "modifications.json"
    save_metadata_json(metadata_dict, json_path)
    print(f"  Saved to: {json_path}")
    print()

    # Step 4: Create FASTA and annotate
    print("Step 4: Creating annotated FASTA files...")

    # Create guide FASTA
    guide_fasta = output_dir / "candidates_guide.fasta"
    with guide_fasta.open("w") as f:
        for candidate in candidates:
            f.write(f">{candidate.id}_guide\n")
            f.write(f"{candidate.guide_sequence}\n")

    # Create passenger FASTA
    passenger_fasta = output_dir / "candidates_passenger.fasta"
    with passenger_fasta.open("w") as f:
        for candidate in candidates:
            f.write(f">{candidate.id}_passenger\n")
            f.write(f"{candidate.passenger_sequence}\n")

    # Annotate both
    guide_annotated = output_dir / "candidates_guide_annotated.fasta"
    passenger_annotated = output_dir / "candidates_passenger_annotated.fasta"

    merge_metadata_into_fasta(guide_fasta, json_path, guide_annotated)
    merge_metadata_into_fasta(passenger_fasta, json_path, passenger_annotated)

    print(f"  Guide (annotated): {guide_annotated}")
    print(f"  Passenger (annotated): {passenger_annotated}")
    print()

    # Step 5: Display results
    print("Step 5: Sample annotated FASTA content...")
    print("-" * 80)
    with guide_annotated.open() as f:
        # Show first record
        lines = f.readlines()
        if len(lines) >= 2:
            print(lines[0].strip())
            print(lines[1].strip())
    print("-" * 80)
    print()

    # Step 6: Create summary
    print("Step 6: Generating summary report...")
    summary = {
        "workflow": "chemical_modification_annotation",
        "candidates_processed": len(candidates),
        "patterns_applied": ["standard_2ome", "minimal_terminal"],
        "output_files": {
            "metadata_json": str(json_path),
            "guide_annotated_fasta": str(guide_annotated),
            "passenger_annotated_fasta": str(passenger_annotated),
        },
        "modification_stats": {
            "total_sequences": len(metadata_dict),
            "confirmed": sum(1 for m in metadata_dict.values() if m.confirmation_status == ConfirmationStatus.CONFIRMED),
            "pending": sum(1 for m in metadata_dict.values() if m.confirmation_status == ConfirmationStatus.PENDING),
        },
    }

    summary_path = output_dir / "summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Summary: {summary_path}")
    print()

    print("=" * 80)
    print("Workflow Complete!")
    print("=" * 80)
    print()
    print("Output files created in:", output_dir)
    print("  - modifications.json          - Chemical modification metadata")
    print("  - candidates_guide_annotated.fasta     - Guide strands with modifications")
    print("  - candidates_passenger_annotated.fasta - Passenger strands with modifications")
    print("  - summary.json                - Workflow summary")
    print()
    print("Next steps:")
    print("  1. Review the annotated FASTA files")
    print("  2. Use metadata JSON for synthesis planning")
    print("  3. Share with collaborators or synthesis vendors")
    print("  4. Update confirmation_status after experimental validation")
    print()


if __name__ == "__main__":
    demonstrate_workflow()
