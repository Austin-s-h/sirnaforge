#!/usr/bin/env python3
"""Example: Working with Chemical Modifications Metadata in siRNAforge

This script demonstrates how to:
1. Create chemical modification annotations
2. Generate FASTA files with metadata
3. Load and parse metadata from FASTA headers
4. Merge metadata from JSON files

For complete documentation, see: docs/modification_annotation_spec.md
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
from sirnaforge.modifications import load_metadata, merge_metadata_into_fasta, save_metadata_json


def example_create_metadata():
    """Example 1: Create metadata programmatically."""
    print("=" * 60)
    print("Example 1: Creating Chemical Modification Metadata")
    print("=" * 60)

    # Create chemical modifications
    mods = [
        ChemicalModification(
            type="2OMe",
            positions=[1, 4, 6, 11, 13, 16, 19],  # Every third base alternating
        ),
        ChemicalModification(
            type="2F",
            positions=[2, 5, 8, 12, 15, 18],  # Different positions
        ),
    ]

    # Create provenance information
    provenance = Provenance(
        source_type=SourceType.PATENT,
        identifier="US10060921B2",
        url="https://patents.google.com/patent/US10060921B2",
    )

    # Create complete strand metadata
    metadata = StrandMetadata(
        id="patisiran_ttr_guide",
        sequence="AUGGAAUACUCUUGGUUAC",
        overhang="dTdT",
        chem_mods=mods,
        provenance=provenance,
        confirmation_status=ConfirmationStatus.CONFIRMED,
        notes="Patisiran (Onpattro) - FDA approved 2018",
    )

    # Create sequence record
    record = SequenceRecord(
        target_gene="TTR",
        strand_role=StrandRole.GUIDE,
        metadata=metadata,
    )

    # Generate FASTA
    fasta = record.to_fasta()
    print("\nGenerated FASTA with metadata:")
    print(fasta)

    return {metadata.id: metadata}


def example_save_and_load_json(metadata_dict):
    """Example 2: Save and load metadata as JSON."""
    print("=" * 60)
    print("Example 2: Saving and Loading Metadata JSON")
    print("=" * 60)

    # Save to JSON
    json_path = Path("/tmp/example_metadata.json")
    save_metadata_json(metadata_dict, json_path)
    print(f"\n✓ Saved metadata to: {json_path}")

    # Load from JSON
    loaded_metadata = load_metadata(json_path)
    print(f"✓ Loaded {len(loaded_metadata)} metadata records")

    for strand_id, data in loaded_metadata.items():
        print(f"\nStrand: {strand_id}")
        print(f"  Sequence: {data['sequence'][:30]}...")
        print(f"  Overhang: {data['overhang']}")
        print(f"  Modifications: {len(data['chem_mods'])}")

    return json_path


def example_merge_metadata():
    """Example 3: Merge metadata into FASTA file."""
    print("=" * 60)
    print("Example 3: Merging Metadata into FASTA")
    print("=" * 60)

    # Create simple FASTA file
    fasta_path = Path("/tmp/example_sequences.fasta")
    with fasta_path.open("w", encoding="utf-8") as f:
        f.write(">patisiran_ttr_guide\n")
        f.write("AUGGAAUACUCUUGGUUAC\n")

    # Create metadata JSON
    json_path = Path("/tmp/example_metadata.json")
    metadata = {
        "patisiran_ttr_guide": {
            "id": "patisiran_ttr_guide",
            "sequence": "AUGGAAUACUCUUGGUUAC",
            "target_gene": "TTR",
            "strand_role": "guide",
            "overhang": "dTdT",
            "chem_mods": [
                {"type": "2OMe", "positions": [1, 4, 6, 11, 13, 16, 19]},
                {"type": "2F", "positions": [2, 5, 8, 12, 15, 18]},
            ],
            "provenance": {
                "source_type": "patent",
                "identifier": "US10060921B2",
                "url": "https://patents.google.com/patent/US10060921B2",
            },
            "confirmation_status": "confirmed",
        }
    }

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    # Merge metadata into FASTA
    output_path = Path("/tmp/example_annotated.fasta")
    count = merge_metadata_into_fasta(fasta_path, json_path, output_path)

    print(f"\n✓ Updated {count} sequences with metadata")
    print(f"✓ Output saved to: {output_path}")

    # Show result
    print("\nAnnotated FASTA content:")
    with output_path.open(encoding="utf-8") as f:
        print(f.read())


def example_common_modifications():
    """Example 4: Common modification patterns."""
    print("=" * 60)
    print("Example 4: Common Chemical Modification Patterns")
    print("=" * 60)

    patterns = {
        "Alternate 2'-OMe": ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]),
        "End-capping PS": ChemicalModification(type="PS", positions=[1, 2, 20, 21]),
        "Seed region 2F": ChemicalModification(type="2F", positions=[2, 3, 4, 5, 6, 7, 8]),
        "Full 2'-OMe": ChemicalModification(type="2OMe", positions=list(range(1, 22))),  # All positions for 21-mer
    }

    print("\nCommon modification patterns:")
    for name, mod in patterns.items():
        print(f"\n  {name}:")
        print(f"    Type: {mod.type}")
        print(f"    Positions: {mod.positions[:5]}{'...' if len(mod.positions) > 5 else ''}")
        print(f"    Header format: {mod.to_header_string()[:50]}{'...' if len(mod.to_header_string()) > 50 else ''}")


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("siRNAforge Chemical Modifications Examples")
    print("=" * 60 + "\n")

    # Run examples
    metadata_dict = example_create_metadata()
    print()

    example_save_and_load_json(metadata_dict)
    print()

    example_merge_metadata()
    print()

    example_common_modifications()
    print()

    print("=" * 60)
    print("✅ All examples completed successfully!")
    print("=" * 60)
    print("\nNext steps:")
    print("  • Review the generated files in /tmp/example_*")
    print("  • Read docs/modification_annotation_spec.md for full API")
    print("  • Try: sirnaforge sequences --help")


if __name__ == "__main__":
    main()
