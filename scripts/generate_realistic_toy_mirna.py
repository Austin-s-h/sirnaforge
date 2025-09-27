#!/usr/bin/env python3
"""
Create a realistic miRNA database with proper seed region matching for siRNA off-target analysis.
Uses a curated set of real human miRNAs and creates variants that will match our test siRNAs.
"""

import random
from pathlib import Path

# Real human miRNA sequences (curated from miRBase)
REAL_HUMAN_MIRNAS = [
    ("hsa-let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU"),
    ("hsa-let-7b-5p", "UGAGGUAGUAGGUUGUGUGGUU"),
    ("hsa-let-7c-5p", "UGAGGUAGUAGGUUGUAUGGUU"),
    ("hsa-miR-1-3p", "UGGAAUGUAAAGAAGUAUGUAU"),
    ("hsa-miR-16-5p", "UAGCAGCACGUAAAUAUUGGCG"),
    ("hsa-miR-21-5p", "UAGCUUAUCAGACUGAUGUUGA"),
    ("hsa-miR-22-3p", "AAGCUGCCAGUUGAAGAACUGU"),
    ("hsa-miR-23a-3p", "AUCACAUUGCCAGGGAUUUCC"),
    ("hsa-miR-27a-3p", "UUCACAGUGGCUAAGUUCCGC"),
    ("hsa-miR-29a-3p", "UAGCACCAUCUGAAAUCGGUUA"),
    ("hsa-miR-92a-3p", "UAUUGCACUUGUCCCGGCCUGU"),
    ("hsa-miR-122-5p", "UGGAGUGUGACAAUGGUGUUUG"),
    ("hsa-miR-124-3p", "UAAGGCACGCGGUGAAUGCC"),
    ("hsa-miR-125b-5p", "UCCCUGAGACCCUAACUUGUGA"),
    ("hsa-miR-143-3p", "UGAGAUGAAGCACUGUAGCUC"),
    ("hsa-miR-145-5p", "GUCCAGUUUUCCCAGGAAUCCCU"),
    ("hsa-miR-155-5p", "UUAAUGCUAAUCGUGAUAGGGGU"),
    ("hsa-miR-200a-3p", "UAACACUGUCUGGUAACGAUGU"),
    ("hsa-miR-200b-3p", "UAAUACUGCCUGGUAAUGAUGA"),
    ("hsa-miR-221-3p", "AGCUACAUUGUCUGCUGGGUUUC"),
]

# Test siRNA sequences (from smoke test data)
TEST_SIRNAS = [
    ("sirna_perfect_match_1", "GTAGTCGATCAGCATCGTAGT"),
    ("sirna_perfect_match_2", "GTAGTCGATCAGGCGTTTCAT"),
    ("sirna_test_mismatch", "GTAGTCGATCAGCTTCGTAGA"),
]


def convert_to_dna(rna_seq: str) -> str:
    """Convert RNA sequence to DNA (U -> T)."""
    return rna_seq.replace("U", "T")


def extract_seed_region(sequence: str, start: int = 1, end: int = 7) -> str:
    """Extract seed region from sequence (1-based indexing)."""
    return sequence[start - 1 : end]


def create_mirna_with_seed_match(sirna_seq: str, base_mirna: tuple) -> str:
    """Create a miRNA that contains the siRNA seed region for matching."""
    mirna_name, mirna_seq = base_mirna

    # Convert RNA miRNA to DNA
    mirna_dna = convert_to_dna(mirna_seq)

    # Get siRNA seed region (positions 2-7, as used in the off-target code)
    sirna_seed = extract_seed_region(sirna_seq, 2, 7)

    # Create a variant of the miRNA that contains the siRNA seed
    # Insert the siRNA seed at position 2-7 of the miRNA
    if len(mirna_dna) >= 7:
        # Replace positions 2-7 with the siRNA seed
        modified_mirna = mirna_dna[0] + sirna_seed + mirna_dna[7:]
    else:
        # If miRNA is too short, just use the seed + some padding
        modified_mirna = "A" + sirna_seed + "AAAAAAAAAAA"

    return modified_mirna[:22]  # Trim to typical miRNA length


def generate_mirna_database() -> None:
    """Generate a realistic miRNA database with seed matches for our test siRNAs."""

    script_dir = Path(__file__).parent
    output_file = script_dir.parent / "tests" / "unit" / "data" / "toy_mirna_db.fasta"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"🧬 Creating miRNA database: {output_file}")

    with output_file.open("w") as f:
        mirna_id = 1

        # 1. Create miRNAs with exact seed matches for our test siRNAs
        print("📍 Creating miRNAs with exact seed matches...")
        for sirna_name, sirna_seq in TEST_SIRNAS:
            # Create 2-3 miRNAs that should match this siRNA's seed region
            for i in range(3):
                base_mirna = REAL_HUMAN_MIRNAS[i % len(REAL_HUMAN_MIRNAS)]
                modified_seq = create_mirna_with_seed_match(sirna_seq, base_mirna)

                f.write(f">test-miR-{mirna_id}-{i + 1} {base_mirna[0]}_seed_match_{sirna_name}\n")
                f.write(f"{modified_seq}\n")

                print(
                    f"  Created match for {sirna_name}: seed={extract_seed_region(sirna_seq, 2, 7)} in {modified_seq}"
                )

            mirna_id += 1

        # 2. Add some real human miRNAs (converted to DNA)
        print("🧬 Adding real human miRNAs...")
        for name, seq in REAL_HUMAN_MIRNAS[:10]:  # Add first 10
            dna_seq = convert_to_dna(seq)
            f.write(f">{name} {name.replace('hsa-', 'Homo sapiens ')}\n")
            f.write(f"{dna_seq}\n")

    file_size = output_file.stat().st_size
    print("✅ miRNA database created!")
    print(f"   File: {output_file}")
    print(f"   Size: {file_size / 1024:.1f} KB")

    # Verify seed regions
    print("\n🔍 Verification - Expected seed matches:")
    for sirna_name, sirna_seq in TEST_SIRNAS:
        seed = extract_seed_region(sirna_seq, 2, 7)
        print(f"  {sirna_name}: seed region = {seed}")


if __name__ == "__main__":
    random.seed(42)  # For reproducible results
    generate_mirna_database()
