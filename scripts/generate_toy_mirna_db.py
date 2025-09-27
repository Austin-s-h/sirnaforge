#!/usr/bin/env python3
"""
Generate a realistic toy miRNA database for testing siRNA off-target analysis.

This creates a smaller FASTA file (~500KB) with miRNA sequences that may have
seed-region matches with our test siRNAs.
"""

import random
from pathlib import Path

# siRNA candidates from smoke test data (DNA versions)
SMOKE_TEST_SIRNAS = [
    "GTAGTCGATCAGCATCGTAGT",  # Best candidate from smoke test
    "GTAGTCGATCAGGCGTTTCAT",  # From short sequence
    "ATCGTAGTCGATCAGGCGTTT",  # Another candidate
    "TAGTCGATCAGCATCGTAGTC",  # Additional candidate
    "GATCAGCATCGTAGTCGATCA",  # More candidates for testing
]


def generate_mirna_sequence(length: int = 22) -> str:
    """Generate a miRNA-like sequence."""
    # miRNAs tend to have specific patterns, but for testing we'll use random
    return "".join(random.choices("ATCG", k=length))


def extract_seed_region(sirna: str) -> str:
    """Extract seed region (positions 2-7) from siRNA for miRNA matching."""
    return sirna[1:7]  # Positions 2-7 (0-indexed: 1-6)


def generate_mirna_with_seed_match(sirna_seed: str) -> str:
    """Generate a miRNA that has complementary seed region to siRNA."""
    # For miRNA seed matching, we need to create sequences that will align
    # The miRNA should contain the siRNA seed region (not complement)
    # because BWA will find alignments between siRNA and miRNA sequences

    # Embed the siRNA seed region in the miRNA
    prefix = generate_mirna_sequence(random.randint(1, 3))
    suffix_len = 22 - len(prefix) - len(sirna_seed)
    suffix = generate_mirna_sequence(max(1, suffix_len))

    # Truncate if too long
    full_seq = prefix + sirna_seed + suffix
    return full_seq[:22] if len(full_seq) > 22 else full_seq


def generate_toy_mirna_db(output_file: Path) -> None:
    """Generate a toy miRNA database."""
    print(f"Generating toy miRNA database: {output_file}")

    mirna_id = 1

    with output_file.open("w") as f:
        # 1. Create miRNAs with seed matches to our test siRNAs
        print("Creating miRNAs with seed region matches...")
        for sirna in SMOKE_TEST_SIRNAS:
            seed = extract_seed_region(sirna)

            # Create exact matches for testing - embed the full siRNA sequence
            mirna_seq = sirna  # Use the siRNA sequence directly as a miRNA
            mirna_name = f"hsa-miR-{mirna_id}a"
            f.write(
                f">{mirna_name} Homo sapiens miR-{mirna_id} exact_match_sirna_{SMOKE_TEST_SIRNAS.index(sirna) + 1}\n"
            )
            f.write(f"{mirna_seq}\n")

            # Create seed-only matches
            for variant in range(2):
                mirna_seq = generate_mirna_with_seed_match(seed)
                mirna_name = f"hsa-miR-{mirna_id}{['b', 'c'][variant]}"

                f.write(
                    f">{mirna_name} Homo sapiens miR-{mirna_id} seed_match_sirna_{SMOKE_TEST_SIRNAS.index(sirna) + 1}\n"
                )
                f.write(f"{mirna_seq}\n")

            mirna_id += 1

        # 2. Add some common human miRNAs (realistic sequences)
        print("Adding common human miRNAs...")
        common_mirnas = [
            ("hsa-let-7a-5p", "TGAGGTAGTAGGTTGTATAGTT"),
            ("hsa-miR-21-5p", "TAGCTTATCAGACTGATGTTGA"),
            ("hsa-miR-155-5p", "TTAATGCTAATCGTGATAGGGGT"),
            ("hsa-miR-122-5p", "TGGAGTGTGACAATGGTGTTTG"),
            ("hsa-miR-16-5p", "TAGCAGCACGTAAATATTGGCG"),
            ("hsa-miR-92a-3p", "TATTGCACTTGTCCCGGCCTGT"),
            ("hsa-miR-25-3p", "CATTGCACTTGGTCTCGGTCTG"),
            ("hsa-miR-143-3p", "TGAGATGAAGCACTGTAGCTC"),
            ("hsa-miR-145-5p", "GTCCAGTTTTCCCAGGAATCCCT"),
            ("hsa-miR-200a-3p", "TAACACTGTCTGGTAACGATGT"),
        ]

        for name, seq in common_mirnas:
            f.write(f">{name} {name.replace('hsa-', 'Homo sapiens ')}\n")
            f.write(f"{seq}\n")

        # 3. Fill with random miRNA-like sequences
        print("Adding random miRNA sequences...")
        for i in range(50):  # Add 50 more random miRNAs
            mirna_seq = generate_mirna_sequence(random.randint(18, 24))
            f.write(f">hsa-miR-{1000 + i}-3p Random miRNA {i + 1}\n")
            f.write(f"{mirna_seq}\n")

    final_size = output_file.stat().st_size
    print("✅ Toy miRNA database generated!")
    print(f"   File: {output_file}")
    print(f"   Size: {final_size / 1024:.1f} KB")
    print(f"   miRNAs: ~{mirna_id + len(common_mirnas) + 50}")


def main() -> None:
    """Main function."""
    # Set random seed for reproducible test data
    random.seed(42)

    # Default to tests directory
    script_dir = Path(__file__).parent
    test_data_dir = script_dir.parent / "tests" / "unit" / "data"
    output_file = test_data_dir / "toy_mirna_db.fasta"

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Generate the database
    generate_toy_mirna_db(output_file)


if __name__ == "__main__":
    main()
