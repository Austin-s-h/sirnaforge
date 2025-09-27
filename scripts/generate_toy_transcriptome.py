#!/usr/bin/env python3
"""
Generate a realistic toy transcriptome database for testing siRNA off-target analysis.

This creates a ~2MB FASTA file with:
1. Perfect matches for smoke test siRNAs
2. 1-3 mismatch sequences for off-target testing
3. Random sequences to reach target size
4. Biologically relevant transcript-like sequences
"""

import random
import sys
from pathlib import Path
from typing import Optional

# siRNA candidates from smoke test data (DNA versions)
SMOKE_TEST_SIRNAS = [
    "GTAGTCGATCAGCATCGTAGT",  # Best candidate from smoke test
    "GTAGTCGATCAGGCGTTTCAT",  # From short sequence
    "ATCGTAGTCGATCAGGCGTTT",  # Another candidate
    "TAGTCGATCAGCATCGTAGTC",  # Additional candidate
    "GATCAGCATCGTAGTCGATCA",  # More candidates for testing
]


def generate_random_dna(length: int) -> str:
    """Generate random DNA sequence."""
    return "".join(random.choices("ATCG", k=length))


def introduce_mismatches(sequence: str, num_mismatches: int) -> str:
    """Introduce specified number of mismatches into a sequence."""
    seq_list = list(sequence)
    positions = random.sample(range(len(sequence)), min(num_mismatches, len(sequence)))

    for pos in positions:
        # Replace with a different nucleotide
        original = seq_list[pos]
        choices = [n for n in "ATCG" if n != original]
        seq_list[pos] = random.choice(choices)

    return "".join(seq_list)


def generate_transcript_like_sequence(length: int, sirna_targets: Optional[list[str]] = None) -> str:
    """Generate a transcript-like sequence with optional embedded siRNA targets."""
    sequence = ""
    remaining_length = length

    # Embed some siRNA targets if provided
    if sirna_targets:
        for _i, sirna in enumerate(sirna_targets[:3]):  # Limit to 3 per transcript
            if remaining_length < len(sirna) + 100:
                break

            # Add some random sequence before the target
            prefix_len = random.randint(50, 200)
            prefix_len = min(prefix_len, remaining_length - len(sirna) - 50)

            sequence += generate_random_dna(prefix_len)
            sequence += sirna
            remaining_length -= prefix_len + len(sirna)

    # Fill remaining length with random DNA that mimics transcript patterns
    while remaining_length > 0:
        # Create blocks with different GC content to mimic real transcripts
        block_size = min(random.randint(50, 300), remaining_length)

        if random.random() < 0.3:  # 30% chance for GC-rich regions
            # GC-rich region (like exons)
            gc_choices = ["G", "C"] * 3 + ["A", "T"]  # 60% GC
            block = "".join(random.choices(gc_choices, k=block_size))
        elif random.random() < 0.5:  # 20% chance for AT-rich regions
            # AT-rich region (like introns)
            at_choices = ["A", "T"] * 3 + ["G", "C"]  # 60% AT
            block = "".join(random.choices(at_choices, k=block_size))
        else:  # 50% chance for balanced regions
            block = generate_random_dna(block_size)

        sequence += block
        remaining_length -= block_size

    return sequence


def generate_toy_transcriptome(output_file: Path, target_size_mb: float = 2.0) -> None:
    """Generate a toy transcriptome database of specified size."""
    target_size = int(target_size_mb * 1024 * 1024)  # Convert MB to bytes
    current_size = 0
    transcript_id = 1

    print(f"Generating toy transcriptome database: {output_file}")
    print(f"Target size: {target_size_mb:.1f} MB")

    with output_file.open("w") as f:
        # 1. Create transcripts with perfect siRNA matches
        print("Creating transcripts with perfect siRNA matches...")
        for i, sirna in enumerate(SMOKE_TEST_SIRNAS):
            transcript_length = random.randint(2000, 8000)
            sequence = generate_transcript_like_sequence(transcript_length, [sirna])

            header = f">ENST{transcript_id:08d}.1|{['TP53', 'BRCA1', 'EGFR', 'MYC', 'PTEN'][i % 5]}|protein_coding|perfect_match_{i + 1}"
            f.write(f"{header}\n")

            # Write sequence in 80 character lines
            for j in range(0, len(sequence), 80):
                f.write(f"{sequence[j : j + 80]}\n")

            current_size += len(header) + len(sequence) + (len(sequence) // 80) + 2
            transcript_id += 1

        # 2. Create transcripts with 1-3 mismatch versions for off-target testing
        print("Creating off-target sequences with 1-3 mismatches...")
        for mismatch_count in [1, 2, 3]:
            for i, sirna in enumerate(SMOKE_TEST_SIRNAS):
                if current_size >= target_size:
                    break

                # Create multiple variants with this mismatch count
                for variant in range(2):  # 2 variants per mismatch level
                    if current_size >= target_size:
                        break

                    modified_sirna = introduce_mismatches(sirna, mismatch_count)
                    transcript_length = random.randint(1500, 6000)
                    sequence = generate_transcript_like_sequence(transcript_length, [modified_sirna])

                    header = f">ENST{transcript_id:08d}.1|{['BRCA2', 'ATM', 'CHEK2', 'PALB2', 'RAD51'][i % 5]}|protein_coding|mismatch_{mismatch_count}_{variant + 1}"
                    f.write(f"{header}\n")

                    for j in range(0, len(sequence), 80):
                        f.write(f"{sequence[j : j + 80]}\n")

                    current_size += len(header) + len(sequence) + (len(sequence) // 80) + 2
                    transcript_id += 1

        # 3. Fill remaining space with random transcript-like sequences
        print("Filling with random transcript-like sequences...")
        while current_size < target_size:
            transcript_length = random.randint(1000, 10000)

            # Occasionally embed random siRNA targets for more realistic testing
            embed_targets = []
            if random.random() < 0.2:  # 20% chance
                embed_targets = random.sample(SMOKE_TEST_SIRNAS, k=random.randint(1, 2))
                embed_targets = [introduce_mismatches(s, random.randint(1, 3)) for s in embed_targets]

            sequence = generate_transcript_like_sequence(transcript_length, embed_targets)

            # Use realistic gene names
            gene_names = [
                "ACTB",
                "GAPDH",
                "TUBB",
                "HIST1H4C",
                "RPS27A",
                "UBC",
                "RPL13A",
                "SDHA",
                "HPRT1",
                "TBP",
                "GUSB",
                "HMBS",
                "TFRC",
                "PGK1",
                "LDHA",
                "ENO1",
                "ALDOA",
                "PKM",
                "GAPDH",
                "ACTG1",
                "CFL1",
                "EEF1A1",
                "RPL32",
                "RPS18",
            ]

            gene_name = random.choice(gene_names)
            biotype = random.choices(
                ["protein_coding", "lncRNA", "miRNA", "processed_transcript"], weights=[0.7, 0.15, 0.1, 0.05]
            )[0]

            header = f">ENST{transcript_id:08d}.1|{gene_name}|{biotype}|random_transcript"
            f.write(f"{header}\n")

            for j in range(0, len(sequence), 80):
                f.write(f"{sequence[j : j + 80]}\n")

            current_size += len(header) + len(sequence) + (len(sequence) // 80) + 2
            transcript_id += 1

            if transcript_id % 100 == 0:
                print(f"  Generated {transcript_id} transcripts, size: {current_size / 1024 / 1024:.2f} MB")

    final_size = output_file.stat().st_size
    print("✅ Toy transcriptome generated!")
    print(f"   File: {output_file}")
    print(f"   Size: {final_size / 1024 / 1024:.2f} MB")
    print(f"   Transcripts: {transcript_id - 1}")
    print(f"   Perfect matches: {len(SMOKE_TEST_SIRNAS)}")
    print(f"   Off-target variants: {len(SMOKE_TEST_SIRNAS) * 3 * 2}")


def main() -> None:
    """Main function."""
    # Set random seed for reproducible test data
    random.seed(42)

    # Determine output file
    if len(sys.argv) > 1:
        output_file = Path(sys.argv[1])
    else:
        # Default to tests directory
        script_dir = Path(__file__).parent
        test_data_dir = script_dir.parent / "tests" / "unit" / "data"
        output_file = test_data_dir / "toy_transcriptome_db.fasta"

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Generate the database
    generate_toy_transcriptome(output_file, target_size_mb=2.0)

    # Also generate simpler toy candidates that match our siRNAs
    candidates_file = output_file.parent / "toy_candidates.fasta"
    print(f"\n📝 Updating toy candidates: {candidates_file}")

    with candidates_file.open("w") as f:
        f.write(">sirna_perfect_match_1\n")
        f.write(f"{SMOKE_TEST_SIRNAS[0]}\n")  # Perfect match candidate
        f.write(">sirna_perfect_match_2\n")
        f.write(f"{SMOKE_TEST_SIRNAS[1]}\n")  # Another perfect match
        f.write(">sirna_test_mismatch\n")
        f.write(f"{introduce_mismatches(SMOKE_TEST_SIRNAS[0], 2)}\n")  # Should find off-targets

    print(f"✅ Updated toy candidates with {len(SMOKE_TEST_SIRNAS)} test sequences")


if __name__ == "__main__":
    main()
