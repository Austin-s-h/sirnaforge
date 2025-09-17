#!/usr/bin/env python3
"""
Standalone validation script for smoke test data.

This script can be used to verify that the smoke test data is valid
without requiring the full pytest framework or siRNAforge installation.
Useful for quick validation in CI environments.
"""

import sys
from pathlib import Path


def validate_fasta_basic(fasta_path: Path) -> bool:
    """Basic FASTA format validation."""
    if not fasta_path.exists():
        print(f"❌ File not found: {fasta_path}")
        return False

    try:
        content = fasta_path.read_text()
        lines = content.strip().split("\n")

        # Basic format checks
        if not content.strip():
            print("❌ File is empty")
            return False

        # Count headers and sequences
        headers = [line for line in lines if line.startswith(">")]
        sequences = [line for line in lines if line and not line.startswith(">")]

        if len(headers) == 0:
            print("❌ No FASTA headers found")
            return False

        if len(sequences) == 0:
            print("❌ No sequences found")
            return False

        print(f"✅ Found {len(headers)} sequences")

        # Validate sequences are RNA-like
        rna_sequences = 0
        for i, seq in enumerate(sequences):
            seq_upper = seq.upper()
            if "U" in seq_upper:
                rna_sequences += 1
            if any(char not in "AUGC" for char in seq_upper):
                print(f"⚠️ Sequence {i + 1} contains non-RNA characters: {seq[:20]}...")

        if rna_sequences > 0:
            print(f"✅ Found {rna_sequences} RNA sequences (containing U)")
        else:
            print("⚠️ No RNA sequences found (none contain U)")

        # Check total size
        total_size = fasta_path.stat().st_size
        print(f"✅ File size: {total_size} bytes (good for smoke tests: {'<1KB' if total_size < 1024 else '>1KB'})")

        return True

    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return False


def main():
    """Main validation function."""
    # Find the smoke test data file
    script_dir = Path(__file__).parent
    test_data_path = script_dir / "data" / "smoke_test.fasta"

    print("🧪 Smoke Test Data Validation")
    print("=" * 40)

    if validate_fasta_basic(test_data_path):
        print("\n✅ Smoke test data is valid!")
        return 0
    print("\n❌ Smoke test data validation failed!")
    return 1


if __name__ == "__main__":
    sys.exit(main())
