#!/usr/bin/env python3
"""
Download and prepare miRNA databases from gold-standard sources.

This script downloads miRNA sequences from:
1. miRBase - The primary miRNA sequence database
2. MirGeneDB - Curated miRNA database with high-confidence annotations
3. TargetScan - miRNA family data for seed region analysis

Creates standardized FASTA files for siRNA off-target analysis.
"""

import gzip
import sys
import urllib.request
from pathlib import Path
from typing import TextIO


class MiRNADatabaseDownloader:
    """Download and process miRNA databases from multiple sources."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def download_file(self, url: str, output_path: Path, decompress_gz: bool = False) -> bool:
        """Download a file from URL."""
        try:
            print(f"📥 Downloading {url}")
            with urllib.request.urlopen(url) as response:
                data = response.read()

            if decompress_gz and url.endswith(".gz"):
                data = gzip.decompress(data)

            with output_path.open("wb") as f:
                f.write(data)

            print(f"✅ Downloaded to {output_path} ({len(data):,} bytes)")
            return True

        except Exception as e:
            print(f"❌ Failed to download {url}: {e}")
            return False

    def download_mirbase(self) -> Path:
        """Download miRBase mature miRNA sequences."""
        print("\n🧬 Downloading miRBase mature sequences...")

        # Updated miRBase mature sequences (all species) - new URL structure
        mirbase_url = "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz"
        mirbase_raw = self.output_dir / "mirbase_mature_raw.fa"

        if not self.download_file(mirbase_url, mirbase_raw, decompress_gz=True):
            raise RuntimeError("Failed to download miRBase data")

        # Process miRBase to extract human miRNAs
        mirbase_human = self.output_dir / "mirbase_human.fa"
        human_count = 0

        print("🔄 Extracting human miRNAs...")
        with mirbase_raw.open("r") as infile, mirbase_human.open("w") as outfile:
            write_seq = False

            for line in infile:
                if line.startswith(">"):
                    # Check if this is a human miRNA
                    write_seq = "hsa-" in line or "Homo sapiens" in line
                    if write_seq:
                        outfile.write(line)
                        human_count += 1
                elif write_seq:
                    outfile.write(line)

        print(f"✅ Extracted {human_count} human miRNAs from miRBase")
        return mirbase_human

    def _process_fasta_file(self, infile: Path, outfile: TextIO, source_name: str, seen_sequences: set[str]) -> int:
        """Process a FASTA file and write unique sequences to outfile."""
        count = 0
        header = None
        for line in infile.open("r"):
            stripped_line = line.strip()
            if not stripped_line:
                continue
            if stripped_line.startswith(">"):
                header = stripped_line
            elif header is not None:
                seq_upper = stripped_line.upper()
                if seq_upper not in seen_sequences:
                    outfile.write(f"{header} [{source_name}]\n{stripped_line}\n")
                    seen_sequences.add(seq_upper)
                    count += 1
                header = None
        return count

    def download_mirgenedb(self) -> Path:
        """Download MirGeneDB high-confidence miRNA sequences."""
        print("\n🎯 Downloading MirGeneDB sequences...")

        # MirGeneDB human sequences
        mirgenedb_url = "https://mirgenedb.org/fasta/9606?mat=1&ha=1"  # Human (9606) mature sequences
        mirgenedb_file = self.output_dir / "mirgenedb_human.fa"

        if not self.download_file(mirgenedb_url, mirgenedb_file):
            raise RuntimeError("Failed to download MirGeneDB data")

        # Count sequences
        seq_count = 0
        with mirgenedb_file.open("r") as f:
            for line in f:
                if line.startswith(">"):
                    seq_count += 1

        print(f"✅ Downloaded {seq_count} sequences from MirGeneDB")
        return mirgenedb_file

    def create_combined_database(self, mirbase_file: Path, mirgenedb_file: Path) -> Path:
        """Combine miRBase and MirGeneDB into a single non-redundant database."""
        print("\n🔄 Creating combined miRNA database...")

        combined_file = self.output_dir / "combined_human_mirnas.fa"
        seen_sequences: set[str] = set()
        total_sequences = 0

        with combined_file.open("w") as outfile:
            # Add miRBase sequences
            if mirbase_file and mirbase_file.exists():
                total_sequences += self._process_fasta_file(mirbase_file, outfile, "miRBase", seen_sequences)

            # Add MirGeneDB sequences
            if mirgenedb_file and mirgenedb_file.exists():
                total_sequences += self._process_fasta_file(mirgenedb_file, outfile, "MirGeneDB", seen_sequences)

        print(f"✅ Combined database created with {total_sequences} unique sequences")
        print(f"   File: {combined_file}")
        print(f"   Size: {combined_file.stat().st_size / 1024:.1f} KB")

        return combined_file

    def extract_seed_regions(self, mirna_file: Path) -> Path:
        """Extract seed regions from miRNAs for efficient matching."""
        print("\n🎯 Extracting seed regions for efficient matching...")

        seed_file = self.output_dir / "mirna_seed_regions.fa"
        seed_count = 0

        with mirna_file.open("r") as infile, seed_file.open("w") as outfile:
            header = None
            for line in infile:
                stripped_line = line.strip()
                if stripped_line.startswith(">"):
                    header = stripped_line
                elif stripped_line and header and len(stripped_line) >= 7:
                    # Extract seed region (positions 2-7, 0-indexed: 1-6)
                    seed_region = stripped_line[1:7].upper()
                    outfile.write(f"{header}_SEED\n{seed_region}\n")
                    seed_count += 1
                    header = None

        print(f"✅ Extracted {seed_count} seed regions")
        return seed_file


def main() -> int:
    """Main function to download and prepare miRNA databases."""
    print("🧬 miRNA Database Downloader for siRNAforge")
    print("=" * 50)

    # Set up output directory
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / "data" / "mirna_databases"

    # Initialize downloader
    downloader = MiRNADatabaseDownloader(output_dir)

    # Download databases
    print("📥 Downloading miRNA databases...")

    mirbase_file = downloader.download_mirbase()
    mirgenedb_file = downloader.download_mirgenedb()

    # Create combined database
    if mirbase_file or mirgenedb_file:
        combined_file = downloader.create_combined_database(mirbase_file, mirgenedb_file)

        # Create seed region database for efficient matching
        seed_file = downloader.extract_seed_regions(combined_file)

        print("\n✅ miRNA databases ready!")
        print(f"📁 Output directory: {output_dir}")
        print(f"🧬 Combined database: {combined_file.name}")
        print(f"🎯 Seed regions: {seed_file.name}")
        print(f"\n💡 Use in siRNAforge with: --mirna-database {combined_file}")

    else:
        print("❌ Failed to download any miRNA databases")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
