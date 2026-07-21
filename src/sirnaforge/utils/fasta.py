"""Shared FASTA parsing helpers used across siRNAforge modules."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path


def _iter_fasta_records(fasta_path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(record_id, sequence)`` tuples from a FASTA file."""
    if not fasta_path.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}")

    current_id: str | None = None
    chunks: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(chunks)

                header = line[1:].strip()
                current_id = header.split()[0] if header else None
                chunks = []
                continue

            if current_id is not None:
                chunks.append(line)

    if current_id is not None:
        yield current_id, "".join(chunks)


def load_fasta_sequences(fasta_path: Path, *, uppercase: bool = True) -> dict[str, str]:
    """Load FASTA records into a ``{record_id: sequence}`` mapping."""
    sequences: dict[str, str] = {}
    for record_id, sequence in _iter_fasta_records(fasta_path):
        sequences[record_id] = sequence.upper() if uppercase else sequence

    if not sequences:
        raise ValueError(f"No sequences found in genome FASTA: {fasta_path}")
    return sequences


def load_fasta_contig_lengths(fasta_path: Path) -> dict[str, int]:
    """Load FASTA record lengths into a ``{record_id: length}`` mapping."""
    lengths: dict[str, int] = {}
    for record_id, sequence in _iter_fasta_records(fasta_path):
        lengths[record_id] = len(sequence)

    if not lengths:
        raise ValueError(f"No FASTA contigs found in {fasta_path}")
    return lengths
