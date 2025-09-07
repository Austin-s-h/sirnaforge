#!/usr/bin/env bash
# Helper script to build BWA index files for a genome FASTA
# Usage: ./build_index.sh <genome_fasta> <out_prefix>
# Example: ./build_index.sh Homo_sapiens.GRCh38.dna.primary_assembly.fa nextflow_pipeline/genomes_indices/human/GRCh38/GRCh38

set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <genome_fasta> <out_prefix>"
  exit 2
fi

GENOME_FASTA=$1
OUT_PREFIX=$2

# Create output directory
OUT_DIR=$(dirname "$OUT_PREFIX")
mkdir -p "$OUT_DIR"

# Build BWA index (bwa-mem2 supports same index format as bwa)
if command -v bwa-mem2 >/dev/null 2>&1; then
  echo "Using bwa-mem2 to index: $GENOME_FASTA -> $OUT_PREFIX"
  bwa-mem2 index -p "$OUT_PREFIX" "$GENOME_FASTA"
elif command -v bwa >/dev/null 2>&1; then
  echo "Using bwa to index: $GENOME_FASTA -> $OUT_PREFIX"
  bwa index -p "$OUT_PREFIX" "$GENOME_FASTA"
else
  echo "Error: bwa or bwa-mem2 not found in PATH. Install bwa-mem2 or bwa." >&2
  exit 3
fi

# Optional: create a .fasta copy next to index for pipeline existence checks
if [ ! -f "${OUT_PREFIX}.fasta" ] && [ -f "$GENOME_FASTA" ]; then
  cp "$GENOME_FASTA" "${OUT_PREFIX}.fasta"
fi

echo "Index built: ${OUT_PREFIX} (and ${OUT_PREFIX}.fasta copied)"
