#!/usr/bin/env bash
set -euo pipefail
# Script to download Ensembl GRCh38 chr22 and extract NF2 gene sequence (or a region)
# Produces candidates.fasta with 21-nt sliding-window sequences suitable as test input

OUT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
WORK_DIR="$OUT_DIR/tmp_test_data"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

CHR=22
# Ensembl FTP FASTA for Homo sapiens GRCh38 chromosome 22 (primary assembly)
# Using Ensembl release current via REST; fallback to FTP URL
FASTA_URL="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"

echo "Downloading chr${CHR} from Ensembl..."
if command -v curl >/dev/null 2>&1; then
    curl -L -o chr${CHR}.fa.gz "$FASTA_URL"
else
    wget -O chr${CHR}.fa.gz "$FASTA_URL"
fi

echo "Decompressing..."
gunzip -f chr${CHR}.fa.gz

# Coordinates for NF2 on GRCh38: chr22:29912274-29926621 (ENSEMBL gene ID ENSG00000131798)
# Note: coordinates may vary by release; we'll try a small region around expected NF2 locus
START=29912000
END=29927000

echo "Extracting region ${CHR}:${START}-${END}..."
awk 'BEGIN{p=0} /^>/{if(p==1) exit} {if(NR==1) {print; next}}' chr${CHR}.fa > /dev/null 2>&1 || true

# Use samtools faidx if available for region extraction, otherwise grep the sequence and slice
if command -v samtools >/dev/null 2>&1; then
    if [ ! -f chr${CHR}.fa.fai ]; then
        samtools faidx chr${CHR}.fa
    fi
    samtools faidx chr${CHR}.fa "${CHR}:${START}-${END}" > nf2_region.fa
else
    # Fallback: extract sequence into single-line and substring
    # Remove header and join lines
    tail -n +2 chr${CHR}.fa | tr -d '\n' > chr${CHR}.singleline
    echo ">NF2_region_${START}_${END}" > nf2_region.fa
    dd if=chr${CHR}.singleline bs=1 skip=$((START-1)) count=$((END-START+1)) 2>/dev/null | tr '[:lower:]' '[:upper:]' >> nf2_region.fa
fi

# Generate 21-nt sliding window candidates from extracted region
CAND_SIZE=21
OUT_FASTA="$OUT_DIR/nextflow_pipeline/candidates.fasta"

echo "Generating ${CAND_SIZE}-nt sliding-window candidates into ${OUT_FASTA}..."
python3 - <<PY
from pathlib import Path
import os
seqfile = Path('nf2_region.fa')
out = Path(os.environ.get('OUT_FASTA'))
text = seqfile.read_text().splitlines()
header = text[0] if text else '>NF2_region'
seq = ''.join(text[1:]).upper()
size = 21
with out.open('w') as fh:
    count = 0
    for i in range(0, max(0, len(seq)-size+1), 5):  # step 5 to reduce number of candidates
        subseq = seq[i:i+size]
        if set(subseq) <= set('ATCG') and len(subseq)==size:
            fh.write(f">candidate_{count}\n{subseq}\n")
            count += 1
print(f"Wrote {count} candidates to {out}")
PY

# Clean up
rm -rf "$WORK_DIR"

echo "Done. Test candidates available at: $OUT_FASTA"
chmod +r "$OUT_FASTA"
