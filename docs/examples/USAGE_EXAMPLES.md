# siRNAforge Usage Examples

## Basic Examples

```bash
# Complete workflow
uv run sirnaforge workflow TP53 --output-dir results

# Design from FASTA
uv run sirnaforge design transcripts.fasta -o results.tsv

# Search for transcripts
uv run sirnaforge search BRCA1 --all --verbose
```
