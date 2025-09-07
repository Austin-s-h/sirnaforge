bwamem2_index - helper Nextflow module

This small Nextflow script builds a BWA‑MEM2 (or BWA) index for a given FASTA.

Usage example (build human GRCh38 index):

```bash
nextflow run nextflow_pipeline/modules/bwamem2_index/main.nf \
  --species human \
  --fasta_url "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
  --outdir nextflow_pipeline/indices
```

Notes and recommendations
- Use `bwa-mem2` in the execution environment for best performance.
- Building the full human genome index requires substantial disk and memory resources.
- After successful build, update `nextflow_pipeline/genomes.yaml` with the index prefix, e.g.:

```yaml
human:
  index_prefix: /absolute/path/to/nextflow_pipeline/indices/human_bwa_index/GRCh38.fa
  container: biocontainers/bwa-mem2:v2.2.1_cv1
```

You can also run the builder inside a container by ensuring the container image provides `bwa-mem2` and running Nextflow with container support (or by setting `params.container` and modifying the process to use Docker).
