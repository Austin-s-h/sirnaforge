Nextflow genome index instructions

This pipeline requires genome indexes (BWA/BWA-mem2) to perform off-target analysis.

Options to provide indexes:

1) Build repository-local indexes (recommended for reproducible runs)
   - Place genome FASTA files somewhere on disk (e.g., downloaded from Ensembl/NCBI).
   - Run the helper script in this repo to build BWA indexes and copy a .fasta sentinel:

     ./nextflow_pipeline/scripts/build_index.sh /path/to/GRCh38.fa nextflow_pipeline/genomes_indices/human/GRCh38/GRCh38

   - The pipeline will look for either an index directory or a <index_prefix>.fasta sentinel file.

2) Use system/global indexes
   - Edit `nextflow_pipeline/genomes.yaml` and set `index_prefix` to absolute paths where your indexes live.
   - Or pass parameters at runtime using `--genome_indices` with comma separated pairs, for example:

     nextflow run nextflow_pipeline/main.nf --input nextflow_pipeline/candidates.fasta \
       --outdir results --genome_indices "human:/abs/path/to/GRCh38,rat:/abs/path/to/Rnor6"

3) Containerized execution
   - `genomes.yaml` includes an optional `container:` image per species. When the pipeline runs the off-target wrapper
     inside a container it will mount the working directory. Ensure your index paths are accessible inside the container
     (absolute paths on the host will be mounted as-is via `-v $PWD:/work`).

Notes and troubleshooting:
 - The Nextflow process checks for either `${index_prefix}.fasta`, `${index_prefix}.fa` or a directory at `${index_prefix}`.
 - If you get "Warning: Index ... not found" the pipeline will still succeed but produce empty results for that species.
 - To test quickly, create a small mock FASTA and build a tiny index using the helper script.

S3 auto-download support
-----------------------
The pipeline now supports automatically downloading index prefixes from S3 when given an `s3://` path.

Usage:
 - Pass an `s3://` URL in `nextflow_pipeline/genomes.yaml` or via `--genome_indices` for the species you want.
 - Enable the download behavior with the Nextflow flag `--download_indexes true` (or set the env var `DOWNLOAD_INDEXES=true`).

Example:

```bash
nextflow run nextflow_pipeline/main.nf \
  --input nextflow_pipeline/candidates.fasta \
  --outdir results \
  --genome_indices "human:s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/" \
  --download_indexes true
```

Notes and requirements:
 - The worker environment must have the AWS CLI available in PATH for downloads to work.
 - Downloads are performed with `aws s3 sync --no-sign-request` so public buckets (like iGenomes) can be accessed without credentials.
 - Downloaded content is placed under `./.downloaded_indexes/<species>` inside the working directory and used for the index checks and downstream wrapper.
 - Large genome downloads will use disk space and network bandwidth; ensure sufficient resources.

Contact: See repository README for full pipeline instructions.
