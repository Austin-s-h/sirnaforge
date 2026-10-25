# GAA allele-aware design workflow

The earlier illustrative GAA command used an unverified transcript/variant pairing and flags that are not part of the current CLI. It is not a validated GAA run. The unqualified `c.428A>G` string does not identify one allele, and this guide does not substitute a guessed GAA variant. A real GAA run needs a pinned human GRCh38 annotation/sequence release and a user-selected variant whose transcript and genomic reference alleles have been checked against those files.

## Prepare pinned local inputs

Obtain the GENCODE GTF and transcript FASTA from the same release and GRCh38 assembly, and retain the release identifier plus SHA-256 checksums with the files. See {ref}`annotation_providers.md <obtain-a-pinned-release>` for official download and format references. Keep a source VCF containing the exact GRCh38 SNV after checking its chromosome, 1-based position, REF, ALT, and transcript mapping against that pinned reference. Do not derive genomic coordinates by guessing from a coding HGVS string.

Configure the files and run the current allele-aware workflow:

```bash
export SIRNAFORGE_GENCODE_GTF=/data/pinned/gencode.annotation.gtf.gz
export SIRNAFORGE_GENCODE_FASTA=/data/pinned/gencode.transcripts.fa.gz
export SIRNAFORGE_GENCODE_ASSEMBLY=GRCh38

uv run sirnaforge workflow GAA \
  --database gencode \
  --snp-file /data/pinned/user_reviewed_GAA_variants.vcf \
  --variant-mode both \
  --min-af 0.0 \
  --skip-off-targets \
  --output-dir results/gaa_variant_review
```

`--min-af 0.0` is explicit here so that a supplied variant with a reported AF below the default 0.01 is retained. When AF is absent from the VCF, the resolver warns and skips AF filtering; when AF is present, the configured threshold applies. `target` and `both` require at least one resolved variant. Only SNVs can enter allele-aware design; indels can resolve but design rejects them. Local annotation must contain exons for the selected transcripts, and the mapped GRCh38 reference base must exactly match the transcript sequence. See [variant_resolution.md](variant_resolution.md) for the full behavior and error conditions.

Review `results/gaa_variant_review/logs/resolved_variants.json`, the workflow manifest, and `logs/candidate_variants.json` before interpreting any candidate. Confirm the chosen variant is retained under the declared AF and ClinVar settings and that its source, assembly, transcript, coordinates, reference allele, and alternate allele are as intended. A `both` run labels reference and alternate overlapping candidates; it does not combine nearby variants into haplotypes. This workflow is a computational design example, not evidence of GAA knockdown, disease modification, therapeutic benefit, dose, or safety.

## Deterministic synthetic software demo

The repository's deterministic variant-resolution demo uses the bundled transcript FASTA and a synthetic/demo VCF, not a GAA variant or biological result. From the repository root, run the `avoid` or `target` example in {ref}`workflows.md <variant-targeting-documented-runs>` using `examples/sample_transcripts.fasta` and `examples/variant_demo.vcf`. These files exercise the current CLI and report format without making claims about a real GAA transcript or variant.

## GAA acceptance status

No GAA candidate output is included here. The prior issue example's `ENST00000362975` transcript and `c.428A>G` allele pairing were not established against a pinned current GAA annotation and GRCh38 reference in this audit. A real GAA example remains pending a reviewed, checksummed local GENCODE/RefSeq release and a fully resolved user-selected SNV. The `GAA` workflow above is a template; it becomes runnable only after those local paths and VCF contents are supplied and validated.
