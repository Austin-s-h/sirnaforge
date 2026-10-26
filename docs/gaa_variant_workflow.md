# GAA allele-aware design workflow

The earlier illustrative GAA command used an unverified transcript/variant pairing and flags that are not part of the current CLI. It is not a validated GAA run. The unqualified `c.428A>G` string does not identify a GAA allele. For a concrete software example, this guide uses the ClinGen expert-panel curated GAA SNV `NM_000152.5(GAA):c.1935C>A` (p.Asp645Glu), corresponding to `NC_000017.11:g.80112922C>A` on GRCh38 (ClinVar Variation ID 4029). This is a public reference variant selected to make the example reproducible; it is not a patient's genotype, treatment recommendation, or evidence that allele-specific siRNA design is clinically useful. See the [ClinVar record](https://www.ncbi.nlm.nih.gov/clinvar/variation/4029/) and [ClinGen Evidence Repository record](https://erepo.clinicalgenome.org/evrepo/ui/classification/e2861e0f-9989-4040-822a-dcdb1255a3b0) for its classification and transcript mappings.

The scientific basis and limits of the scoring and off-target screens are summarized in [scientific_sources.md](scientific_sources.md); those publications do not validate a GAA-specific candidate or predict its clinical effect.

## Prepare pinned local inputs

Obtain the GENCODE GTF and transcript FASTA from the same release and GRCh38 assembly, and retain the release identifier plus SHA-256 checksums with the files. See {ref}`annotation_providers.md <obtain-a-pinned-release>` for official download and format references. The following minimal VCF records the exact GRCh38 allele above and deliberately omits population-frequency fields. ClinGen maps this variant to GENCODE transcript `ENST00000302262.8`; check that the selected release contains a matching transcript version in both annotation and FASTA and that its sequence agrees with the VCF REF allele. The workflow searches compatible GAA transcripts; it does not convert or select the RefSeq accession `NM_000152.5` by itself.

```text
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr17>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr17	80112922	ClinVar_4029	C	A	.	PASS	.
```

Save this as `/data/pinned/gaa_c1935c_a.vcf`, then run the current allele-aware workflow:

```bash
export SIRNAFORGE_GENCODE_GTF=/data/pinned/gencode.annotation.gtf.gz
export SIRNAFORGE_GENCODE_FASTA=/data/pinned/gencode.transcripts.fa.gz
export SIRNAFORGE_GENCODE_ASSEMBLY=GRCh38

uv run sirnaforge workflow GAA \
  --database gencode \
  --snp-file /data/pinned/gaa_c1935c_a.vcf \
  --variant-mode both \
  --skip-off-targets \
  --output-dir results/gaa_variant_review
```

Because this example VCF omits allele frequency, the resolver warns and skips AF filtering. When AF is present, the configured `--min-af` threshold applies. `target` and `both` require at least one resolved variant. Only SNVs can enter allele-aware design; indels can resolve but design rejects them. Local annotation must contain exons for the selected transcripts, and the mapped GRCh38 reference base must exactly match the transcript sequence. See [variant_resolution.md](variant_resolution.md) for the full behavior and error conditions.

Review `results/gaa_variant_review/logs/resolved_variants.json`, the workflow manifest, and `logs/candidate_variants.json` before interpreting any candidate. Confirm the chosen variant is retained under the declared AF and ClinVar settings and that its source, assembly, transcript, coordinates, reference allele, and alternate allele are as intended. A `both` run labels reference and alternate overlapping candidates; it does not combine nearby variants into haplotypes. This workflow is a computational design example, not evidence of GAA knockdown, disease modification, therapeutic benefit, dose, or safety.

## Deterministic synthetic software demo

The repository's deterministic variant-resolution demo uses the bundled transcript FASTA and a synthetic/demo VCF, not a GAA variant or biological result. From the repository root, run the `avoid` or `target` example in {ref}`workflows.md <variant-targeting-documented-runs>` using `examples/sample_transcripts.fasta` and `examples/variant_demo.vcf`. These files exercise the current CLI and report format without making claims about a real GAA transcript or variant.

## GAA acceptance status

No GAA candidate output is included here because the repository does not bundle a pinned GENCODE/RefSeq release and transcript FASTA. The example now specifies a sourced allele and coordinates, replacing the incorrect prior `c.428A>G` pairing. It is reproducible once matching local GRCh38 inputs are obtained and checksummed, and the selected GAA transcript's presence and sequence are confirmed. Any generated candidates remain computational hypotheses and require independent biological testing.
