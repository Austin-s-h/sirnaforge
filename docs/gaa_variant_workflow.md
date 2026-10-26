# GAA allele-aware design workflow

The earlier illustrative GAA command used an unverified transcript/variant pairing and flags that are not part of the current CLI. It is not a validated GAA run. The unqualified `c.428A>G` string does not identify a GAA allele. For a concrete software example, this guide uses the ClinGen expert-panel curated GAA SNV `NM_000152.5(GAA):c.1935C>A` (p.Asp645Glu), corresponding to `NC_000017.11:g.80112922C>A` on GRCh38 (ClinVar Variation ID 4029). This is a public reference variant selected to make the example reproducible; it is not a patient's genotype, treatment recommendation, or evidence that allele-specific siRNA design is clinically useful. See the [ClinVar record](https://www.ncbi.nlm.nih.gov/clinvar/variation/4029/) and [ClinGen Evidence Repository record](https://erepo.clinicalgenome.org/evrepo/ui/classification/e2861e0f-9989-4040-822a-dcdb1255a3b0) for its classification and transcript mappings.

The scientific basis and limits of the scoring and off-target screens are summarized in [scientific_sources.md](scientific_sources.md); those publications do not validate a GAA-specific candidate or predict its clinical effect.

## Prepare pinned local inputs

This example was run with GENCODE human release 49 (GRCh38.p14), using the [official release files](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/) below. The SHA-256 hashes were computed from the downloaded compressed files and the MD5 values match GENCODE's `MD5SUMS` manifest.

| Input | Official MD5 | SHA-256 |
| --- | --- | --- |
| `gencode.v49.annotation.gtf.gz` | `0ef4a024ea2d35b1b88c12447b0b70b9` | `d6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4` |
| `gencode.v49.transcripts.fa.gz` | `ce72162016c56f64918a7dde76535ddb` | `08ee7ee86514061de6fd62dc643c2fe962fbc393080589ebc491bd11e22b119e` |

See {ref}`annotation_providers.md <obtain-a-pinned-release>` for GENCODE format details. The following minimal VCF records the exact GRCh38 allele above and deliberately omits population-frequency fields. ClinGen maps this variant to GENCODE transcript `ENST00000302262.8`; confirm that the matching version is present in both files and that its sequence agrees with the VCF REF allele. The workflow searches compatible GAA transcripts; it does not convert or select the RefSeq accession `NM_000152.5` by itself.

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

The command completed against the pinned release 49 files above, resolving the one VCF SNV and retrieving 14 protein-coding GAA transcripts. It generated 11,120 general siRNA candidates but no candidate overlapping the variant under the default filters (`candidate_variants.json` reported zero variant-annotated candidates). No candidate output is included because this run did not produce an allele-specific guide. The result demonstrates input and workflow behavior only; it is not evidence of knockdown or clinical utility. The full GENCODE GTF is scanned during transcript lookup and annotation enrichment, so this local run took 8m47s; combine annotation requests and plan for the scan time when using full releases.
