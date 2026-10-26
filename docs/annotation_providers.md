# Transcript annotation providers

Structural annotations can come from Ensembl REST or local GENCODE and NCBI RefSeq files.
Local providers work without Ensembl access and keep their original transcript/gene IDs,
1-based inclusive coordinates, exon/CDS structure, and source identity.

Gene/sequence retrieval and structural enrichment are separate interfaces. RefSeq gene/sequence
search uses NCBI E-utilities and the specific gene-to-RefSeq-RNA links, excluding genomic records.
Requests are paced per client and retry bounded HTTP 429 responses; see the
[NCBI usage guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/). GENCODE gene search now reads a local annotation plus its
matching transcript FASTA. Neither provider changes an Ensembl accession into a RefSeq accession.

(obtain-a-pinned-release)=
## Obtain a pinned release

For GENCODE, obtain the annotation GTF and transcript FASTA from the same release and assembly
using the [official human release catalog](https://www.gencodegenes.org/human/).
GENCODE FASTA headers contain pipe-delimited, versioned transcript identifiers.
The parser preserves repeated GTF tags; explicit `Ensembl_canonical` or `MANE_Select` tags establish
canonical status. See [GENCODE's format documentation](https://www.gencodegenes.org/pages/data_format.html).

For RefSeq, obtain a genomic GFF3 or GTF from a pinned NCBI annotation release. The
[NCBI GFF3 specification](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)
describes the stable transcript accessions, GeneID cross-references, and parent relationships.
NCBI sequence-region labels such as `NC_000017.11` are preserved; `chr17` is not treated as an
alias for a contig accession.

Plain files and gzip files are supported. Retain downloaded release metadata beside the files.
GENCODE gene search currently targets human; local structural annotation accepts an explicitly
declared species. It does not infer species or assembly from a filename.

## Query structural annotations directly

```bash
uv run sirnaforge annotate ENST00000269305 \
  --provider gencode \
  --annotation-file /data/gencode.annotation.gtf.gz \
  --assembly GRCh38 \
  --species human \
  --output annotations.json

uv run sirnaforge annotate NM_000546 \
  --provider refseq \
  --annotation-file /data/refseq_genomic.gff.gz \
  --assembly GRCh38.p14 \
  --output refseq_annotations.json

uv run sirnaforge annotate NC_000017.11:7661779-7687550 \
  --provider refseq \
  --annotation-file /data/refseq_genomic.gff.gz \
  --assembly GRCh38.p14 \
  --regions
```

Use the assembly actually associated with the chosen file; the examples are not release defaults.
NCBI header assembly names/accessions are checked when present. Queries against a different
configured species or assembly fail. Files without such headers rely on the explicit declaration.

Output JSON includes the source URI, SHA-256 of the exact file bytes, assembly, species,
resolved transcript models, and unresolved queries. Exit code 1 means some queries were unresolved;
exit code 2 means invalid configuration or file content. Exons and CDS intervals retain their
individual spans and are sorted in transcript order, including negative strands.
Versionless IDs match a source accession's version; explicit versions require an exact match.
No liftover, cross-provider ID conversion, or RefSeq contig-name guessing is performed.

## Use GENCODE for gene search and design

```bash
export SIRNAFORGE_GENCODE_GTF=/data/gencode.annotation.gtf.gz
export SIRNAFORGE_GENCODE_FASTA=/data/gencode.transcripts.fa.gz
export SIRNAFORGE_GENCODE_ASSEMBLY=GRCh38

uv run sirnaforge search TP53 --database gencode --no-fallback
uv run sirnaforge workflow TP53 --database gencode --skip-off-targets \
  --output-dir gencode_tp53_design
```

The selected database is now tried first when search fallback is enabled. `--no-fallback` makes
a GENCODE search fully local. Workflow gene retrieval uses the requested database directly.
Explicit-version GTF and FASTA accessions must agree; missing sequences cause an error rather
than returning a partially populated design input. A symbol mapping to multiple genes is refused;
use a stable gene ID to disambiguate. Missing canonical tags are not replaced with a guessed
longest transcript. `get_sequence` supports cDNA from the transcript FASTA, not protein or CDS
extraction.

`--skip-off-targets` makes this example a design-only run. A full screen still needs separately
configured transcriptome and miRNA references. Switching gene/annotation provider does not replace
the screening references, orthology service, variant services, or every other Ensembl dependency.

## Select workflow structural enrichment

GENCODE gene workflows automatically use their configured local annotation. For another gene
retrieval database or a FASTA input, select structural enrichment explicitly:

```bash
export SIRNAFORGE_ANNOTATION_PROVIDER=refseq
export SIRNAFORGE_REFSEQ_ANNOTATION=/data/refseq_genomic.gff.gz
export SIRNAFORGE_REFSEQ_ASSEMBLY=GRCh38.p14
export SIRNAFORGE_REFSEQ_SPECIES=human

uv run sirnaforge workflow TP53 --database refseq --skip-off-targets \
  --output-dir refseq_tp53_design
```

This example retrieves sequences from NCBI, then enriches matching RefSeq accessions locally.
For a fully offline RefSeq run, provide a selected local transcript FASTA via `--input-fasta`.
The annotation file must describe those same accessions and assembly.

| Setting | Meaning |
| --- | --- |
| `SIRNAFORGE_ANNOTATION_PROVIDER` | `ensembl`, `gencode`, `refseq`, or `none`. Default is Ensembl except for GENCODE gene workflows. |
| `SIRNAFORGE_GENCODE_GTF` | Local GENCODE structural annotation; `SIRNAFORGE_GENCODE_ANNOTATION` also selects an annotation file for structural enrichment. |
| `SIRNAFORGE_GENCODE_FASTA` | Required for GENCODE sequence retrieval, optional for metadata-only queries. |
| `SIRNAFORGE_GENCODE_ASSEMBLY` | Explicit assembly for GENCODE files. |
| `SIRNAFORGE_GENCODE_SPECIES` | Structural annotation species, default human. Gene retrieval remains human. |
| `SIRNAFORGE_REFSEQ_ANNOTATION` | Local NCBI GTF/GFF3; `SIRNAFORGE_REFSEQ_GTF` is an alias. |
| `SIRNAFORGE_REFSEQ_ASSEMBLY` | Explicit assembly for the NCBI annotation. |
| `SIRNAFORGE_REFSEQ_SPECIES` | Structural annotation species, default human. |

Set `SIRNAFORGE_ANNOTATION_PROVIDER=none` to disable optional enrichment.
No automatic cross-provider structural fallback is attempted: mixing assemblies or transcript
namespaces would make annotations unreliable. A missing local accession is recorded as unresolved.
The workflow summary and manifest identify the selected provider and file checksum only when
annotation was actually attempted. Optional enrichment is limited to the first ten transcripts;
use `sirnaforge annotate` for a full requested set.

## Python API

```python
from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.base import GencodeClient
from sirnaforge.data.local_annotation import LocalTranscriptModelClient

client = GencodeClient(
    gtf_path="/data/gencode.annotation.gtf.gz",
    fasta_path="/data/gencode.transcripts.fa.gz",
    assembly="GRCh38",
)
gene, transcripts = await client.search_gene("TP53")

annotations = LocalTranscriptModelClient(
    "/data/refseq_genomic.gff.gz",
    provider="refseq",
    species="human",
    assembly="GRCh38.p14",
)
bundle = await annotations.fetch_by_ids(
    ["NM_000546"], species="human",
    reference=ReferenceChoice.explicit("GRCh38.p14"),
)
```

Queries stream annotation files in up to three passes, retaining gene records and selected models.
They do not load every exon in a full release into memory, but each uncached query scans the source.
Large files take time to scan; combine IDs into one request. File reading runs off the async event
loop. No new package dependency or background download is required.
