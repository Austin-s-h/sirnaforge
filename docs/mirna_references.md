# Pinned local miRNA references

Use a local mature-FASTA release when miRBase is unavailable or when a run must use
fixed reference bytes. This works with both `workflow` and `offtarget`, including
Nextflow tasks in containers. It does not require the miRBase service to recover.

Create a manifest beside your reference file:

```json
{
  "schema_version": 1,
  "references": [
    {
      "source": "mirbase",
      "species": "human",
      "release": "YOUR_VERIFIED_RELEASE",
      "sha256": "REPLACE_WITH_SHA256_OF_THE_FILE",
      "path": "mature.fa"
    }
  ]
}
```

Calculate the digest with `sha256sum mature.fa`. Record the actual upstream release;
the release field is an explicit declaration, not a version discovered from FASTA
headers. An optional `original_url` records where the release originated. Paths are
resolved relative to the manifest. Supported pinned sources are `mirbase` and
`mirgenedb`; use canonical species names or their supported aliases. miRBase
records need their original species prefixes, such as `hsa-` for human.

```bash
sirnaforge workflow TP53 --mirna-db mirbase --mirna-species human \
  --mirna-reference-manifest /refs/references.json --output-dir results/TP53

sirnaforge offtarget --input-candidates-fasta guides.fa \
  --mirna-db mirbase --mirna-species human \
  --mirna-reference-manifest /refs/references.json --output-dir results/screen
```

These commands only make the selected miRNA reference local. Gene retrieval,
transcriptome references, orthology and other configured services may still need
their own local inputs for a fully offline run.

For `--mirna-db mirgenedb+mirbase`, supply entries for both sources and every
requested species to avoid live acquisition entirely. Missing entries retain the
normal acquisition path for that source/species. Invalid explicit pins fail; they
never fall back to another release or silently reduce a union to one source.

The workflow stages local FASTAs for Nextflow. The task verifies the source file's
SHA-256, validates mature sequences and filters species before scanning. Pinned
references are immutable and do not expire under the live-download cache TTL.
Failed refreshes preserve previous cache bytes and metadata; a failed refresh still
returns a failure rather than silently claiming success with an older reference.

Live acquisition has bounded attempts and socket timeouts. Python callers can set
`MiRNADatabaseManager(download_timeout=30, download_attempts=2,
retry_backoff_seconds=0.25)`. This is a per-socket timeout, not a strict total
wall-clock deadline for a server that continually trickles bytes.

The published `off_target/results/mirna/mirna_reference_evidence.json` records actual
task acquisition, release/checksums, species, screening outcome and submitted guide
membership. The manifest carries this evidence into the existing HTML report.
Unavailable or unproven screening is distinct from a completed zero-hit result.
No archive or mirror is trusted automatically; obtain and verify the release and
its redistribution conditions before provisioning it.
