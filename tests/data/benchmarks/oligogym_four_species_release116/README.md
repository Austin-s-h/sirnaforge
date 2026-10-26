# Release 116 paired-target mappings

These tables preserve 4,113 human assay observations. Human matches are native
sequence mappings; mouse, rat, and macaque matches are sequence projections, with
no efficacy label transfer or claim of measured activity in those species.
Unsupported duplex geometries remain review-only.

The CSV artifacts are stored as deterministic gzip files to keep the checkout
small. The human mapping gzip stream is split into numbered parts below the
repository per-file size limit; concatenate the parts in lexical order before
decompression. Manifests record the complete gzip-stream SHA-256 as `sha256`
and the decompressed CSV SHA-256 as `uncompressed_sha256`, together with Ensembl
assembly/release, exact transcript sites, and mapping status. Reference acquisition records are
preserved in `references/`; their absolute FASTA paths describe the generation
machine, not a dependency on an installed file at that location.

On 2026-10-25, all mapped sites were independently checked against the SHA-256
verified full release FASTAs: human 66,204; macaque 5,150; mouse 18,624; rat 2,782.
Each site matched both its declared transcript interval and its oriented query.
This verifies sequence mapping only, not biological efficacy or scoring defaults.

To regenerate, obtain the OligoGym extracts using the repository's
`scripts/fetch_public_benchmarks.py` workflow and compare their checksums to the
root manifest. For each of human, mouse, rat, and macaque, run:

```sh
uv run python scripts/fetch_native_reference_release.py \
  --release 116 --species human --output work/release116
```

Repeat with the other species, then generate the tables:

```sh
uv run python scripts/map_oligo_benchmarks_to_ensembl.py \
  --source-dir /path/to/oligogym --output-dir work/mapped116 \
  --species-reference-manifest work/release116/human_manifest.json \
  --species-reference-manifest work/release116/mouse_manifest.json \
  --species-reference-manifest work/release116/rat_manifest.json \
  --species-reference-manifest work/release116/macaque_manifest.json

uv run python scripts/package_four_species_mappings.py \
  --artifact-dir work/mapped116
```

The reference FASTAs and upstream assay extracts are not bundled here. Existing
upstream provenance and redistribution limitations still apply; these generated
mapping tables do not establish an independent held-out validation cohort.

Decompress `species/human/observations.csv.gz` before passing it to the native
benchmark importer. Reconstruct the human mapping table by concatenating the parts:
`cat species/human/mappings.csv.gz.part-* | gzip -dc > work/human-mappings.csv`.
The same command works on an unsplit mapping file by passing that single `.gz` file.
