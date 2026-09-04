# ZFN Nextflow Bridge: Sanity Validation

> ## ⚠️ EXPERIMENTAL — this check is not authoritative
>
> The ZFN arm ships **experimental** in 0.6.0 with known unfixed defects tracked in [#82](https://github.com/Austin-s-h/sirnaforge/issues/82). **Do not use ZFN output for any decision without independent
> validation.**
>
> Two caveats apply to the run recorded below. First, the fixtures were built around the published
> CCR5 half-site pair, which does not match its own on-target site under the default strand-pairing
> rule, so the site counts are not evidence about real genomic behaviour. Second, the Nextflow ZFN
> arm exercised here is not reachable from either CLI entry point — the workflow always runs the
> in-process path — and the aggregator emits a narrower CSV schema and a different ranking than the
> direct route.

Date: 2026-03-01

This note records a focused sanity check of the ZFN Nextflow bridge refactor where Nextflow orchestrates processes and Python (via CLI entrypoints) owns the higher-order logic.

## Scope

Validated hidden internal CLI commands used by Nextflow modules:

- `sirnaforge _internal zfn-make-shards`
- `sirnaforge _internal zfn-search-shard`
- `sirnaforge _internal zfn-aggregate-shards`

All commands were executed with `uv run`.

## uv invocation guidance

- From repository root: `uv run sirnaforge ...`
- From outside the repository (for example, staged Nextflow work dirs under `/tmp`):
  `uv run --project /home/hovland/sirnaforge/sirnaforge sirnaforge ...`

Without `--project` outside the workspace, `uv` cannot resolve the `sirnaforge` entrypoint.

## Fixtures

Temporary validation artifacts were generated under `/tmp/sf_zfn_bridge_check`:

- `multi.fa` (2 contigs, both containing one canonical ZFN site)
- `single.fa` (1 contig, one canonical ZFN site)

Settings used:

- left half-site: `GTCATCCTCATC`
- right half-site: `AAACTGCAAAAG`
- spacer lengths: `5,6`
- `max_mismatches=2`
- `shard_chunk_mb=0.00004` (nominal 40 bp)
- `shard_overlap_bp=0` (auto-raised by safety logic)

## Results

### 1) Manifest generation

- Multi-contig (`multi_shards.tsv`): **6 shards** across 2 contigs.
- Single-contig (`single_shards.tsv`): **1 shard** (`chr3:1-121`).

This confirms default behavior: sharding remains enabled generically, but chunking is not applied when the selected input resolves to only one contig.

### 2) Shard-level search

- Generated shard-scoped outputs for all 6 multi-contig manifest rows:
  - `zfn_offtarget_sites_<shard_id>.csv`
  - `zfn_candidate_summary_<shard_id>.json`

### 3) Aggregation

Running aggregate on all shard CSVs produced:

- `zfn_candidate_summary.json` summary:
  - `off_target_sites: 2`
  - `shards: 6`
- `zfn_offtarget_sites.csv` rows: **2** (one retained site per contig)

This confirms global deduplication and deterministic merge behavior over overlapping shards.

## Unit verification

Targeted suite after refactor:

```bash
uv run python -m pytest -n 0 tests/unit/test_zfn_search_and_design.py tests/unit/test_zfn_workflow_integration.py
```

Result: **51 passed**.

## Notes for maintainers

- The Nextflow ZFN modules should remain thin wrappers over the internal CLI commands.
- Keep domain logic and typing in Python modules (`sirnaforge.zfn.*`), and use Nextflow for channel/process orchestration only.
