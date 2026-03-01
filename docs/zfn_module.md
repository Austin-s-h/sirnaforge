# ZFN Module Guide

This guide is the user-facing entry point for ZFN design and off-target analysis in siRNAforge.

## What the ZFN module does

- Scores a user-provided left/right half-site pair.
- Searches a genomic FASTA for candidate off-target pairings.
- Ranks results and writes candidate + site reports.
- Supports scalable shard-based execution for large references.

## CLI usage

### Required inputs for `--design-mode zfn`

When using ZFN mode, both half-sites are required:

- `--zfn-left-half-site`
- `--zfn-right-half-site`

The workflow exits with a validation error if either flag is missing.

### ZFN workflow mode

```bash
sirnaforge workflow CCR5_ZFN_RUN \
  --design-mode zfn \
  --zfn-left-half-site GTCATCCTCATC \
  --zfn-right-half-site AAACTGCAAAAG \
  --zfn-search-space ensembl_human_hg38_primary \
  --zfn-spacer-lengths 5,6 \
  --zfn-max-mismatches 2 \
  --zfn-algorithm zfn_v2 \
  --output-dir zfn_output
```

Outputs are written under:

- `zfn_output/sirnaforge/zfn_candidate_summary.json`
- `zfn_output/sirnaforge/zfn_offtarget_sites.csv`
- `zfn_output/logs/workflow_summary.json`

### Common ZFN options

- `--zfn-algorithm`: `homology`, `conserved_g`, `zfn_v2`
- `--zfn-dimer-mode`: `heterodimer_only` (default) or `include_homodimers`
- `--zfn-spacer-lengths`: comma-separated allowed spacers (for example `5,6`)
- `--zfn-max-mismatches`: mismatch cap per half-site during exhaustive search
- `--zfn-annotation`: optional GTF/GFF annotation for region classification

For larger genomes or advanced runtime control, see the workflow-level environment toggles documented in [Workflows](workflows.md).

## Python API quickstart

```python
from sirnaforge.models.zfn import ZFNDesignParameters
from sirnaforge.zfn.design import ZFNDesigner

params = ZFNDesignParameters(
    left_half_site="GTCATCCTCATC",
    right_half_site="AAACTGCAAAAG",
    search_space_reference="ensembl_human_hg38_primary",
)

result = ZFNDesigner().evaluate_pair(params=params)
print(result.get_summary())
```

## API map for ZFN users

See these sections in [API Reference](api_reference.rst):

- **ZFN Models** (`sirnaforge.models.zfn`) for typed constraints and result structures.
- **ZFN Design** (`sirnaforge.zfn.design`) for pair evaluation orchestration.
- **ZFN Search** (`sirnaforge.zfn.search`) for exhaustive/site search and sharding behavior.
- **ZFN Ranking** (`sirnaforge.zfn.rank`) for deterministic ranking utilities.
- **Workflow Orchestration** (`sirnaforge.workflow`) for `run_sirna_workflow(...)` with `design_mode="zfn"`.

## Related documentation

- [Workflows](workflows.md)
- [CLI Reference](cli_reference.md)
- [ZFN Ranking & Scoring](zfn_ranking.md)
