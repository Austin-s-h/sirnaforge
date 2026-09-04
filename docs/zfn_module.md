# ZFN Module Guide

> ## ⚠️ EXPERIMENTAL — results are not decision-grade
>
> The ZFN arm ships **experimental** in 0.6.0 with known unfixed defects, tracked in the ZFN
> experimental-status issue. **Do not use ZFN output for any decision without independent
> validation.** The known defects affect:
>
> - **Half-site orientation.** With the default `require_opposite_strands=True`, the CCR5
>   half-site pair used in the examples below (`GTCATCCTCATC` / `AAACTGCAAAAG`) does not match
>   its own on-target site, because both sequences occur on the hg38 plus strand. The
>   undocumented convention is that `right_half_site` must be given as the reverse complement
>   of the published text (`CTTTTGCAGTTT`), and there is no CLI escape hatch.
> - **FokI seed region and polarity weighting** on the right half-site, which guard the end of
>   the half-site farthest from FokI rather than the nearest.
> - **Off-target region classification**, which can report a site inside a large containing gene
>   as `intergenic` and therefore undercount exonic/promoter hits used by the pass/fail filters.
>
> Because of the half-site convention issue, the ZFN validation runs recorded in this
> documentation set are **not authoritative** and must be re-run once the convention is fixed.

This guide is the user-facing entry point for ZFN design and off-target analysis in siRNAforge.

## Current implementation status (PR #61: dev → master)

### Completed

- ✅ ZFN support is available through the dedicated `sirnaforge zfn ...` command.
- ✅ Typed models and JSON-serializable workflow artifacts are implemented (`sirnaforge.models.zfn` and `zfn_candidate_summary.json` output).
- ✅ Output/reporting follows existing workflow conventions, including:
  - `sirnaforge/zfn_candidate_summary.json`
  - `sirnaforge/zfn_offtarget_sites.csv`
  - `logs/workflow_summary.json`
- ✅ Workflow provenance is captured in `workflow_summary.json`.
- ⚠️ Nextflow-backed ZFN execution path is present (including internal shard commands), but it is
  **not reachable** from either CLI entry point: every run uses the in-process path. See the ZFN
  experimental-status issue.
- ✅ Automated tests exist for ZFN workflow, search/design behavior, and integration paths.
- ✅ User/developer documentation for ZFN workflow usage is present.

### Not implemented in this PR branch

- ❌ Distinct `repeat_motif` mode is not implemented; current implementation focuses on nuclease-style half-site pair workflow.

## What the ZFN module does

- Scores a user-provided left/right half-site pair.
- Searches a genomic FASTA for candidate off-target pairings.
- Ranks results and writes candidate + site reports.
- Supports scalable shard-based execution for large references.

## CLI usage

### Required inputs for `sirnaforge zfn`

When using the ZFN command, both half-sites are required:

- `--zfn-left-half-site`
- `--zfn-right-half-site`

The workflow exits with a validation error if either flag is missing.

### ZFN command

```bash
sirnaforge zfn \
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

`zfn_candidate_summary.json` uses a stable envelope with:

- `schema_version` (currently `zfn_candidate_summary.v1`)
- `search_contract` (canonical half-sites, spacer list, mismatch budgets, dimer mode, algorithm)
- `candidates` and `summary`

`workflow_summary.json` records ZFN provenance fields including `algorithm`, `dimer_mode`, `spacer_lengths`, mismatch thresholds, and resolved search-space input (`search_space_reference`/`search_space_fasta`).

### Common ZFN options

- `--zfn-algorithm`: `homology`, `conserved_g`, `zfn_v2`
- `--zfn-dimer-mode`: `heterodimer_only` (default) or `include_homodimers`
- `--zfn-spacer-lengths`: comma-separated allowed spacers (for example `5,6`)
- `--zfn-max-mismatches`: mismatch cap per half-site during exhaustive search
- `--zfn-annotation`: optional GTF/GFF annotation for region classification

### Search-space references

`--zfn-search-space` accepts a local FASTA path, an HTTP(S) URL, or one of the built-in
Ensembl genome keys (auto-downloaded and cached on first use):

| Key                                    | Species        | Assembly | Ensembl DNA file |
| -------------------------------------- | -------------- | -------- | ---------------- |
| `ensembl_human_hg38_primary` (default) | human          | GRCh38   | primary assembly |
| `ensembl_mouse_grcm39_primary`         | mouse          | GRCm39   | primary assembly |
| `ensembl_rat_grcr8_toplevel`           | rat            | GRCr8    | toplevel         |
| `ensembl_macaque_mmul10_toplevel`      | rhesus macaque | Mmul_10  | toplevel         |

Rat and macaque use Ensembl's `toplevel` DNA file because those assemblies do not publish a
separate `primary_assembly` FASTA. Genomic DNA references are large (~700 MB–3 GB); the first
run downloads and caches them under `~/.cache/sirnaforge/genomes/`.

### Backend choice

The backend tuning work settled on this operational order:

1. `pyahocorasick` for the practical first run on large genomic references
2. `fm_index` only for repeated persisted-index reuse workflows (experimental on large references)
3. `exhaustive_python` as the strict baseline and fallback path

See [Developer Documentation](developer/index.rst) for the backend tuning summary and the hg38 primary
validation commands used to reach that conclusion. Those runs used the published CCR5 half-site pair
directly, so they matched no on-target site; treat the backend ordering as a runtime observation only,
not as a correctness result.

For larger genomes or advanced runtime control, see the workflow-level environment toggles documented in [Workflows](workflows.md).

Legacy note: `sirnaforge workflow --design-mode zfn` remains a compatibility path, but `sirnaforge zfn` is the primary entrypoint.

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
