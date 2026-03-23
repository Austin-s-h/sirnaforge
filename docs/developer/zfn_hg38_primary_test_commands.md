# ZFN Hg38 Primary Validation Commands

These commands are the manual validation reference for real hg38 primary ZFN runs. They are the heavyweight reruns behind the backend conclusions summarized in [zfn_backend_tuning.md](zfn_backend_tuning.md).

They use:

- the real hg38 primary reference: `ensembl_human_hg38_primary`
- the Ensembl human GTF annotation
- the CCR5 ZFN half-sites already used elsewhere in this repo
- the `workflow --design-mode zfn` entrypoint

## Annotation URL



## Recommended output root

Use a scratch directory for manual runs so the notebook tree stays focused on curated artifacts.

```bash
mkdir -p /tmp/sirnaforge-zfn-real-runs
```

## Measured Backend Order

The backend order from the tuning work was:

1. `pyahocorasick`
2. `fm_index`
3. `exhaustive_python`

That order reflects the current qualitative findings:

- `pyahocorasick` is the fastest practical backend across the synthetic sweep and the real hg38 chr3 slice sweep
- `fm_index` is the alternate backend when persisted-index reuse is part of the workflow, but its preparation cost is not the best default for one-off runs
- `exhaustive_python` remains the strict baseline rather than the practical first pass

Use this runbook to revalidate those conclusions against the current code and dependency set.

## Real Run on hg38 Primary with `pyahocorasick`

```bash
export HG38_ANNOTATION_URL="https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"

uv run sirnaforge workflow CCR5_ZFN_HG38_PRIMARY_PYAHO \
  --design-mode zfn \
  --zfn-left-half-site GTCATCCTCATC \
  --zfn-right-half-site AAACTGCAAAAG \
  --zfn-search-space ensembl_human_hg38_primary \
  --zfn-search-backend pyahocorasick \
  --zfn-spacer-lengths 5,6 \
  --zfn-max-mismatches 2 \
  --zfn-algorithm zfn_v2 \
  --zfn-annotation "$HG38_ANNOTATION_URL" \
  --cores 12 \
  --output-dir notebooks/zfn_experiment_runs/CCR5_ZFN_HG38_PRIMARY_PYAHO

```

## Real Run on hg38 Primary with `exhaustive_python`

```bash
uv run sirnaforge workflow CCR5_ZFN_HG38_PRIMARY_EXHAUSTIVE \
  --design-mode zfn \
  --zfn-left-half-site GTCATCCTCATC \
  --zfn-right-half-site AAACTGCAAAAG \
  --zfn-search-space ensembl_human_hg38_primary \
  --zfn-search-backend exhaustive_python \
  --zfn-spacer-lengths 5,6 \
  --zfn-max-mismatches 2 \
  --zfn-algorithm zfn_v2 \
  --zfn-annotation "$HG38_ANNOTATION_URL" \
  --cores 12 \
  --output-dir /tmp/sirnaforge-zfn-real-runs/CCR5_ZFN_HG38_PRIMARY_EXHAUSTIVE
```

## Expected Artifacts

Each run writes these files under its output directory:

- `sirnaforge/zfn_candidate_summary.json`
- `sirnaforge/zfn_offtarget_sites.csv`
- `logs/workflow_summary.json`
- `logs/sirnaforge.log`

## Quick Checks

Inspect the summary JSON:

```bash
jq '{search_backend, search_space_reference, annotation_source, predicted_sites_total, predicted_sites_exonic, total_workflow_time_s}' \
  /tmp/sirnaforge-zfn-real-runs/CCR5_ZFN_HG38_PRIMARY_PYAHO/logs/workflow_summary.json
```

Compare row counts across backends:

```text
python - <<PY
import csv
from pathlib import Path

base = Path('/tmp/sirnaforge-zfn-real-runs')
runs = [
    'CCR5_ZFN_HG38_PRIMARY_PYAHO',
    'CCR5_ZFN_HG38_PRIMARY_EXHAUSTIVE',
]
for run in runs:
    csv_path = base / run / 'sirnaforge' / 'zfn_offtarget_sites.csv'
    with csv_path.open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    print(f'{run}: {len(rows)} rows')
PY
```
