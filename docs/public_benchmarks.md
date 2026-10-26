# Public benchmark catalog and regeneration guide

> **ZFN is EXPERIMENTAL** and requires independent validation before use in decisions.
> Known defects are tracked in [#82](https://github.com/Austin-s-h/sirnaforge/issues/82).

This is the entry point for public benchmarking in siRNAforge. Keep provenance, curated extracts,
regression fixtures, scripts, and compact historical findings in Git. Fetch upstream downloads and
regenerate large references, candidate tables, reports, and figures locally.

## Catalog

| Collection | Retained material | Purpose and limits |
| --- | --- | --- |
| Huesken subset | `tests/unit/data/sirna_efficacy_subset.csv`, 180 records | Scoring regression fixture with a declared development/held-out split; values not verified against the non-open primary paper. |
| OligoGym synthetic adapters | `tests/data/benchmarks/oligogym/`, 4,113 observations | Design enumeration and filter regression in synthetic flanks; no native accessibility evidence. |
| OligoGym native design | `tests/data/benchmarks/oligogym_native_design/` mappings, inputs, trace tables, orthology, and manifests | Experimental native transcript-site mapping for three runnable panels; representative four-species screening. |
| 2024/2026 ML-review collections | `tests/data/benchmarks/ml_reviews/`, eight collections, 8,287 assay records, 4,091 sequence clusters | Partial recovery proposal with heterogeneous labels, chemistry, coverage gaps, and overlap audit; not independent validation. |
| CCR5 ZFN / PROGNOS extract | `tests/unit/data/zfn/` | Visible supplementary rows for coordinate, annotation, and ranking regression. Published PROGNOS recovery is not siRNAforge performance. |
| ASO Atlas 2.0 | Pinned download recipe | ASO inhibition, dose-response, hepatotoxicity, and neurotoxicity review; not integrated siRNA efficacy validation. |
| TxBench-Oligo Discovery | Pinned download recipe | Public prompts, trajectories, paper, and taxonomy; underlying data, graders, and most evaluations withheld. |

Detailed contracts and attribution remain in the [artifact guide](benchmark_artifacts.md),
[ML-review README](https://github.com/Austin-s-h/sirnaforge/blob/main/tests/data/benchmarks/ml_reviews/README.md),
[Huesken README](https://github.com/Austin-s-h/sirnaforge/blob/main/tests/unit/data/README.md), and [CCR5 guide](ccr5_zfn_benchmark.md).
Preserve the curated ML-review source extracts: upstream extraction scripts are not all vendored.
Checked-in fixtures stay in Git so tests work offline.

## Download public sources

The [download manifest](../tests/data/benchmarks/public_sources.json) records 88 file URLs,
pinned revisions, SHA-256 values, and sizes. All 88 pinned files were fetched and verified during
the October 17, 2026 cleanup. Two upstream helper scripts had only local formatting changes;
the manifest uses upstream bytes.

| Public source | Reviewed snapshot | Notes |
| --- | --- | --- |
| [OligoGym](https://huggingface.co/datasets/CollageBio/oligogym), [curation code](https://github.com/Roche/OligoGym) | 12 compressed CSVs, 337,593 rows | Mixed siRNA, ASO, and shRNA endpoints; metadata and source citations included in downloads. |
| [ASO Atlas 2.0](https://huggingface.co/datasets/barneyhill/aso-atlas-2) | Four assay parquet tables plus endpoint-specific folds | 174,459 inhibition, 99,821 dose-response, 4,341 hepatotoxicity, 4,900 neurotoxicity rows. Preserve patent-grouped folds. |
| [TxBench public repository](https://github.com/latchbio/txbench-oligo-discovery) | Public snapshot | Four example prompts, released trajectories, paper/source, and notices. Preserve its canary and training-use restrictions. |

```bash
uv run python scripts/fetch_public_benchmarks.py --dataset oligogym
uv run python scripts/fetch_public_benchmarks.py --dataset aso_atlas_2
uv run python scripts/fetch_public_benchmarks.py --dataset txbench_oligo_discovery
# Omit --dataset to fetch everything. Verify selected local files with --verify-only.
```

Files go to ignored `tests/data/external/`; `--output-dir` selects another location.
Valid local files are reused; checksum failures prevent replacement. Original dataset terms and
attribution apply; the software license does not relicense third-party data.

| OligoGym dataset | Observations | Interpretation |
| --- | ---: | --- |
| `ichihara_2007_1` | 2,431 | Huesken-derived; overlaps the existing efficacy basis. |
| `ichihara_2007_2` | 419 | Mixed-study siRNA activity; 19-nt paired-core native input. |
| `martinelli_2023_1` | 907 | Modified siRNA; nominal 21-nt input, target identity unqualified. |
| `shmushkovich_2018_1` | 356 | Asymmetric 15/20-nt architecture; native review only. |

Other downloads cover ASO activity/toxicity/immunomodulation and 291,551 shRNA observations.
Do not pool these endpoints into siRNA efficacy labels. Martinelli's nominal strand is not proof
of a fully complementary 21-nt duplex; the artifact guide records its mixed geometry.

## Rebuild review and synthetic fixtures

Use isolated output directories for comparison before replacing checked-in fixtures:

```bash
uv run python scripts/prepare_oligo_benchmarks.py \
  --output-dir benchmark_results/rebuilt_oligogym

uv run python scripts/build_review_benchmarks.py \
  --source-dir tests/data/benchmarks/ml_reviews \
  --existing-records-csv tests/data/benchmarks/oligogym/records.csv \
  --out-dir benchmark_results/rebuilt_ml_reviews
```

Synthetic inputs use `70 × A + reverse-complement(guide) + 70 × A` with site start 71.
Accessibility scores describe fabricated context. Source measurements and label direction remain
in the tables. Review collections preserve dose/replicate identity, censored and categorical
endpoints, and chemistry. Collection names and sequence clusters do not establish independent
assays. The ML-review table is not a registered artifact panel.

## Rebuild native design and references

Obtain uncompressed human, mouse, rat, and macaque Ensembl cDNA FASTAs.
[Historical reference metadata](_data/ensembl_reference_sources_historical.json) retains public
URLs, download dates, and cache checksums. These downloads used mutable `current_fasta` URLs,
not pinned releases. New downloads permit a rerun, but do not guarantee exact historical bytes.
Pin an Ensembl release and record fresh checksums for new experiments.

Assuming reference files at `/ref/human.fa`, `/ref/mouse.fa`, `/ref/rat.fa`, and `/ref/macaque.fa`:

```bash
uv run python scripts/prepare_native_design_benchmarks.py \
  --cdna-fasta /ref/human.fa \
  --output-dir benchmark_results/native_design
```

Established native mapping coverage: 216/2,431 Ichihara-1 observations, 65/419 Ichihara-2,
35/907 Martinelli. Shmushkovich remained review-only. A Martinelli strand match does not confirm
its target gene. Missing mappings remain explicit.

The established runner uses inputs under `tests/data/benchmarks/oligogym_native_design`.
That is the preparer's default and contains tracked fixtures; using it deliberately rewrites
them. Review an isolated build first.

```bash
for species in human mouse rat macaque; do
  uv run python scripts/prepare_representative_ensembl_transcriptomes.py \
  --input-fasta "/ref/${species}.fa" \
  --mapping-path tests/data/benchmarks/oligogym_native_design/mappings.csv \
  --output-fasta "tests/data/benchmarks/oligogym_native_design/references/representative/${species}_protein_coding_gene_representatives.fa" \
  --manifest "benchmark_results/${species}_representative_manifest.json"
done
```

Selection retains all protein-coding isoforms for mapped gene IDs, and the longest protein-coding
transcript per other gene. Human IDs do not preserve orthologous isoforms in other species:
those species get longest-transcript representatives unless their own IDs are supplied.
This is reduced isoform coverage, not complete-transcriptome validation.

The historical human reference contained 13 duplicate transcript IDs. The preparer now emits
each transcript once; new checksums must remain separate from historical manifests.

### Orthology and workflow

Offline mappings for the established fixtures are tracked. To refresh an Ichihara mapping:

```bash
uv run python scripts/build_panel_ortholog_mapping.py \
  --mapping-path tests/data/benchmarks/oligogym_native_design/mappings.csv \
  --dataset ichihara_2007_1 \
  --query-key benchmark_oligogym_ichihara_2007_1 \
  --output benchmark_results/ichihara_2007_1_orthology.json \
  --metadata benchmark_results/ichihara_2007_1_orthology_metadata.json
```

Repeat for panel 2. Martinelli's tracked mapping is explicitly empty because target identity is
unqualified. Synthetic panel identifiers are not valid genes for online orthology lookup.

Inside the benchmark container with the checkout mounted at `/workspace`:

```bash
bash scripts/run_oligogym_preproduction_panels.sh
```

The runner uses tracked inputs, generated references, explicit offline orthology, and sequential
human/mouse/rat/macaque screening. Settings: qualified mode, GC 0–100, active homopolymer filter,
no modifications, no top-N cap. It uses checkout scoring defaults; record commit and profiles.
Outputs go to ignored `workflow_results/oligogym_native_design_preproduction_v1/`.
Verify a final workflow summary for every panel before claiming completion.

These multi-gene panel FASTAs define a pooled intended target set: the workflow considers all
input genes on-target. Post-screen associations from this runner are therefore panel-level;
gene-specific specificity requires separate per-gene runs. Check every requested species/channel's
evidence envelope, including its submitted/processed guide counts, rather than relying on exit zero.

## Historical scoring findings

The [historical analysis JSON](_data/oligogym_scoring_analysis_historical.json) preserves the
September 13, 2026 investigation. Its paths are historical provenance strings; original run
directories are absent. Exact baseline runs cannot be regenerated from the remaining files.

| Panel | Observations used | Linked candidates | Historical Spearman rho | Alternative rho |
| --- | ---: | ---: | ---: | ---: |
| Ichihara 1 | 186 | 5,205 | 0.265 | 0.296 |
| Ichihara 2 | 55 | 3,270 | 0.228 | 0.259 |
| Martinelli | 35 | 1,585 | −0.317 | −0.223 |

Median aggregation collapses multiple transcript-site candidates per observation.
Historical weights were 0.25/0.30/0.25/0.20; the constrained comparison used
0.25/0.10/0.30/0.35, ordered off-target/accessibility/asymmetry/GC. Historical rho bootstrap
intervals were 0.132–0.400, −0.041–0.472, and −0.670–0.085 respectively.

These are exploratory in-sample associations, not independently validated efficacy predictions.
GC and asymmetry showed the strongest repeatable activity associations. Off-target allocation
should not be fitted to activity labels; it had no variation in two panels.
Full-expanded-isoform comparisons were unavailable.

Historical work notes reported 736,204 candidates and 66,683 unique screened guides over three
baseline runs. Later preproduction outputs contained final summaries for both Ichihara panels
but no final Martinelli summary; they were not three verified completed reruns.

### Recompute analysis and figures

The runner emits `ichihara_2007_1_len19`, `ichihara_2007_2_len19`, and
`martinelli_2023_1_len21` directories, matching the linker's expectations:

```bash
uv run python scripts/link_native_design_candidates.py \
  --workflow-root workflow_results/oligogym_native_design_preproduction_v1
```

This rewrites trace tables beside tracked native fixtures; review before committing.
For analysis, supply the actual weights used to create source CSV score components.
This example applies only to source runs with weights 0.25/0.10/0.30/0.35:

```bash
uv run python scripts/analyze_oligogym_scoring.py \
  --workflow-root workflow_results/oligogym_native_design_preproduction_v1 \
  --run ichihara_2007_1=ichihara_2007_1_len19 \
  --run ichihara_2007_2=ichihara_2007_2_len19 \
  --run martinelli_2023_1=martinelli_2023_1_len21 \
  --component-weights 0.25 0.10 0.30 0.35 \
  --output benchmark_results/oligogym_scoring_analysis.json
```

The analysis records normalization weights. Do not apply historical normalization to newer profiles.
Published missing `composite_score` values remain missing; available components do not substitute
for a score the workflow did not publish. The `current` comparison uses `--component-weights`.
The exploratory weight search may take time. The renderer was written for historical comparison
labels; review captions against source profiles for new experiments.

Recreate historical figures without workflow outputs:

```bash
uv run python scripts/render_oligogym_scoring_report.py \
  --analysis docs/_data/oligogym_scoring_analysis_historical.json
```

## Disposable files

`benchmark_results/`, `workflow_results/`, `tests/data/external/`, both
`oligogym_native_smoke*` experiment directories, representative reference FASTAs, and generated
`docs/_static/oligogym_scoring/` figures are ignored. Remove them after preserving any new findings
and provenance. Keep tracked regression fixtures, curated source extracts, and historical summaries.
