# CCR5 ZFN Benchmark (PROGNOS 2014 Extract)

> ## ⚠️ EXPERIMENTAL — literature evidence is not a Sirnaforge validation panel
>
> The implementation defects tracked in [#82](https://github.com/Austin-s-h/sirnaforge/issues/82)
> have regression coverage. ZFN remains experimental because these software checks do not establish
> nuclease activity, biological specificity, or ranking performance against an independent panel.
> **Require independent validation before using ZFN output for decisions.**
>
> This page preserves the PROGNOS visible rows and their hg19 coordinates. The published CCR5
> on-target row is chr3:46,414,544 on hg19; do not use that coordinate as GRCh38. The exact target
> sequence occurs at chr3:46,373,053–46,373,081 in the GRCh38 chr3 reference slice at
> `tests/unit/data/zfn/ccr5_grch38_on_target.fa`. With the default
> `require_opposite_strands=True`, the CLI right query is the **reverse complement** of the
> publication's `(−)` text. The pair as printed in the source table therefore finds no inward-facing
> site, while the model-oriented query finds the GRCh38 on-target:
>
> ```bash
> # published text, as tabulated below -> 0 sites
> sirnaforge zfn --zfn-left-half-site GTCATCCTCATC --zfn-right-half-site AAACTGCAAAAG ...
> # model-oriented query; GRCh38 chr3:46373053 -> 1 exact site
> sirnaforge zfn --zfn-left-half-site GTCATCCTCATC --zfn-right-half-site CTTTTGCAGTTT ...
> ```
>
> Apply the same reverse complement to every `(−) half-site` value in the fixtures before feeding it
> to the tool. See [ZFN Module Guide](zfn_module.md).

This page documents the checked-in CCR5 ZFN benchmark fixtures derived from visible rows in Fine et al. (2014) supplementary extracts.

The benchmark is intentionally maintained as **cleaned test data + annotation metadata**, not as a loader-validation target.

## Included Fixtures

- `tests/unit/data/zfn/ccr5_s10_visible_rows.csv`
- `tests/unit/data/zfn/ccr5_s11_homology_visible_rows.csv`
- `tests/unit/data/zfn/ccr5_benchmark_annotations.json`

These fixtures preserve extraction ambiguity markers (`?`, `Sequencing Failure`) exactly as source-visible text.

## S10 Coverage (Visible Rows)

`ccr5_s10_visible_rows.csv` includes 17 visible rows with these fields:

- `Closest gene`, `Match type`, `hg19 coordinate`
- `(+) half-site`, `(−) half-site`
- empty/active indel and read columns
- active mutation frequency and one-sided `p-value`
- extraction notes

High-value rows captured in tests:

- `CCR5` on-target control (`chr3:46414544`, `p = 2.7E-33`)
- `CSNK1G3` novel validated off-target (`chr5:123393701`, active indels `17`, mutation frequency `0.086%`, `p = 0.000019`)
- `ZNF587` sequencing-failure row handling

## S11 Coverage (Homology Ranking)

`ccr5_s11_homology_visible_rows.csv` includes 32 visible rows (rank 2 through 33) with:

- homology rank and mismatch split (`T`, `+`, `−`)
- interrogation provenance
- closest gene, match type, coordinate, and half-sites

Tests assert:

- rank range integrity (2..33)
- presence of `CSNK1G3` at rank 33
- stable row-level values for deterministic regression use

## Annotation Model Fixture

`ccr5_benchmark_annotations.json` captures reusable benchmark context:

- nuclease/platform assumptions (CCR5-224, 4-finger monomers, heterodimeric FokI)
- performance summary (12 prior validated off-targets, 10/12 recovery at 3×) — **this is PROGNOS's
  own published figure for PROGNOS, quoted from Fine et al. 2014, not a siRNAforge measurement.**
  The JSON field is named `performance.prognos_recovery_at_3x` for that reason, and
  `ccr5_study_facts.json` records it as a target to regress against "once golden-site
  coordinates/sequences are available". Sirnaforge has not been scored against this experimental
  panel. The GRCh38 primary search confirms the exact CCR5 on-target and reports two-mismatch
  candidate sites, but that is not a recovery comparison against the PROGNOS panel; its rows still
  need verified sequence/coordinate matching and independent review before they can validate this
  tool.
- tie-break policy metadata (Exon > Promoter > Intron > Intergenic, then chromosomal location) —
  nested-feature classification has regression coverage, but annotation accuracy still needs
  independent real-genome review.
- fixture scope and ambiguity markers (`?`, `Sequencing Failure`, `N/A`)

## Notes and Limits

These fixtures are intentionally restricted to rows that are legible in extracted text. Full coordinate/sequence universes, complete rank matrices, and all validation statistics still require direct ingestion of original supplementary tables (S1, S10, S11, S14, S15).
