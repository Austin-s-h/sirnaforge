# CCR5 ZFN Benchmark (PROGNOS 2014 Extract)

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
- performance summary (12 prior validated off-targets, 10/12 recovery at 3×)
- tie-break policy metadata (Exon > Promoter > Intron > Intergenic, then chromosomal location)
- fixture scope and ambiguity markers (`?`, `Sequencing Failure`, `N/A`)

## Notes and Limits

These fixtures are intentionally restricted to rows that are legible in extracted text. Full coordinate/sequence universes, complete rank matrices, and all validation statistics still require direct ingestion of original supplementary tables (S1, S10, S11, S14, S15).
