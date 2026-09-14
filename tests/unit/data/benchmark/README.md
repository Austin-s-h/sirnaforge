# Benchmark panel registry fixtures

Provenance, citations and checksums for the files backing `tests/unit/test_benchmark_panels.py` and
`src/sirnaforge/benchmark/panels.py`. Issue #109.

## The two verified facts this directory works around

1. **The issue's cited source does not exist.** #109 cites `docs/prd_benchmark_artifacts_and_variable
_length.md` as its specification. That file is not in the working tree and not in any branch's
   history of this repository. The issue body is the entire specification; nothing here or in
   `panels.py` claims to summarise a document that was read.
2. **Only one of the five named panels ships real bytes.** #109 and #110 name five benchmark panels:
   Huesken, Ichihara, Martinelli, Shmushkovich, OligoGym. ⚠️ **Ichihara, Martinelli, Shmushkovich and
   OligoGym are not vendored anywhere in this repository.** No bytes, no redistribution, no prior
   reference to any of the four exists in this repository's history. The only real measured data
   `panels.py` exercises is `../sirna_efficacy_subset.csv` — the 180-row Huesken redistribution
   documented in `../README.md`, reused here rather than duplicated.

   Every fixture in this directory is therefore **synthetic**, and is named `synthetic_<architecture>
_*` — deliberately _not_ named after any of the four missing panels — so a passing test here is
   never mistaken for evidence over real panel data. `PanelDescriptor` entries for the four missing
   panels exist in the registry (`data_present=False`) so the CLI and manifest can name them
   honestly; their `citation` fields say directly that no primary source has been independently
   verified in this repository, rather than inventing a DOI or PMID this repository cannot back.

## Fixtures

### `synthetic_paired_core_with_overhang.csv`

Five synthetic Tuschl-rule-style records: a 21 nt guide, a 19 nt paired core, a measured 2 nt 3'
overhang — the same duplex shape as the vendored Huesken subset. `SYN0005`'s guide is 18 nt,
deliberately shorter than any fixed-length core #109 supports, to exercise the "too short to slice"
incompatible branch. No accession is a real GenBank/RefSeq id; all four `SYN000N` prefixes are
invented for this fixture and carry no biological meaning.

### `synthetic_fully_complementary.csv`

Four synthetic records: guide and passenger paired over their full, equal 21 nt length (blunt, no
overhang). `FC0001`–`FC0003` have a guide and an exact reverse-complement passenger of equal length
(compatible). `FC0004`'s passenger is 20 nt against a 21 nt guide — a deliberate length mismatch, to
exercise the "declared fully complementary but the strands differ in length" incompatible branch.

### `synthetic_asymmetric_15_20.csv`

Two synthetic records reproducing the asymmetric 15/20 duplex geometry #110 attributes to
Shmushkovich (a 20 nt guide paired against a 15 nt passenger). Every asymmetric record is
incompatible with #109's fixed-length interface by construction (see `panels.py::_check_compatibility`
and #110's own scope) — this fixture exists to prove the incompatibility is reported honestly, with
both lengths named in the reason, and that the guide is never truncated or padded to a symmetric
19–23 nt length to force it through.

### `synthetic_context.fa`

Two short synthetic "transcript context" records, each a guide's (or a fully-complementary guide's)
reverse complement embedded in unrelated flanking sequence. Not consumed by `panels.py` itself —
provided as a placeholder for a future `--panel-transcripts` fixture (the `design_context_source
=panel_transcript` path #109's artifact schema declares fields for but does not populate; see
`derive_observation`'s docstring and #110). Purely synthetic; no accession, no real transcript, no
biological claim.

## What these fixtures do and do not prove

They prove that `derive_observation` and `PanelDescriptor` correctly implement the three declared
architectures' slicing and compatibility rules, including the boundary cases (too-short guide,
strand-length mismatch, asymmetric geometry) that the vendored 180-row Huesken subset alone does not
exercise (it contains only compatible 21 nt paired-core-with-overhang records).

They do **not** prove anything about Ichihara, Martinelli, Shmushkovich or OligoGym as real assays —
no bytes from any of the four are present here, and a synthetic record's efficacy value is an
arbitrary placeholder, not a measurement. They also do not prove the paired-core/fully-complementary
split rule generalises beyond the three architectures #109 names: a fourth architecture would need
its own descriptor, its own compatibility branch, and its own fixture.

## Checksums

SHA-256 of every fixture this file describes, so a number or a test assertion quoted elsewhere can
be traced to the exact bytes that produced it.

| File                                      | Bytes | SHA-256                                                            |
| ----------------------------------------- | ----: | ------------------------------------------------------------------ |
| `synthetic_paired_core_with_overhang.csv` |   205 | `4e77463d0f5a3f3d2ae9793bd8af735e20a013a723178e461903063d1b4389a3` |
| `synthetic_fully_complementary.csv`       |   275 | `c6cb2f4a8b8173e4d6592c01e39f2116ee46c73383beca48a9fc3355d80bf07a` |
| `synthetic_asymmetric_15_20.csv`          |   155 | `07958b60dcab55004e002a6d45aedeb4b22cc150aa867602dc4224f98110f3f1` |
| `synthetic_context.fa`                    |   384 | `550b05b8d03443835d3afe5f24f935659f0bfd6e8ddadba3ceeebe801013ed14` |

Regenerate with `shasum -a 256 <file>` from this directory; every value above is over the tracked
file exactly as committed, no line-ending or encoding transform.

## Predeclared split

`panels.py::predeclared_split` implements the identical accession-parity rule
`scripts/validate_scoring_profiles.py::split_of` pins for the Huesken panel (`sha256(accession)[0] %
2 == 0 -> development`), under its own declared id, `sha256_accession_parity_v1`
(`panels.py::PREDECLARED_SPLIT_RULE_ID`). It is deliberate duplication with a named authority: that
script is a standalone entry point, not an importable package member. `test_benchmark_panels.py`
pins the function against the exact `development`/`held_out` accession lists published in
`../README.md`, so the duplication cannot silently drift from the one place those lists are
authoritative. None of the synthetic fixtures in this directory declare a split — only
`huesken_subset` and `huesken_full` do, and only `huesken_subset` has rows to split.
