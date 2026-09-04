# ZFN Backend Tuning Summary

> ## ⚠️ EXPERIMENTAL — these runs are not authoritative
>
> The ZFN arm ships **experimental** in 0.6.0 with known unfixed defects tracked in [#82](https://github.com/Austin-s-h/sirnaforge/issues/82). **Do not use ZFN output for any decision without independent
> validation.**
>
> Every run summarised below used the published CCR5 half-site pair as written. Under the default
> `require_opposite_strands=True` that pair matches no site at all — including its own on-target
> site — so the conclusions here describe runtime behaviour only and must be re-measured once the
> half-site convention is fixed.

This note captures the engineering outcome from the ZFN backend tuning work so the rationale no longer lives only in notebooks.

## Question Answered

Which in-process backend should be the practical first choice for ZFN off-target discovery on large genomic references, and when should `fm_index` be preferred instead?

## Benchmark Surfaces Used

The tuning work compared the three in-process backends across the benchmark surfaces that mattered for runtime decisions:

- synthetic fixed-size fragmentation sweeps for controlled backend comparison
- real hg38 chr3 slices to expose fragmentation and contig-shape effects on realistic Ensembl sequence
- full hg38 primary manual runs using the CCR5 ZFN half-sites and Ensembl annotation

The intent was not just raw parity, but the operational question of which backend remains fastest and easiest to use when the search space is large enough that preparation costs and sharding behavior matter.

## Decision

The tuning conclusion is:

1. `pyahocorasick` is the operational default for direct ZFN runs.
2. `fm_index` is the alternate backend only when the same search space will be reused enough to amortize persisted-index construction; on large references it should be treated as experimental.
3. `exhaustive_python` remains the strict baseline and fallback implementation, not the recommended first-pass runtime choice.

## Why `pyahocorasick` Won

`pyahocorasick` was selected as the practical default because it had the best overall tradeoff across both synthetic and real chr3 slices:

- fastest practical time-to-first-result across the measured sweeps
- no persisted search-space bundle required before use
- less sensitivity to fragmented search-space inputs than `fm_index`
- a simpler operational story for one-off and exploratory hg38 runs

That combination matters more than peak theoretical throughput. The default path has to be the backend that is fastest to use correctly, not only the backend that can win once more setup has already been paid.

## Why `fm_index` Is the Alternate

`fm_index` stays important, but as a deliberate alternate path rather than the first recommendation:

- it benefits most when one indexed reference is reused across many runs
- its dominant extra cost is backend preparation, especially on full-primary references
- it prefers consolidated inputs and is less attractive when the workload is repeatedly fragmented or rebuilt
- it now uses a more conservative sharding profile because it can over-fragment and over-allocate more easily on full-primary genomes

In practice, use `fm_index` when you are doing repeated searches against the same hg38 reference and can reuse a persisted bundle. Use `pyahocorasick` when you want the best first-run experience.

Recent full-primary runs also showed memory pressure can become severe with `fm_index` on large references, so release-facing guidance now treats that path as experimental until memory behavior is improved.

## Why `exhaustive_python` Did Not Stay First

`exhaustive_python` remains valuable for correctness work:

- backend-parity checks
- dependency-light fallback behavior
- debugging when accelerated backend behavior needs a baseline comparison

It did not remain the recommended first choice because the tuned runs showed that the accelerated backends provide a better runtime profile on the references that motivated this work.

## Resulting Runtime Policy

The tuned runtime policy is:

- use `pyahocorasick` for the default operational path
- use `fm_index` for repeated runs on a reused indexed reference
- keep `exhaustive_python` as the baseline and regression comparator

The sharding defaults follow that same policy. `fm_index` uses a more conservative autotuned shard profile, while `pyahocorasick` and `exhaustive_python` keep the higher-concurrency profile.

## Durable Artifacts

The tuning outcome is preserved in these durable places:

- code-level sharding behavior in `src/sirnaforge/cli.py`
- typed sharding defaults in `src/sirnaforge/models/zfn.py`
- backend parity coverage in `tests/integration/` and `tests/container/`
- real-genome manual validation commands in [zfn_hg38_primary_test_commands.md](zfn_hg38_primary_test_commands.md)

## Manual Validation Scope

The manual validation commands remain important for the cases that are too heavy for normal fast test runs:

- full hg38 primary runs
- real Ensembl chr3 validation slices
- backend-to-backend row-count and runtime inspection on cached real references

Use [zfn_hg38_primary_test_commands.md](zfn_hg38_primary_test_commands.md) when you need to rerun that real-reference validation.
