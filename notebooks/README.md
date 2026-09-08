# Notebook Layout

This directory should stay small and purpose-driven.

## Active Notebooks

- `zfn_backend_runtime_comparison.ipynb`: canonical ZFN backend benchmark and fragmentation analysis

> **⚠️ EXPERIMENTAL — the recorded ZFN runs are timing evidence, not correctness evidence.** The ZFN
> arm ships experimental in 0.6.0 with known unfixed defects tracked in
> [#82](https://github.com/Austin-s-h/sirnaforge/issues/82). `zfn_backend_runtime_comparison.ipynb`
> pins the published CCR5 pair `GTCATCCTCATC` / `AAACTGCAAAAG`, which matches no genomic site under
> the default strand-pairing rule; the notebook only finds sites because it builds its own fixture
> with `canonical_site()`, which reverse-complements the right half-site for it. Its parity and
> timing conclusions stand; **do not read its site counts as genome-scale behaviour**, and do not
> cite it as ZFN validation — ZFN output needs independent validation before it informs any
> decision. Anything derived from it needs re-measuring once #82 is resolved.

The two notebooks previously listed here — `zfn_default_behavior_bottleneck_profile.ipynb` and
`zfn_dev_image_experiment_series.ipynb` — are **not in the repository** (untracked and absent, not
merely archived). The entries were stale; if either is restored it needs the caveat above.

## Historical Material

- `archive/`: notebooks retained for provenance after their key conclusions have been extracted into tests, docs, or code defaults

## Run Artifacts

- `zfn_experiment_runs/`: saved benchmark artifacts that are intentionally kept with the repo

## Rule of Thumb

If a notebook conclusion affects shipped defaults, developer workflow, or expected correctness, extract that knowledge into tests or docs before letting the notebook become the only source of truth.
