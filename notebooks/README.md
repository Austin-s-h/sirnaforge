# Notebook Layout

This directory should stay small and purpose-driven.

## Active Notebooks

- `zfn_backend_runtime_comparison.ipynb`: canonical ZFN backend benchmark and fragmentation analysis
- `zfn_default_behavior_bottleneck_profile.ipynb`: focused search-path bottleneck analysis that still contains unique performance rationale
- `zfn_dev_image_experiment_series.ipynb`: staged container-validation narrative

## Historical Material

- `archive/`: notebooks retained for provenance after their key conclusions have been extracted into tests, docs, or code defaults

## Run Artifacts

- `zfn_experiment_runs/`: saved benchmark artifacts that are intentionally kept with the repo

## Rule of Thumb

If a notebook conclusion affects shipped defaults, developer workflow, or expected correctness, extract that knowledge into tests or docs before letting the notebook become the only source of truth.
