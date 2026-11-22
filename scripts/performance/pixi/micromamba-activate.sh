#!/usr/bin/env bash
# Nextflow sources this file for a Conda-backed process environment.
eval "$("${CONDA_PREFIX:?CONDA_PREFIX is required}/bin/micromamba" shell hook --shell bash)"
micromamba activate "$@"
