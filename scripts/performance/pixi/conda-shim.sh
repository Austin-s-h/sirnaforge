#!/usr/bin/env bash
# Conda compatibility entry point for the disposable Pixi review image.
if [[ "${1:-}" == "info" && "${2:-}" == "--json" ]]; then
    printf '{\n  "conda_prefix": "%s",\n  "channels": []\n}\n' "${CONDA_PREFIX:?CONDA_PREFIX is required}"
else
    exec "${CONDA_PREFIX:?CONDA_PREFIX is required}/bin/micromamba" "$@"
fi
