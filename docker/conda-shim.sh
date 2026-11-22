#!/usr/bin/env bash
# Nextflow 25.10 expects `conda info --json` to report conda_prefix and sources
# that prefix's bin/activate script for every Conda-backed process. Micromamba's
# JSON uses a different schema, so expose the one field Nextflow needs while
# forwarding environment operations to the installed micromamba binary.
if [[ "${1:-}" == "info" && "${2:-}" == "--json" ]]; then
    printf '{\n  "conda_prefix": "/opt/conda",\n  "channels": []\n}\n'
else
    exec /usr/local/bin/micromamba "$@"
fi
