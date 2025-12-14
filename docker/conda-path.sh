#!/bin/bash
# Ensure conda toolchain is on PATH for login shells
# This script is sourced by /etc/profile.d/ during bash login initialization

# Prepend conda binaries to PATH if not already present
# Check for /opt/conda/bin since both paths are added as a pair
if [[ ":$PATH:" != *":/opt/conda/bin:"* && ":$PATH:" != *":/opt/conda/condabin:"* ]]; then
    export PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
fi
