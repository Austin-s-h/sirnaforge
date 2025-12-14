#!/bin/bash
# Ensure conda toolchain is on PATH for login shells
# This script is sourced by /etc/profile.d/ during bash login initialization

# Prepend conda binaries to PATH if not already present
if [[ ":$PATH:" != *":/opt/conda/bin:"* ]]; then
    export PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
fi
