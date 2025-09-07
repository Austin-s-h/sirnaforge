#!/usr/bin/env python3
"""Resolve genome index path for a species using genomes.yaml or a comma-separated --genome_indices string.

Prints the resolved index prefix to stdout (or prints nothing if not found).

Usage:
    python3 resolve_index.py <species> [genome_config]

It will also consult the environment variable NXF_PARAMS_GENOME_INDICES if present.
This script prefers PyYAML when available but falls back to a minimal parser so it
can run in lightweight containers that don't include PyYAML.
"""

import os
import sys
from pathlib import Path

try:
    import yaml  # type: ignore
except ImportError:
    yaml = None


def parse_simple_yaml(path):
    # Minimal YAML fallback: parse top-level key: value pairs where value is a string
    data = {}
    try:
        with Path(path).open() as fh:
            for file_line in fh:
                line = file_line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    key, val = line.split(":", 1)
                    key = key.strip()
                    val = val.strip().strip("\"'")
                    data[key] = val
    except Exception:
        pass
    return data


def main():  # noqa: PLR0912
    if len(sys.argv) < 2:
        print("", end="")
        return 0

    sp = sys.argv[1]
    cfg = None
    if len(sys.argv) >= 3 and sys.argv[2]:
        cfg = sys.argv[2]

    cfg = cfg or os.environ.get("NXF_PARAMS_GENOME_CONFIG") or "nextflow_pipeline/genomes.yaml"
    gindices = os.environ.get("NXF_PARAMS_GENOME_INDICES")

    idx = None

    # Try PyYAML if available

    if Path(cfg).exists():
        try:
            if yaml:
                with Path(cfg).open() as fh:
                    g = yaml.safe_load(fh) or {}
            else:
                g = parse_simple_yaml(cfg)

            # g may be a mapping where species key maps to either a string or a dict
            if isinstance(g, dict) and sp in g:
                val = g[sp]
                if isinstance(val, dict) and val.get("index_prefix"):
                    idx = val.get("index_prefix")
                elif isinstance(val, str):
                    idx = val
        except Exception:
            idx = None

    # Fallback to NXF_PARAMS_GENOME_INDICES env (comma-separated list like "human:/path,rat:/path")
    if not idx and gindices:
        for entry in gindices.split(","):
            parts = entry.split(":", 1)
            if len(parts) == 2 and parts[0].strip() == sp:
                idx = parts[1].strip()
                break

    if idx:
        sys.stdout.write(idx)
        return None
    return None


if __name__ == "__main__":
    try:
        main()
    except Exception:
        # Never raise an exception to the caller; Nextflow will treat missing output as empty
        sys.exit(0)
