"""Token parsing with no sirnaforge imports at all.

A leaf on purpose: the CLI, the workflow, the Nextflow entry points and the ZFN bridge all split
comma-separated values, and the bridge runs inside a Nextflow process where importing
``utils.cli_inputs`` (which reaches the miRNA manager and the reference policy) would be too heavy.
"""

from __future__ import annotations

__all__ = ["parse_csv"]


def parse_csv(value: str) -> list[str]:
    """Split a comma-separated string into normalized non-empty tokens."""
    return [token.strip() for token in value.split(",") if token.strip()]
