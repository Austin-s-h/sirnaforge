"""Ensembl identifier normalisation, shared by every module that compares IDs.

A leaf on purpose: hit classification, the transcript index, the orthologue resolver and the
workflow all need this and must not import each other to get it.
"""

from __future__ import annotations

import re

__all__ = ["strip_version"]

_VERSION_SUFFIX = re.compile(r"\.\d+$")


def strip_version(identifier: str) -> str:
    """Drop an Ensembl version suffix (e.g. ``.9``) so IDs compare equal across releases."""
    return _VERSION_SUFFIX.sub("", identifier.strip())
