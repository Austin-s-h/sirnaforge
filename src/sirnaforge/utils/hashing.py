"""One streaming file digest for every manifest, envelope and provenance record.

Plain hex, no ``sha256:`` prefix, matching ``source_sha256``'s convention in the benchmark
manifests. Chunked so a 1.5 GB reference costs a buffer rather than memory.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

__all__ = ["file_sha256"]

_CHUNK_BYTES = 1024 * 1024


def file_sha256(path: str | Path) -> str:
    """SHA-256 of a file, for integrity (non-security) attestation."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()
