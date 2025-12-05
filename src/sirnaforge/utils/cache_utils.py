"""Shared cache helpers for siRNAforge.

Centralizes cache path resolution so every subsystem honors the same
SIRNAFORGE_CACHE_DIR / XDG cache layout.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

_CACHE_ROOT_SENTINEL = object()


def resolve_cache_subdir(subdir: str, *, override: str | os.PathLike[str] | None = None) -> Path:
    """Resolve a writable cache directory for the requested subdir.

    The lookup order matches ReferenceManager: explicit override, env vars,
    XDG, $HOME/.cache, workspace-local fallback, then temp dir.
    """
    if not subdir:
        raise ValueError("subdir must be provided")

    candidates: list[Path] = []
    if override is not None:
        candidates.append(Path(override))
    else:
        env_override = os.getenv("SIRNAFORGE_CACHE_DIR")
        if env_override:
            candidates.append(Path(env_override) / subdir)

        xdg_cache = os.getenv("XDG_CACHE_HOME")
        if xdg_cache:
            candidates.append(Path(xdg_cache) / "sirnaforge" / subdir)

        candidates.append(Path.home() / ".cache" / "sirnaforge" / subdir)
        candidates.append(Path.cwd() / ".sirnaforge_cache" / subdir)
        candidates.append(Path(tempfile.gettempdir()) / "sirnaforge" / subdir)

    first_error: Exception | None = None
    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            return candidate
        except PermissionError as exc:
            if first_error is None:
                first_error = exc
            continue

    if override is not None and first_error is not None:
        raise first_error
    raise RuntimeError(f"Cannot create cache directory for '{subdir}'")


def stable_cache_key(payload: Mapping[str, Any]) -> str:
    """Create a deterministic sha256 digest for arbitrary JSON-serializable data."""
    normalized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()[:16]
