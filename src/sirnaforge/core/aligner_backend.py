"""Backend definitions for transcriptome alignment and index handling."""

from __future__ import annotations

import re
import subprocess  # nosec B404
from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

from sirnaforge.utils.subprocess_utils import _get_executable_path

AlignerName = Literal["bwa-mem2", "bwa-mem3"]


@dataclass(frozen=True, slots=True)
class AlignerBackend:
    """Executable and index contract for one supported BWA implementation."""

    name: AlignerName
    version: str
    index_extensions: tuple[str, ...]


def normalize_aligner(value: str) -> AlignerName:
    """Validate a public backend selector."""
    if value not in ("bwa-mem2", "bwa-mem3"):
        raise ValueError(f"Unsupported aligner {value!r}; choose 'bwa-mem2' or 'bwa-mem3'")
    return value  # type: ignore[return-value]


@lru_cache(maxsize=2)
def get_aligner_version(name: str) -> str | None:
    """Read the selected executable's version banner, or None when it is unavailable."""
    backend = normalize_aligner(name)
    executable = _get_executable_path(backend)
    if executable is None:
        return None
    try:
        completed = subprocess.run(  # nosec B603
            [executable, "version"], capture_output=True, text=True, timeout=10, check=True
        )
    except (OSError, subprocess.SubprocessError):
        return None
    output = (completed.stdout or completed.stderr).strip()
    match = re.search(r"(?m)^(\d+(?:\.\d+)+(?:[-+][\w.-]+)?)$", output)
    return match.group(1) if match else (output.splitlines()[0] if output else None)


def aligner_identity(name: str) -> dict[str, str | None]:
    """Stable backend identity persisted with an index stamp."""
    backend = normalize_aligner(name)
    return {"name": backend, "version": get_aligner_version(backend)}
