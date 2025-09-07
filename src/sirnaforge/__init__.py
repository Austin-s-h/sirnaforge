"""siRNAforge - Comprehensive siRNA design toolkit for gene silencing. Comprehensive gene silencing design and analysis.

This module exposes package metadata (author/email/version) in a single place.
The version is resolved from installed package metadata (importlib.metadata). When
running from a source checkout (not installed), it falls back to reading
`pyproject.toml` if available, otherwise uses a conservative placeholder.
"""

from __future__ import annotations

from importlib import metadata
from importlib.metadata import PackageNotFoundError
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except Exception:  # pragma: no cover - fallback handled at runtime
    tomllib = None

__author__ = "Austin S. Hovland."
__email__ = "your-email@example.com"


def _read_version_from_pyproject() -> str | None:
    """Attempt to read the project version from pyproject.toml in the repo root.

    Returns the version string or None if it cannot be determined.
    """
    try:
        if tomllib is None:
            return None
        root = Path(__file__).resolve().parents[1]
        pyproject = root.joinpath("pyproject.toml")
        if not pyproject.exists():
            return None
        data = tomllib.loads(pyproject.read_text(encoding="utf8"))
        version = data.get("project", {}).get("version")
        return str(version) if version is not None else None
    except Exception:
        return None


def _get_version() -> str:
    """Resolve the package version from metadata or pyproject fallback.

    Priority:
    1. importlib.metadata.version("sirnaforge") when installed
    2. pyproject.toml project.version when available (dev checkout)
    3. fallback to '0.0.0+unknown'
    """
    try:
        return metadata.version("sirnaforge")
    except PackageNotFoundError:
        v = _read_version_from_pyproject()
        if v:
            return v
        return "0.0.0+unknown"


__version__ = _get_version()

# Core imports will be added as modules are implemented
# from .core.design import SiRNADesigner
# from .core.scoring import ScoringEngine
# from .core.filters import FilterEngine
# from .models.sirna import SiRNACandidate, DesignParameters

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    # "SiRNADesigner",
    # "ScoringEngine",
    # "FilterEngine",
    # "SiRNACandidate",
    # "DesignParameters",
]
