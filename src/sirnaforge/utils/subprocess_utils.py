"""Shared preconditions for the repository's ``subprocess`` call sites.

The absolute-path check is what the ``# nosec B603`` annotations on those calls rest on, so it lives
in one place and every caller imports it rather than restating it.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from sirnaforge.utils.logging_utils import get_logger

__all__ = ["_get_executable_path", "_validate_command_args"]

logger = get_logger(__name__)


def _get_executable_path(tool_name: str) -> str | None:
    """Get the full path to an executable, ensuring it exists."""
    path = shutil.which(tool_name)
    if path is None:
        logger.warning(f"Tool '{tool_name}' not found in PATH")
    return path


def _validate_command_args(cmd: list[str]) -> None:
    """Validate command arguments for subprocess execution."""
    if not cmd:
        raise ValueError("Command list cannot be empty")

    executable = cmd[0]
    if not executable:
        raise ValueError("Executable path cannot be empty")

    # Ensure we have an absolute path to the executable
    if not Path(executable).is_absolute():
        raise ValueError(f"Executable must be an absolute path: {executable}")
