"""Logging utilities for siRNAforge toolkit.

This module provides a single point to configure logging for both console
and an optional centralized log file. Call `configure_logging` once at
application startup (CLI entrypoint) to enable file logging. Individual
modules should use `get_logger(__name__)` to obtain a configured logger.
"""

import logging
import os
import sys
from contextlib import suppress
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional

DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


def configure_logging(level: Optional[str] = None, log_file: Optional[str] = None) -> None:
    """Configure root logger: console + optional rotating file handler.

    Args:
        level: Logging level name (DEBUG/INFO/WARNING/ERROR). If None, uses INFO.
        log_file: Path to central log file. If None, will read env var
            SIRNAFORGE_LOG_FILE. If still None, no file handler is added.
    """
    root = logging.getLogger()
    # Avoid re-configuring if already configured
    if root.handlers:
        return

    env_level = os.getenv("SIRNAFORGE_LOG_LEVEL")
    lvl_name = (level or env_level or "INFO").upper()
    lvl = getattr(logging, lvl_name, logging.INFO)
    root.setLevel(lvl)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))
    root.addHandler(console_handler)

    logfile = log_file or os.getenv("SIRNAFORGE_LOG_FILE")
    if logfile:
        # Ensure directory exists
        # Ensure directory exists (ignore errors)
        with suppress(Exception):
            parent = Path(logfile).parent
            if parent and not parent.exists():
                parent.mkdir(parents=True, exist_ok=True)

        file_handler = RotatingFileHandler(logfile, maxBytes=10 * 1024 * 1024, backupCount=5)
        file_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))
        root.addHandler(file_handler)


def get_logger(name: str, level: Optional[str] = None) -> logging.Logger:
    """Get a logger with standard configuration.

    This will return a child logger of the root logger configured by
    `configure_logging`. For scripts that don't call `configure_logging`,
    get_logger will still set a console handler on first use.
    """
    logger = logging.getLogger(name)

    if not logging.getLogger().handlers:
        # Lazy configure console-only if the root hasn't been configured
        configure_logging(level=level)

    # Allow per-logger override
    if level:
        logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    return logger
