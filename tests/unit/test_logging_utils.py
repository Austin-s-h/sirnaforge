"""Unit tests for logging utility configuration behavior."""

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

from sirnaforge.utils.logging_utils import configure_logging


def _detach_root_handlers() -> list[logging.Handler]:
    """Detach and return current root handlers for test isolation."""
    root = logging.getLogger()
    handlers = list(root.handlers)
    for handler in handlers:
        root.removeHandler(handler)
    return handlers


def _restore_root_handlers(handlers: list[logging.Handler]) -> None:
    """Restore previously detached root handlers."""
    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
        handler.close()
    for handler in handlers:
        root.addHandler(handler)


def test_configure_logging_adds_file_handler_when_root_already_has_handlers(tmp_path: Path) -> None:
    """Repeated configure calls should still attach requested rotating file handlers."""
    previous_handlers = _detach_root_handlers()
    root = logging.getLogger()

    try:
        root.addHandler(logging.StreamHandler())
        log_file = tmp_path / "sirnaforge.log"

        configure_logging(level="INFO", log_file=str(log_file))
        logging.getLogger("sirnaforge.test").info("log-line-from-test")

        assert log_file.exists()
        assert "log-line-from-test" in log_file.read_text(encoding="utf-8")

        configure_logging(level="INFO", log_file=str(log_file))
        file_handlers = [
            handler
            for handler in logging.getLogger().handlers
            if isinstance(handler, RotatingFileHandler) and Path(handler.baseFilename).resolve() == log_file.resolve()
        ]
        assert len(file_handlers) == 1
    finally:
        _restore_root_handlers(previous_handlers)
