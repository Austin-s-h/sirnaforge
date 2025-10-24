"""Utilities for resolving user-provided input resources.

Supports downloading transcript FASTA files from remote locations and
normalises them into local paths that the workflow can consume.
"""

from __future__ import annotations

import hashlib
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from urllib.parse import ParseResult, urlparse

import httpx

from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class InputSource:
    """Normalized representation of a workflow input resource."""

    original: str
    local_path: Path
    source_type: str
    downloaded: bool
    size_bytes: int
    sha256: str | None = None

    @property
    def stem(self) -> str:
        """Return the filesystem stem for the local resource."""
        return self.local_path.stem


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ensure_destination(destination_root: Path) -> Path:
    destination_root.mkdir(parents=True, exist_ok=True)
    return destination_root


def _build_input_source(
    *,
    original: str,
    path: Path,
    source_type: str,
    downloaded: bool,
) -> InputSource:
    size_bytes = path.stat().st_size
    sha_value = _hash_file(path)
    return InputSource(
        original=original,
        local_path=path.resolve(),
        source_type=source_type,
        downloaded=downloaded,
        size_bytes=size_bytes,
        sha256=sha_value,
    )


def _handle_local_path(input_location: str) -> InputSource:
    local_path = Path(input_location).expanduser().resolve()
    if not local_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {local_path}")
    return _build_input_source(
        original=input_location,
        path=local_path,
        source_type="local",
        downloaded=False,
    )


def _handle_file_uri(parsed: ParseResult) -> InputSource:
    target = Path(parsed.path).expanduser()
    if parsed.netloc:
        target = Path(f"/{parsed.netloc}{parsed.path}").expanduser()
    target = target.resolve()
    if not target.exists():
        raise FileNotFoundError(f"Input FASTA not found: {target}")
    return _build_input_source(
        original=parsed.geturl(),
        path=target,
        source_type="file",
        downloaded=False,
    )


def _handle_http_like(
    *,
    input_location: str,
    parsed: ParseResult,
    destination_root: Path,
    timeout: float,
    downloader: Callable[[str, Path, float], None],
) -> InputSource:
    destination = _ensure_destination(destination_root)
    filename = Path(parsed.path).name or "remote_input.fasta"
    target_path = destination / filename

    if target_path.exists():
        logger.info("Reusing previously downloaded FASTA: %s", target_path)
    else:
        downloader(input_location, target_path, timeout)

    return _build_input_source(
        original=input_location,
        path=target_path,
        source_type=parsed.scheme.lower(),
        downloaded=True,
    )


def _download_http(url: str, target_path: Path, timeout: float) -> None:
    logger.info("Downloading FASTA from %s", url)
    with httpx.stream("GET", url, timeout=timeout, follow_redirects=True) as response:
        response.raise_for_status()
        with target_path.open("wb") as handle:
            for chunk in response.iter_bytes():
                handle.write(chunk)


def _download_ftp(url: str, target_path: Path, timeout: float) -> None:  # noqa: ARG001
    logger.info("Downloading FASTA from FTP: %s", url)
    with urllib.request.urlopen(url) as response, target_path.open("wb") as handle:  # noqa: S310
        shutil.copyfileobj(response, handle)


def resolve_input_source(input_location: str, destination_root: Path, *, timeout: float = 30.0) -> InputSource:
    """Resolve a workflow input location into a local path.

    Args:
        input_location: Raw string provided by the user (path or URI).
        destination_root: Directory where downloaded inputs should be stored.
        timeout: Timeout for remote downloads in seconds.

    Returns:
        InputSource describing the normalized local resource.

    Raises:
        FileNotFoundError: If a local file doesn't exist.
        ValueError: If the URI scheme is unsupported.
        httpx.HTTPStatusError: If the remote download fails with non-2xx status.
    """

    if not input_location:
        raise ValueError("input_location must be a non-empty string")

    parsed = urlparse(input_location)

    if not parsed.scheme:
        return _handle_local_path(input_location)

    scheme = parsed.scheme.lower()
    if scheme == "file":
        return _handle_file_uri(parsed)
    if scheme in {"http", "https"}:
        return _handle_http_like(
            input_location=input_location,
            parsed=parsed,
            destination_root=destination_root,
            timeout=timeout,
            downloader=_download_http,
        )
    if scheme == "ftp":
        return _handle_http_like(
            input_location=input_location,
            parsed=parsed,
            destination_root=destination_root,
            timeout=timeout,
            downloader=_download_ftp,
        )

    raise ValueError(f"Unsupported input source scheme: {scheme}")


__all__ = ["InputSource", "resolve_input_source"]
