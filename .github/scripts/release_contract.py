#!/usr/bin/env python3
"""Verify release identity before a tag can publish artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import urllib.error
import urllib.request
import zipfile
from email.parser import Parser
from pathlib import Path

import tomllib

_VERSION = re.compile(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:rc([1-9]\d*))?$")


def parse_version(version: str) -> tuple[int, int, int, int | None]:
    """Accept only versions shared by SemVer's core and PyPI's rc syntax."""
    match = _VERSION.fullmatch(version)
    if not match:
        raise ValueError(f"not a supported SemVer/PEP 440 release: {version!r}")
    major, minor, patch, rc = match.groups()
    return int(major), int(minor), int(patch), int(rc) if rc else None


def validate_identity(version: str, project_version: str, changelog: str, *, tagged: bool) -> None:
    """Require one package version and a release-specific changelog entry."""
    parse_version(version)
    if version != project_version:
        raise ValueError(f"tag v{version} disagrees with package version {project_version}")
    heading = re.compile(rf"^## \[{re.escape(version)}\] - (.+)$", re.MULTILINE)
    match = heading.search(changelog)
    if not match:
        raise ValueError(f"missing changelog entry for {version}")
    if tagged and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", match.group(1)):
        raise ValueError(f"release requires a dated changelog entry for {version}")


def validate_wheel(wheel: Path, version: str) -> None:
    """Check wheel metadata rather than trusting its filename."""
    with zipfile.ZipFile(wheel) as archive:
        metadata_files = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(metadata_files) != 1:
            raise ValueError("wheel must contain exactly one dist-info/METADATA")
        metadata = Parser().parsestr(archive.read(metadata_files[0]).decode("utf-8"))
    if metadata.get("Name", "").lower() != "sirnaforge" or metadata.get("Version") != version:
        raise ValueError(f"wheel metadata disagrees with release {version}")


def validate_published_wheel(wheel: Path, published: dict) -> None:
    """Require PyPI to host the tested wheel bytes."""
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    matches = [item for item in published["urls"] if item["filename"] == wheel.name]
    if len(matches) != 1 or matches[0]["digests"]["sha256"] != digest:
        raise ValueError(f"PyPI wheel digest disagrees with tested artifact: {wheel.name}")


def assert_unpublished(version: str) -> None:
    """Reject immutable PyPI versions already in use."""
    request = urllib.request.Request(f"https://pypi.org/pypi/sirnaforge/{version}/json")
    try:
        with urllib.request.urlopen(request, timeout=15):
            pass
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return
        raise
    raise ValueError(f"PyPI already has sirnaforge {version}; versions cannot be republished")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def main() -> None:
    """Run the requested release identity or artifact preflight."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("candidate", "tag", "wheel", "published"))
    parser.add_argument("--version", required=True)
    parser.add_argument("--wheel", type=Path)
    args = parser.parse_args()
    version = args.version
    if args.mode == "wheel":
        if not args.wheel:
            parser.error("--wheel is required in wheel mode")
        validate_wheel(args.wheel, version)
        return
    if args.mode == "published":
        if not args.wheel:
            parser.error("--wheel is required in published mode")
        with urllib.request.urlopen(f"https://pypi.org/pypi/sirnaforge/{version}/json", timeout=15) as response:
            published = json.load(response)
        validate_published_wheel(args.wheel, published)
        print(f"PyPI wheel matches tested bytes: {args.wheel.name}")
        return

    project_version = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    changelog = Path("CHANGELOG.md").read_text(encoding="utf-8")
    validate_identity(version, project_version, changelog, tagged=args.mode == "tag")
    if args.mode == "tag":
        tag = f"v{version}"
        if _git("rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}") != _git("rev-parse", "HEAD"):
            raise ValueError("tag target differs from checked-out commit")
        if subprocess.run(["git", "merge-base", "--is-ancestor", "HEAD", "origin/master"], check=False).returncode:
            raise ValueError("release tag must point to a commit on master")
        assert_unpublished(version)
    print(f"Release {args.mode} identity verified: {version}")


if __name__ == "__main__":
    main()
