"""Validate installed wheel requirements against a Pixi prefix."""

from __future__ import annotations

import json
import sys
from importlib.metadata import PackageNotFoundError, requires, version
from pathlib import Path
from shutil import which

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

PREFIX = Path(sys.prefix)


def conda_version(name: str) -> str:
    """Read a unique Conda package record."""
    records = [
        record
        for path in (PREFIX / "conda-meta").glob(f"{name}-*.json")
        if (record := json.loads(path.read_text()))["name"] == name
    ]
    if len(records) != 1:
        msg = f"Expected one Conda provider for {name}, found {len(records)}"
        raise ValueError(msg)
    return str(records[0]["version"])


def verify_application() -> None:
    """Require the installed closure to satisfy wheel metadata."""
    for entry in requires("sirnaforge") or []:
        req = Requirement(entry)
        if req.marker and not req.marker.evaluate():
            continue
        name = canonicalize_name(req.name)
        # These Conda packages provide the runtime without Python dist-info.
        if name in {"nextflow", "viennarna"}:
            try:
                version(name)
            except PackageNotFoundError:
                pass
            else:
                msg = f"Unexpected second provider for {name}"
                raise ValueError(msg)
            actual = conda_version(name)
        else:
            actual = version(name)
        if actual not in req.specifier:
            msg = f"{req} is not satisfied by installed version {actual}"
            raise ValueError(msg)

    for tool, package in {
        "nextflow": "nextflow",
        "java": "openjdk",
        "bwa-mem2": "bwa-mem2",
        "bwa-mem3": "bwa-mem3",
        "samtools": "samtools",
        "RNAfold": "viennarna",
    }.items():
        executable = which(tool)
        actual = conda_version(package)
        if executable is None or not Path(executable).resolve().is_relative_to(PREFIX):
            msg = f"{tool} does not resolve inside {PREFIX}"
            raise ValueError(msg)
        print(f"{tool}: {executable} ({package} {actual})")
    jars = list((PREFIX / "share/nextflow/dist").glob("*/*-one.jar"))
    if len(jars) != 1 or conda_version("nextflow") not in jars[0].name:
        msg = f"Expected the JAR paired with the Conda Nextflow launcher, found {jars}"
        raise ValueError(msg)
    print(f"Nextflow JAR: {jars[0]}")


def verify_development_constraints(root: Path) -> None:
    """Check agreement between Pixi and portable development metadata."""
    project = tomllib.loads((root / "pyproject.toml").read_text())
    manifest = tomllib.loads((root / "pixi.toml").read_text())
    features = manifest["feature"]
    for entry in project["dependency-groups"]["dev"]:
        req = Requirement(entry)
        spec = next(
            (
                features[group]["pypi-dependencies"][req.name]
                for group in ("dev", "docs")
                if req.name in features[group]["pypi-dependencies"]
            ),
            None,
        )
        version_spec, extras = (
            (spec.get("version"), set(spec.get("extras", []))) if isinstance(spec, dict) else (spec, set())
        )
        if version_spec != str(req.specifier) or extras != req.extras:
            msg = f"Pixi dev/docs constraints differ from {entry}"
            raise ValueError(msg)

    supported = SpecifierSet(project["project"]["requires-python"])
    current = ".".join(map(str, sys.version_info[:3]))
    if current not in supported:
        msg = f"Locked Python {current} is outside {supported}"
        raise ValueError(msg)


if __name__ == "__main__":
    verify_application()
    if len(sys.argv) == 1:
        verify_development_constraints(Path(__file__).resolve().parents[1])
