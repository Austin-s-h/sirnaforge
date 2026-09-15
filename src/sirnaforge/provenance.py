"""Build and artifact identity for ``manifest.json``: what produced this run, or why we cannot say.

The manifest recorded ``tool_version: "0.7.1"`` -- a release string, not a build identity -- and no
digest at all for ``report.html``, the artifact everyone actually reads. This module supplies both
under one rule: **never assert an identity that cannot be verified.** Anything unverifiable is a
``{value, state, reason}`` triple whose ``state`` names the reason-class, which is the shape
:meth:`~sirnaforge.config.reference_policy.ReferenceChoice.to_metadata` already emits. The bare
string ``"unknown"`` is forbidden, and ``runner.py``'s synthesised ``nogit-<mtime_ns>`` fingerprint
may never surface under a key that reads like a revision. That is #100's rule for screening
evidence, applied to the build.

Imports stay stdlib-only apart from ``__version__`` and one leaf cache util, so the manifest builder
can reach this module without pulling ``sirnaforge.data`` into a fresh import cycle.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import time
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

from sirnaforge import __version__
from sirnaforge.utils.cache_utils import read_artifact_stamp

PROVENANCE_SCHEMA_VERSION = "1.0"
REPORT_MANIFEST_SCHEMA_VERSION = "1.0"

#: The aligner the screen is built around. Named as a constant because its *version* is a declared
#: non-availability everywhere it appears, and the two must not drift apart.
ALIGNER_NAME = "bwa-mem2"

#: Digests recorded by ``cache_utils`` are md5, while ``manifest.files.*`` is sha256. Republishing
#: one beside the other without naming the algorithm invites a false mismatch.
CACHE_DIGEST_ALGORITHM = "md5"

#: Env names the pipeline config already uses to detect a container (pipeline/nextflow/config.py).
_CONTAINER_ENV_MARKERS = ("SIRNAFORGE_IN_CONTAINER", "CONTAINER")

_GIT_TIMEOUT_SECONDS = 5

#: When this process started, so a cache stamp written afterwards can be attributed to this run
#: rather than published as a bare ``built_by_this_run: false``.
_PROCESS_START = datetime.now()

_DESCRIBE_DISTANCE = re.compile(r"-(\d+)-g[0-9a-f]+")
_SHA1_HEX = re.compile(r"[0-9a-f]{40}")


def present(value: Any, state: str, reason: str) -> dict[str, Any]:
    """A verified fact: the value, the authority that supplied it, and how it was obtained."""
    return {"value": value, "state": state, "reason": reason}


def absent(state: str, reason: str) -> dict[str, Any]:
    """A named non-availability. ``state`` is the reason-class; never the bare string ``unknown``."""
    return {"value": None, "state": state, "reason": reason}


def file_sha256(path: Path) -> str:
    """SHA-256 of a file, for integrity (non-security) attestation."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(*args: str, cwd: Path) -> str:
    """Run one git command in ``cwd`` and return its stdout, raising on anything else.

    The single seam for every git invocation here, so a test can observe both hazards that are
    otherwise invisible: a binary that cannot be run, and a repository that vouches for somebody
    else's commit.
    """
    result = subprocess.run(  # nosec B603
        ["git", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=True,
        timeout=_GIT_TIMEOUT_SECONDS,
    )
    return result.stdout.strip()


def _package_repo_root() -> Path | None:
    """Nearest ancestor of this module carrying a ``.git``, or None when the package is not in one."""
    start = Path(__file__).resolve().parent
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    return None


def _vcs_unavailable(state: str, reason: str) -> dict[str, Any]:
    """Every VCS field as one named non-availability.

    ``dirty`` goes null rather than False on purpose: False is a positive claim that the tree was
    clean, and a `git status` that never ran cannot make it.
    """
    return {
        "commit": absent(state, reason),
        "describe": absent(state, reason),
        "dirty": absent(state, reason),
        "repo_verified": False,
    }


def _vcs_identity() -> dict[str, Any]:
    """Resolve the commit, describe string and cleanliness of the repo that tracks this package."""
    root = _package_repo_root()
    if root is None:
        return _vcs_unavailable(
            "not_a_git_checkout",
            "no ancestor of the installed package carries a .git, so this build has no revision to name",
        )

    package_init = Path(__file__).resolve().parent / "__init__.py"
    try:
        # The guard runner.py lacks: three ancestor directories of the installed package carry .git,
        # so walking up until one appears will happily attribute a foreign commit to this build.
        _git("ls-files", "--error-unmatch", str(package_init), cwd=root)
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return _vcs_unavailable(
            "git_binary_missing",
            f"git could not be run to identify this build ({type(exc).__name__}), so no revision was read",
        )
    except subprocess.CalledProcessError:
        return _vcs_unavailable(
            "repo_does_not_track_package",
            f"the nearest repository at {root} does not track {package_init.name}, so it cannot vouch for this build",
        )

    try:
        commit = _git("rev-parse", "HEAD", cwd=root)
        describe = _git("describe", "--tags", "--always", "--dirty", cwd=root)
        # --untracked-files=no: an untracked scratch file is not a modification of what git tracks.
        porcelain = _git("status", "--porcelain", "--untracked-files=no", cwd=root)
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return _vcs_unavailable(
            "git_binary_missing",
            f"git could not be run to identify this build ({type(exc).__name__}), so no revision was read",
        )
    except subprocess.CalledProcessError as exc:
        return _vcs_unavailable(
            "git_command_failed",
            f"the repository at {root} was verified but git exited {exc.returncode}, so no revision was read",
        )

    dirty = bool(porcelain)
    return {
        "commit": present(commit, "git_worktree", "resolved from the repo that tracks src/sirnaforge/__init__.py"),
        # Never falls back to __version__: disagreement between the two is the signal.
        "describe": present(describe, "git_worktree", "git describe --tags --always --dirty"),
        "dirty": present(
            dirty,
            "git_status_porcelain",
            "git status reported modified tracked files" if dirty else "working tree had no modified tracked files",
        ),
        "repo_verified": True,
    }


def _container_identity() -> dict[str, Any]:
    """Whether this is a container run, and the digest it cannot claim.

    A tag is mutable, so it is never substituted for a digest: OCI labels are unreadable from inside
    a running container and no build-arg carries the digest, which makes both an explicit absence.
    """
    detection = "none"
    for name in _CONTAINER_ENV_MARKERS:
        if os.getenv(name):
            detection = f"env:{name}"
            break
    if detection == "none" and Path("/.dockerenv").exists():
        detection = "dockerenv"

    build_version = os.getenv("BUILD_VERSION")
    return {
        "in_container": detection != "none",
        "detection": detection,
        "image_digest": absent(
            "not_injected",
            "no build-arg carries the digest and OCI labels are unreadable from inside a running container",
        ),
        "image_ref": absent("not_injected", "a tag is mutable and is never substituted for a digest"),
        "build_version": present(
            build_version,
            "env:BUILD_VERSION",
            "recorded verbatim; the literal 0.0.0+local means the image came from docker-dev.yml, "
            "which passes no VERSION build-arg",
        )
        if build_version
        else absent("not_injected", "no BUILD_VERSION build-arg reached this environment, so the image is unnamed"),
    }


def _normalize_distribution_name(name: str) -> str:
    """PEP 503 normalisation, so a declared name matches its installed distribution."""
    return re.sub(r"[-_.]+", "-", name).lower()


def _declared_dependencies() -> list[str]:
    """Runtime dependency names as the package metadata declares them.

    Derived from ``importlib.metadata.requires`` rather than hand-listed: a hand-list is the same
    drift class as the ``reported_not_scored`` literal, and a new dependency would escape the record.
    Extras are excluded because they are not what a default run installs.
    """
    try:
        requirements = importlib.metadata.requires("sirnaforge") or []
    except importlib.metadata.PackageNotFoundError:
        return []
    declared: list[str] = []
    for requirement in requirements:
        specifier, _, marker = requirement.partition(";")
        if "extra" in marker:
            continue
        name = re.split(r"[\[<>=!~ (]", specifier.strip(), maxsplit=1)[0].strip()
        if name:
            declared.append(_normalize_distribution_name(name))
    return sorted(set(declared))


def _dependency_identity() -> dict[str, Any]:
    """Digest the whole installed environment, and name each declared dependency's resolved version.

    ``uv.lock`` is not shipped in the wheel, so the lock digest is an explicit absence rather than a
    field a container run silently drops. The inventory digest plus the count is what distinguishes
    a four-package environment from a two-hundred-package one.
    """
    installed: dict[str, str] = {}
    for distribution in importlib.metadata.distributions():
        name = distribution.metadata["Name"]
        if name:
            installed[_normalize_distribution_name(name)] = distribution.version or ""

    inventory = "\n".join(f"{name}=={version}" for name, version in sorted(installed.items()))
    return {
        "lock_sha256": absent(
            "not_packaged",
            "uv.lock is not shipped in the wheel: pyproject packages only src/sirnaforge",
        ),
        "lock_revision": absent(
            "not_packaged",
            "uv.lock is not shipped in the wheel, so its revision cannot be read at runtime",
        ),
        "distributions_sha256": hashlib.sha256(inventory.encode("utf-8")).hexdigest(),
        "distribution_count": len(installed),
        # A declared-but-uninstalled dependency maps to null, never to omission.
        "declared": {name: installed.get(name) for name in _declared_dependencies()},
    }


def _runtime_identity() -> dict[str, Any]:
    """The interpreter, the ViennaRNA that drew the report's structures, and the aligner's absence."""
    try:
        import RNA  # noqa: PLC0415

        viennarna = present(
            str(RNA.__version__),
            "python_bindings",
            "RNA.__version__ of the module this process imported",
        )
    except Exception:
        viennarna = absent(
            "python_bindings_unavailable",
            "the ViennaRNA Python bindings are not importable here, so no version can be named",
        )
    # Never quote `RNAfold --version`: the image carries two ViennaRNA installs and no code execs the
    # binary, so its banner would describe something this run did not use.
    viennarna["binary_probed"] = False

    return {
        "python": {"version": platform.python_version(), "implementation": platform.python_implementation()},
        "viennarna": viennarna,
        # path-present must never read as version-known: a conda pin is a constraint, not an observation.
        "aligner": {
            "name": ALIGNER_NAME,
            "path": shutil.which(ALIGNER_NAME),
            "value": None,
            "state": "not_captured",
            "reason": (
                f"no code path captures a {ALIGNER_NAME} version; the pipeline's versions.yml recorded the "
                "launcher banner plus 'not available'"
            ),
        },
        "probed_not_used": {"samtools": {"invoked": False}},
    }


def _release_verdict(vcs: dict[str, Any], container: dict[str, Any]) -> tuple[bool, str, str]:
    """Decide ``release``/``kind``/``detail``, failing closed on every undecidable case.

    A run is a release build only if a repository that verifiably tracks this package reports a clean
    tree at a tag equal to ``v{__version__}``. Everything else -- including everything undecidable --
    is false, with a kind and a one-sentence detail saying which counts it failed on.
    """
    expected_tag = f"v{__version__}"

    if not vcs["repo_verified"]:
        state = str(vcs["commit"]["state"])
        if container["in_container"]:
            kind = "dev_container"
        elif state in {"not_a_git_checkout", "repo_does_not_track_package"}:
            kind = "installed_no_vcs"
        else:
            kind = "unknown_environment"
        detail = (
            f"not a release build: {vcs['commit']['reason']}, and no build-time commit or image digest was "
            f"injected to stand in for it, so __version__ {__version__} is unsupported by any revision"
        )
        return False, kind, detail

    describe = str(vcs["describe"]["value"])
    dirty = bool(vcs["dirty"]["value"])
    if describe == expected_tag and not dirty:
        return (
            True,
            "release_tag",
            f"release build: git describe is {describe}, equal to v{__version__}, on a clean tree",
        )

    shortfalls: list[str] = []
    if describe != expected_tag:
        distance = _DESCRIBE_DISTANCE.search(describe)
        if distance:
            shortfalls.append(f"it is {distance.group(1)} commits past the tag it describes")
        shortfalls.append(f"no tag on this branch equals {expected_tag}")
    if dirty:
        shortfalls.append("the working tree has modified tracked files")
    detail = f"dev build: git describe is {describe} against __version__ {__version__} -- " + "; ".join(shortfalls)
    return False, "dev_worktree", detail


@lru_cache(maxsize=1)
def build_identity() -> dict[str, Any]:
    """The manifest's ``build`` block: what produced this run, and what it cannot claim.

    Cached for the process because it shells out to git; ``build_identity.cache_clear()`` is the way
    a test that changes the environment gets a fresh answer.
    """
    vcs = _vcs_identity()
    container = _container_identity()
    release, kind, detail = _release_verdict(vcs, container)
    return {
        "release": release,
        "kind": kind,
        "detail": detail,
        # Repeated from the manifest's top level so `build` reads standalone.
        "tool_version": __version__,
        "vcs": vcs,
        "container": container,
        "dependencies": _dependency_identity(),
        "runtime": _runtime_identity(),
    }


def pipeline_revision_identity(recorded: str | None) -> dict[str, Any]:
    """Translate ``runner.py``'s pipeline revision, which is not always a revision.

    ``_detect_pipeline_revision`` falls back to ``nogit-<mtime_ns>``, and a filesystem mtime under a
    key called ``commit`` is a fabricated identity. It is republished only under a key that says
    fingerprint.
    """
    if recorded and _SHA1_HEX.fullmatch(recorded):
        return {
            "commit": present(
                recorded,
                "git_worktree",
                "git rev-parse HEAD of the repository holding the pipeline directory",
            )
        }

    if recorded and recorded.startswith("nogit-"):
        return {
            "commit": absent(
                "not_a_git_checkout",
                ".git never enters the image: the Dockerfile copies only pyproject.toml, uv.lock, README, "
                "LICENSE and src/",
            ),
            "workflow_dir_fingerprint": f"mtime-{recorded.removeprefix('nogit-')}",
        }

    return {
        "commit": absent(
            "pipeline_revision_undetermined",
            "neither a git revision nor a directory fingerprint could be read for the pipeline directory",
        )
    }


def index_attestation(index_prefix: Path, *, reference_fasta_name: str | None = None) -> dict[str, Any]:
    """Whether the index bytes this run aligned against can be shown to be the bytes someone built.

    bwa-mem2 aligns against the index files alone, so an index nobody stamped is unattested bytes and
    is named as such. ``built_from_digest`` is what proves index and reference are the same release;
    ``members`` is never collapsed into one "index digest", because a multi-file artifact proves its
    own completeness only member by member, and ``digest_mode`` is never dropped -- every member of a
    real human index exceeds the full-digest budget, so a sampled digest must not read as a whole-file
    one.
    """
    aligner_version = absent(
        "not_captured",
        f"no code path captures a {ALIGNER_NAME} version; see build.runtime.aligner",
    )
    stamp = read_artifact_stamp(index_prefix)
    if stamp is None:
        return {
            "attested": False,
            "state": "no_stamp_for_this_index",
            "reason": (
                "this index prefix has no .sirnaforge-cache.json sidecar, so the index bytes this run aligned "
                "against are unattested"
            ),
            "aligner": ALIGNER_NAME,
            "built_by_this_run": False,
            "built_from_digest": None,
            "aligner_version": aligner_version,
            "members": [],
        }

    inputs: dict[str, Any] = dict(stamp.get("inputs") or {})
    built_from = _stamped_input_digest(inputs, reference_fasta_name)
    members = [
        {
            "name": name,
            "size": fingerprint.get("size"),
            "digest": fingerprint.get("digest"),
            "digest_algorithm": CACHE_DIGEST_ALGORITHM,
            "digest_mode": fingerprint.get("digest_mode"),
        }
        for name, fingerprint in sorted((stamp.get("outputs") or {}).items())
        if isinstance(fingerprint, dict)
    ]
    return {
        "attested": True,
        "state": "artifact_stamp",
        "reason": "the cache stamp beside this prefix records the inputs it was built from and its own bytes",
        "aligner": ALIGNER_NAME,
        "built_by_this_run": _stamped_by_this_process(stamp.get("stamped_at")),
        "built_from_digest": built_from,
        "aligner_version": aligner_version,
        "members": members,
    }


def _stamped_input_digest(inputs: dict[str, Any], reference_fasta_name: str | None) -> str | None:
    """The stamp's digest for the reference FASTA, which is what binds index to reference."""
    if reference_fasta_name:
        for key, digest in inputs.items():
            if str(key).endswith(f":{reference_fasta_name}"):
                return str(digest)
    if len(inputs) == 1:
        return str(next(iter(inputs.values())))
    return None


def _stamped_by_this_process(stamped_at: Any) -> bool | None:
    """Whether the stamp was written after this process started; None when it cannot be read."""
    try:
        return datetime.fromisoformat(str(stamped_at)) >= _PROCESS_START
    except (TypeError, ValueError):
        return None


def report_html_artifact(*, rendered: bool, path: str = "report.html") -> dict[str, Any]:
    """How ``report.html`` is attested, at exactly the point a reader looks for its digest.

    Not a second write pass: ``build_payload`` reads manifest.json and the report embeds its content,
    so rewriting the manifest afterwards would leave the artifact everyone reads quoting a manifest
    that no longer exists on disk. The sidecar attests both instead. When no report was rendered at
    all, "not attested" is the wrong word -- nothing was expected.
    """
    if not rendered:
        return {
            "path": None,
            "attested": False,
            "state": "not_rendered_by_this_entry_point",
            "reason": "this entry point publishes candidate tables and a manifest but renders no HTML report",
        }
    return {
        "path": path,
        "attested": False,
        "state": "attested_by_sidecar",
        "attestation_artifact": "report_manifest.json",
        "reason": (
            "report.html is rendered after this manifest is written and embeds this manifest's content, so its "
            "digest is published in report_manifest.json rather than by rewriting the bytes the report quoted"
        ),
    }


def write_report_manifest(report_html: Path, manifest_json: Path, out: Path) -> Path:
    """Attest ``report.html`` and the manifest it quotes, in a sidecar rather than by a rewrite.

    Both digests in one file, so the pair is mutually verifiable with no circularity: a mismatch on
    ``manifest_json`` means the manifest was rewritten after the report rendered. Deliberately the
    minimal form of the ``report_manifest`` artifact #111 proposes, so #111 can add fields rather
    than undo a name.
    """
    payload = {
        "schema_version": REPORT_MANIFEST_SCHEMA_VERSION,
        "written_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "report_html": _attest_file(report_html),
        "manifest_json": _attest_file(manifest_json),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        json.dump(payload, handle, indent=2)
    return out


def _attest_file(path: Path) -> dict[str, Any]:
    """Path, size and sha256 of one artifact, or a named absence when it is not there."""
    if not path.exists():
        return {"path": path.name, "exists": False, "size_bytes": None, "sha256": None}
    return {
        "path": path.name,
        "exists": True,
        "size_bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }
