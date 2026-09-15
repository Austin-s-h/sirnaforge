"""Build and artifact identity for ``manifest.json``: what produced this run, or why we cannot say.

The manifest recorded ``tool_version: "0.7.1"`` -- a release string, not a build identity -- and no
digest at all for ``report.html``, the artifact everyone actually reads. This module supplies both
under one rule: **never assert an identity that cannot be verified.** Anything unverifiable is a
``{value, state, reason}`` triple whose ``state`` names the reason-class, which is the shape
:meth:`~sirnaforge.config.reference_policy.ReferenceChoice.to_metadata` already emits. The bare
string ``"unknown"`` is forbidden, ``runner.py``'s synthesised ``nogit-<mtime_ns>`` fingerprint may
never surface under a key that reads like a revision, and neither may a HEAD read from a repository
that cannot be shown to hold what it is being asked to name -- the package for ``build.vcs``, the
pipeline for ``build.pipeline_revision``. That is #100's rule for screening evidence, applied to the
build.

Two corollaries, each learned from a field here that broke the rule while stating it. An **inference**
is published with the comparison it rests on, never as the thing it approximates: a cache stamp dated
after this process started bounds a build to the process, not to a run. And an absence is **tested for
rather than declared**: ``uv.lock`` really is missing from the wheel, but saying so permanently made a
dev worktree that has the lock at its verified root report it as unobtainable.

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
from collections.abc import Container, Iterable
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

#: Metrics the report renders per candidate that no weight vector can score, because they carry no
#: ``TermRecord``, mapped to the manifest path that does define each one. They belong in
#: ``reported_not_scored`` all the same: it is a claim about what is reported, not about the scoring
#: vocabulary, and deriving the list from ``TERM_REGISTRY`` alone made ``paired_fraction`` vanish from
#: both manifest lists. Add a term here rather than to a literal, and ``reported_term_partition``
#: keeps it accounted for.
REPORTED_TERMS_WITHOUT_TERM_RECORD: dict[str, str] = {
    "paired_fraction": "design_parameters.filters.max_paired_fraction",
}

#: Env names the pipeline config already uses to detect a container (pipeline/nextflow/config.py).
_CONTAINER_ENV_MARKERS = ("SIRNAFORGE_IN_CONTAINER", "CONTAINER")

_GIT_TIMEOUT_SECONDS = 5

#: The packaged pipeline directory ``runner.py`` fingerprints, and the file a repository must track to
#: be allowed to name its revision. Same guard as the build block's: a recorded HEAD is only the
#: pipeline's revision if the repository reporting it demonstrably holds the pipeline.
_PIPELINE_WORKFLOW_DIR = Path(__file__).resolve().parent / "pipeline" / "nextflow" / "workflows"
_PIPELINE_TRACKED_FILE = "main.nf"

#: When this process started, so a cache stamp can be dated against it. It bounds a build to the
#: *process*, never to one run: a process may execute more than one workflow, which is why the field
#: it feeds is named ``built_in_this_process`` and states its own basis.
_PROCESS_START = datetime.now()

#: ``uv.lock``, present at a verified repo root and absent from the wheel. The absence is conditional,
#: so it is tested for rather than declared permanent.
_LOCK_FILENAME = "uv.lock"

_DESCRIBE_DISTANCE = re.compile(r"-(\d+)-g[0-9a-f]+")
_SHA1_HEX = re.compile(r"[0-9a-f]{40}")
_LOCK_REVISION = re.compile(r"^revision\s*=\s*(\d+)\s*$", re.MULTILINE)

#: Ensembl writes the assembly into its own FASTA headers, e.g.
#: ``>ENST00000632684.1 cdna chromosome:GRCh38:14:22438547:22438554:1``. Reading it there consults the
#: bytes; a table keyed on the source *name* only restates what the source was called.
_HEADER_ASSEMBLY = re.compile(r"\b(?:chromosome|scaffold|primary_assembly|contig|supercontig):([A-Za-z0-9_.\-]+):")


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


def fasta_header_assembly(fasta: Path | None, *, max_records: int = 20) -> str | None:
    """The assembly token these bytes name in their own Ensembl headers, or None.

    The artifact, not a table: ``ENSEMBL_ASSEMBLIES`` is keyed on the source *name*, so it says what
    the reference was called and nothing about what was downloaded. Only the first ``max_records``
    headers are read -- a cDNA file is single-species and the first line of one is a header, so this
    costs a few KB of a 1.5 GB file -- and a disagreement inside the sample yields None rather than a
    coin toss. Returns None for an unreadable, compressed or absent path, which is the ordinary case
    when the manifest is written on the host for a container run's paths.
    """
    if fasta is None:
        return None
    seen: set[str] = set()
    try:
        with fasta.open() as handle:
            records = 0
            for line in handle:
                if not line.startswith(">"):
                    continue
                records += 1
                if records > max_records:
                    break
                match = _HEADER_ASSEMBLY.search(line)
                if match:
                    seen.add(match.group(1))
    except (OSError, UnicodeDecodeError):
        return None
    return seen.pop() if len(seen) == 1 else None


def assembly_identity(*, tabled: str | None, url: str | None, fasta: Path | None) -> dict[str, Any]:
    """Which genome build these bytes are against, preferring the artifact over the source table.

    ``tabled`` is ``ENSEMBL_ASSEMBLIES``' answer for the source *name*, which is an association and
    not an observation -- it was published as a plain fact in the same dict that admits the release
    floated and is unpinned. So the bytes are asked first (their headers carry the assembly), the URL
    second (it contains the assembly for every bundled source), and only then does the table speak,
    under a state that says nothing checked it. A header that contradicts the table is published as
    the header's value: the disagreement is the signal, and the bytes are what was screened.
    """
    observed = fasta_header_assembly(fasta)
    if observed and tabled and observed != tabled:
        return present(
            observed,
            "fasta_header_conflicts_with_source_table",
            f"the cached FASTA's own Ensembl headers name {observed} while the bundled source table names "
            f"{tabled} for this source; the bytes screened are what is published here",
        )
    if observed:
        return present(
            observed,
            "fasta_header",
            "read from the cached FASTA's own Ensembl headers, so this describes the bytes rather than the "
            "name the source was requested under",
        )
    if tabled and url and tabled in url:
        return present(
            tabled,
            "source_table_corroborated_by_url",
            f"the bundled source table names {tabled} for this source and the recorded download URL contains "
            "it; no assembly was readable from the FASTA's own headers here",
        )
    if tabled:
        return present(
            tabled,
            "source_table_unchecked",
            f"the bundled source table names {tabled} for this source name only: no assembly was readable from "
            "the bytes and none appears in the recorded URL, so this is the table's association and not an "
            "observation",
        )
    return absent(
        "not_a_bundled_source",
        "this reference is not a bundled Ensembl source and its headers name no assembly, so no genome build "
        "can be named for the bytes that were screened",
    )


def digest_scope(size_scope: str | None, *, filtered: bool) -> str:
    """What the recorded digest covers, derived from the same evidence as ``size_scope``.

    It was the literal ``decompressed_fasta_full`` beside a derived ``size_scope``, so the two labels
    contradicted each other for any uncompressed reference and the digest said "full" for a filtered
    subset -- while the digest is in fact computed over whatever file the cache entry points at.
    """
    if not size_scope or size_scope == "not_recorded":
        return "not_recorded"
    return f"{size_scope}_filtered_subset" if filtered else f"{size_scope}_whole_file"


def reported_term_partition(registry_terms: Iterable[str], scored_terms: Container[str]) -> dict[str, Any]:
    """Split every reported term into the scored and the reported-but-unscored, over one total universe.

    Derived, never hand-listed: whether a term is scored is a per-run property, and promoting one is
    invisible to a literal (#96). But the scoring registry is not the whole universe of what a reader
    sees, and deriving the lists from ``TERM_REGISTRY`` alone dropped ``paired_fraction`` out of *both*
    of them -- a metric the report renders beside the structure it came from, gating EXCESS_PAIRING,
    reported and unscored and recorded nowhere. ``reported_not_scored`` is a claim about what is
    reported, so a term qualifies whether or not the scoring vocabulary has a record for it.

    The two lists are a partition of that universe by construction, so nothing can fall out of both
    again. Registry order first, so both lists are deterministic.
    """
    universe = list(registry_terms)
    universe += [term for term in REPORTED_TERMS_WITHOUT_TERM_RECORD if term not in universe]
    unscored = [term for term in universe if term not in scored_terms]
    return {
        "scored_terms": [term for term in universe if term in scored_terms],
        "reported_not_scored": unscored,
        "reported_not_scored_source": (
            "sirnaforge.models.scoring_profile.TERM_REGISTRY plus "
            "sirnaforge.provenance.REPORTED_TERMS_WITHOUT_TERM_RECORD, minus the union of this run's "
            "weight-vector terms"
        ),
        # A term with no TermRecord carries its definition elsewhere, so name where: a term listed here
        # and defined nowhere in the manifest is a dead claim, and this is what a test can falsify.
        "reported_not_scored_definitions": {
            term: location for term, location in REPORTED_TERMS_WITHOUT_TERM_RECORD.items() if term in unscored
        },
    }


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


def _repo_root_for(start: Path) -> Path | None:
    """Nearest ancestor of ``start`` carrying a ``.git``, or None when it is not inside one.

    Discovery only. On its own it proves nothing about *what* that repository tracks, which is why
    every caller pairs it with a ``git ls-files --error-unmatch`` check.
    """
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    return None


def _package_repo_root() -> Path | None:
    """Nearest ancestor of this module carrying a ``.git``, or None when the package is not in one."""
    return _repo_root_for(Path(__file__).resolve().parent)


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


def _lock_identity(lock_root: Path | None) -> dict[str, dict[str, Any]]:
    """The resolution lock's digest and revision when it is on disk, or why it is not.

    Both were declared *permanently* absent as "uv.lock is not shipped in the wheel", which is true of
    a wheel or container run and false of the dev worktree the manifest under review came from: the
    lock sits at the repo root this module already verified, readable. A permanent excuse for a
    conditional absence is its own fabrication, so the file is tested for. ``lock_root`` is None
    whenever no repository could be *verified* to hold this package -- an unverified root's lock is
    some other project's resolution.
    """
    if lock_root is None:
        return {
            "lock_sha256": absent(
                "no_verified_repo_to_read_a_lock_from",
                f"{_LOCK_FILENAME} is not shipped in the wheel (pyproject packages only src/sirnaforge) and no "
                "repository was verified to hold this package, so no lock may be attributed to this build",
            ),
            "lock_revision": absent(
                "no_verified_repo_to_read_a_lock_from",
                f"no verified repository root was available to read {_LOCK_FILENAME} from, so its revision "
                "cannot be named",
            ),
        }

    lock = lock_root / _LOCK_FILENAME
    if not lock.is_file():
        return {
            "lock_sha256": absent(
                "not_packaged",
                f"{_LOCK_FILENAME} is absent from {lock_root}: it is not shipped in the wheel, because pyproject "
                "packages only src/sirnaforge",
            ),
            "lock_revision": absent(
                "not_packaged",
                f"{_LOCK_FILENAME} is not on disk beside this build, so its revision cannot be read at runtime",
            ),
        }

    try:
        digest = file_sha256(lock)
        revision = _LOCK_REVISION.search(lock.read_text(encoding="utf-8", errors="replace"))
    except OSError as exc:
        unreadable = absent(
            "lock_unreadable",
            f"{lock} exists but could not be read ({type(exc).__name__}), so the resolution behind this "
            "environment is not identified",
        )
        return {"lock_sha256": unreadable, "lock_revision": dict(unreadable)}

    return {
        "lock_sha256": present(digest, "repo_worktree_lock", f"sha256 of {lock}, at the verified repo root"),
        "lock_revision": present(
            revision.group(1),
            "repo_worktree_lock",
            f"the `revision` key of {lock}, which is uv's own lock-format revision",
        )
        if revision
        else absent(
            "revision_not_declared_in_lock",
            f"{lock} declares no `revision` key, so the lock format it was written in cannot be named",
        ),
    }


def _dependency_identity(lock_root: Path | None) -> dict[str, Any]:
    """Digest the whole installed environment, and name each declared dependency's resolved version.

    The inventory digest plus the count is what distinguishes a four-package environment from a
    two-hundred-package one; ``_lock_identity`` says whether the resolution that produced it can be
    named as well.
    """
    installed: dict[str, str] = {}
    for distribution in importlib.metadata.distributions():
        name = distribution.metadata["Name"]
        if name:
            installed[_normalize_distribution_name(name)] = distribution.version or ""

    inventory = "\n".join(f"{name}=={version}" for name, version in sorted(installed.items()))
    return {
        **_lock_identity(lock_root),
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
    # Only a repository verified to track this package may supply its lock, for the same reason it is
    # the only one allowed to supply its commit.
    lock_root = _package_repo_root() if vcs["repo_verified"] else None
    return {
        "release": release,
        "kind": kind,
        "detail": detail,
        # Repeated from the manifest's top level so `build` reads standalone.
        "tool_version": __version__,
        "vcs": vcs,
        "container": container,
        "dependencies": _dependency_identity(lock_root),
        "runtime": _runtime_identity(),
    }


def _pipeline_repo_verification(recorded: str, workflow_dir: Path) -> tuple[bool, str, str]:
    """Whether a repository verified to hold the pipeline is the one at the recorded revision.

    Returns ``(verified, state, detail)``. The same three failure modes ``_vcs_identity`` handles, plus
    one it cannot have: the recorded HEAD may belong to a repository that does track the pipeline but
    has since moved, and a stale SHA is no more this run's revision than a foreign one.
    """
    root = _repo_root_for(workflow_dir)
    if root is None:
        return (
            False,
            "pipeline_not_a_git_checkout",
            f"no ancestor of {workflow_dir} carries a .git, so the repository that reported this revision "
            "cannot be found and its claim on the pipeline cannot be checked",
        )

    tracked = workflow_dir / _PIPELINE_TRACKED_FILE
    try:
        _git("ls-files", "--error-unmatch", str(tracked), cwd=root)
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return (
            False,
            "git_binary_missing",
            f"git could not be run to check which repository holds the pipeline ({type(exc).__name__}), so the "
            "recorded revision cannot be attributed to it",
        )
    except subprocess.CalledProcessError:
        return (
            False,
            "repo_does_not_track_pipeline",
            f"the nearest repository at {root} does not track {tracked.name}, so it cannot vouch for the "
            "pipeline this run executed",
        )

    try:
        head = _git("rev-parse", "HEAD", cwd=root)
    except (FileNotFoundError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
        return (
            False,
            "git_command_failed",
            f"the repository at {root} tracks {tracked.name} but its HEAD could not be read "
            f"({type(exc).__name__}), so the recorded revision was not confirmed",
        )

    if head != recorded:
        return (
            False,
            "recorded_revision_is_not_the_pipeline_repo_head",
            f"the repository at {root} that does track {tracked.name} is at {head}, not the recorded revision, "
            "so the recorded value came from elsewhere or from an earlier checkout",
        )

    return (
        True,
        "git_worktree",
        f"git rev-parse HEAD of {root}, verified by git ls-files to track {tracked.name}",
    )


def pipeline_revision_identity(recorded: str | None, *, workflow_dir: Path | None = None) -> dict[str, Any]:
    """Translate ``runner.py``'s pipeline revision, which is neither always a revision nor ever verified.

    Two hazards, not one. ``_detect_pipeline_revision`` falls back to ``nogit-<mtime_ns>``, and a
    filesystem mtime under a key called ``commit`` is a fabricated identity, so it is republished only
    under a key that says fingerprint. But it also walks up from the packaged workflow directory to the
    first ancestor ``.git`` and returns *that* repository's HEAD -- the foreign-commit hazard
    ``_vcs_identity`` guards against with ``repo_verified``, and three ancestors of an installed package
    carry ``.git``. A well-formed SHA1 is therefore only published as ``commit`` when a repository
    verified to track ``pipeline/nextflow/workflows/main.nf`` is at that HEAD. Otherwise it is some
    repository's HEAD and appears under ``ancestor_repo_head``, which is what it is.
    """
    if recorded and _SHA1_HEX.fullmatch(recorded):
        verified, state, detail = _pipeline_repo_verification(recorded, workflow_dir or _PIPELINE_WORKFLOW_DIR)
        if verified:
            return {"commit": present(recorded, state, detail), "repo_verified": True}
        return {
            "commit": absent(state, detail),
            # The value is real; its attribution to this pipeline is not. Never under `commit`.
            "ancestor_repo_head": recorded,
            "repo_verified": False,
        }

    if recorded and recorded.startswith("nogit-"):
        return {
            "commit": absent(
                "not_a_git_checkout",
                ".git never enters the image: the Dockerfile copies only pyproject.toml, uv.lock, README, "
                "LICENSE and src/",
            ),
            "workflow_dir_fingerprint": f"mtime-{recorded.removeprefix('nogit-')}",
            "repo_verified": False,
        }

    return {
        "commit": absent(
            "pipeline_revision_undetermined",
            "neither a git revision nor a directory fingerprint could be read for the pipeline directory",
        ),
        "repo_verified": False,
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
            "built_in_this_process": absent(
                "no_stamp_for_this_index",
                "no cache stamp beside this prefix records when the index was built, so nothing about it can be "
                "attributed to this process either way",
            ),
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
        "built_in_this_process": _stamped_in_this_process(stamp.get("stamped_at")),
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


def _stamped_in_this_process(stamped_at: Any) -> dict[str, Any]:
    """Date the stamp against this process, and say that that is all the comparison establishes.

    A timestamp inequality is not authorship. It was published as ``built_by_this_run``, a positive
    claim no comparison against a module-import time can support -- one process may run more than one
    workflow -- and an unparseable timestamp fell out as a bare ``null`` with no reason-class. So the
    field is named for the process, and the reason states the inequality it rests on.
    """
    try:
        stamp_time = datetime.fromisoformat(str(stamped_at))
    except (TypeError, ValueError):
        return absent(
            "stamp_timestamp_unreadable",
            f"the cache stamp's stamped_at ({stamped_at!r}) is not an ISO timestamp, so when these index bytes "
            "were built cannot be established at all",
        )
    if stamp_time >= _PROCESS_START:
        return present(
            True,
            "stamp_dated_after_process_start",
            f"the stamp is dated {stamp_time.isoformat()}, at or after this process started "
            f"({_PROCESS_START.isoformat()}); a process may run more than one workflow, so this bounds the build "
            "to the process and not to this run",
        )
    return present(
        False,
        "stamp_dated_before_process_start",
        f"the stamp is dated {stamp_time.isoformat()}, before this process started "
        f"({_PROCESS_START.isoformat()}), so these index bytes were built by an earlier run and reused",
    )


def report_html_artifact(*, expected: bool, path: str = "report.html") -> dict[str, Any]:
    """How ``report.html`` is attested, at exactly the point a reader looks for its digest.

    Not a second write pass: ``build_payload`` reads manifest.json and the report embeds its content,
    so rewriting the manifest afterwards would leave the artifact everyone reads quoting a manifest
    that no longer exists on disk. The sidecar attests both instead. When no report was expected at
    all, "not attested" is the wrong word -- nothing was asked for.

    ``expected`` is the caller's *intent*, which is all this manifest can know: it is written before
    the render runs. So the claim is stated as pending and the outcome is the sidecar's to record --
    an earlier version keyed this on the caller merely passing a path, which left manifest.json
    asserting a report.html that a failed render never wrote.
    """
    if not expected:
        return {
            "path": None,
            "attested": False,
            "state": "not_rendered_by_this_entry_point",
            "reason": "this entry point publishes candidate tables and a manifest but renders no HTML report",
        }
    return {
        "path": path,
        "attested": False,
        "state": "pending_sidecar_attestation",
        "attestation_artifact": "report_manifest.json",
        # A named absence, not False: whether the render succeeded is unknowable here, and the sidecar
        # is written either way, so an absent sidecar leaves this claim unresolved rather than met.
        "rendered": absent(
            "outcome_follows_this_manifest",
            "the render runs after this manifest is written, so whether report.html exists is recorded in "
            "report_manifest.json; this manifest asserts only that one was asked for",
        ),
        "reason": (
            "report.html is rendered after this manifest is written and embeds this manifest's content, so its "
            "digest is published in report_manifest.json rather than by rewriting the bytes the report quoted"
        ),
    }


def write_report_manifest(
    report_html: Path,
    manifest_json: Path,
    out: Path,
    *,
    render_error: str | None = None,
) -> Path:
    """Attest ``report.html`` and the manifest it quotes, in a sidecar rather than by a rewrite.

    Both digests in one file, so the pair is mutually verifiable with no circularity: a mismatch on
    ``manifest_json`` means the manifest was rewritten after the report rendered. Deliberately the
    minimal form of the ``report_manifest`` artifact #111 proposes, so #111 can add fields rather
    than undo a name.

    Written whether or not the render succeeded, because this is the only place the outcome can be
    stated: the manifest claims the report as pending, so a sidecar withheld on failure would leave
    that claim permanently unresolved. ``exists: false`` plus ``render_error`` is the record of a
    render that raised; whatever bytes a failed render did leave -- partial, or an earlier run's -- are
    digested as what is actually on disk beside the error, never asserted as this run's report.
    """
    report = _attest_file(report_html)
    # The render's outcome, not the caller's intent. Null means no exception was reported -- `exists`
    # is what says whether report.html is there.
    report["render_error"] = render_error
    payload = {
        "schema_version": REPORT_MANIFEST_SCHEMA_VERSION,
        "written_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "report_html": report,
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
