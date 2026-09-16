"""Identifiers that must not appear in this public repository.

Comments, docstrings and fixture values in this tree once justified a threshold or a design decision by
citing measurements from an internal run, and named that run's subject. The measurements are the
evidence and stay; the identifier does not, because a public reader learns from it something the
repository is not entitled to publish. The established phrasing is "one internal run" / "one internal
reference set".

**The forbidden terms are stored as digests, not as text.** A guard that spells out what it forbids
does not remove the disclosure -- it concentrates it, and adds a signpost saying which word mattered.
That is strictly worse than the buried comment it replaced, so this file holds
:data:`FORBIDDEN_DIGESTS` and never the pre-images.

What that buys and what it does not, stated plainly so nobody over-trusts it:

* It buys removal of the plaintext and of the signpost. Reading this file tells you that some token is
  forbidden; it does not tell you which, or what kind of thing it names.
* It does **not** buy secrecy. These are unsalted SHA-256 digests of short tokens, so anyone with a
  candidate word list can confirm a guess offline. Salting would defeat the point, because the salt
  would have to live here too.

So this is a **regression guard**, not a secret store: its job is to stop a term that has already been
removed from coming back. The control that keeps a new one out in the first place is review.

Detection is per token, not per substring -- the price of holding digests rather than text, since a
digest cannot be matched against the inside of a word. A line is lower-cased and split on non-alphanumeric
characters, so a reference is caught in prose, in ``snake_case`` and in hyphenated identifiers alike, but a
term welded into a longer word with no boundary would not be. A leak is written as a word, so that is the
right granularity; it is a limitation rather than a hole, and it is recorded here rather than discovered
later.

``SIRNAFORGE_FORBIDDEN_DIGESTS`` extends the set with a comma-separated list of hex digests, so CI can
police more terms than a public file should enumerate even in digest form.

**Scope is the whole repository, not two trees of ``*.py``.** The first version of this file scanned
``rglob("*.py")`` under ``src`` and ``tests``, and an audit defeated it three separate ways, each
reproduced for real:

* non-``.py`` files *inside* the trees it claimed to scan -- a comment in the Nextflow config, a line
  in ``tests/data/benchmarks/README.md``, and a planted ``.txt`` next to it. A Nextflow module comment
  is a very plausible home for exactly this kind of run-calibration note.
* trees it never looked at -- ``CHANGELOG.md``, ``docs/``, ``README.md``, ``scripts/``, ``.github/``.
* ``.ipynb``, where a notebook cell is plain JSON text and reads like prose to everyone but a glob.

So discovery *excludes* rather than listing what to include: a new file kind, a new top-level directory
or a new prose format is policed by default and has to be argued out of scope. Binary files are
recognised by content, not by an extension list, for the same reason.
:func:`test_scan_sees_the_file_kinds_that_defeated_the_narrow_glob` pins all three defeats, using a
sentinel term of its own so the fixture does not need the real one either.

**What is scanned is what git tracks, not what happens to be on disk.** The concern is what the public
repository contains, and those are different sets: walking the filesystem flagged 29 sites in
gitignored local run output -- ``workflow_test_debug_*/manifest.json`` and a scratch contract under
``.dev/`` -- none of which is published. A guard that fails on any developer's leftover run output is a
guard that gets deleted, and it would still have said nothing about the published tree. So the file
list comes from ``git ls-files``, which also makes almost every path exclusion unnecessary: build
artefacts, caches, worktrees and Nextflow state are untracked already. ``PRUNED_DIRS`` survives for the
filesystem fallback used when git is unavailable, and for the planted-tree fixture.

Cost, measured on this tree: 359 tracked files / ~23 MB, under 200 ms warm. That is fast-tier work, so
nothing is excluded for size -- including the ``tests/unit/data/baseline_0_7_1`` ``.tsv`` baselines and
the multi-megabyte benchmark CSVs under ``tests/data``.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

#: SHA-256 of each forbidden token, lower-cased. Pre-images are deliberately absent; see the module
#: docstring for why, and for what this does and does not protect.
FORBIDDEN_DIGESTS: frozenset[str] = frozenset(
    {
        "dda32cf24f068ffc50751d8ab393323a3744d97f854bd4e29fd002429d37be72",
        "812d73c14960009f6440a8ef28ca7af79adea5e299fe13f068527616bd2b66b0",
    }
)

#: Extra digests for CI, comma-separated hex. Lets a private runner police terms a public file should
#: not enumerate even as digests.
_ENV_DIGESTS = "SIRNAFORGE_FORBIDDEN_DIGESTS"

#: Tokens are alphanumeric runs. Splitting here is what makes a digest usable at all: there is nothing
#: to hash until the line is cut into candidate words.
_TOKEN = re.compile(r"[a-z0-9]+")

#: Directories never walked, each with the reason it is not authored, published prose. This list is
#: the whole scope decision: everything not named here is scanned.
PRUNED_DIRS: dict[str, str] = {
    ".git": "git internals -- published history carries the pre-scrub blobs and only a rewrite can fix that",
    ".claude": "gitignored agent worktrees; each holds its own checkout, which would be scanned twice",
    ".venv": "installed third-party packages, not this repository's prose",
    "_build": "docs/_build is a Sphinx artefact -- a rendered copy of docs/, which is scanned",
    "__pycache__": "compiled bytecode",
    ".mypy_cache": "tool cache",
    ".ruff_cache": "tool cache",
    ".pytest_cache": "tool cache",
    ".pytest_tmp": "tool scratch",
    ".nextflow": "Nextflow run state",
    "work": "Nextflow work dir (gitignored)",
    "node_modules": "vendored JS dependencies",
    "htmlcov": "coverage report artefact",
    "dist": "build artefact",
}

#: Files that must turn up in discovery. Without them a broken walk would scan nothing and pass
#: silently -- and these are precisely the kinds the ``*.py`` glob could not see.
REQUIRED_DISCOVERED_FILES = frozenset(
    {
        "README.md",
        "CHANGELOG.md",
        "Makefile",
        "pyproject.toml",
        "docs/index.rst",
        "docs/developer/zfn_nextflow_bridge_validation.md",
        "notebooks/zfn_backend_runtime_comparison.ipynb",
        "scripts/test_nextflow_integration.sh",
        ".github/workflows/ci.yml",
        "src/sirnaforge/reporting/payload.py",
        "src/sirnaforge/pipeline/nextflow/workflows/nextflow.config",
        "src/sirnaforge/pipeline/nextflow/workflows/modules/local/offtarget_analysis.nf",
        "tests/data/benchmarks/README.md",
        "tests/unit/data/baseline_0_7_1/offtarget_hits_sample.tsv",
    }
)
MIN_DISCOVERED_FILES = 250


def forbidden_digests() -> frozenset[str]:
    """The digests in force: the declared set plus anything the environment adds."""
    extra = {value.strip().lower() for value in os.environ.get(_ENV_DIGESTS, "").split(",") if value.strip()}
    return FORBIDDEN_DIGESTS | extra


def _is_binary(raw: bytes) -> bool:
    """A NUL byte in the first 8 KB. Decided by content so no novel extension can hide prose."""
    return b"\x00" in raw[:8192]


def _line_offends(line: str, digests: frozenset[str]) -> bool:
    """True when any alphanumeric token on the line hashes to a forbidden digest."""
    return any(hashlib.sha256(token.encode()).hexdigest() in digests for token in _TOKEN.findall(line.lower()))


def _tracked_paths(root: Path) -> list[Path] | None:
    """Every file git tracks, or None when git cannot answer.

    ``git ls-files`` is the right question because the concern is what the repository publishes. It
    also excludes build artefacts, caches, agent worktrees and Nextflow state for free, since none of
    them is tracked.
    """
    try:
        listing = subprocess.run(  # noqa: S603 - fixed argv, no shell
            ["git", "-C", str(root), "ls-files", "-z"],
            capture_output=True,
            check=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    names = [name for name in listing.stdout.decode("utf-8", errors="replace").split("\0") if name]
    return [root / name for name in names]


def _walked_paths(root: Path) -> list[Path]:
    """Fallback discovery for a tree git does not know about, and for the planted-tree fixture."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in PRUNED_DIRS and not d.endswith(".egg-info")]
        found.extend(Path(dirpath) / name for name in sorted(filenames))
    return found


def _scan(root: Path, digests: frozenset[str], *, tracked_only: bool = True) -> tuple[list[str], list[str]]:
    """Return ``(offenders, scanned)`` for one tree: ``path:line`` hits, and every file read.

    Symlinks are skipped: they either point back inside the tree (already walked) or outside it.
    """
    paths = (_tracked_paths(root) if tracked_only else None) or _walked_paths(root)
    offenders: list[str] = []
    scanned: list[str] = []
    for path in paths:
        if path.is_symlink() or not path.is_file():
            continue
        raw = path.read_bytes()
        if _is_binary(raw):
            continue
        rel = path.relative_to(root).as_posix()
        scanned.append(rel)
        for n, line in enumerate(raw.decode("utf-8", errors="replace").splitlines(), start=1):
            if _line_offends(line, digests):
                offenders.append(f"{rel}:{n}")
    return offenders, sorted(scanned)


_OFFENDERS, _SCANNED = _scan(REPO_ROOT, forbidden_digests())


@pytest.mark.unit
def test_no_forbidden_identifier_anywhere_in_the_repository() -> None:
    """No text file in the repository carries a forbidden identifier, in any case.

    Case-insensitive by construction -- tokens are lower-cased before hashing -- because both an
    upper-case and a title-case spelling were present when this was first cleaned up.
    """
    assert not _OFFENDERS, (
        f"{len(_OFFENDERS)} site(s) carry a forbidden identifier. Keep the measurement, drop the "
        f"identifier -- say 'one internal run' instead: {_OFFENDERS}"
    )


@pytest.mark.unit
def test_this_file_does_not_carry_the_pre_images() -> None:
    """The guard must not be the disclosure.

    A denylist written in plaintext concentrates what it set out to remove and signposts which word
    mattered. So the module is held to its own rule: every digest in force must fail to match anything
    in this file, including the digest strings themselves.
    """
    offenders = [
        n
        for n, line in enumerate(Path(__file__).read_text(encoding="utf-8").splitlines(), start=1)
        if _line_offends(line, forbidden_digests())
    ]

    assert not offenders, f"the guard spells out a term it forbids, on line(s) {offenders}"


@pytest.mark.unit
def test_discovery_walked_the_whole_tree() -> None:
    """Guard the walk itself: a bad prune or a stale root would make the scan above vacuous."""
    missing = REQUIRED_DISCOVERED_FILES - set(_SCANNED)
    assert not missing, f"the scan never read these files, so it cannot police them: {sorted(missing)}"
    assert len(_SCANNED) >= MIN_DISCOVERED_FILES, f"only {len(_SCANNED)} files scanned; the walk is pruning too much"


@pytest.mark.unit
@pytest.mark.parametrize("pruned_dir", sorted(PRUNED_DIRS))
def test_pruned_directories_state_a_reason(pruned_dir: str) -> None:
    """An exclusion without a reason is how the scope narrowed silently the first time."""
    assert PRUNED_DIRS[pruned_dir].strip(), f"{pruned_dir} is excluded with no stated reason"


@pytest.mark.unit
def test_the_environment_can_add_a_digest(monkeypatch: pytest.MonkeyPatch) -> None:
    """CI must be able to police a term a public file should not enumerate even as a digest."""
    added = hashlib.sha256(b"someothertoken").hexdigest()
    monkeypatch.setenv(_ENV_DIGESTS, f" {added.upper()} , ")

    assert added in forbidden_digests()
    assert forbidden_digests() >= FORBIDDEN_DIGESTS, "the declared set must survive the extension"


#: A term this repository does not forbid, used only to exercise the machinery. The fixture below needs
#: a pre-image, and taking a real one would put it back in the file.
_SENTINEL = "canaryidentifier"
_SENTINEL_DIGESTS = frozenset({hashlib.sha256(_SENTINEL.encode()).hexdigest()})


@pytest.mark.unit
def test_scan_sees_the_file_kinds_that_defeated_the_narrow_glob(tmp_path: Path) -> None:
    """The three audited defeats, plus the exclusions, asserted against a planted tree.

    ``.config``, ``CHANGELOG.md`` and ``.ipynb`` each hid a reference from the ``*.py`` glob. Planted
    with the sentinel rather than a real term, so the fixture proves the walk and the matching without
    the file needing a pre-image of anything it forbids.
    """
    leak = f"calibrated on the {_SENTINEL.upper()} run"
    (tmp_path / "src" / "workflows").mkdir(parents=True)
    (tmp_path / "src" / "workflows" / "nextflow.config").write_text(f"// {leak}\n", encoding="utf-8")
    (tmp_path / "CHANGELOG.md").write_text(f"- {leak}\n", encoding="utf-8")
    (tmp_path / "notebook.ipynb").write_text(f'{{"cells": [{{"source": ["# {leak}"]}}]}}\n', encoding="utf-8")
    (tmp_path / "leak_note.txt").write_text(f"{leak}\n", encoding="utf-8")
    (tmp_path / "clean.py").write_text("# one internal run\n", encoding="utf-8")

    # Exclusions must hold even when the term is present: a pruned dir and a binary file.
    (tmp_path / "_build").mkdir()
    (tmp_path / "_build" / "index.html").write_text(leak, encoding="utf-8")
    (tmp_path / "index.bin").write_bytes(b"\x00\x01" + leak.encode())

    offenders, scanned = _scan(tmp_path, _SENTINEL_DIGESTS, tracked_only=False)

    assert sorted(offenders) == [
        "CHANGELOG.md:1",
        "leak_note.txt:1",
        "notebook.ipynb:1",
        "src/workflows/nextflow.config:1",
    ]
    assert "clean.py" in scanned
    assert "index.bin" not in scanned, "binary file was decoded as prose"
    assert not [p for p in scanned if p.startswith("_build/")], "pruned directory was walked"


@pytest.mark.unit
@pytest.mark.parametrize(
    "line",
    [
        "calibrated on the CANARYIDENTIFIER run",
        "see canaryidentifier_run for the numbers",
        "threshold from canaryidentifier-002",
        "# CanaryIdentifier",
        '{"source": ["taken from the canaryidentifier run"]}',
    ],
)
def test_a_term_is_found_in_prose_snake_case_hyphenation_and_json(line: str) -> None:
    """Token splitting is what makes a digest usable, so the shapes it must cut are pinned."""
    assert _line_offends(line, _SENTINEL_DIGESTS)


@pytest.mark.unit
def test_a_term_welded_into_a_longer_word_is_the_documented_limitation() -> None:
    """Recorded rather than discovered later: per-token matching cannot see inside a word.

    A digest cannot be matched against a substring, which is the price of not storing the pre-image.
    A leak is written as a word, so this is a limitation and not a hole -- but it is asserted so the
    boundary is a decision rather than an accident.
    """
    assert not _line_offends(f"xx{_SENTINEL}yy", _SENTINEL_DIGESTS)
    assert _line_offends(f"xx {_SENTINEL} yy", _SENTINEL_DIGESTS)
