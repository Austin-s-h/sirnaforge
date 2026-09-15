"""This repository is public, so an internal programme target may not name itself in it.

Twenty-four comments, docstrings and fixture values justified a threshold or a design decision by
citing measurements from an internal drug-programme run and naming that programme's target. The
measurements are the evidence and stay; the identifier does not, because a public reader learns the
programme from it. The established phrasing is "one internal run" / "one internal reference set".

The token is assembled at runtime rather than written out, so a plain ``grep`` over the tree stays
clean and this guard is not its own only hit.

**Scope is the whole repository, not two trees of ``*.py``.** The first version of this file scanned
``rglob("*.py")`` under ``src`` and ``tests``, and an audit defeated it three separate ways, each
reproduced for real:

* non-``.py`` files *inside* the trees it claimed to scan -- a comment in the Nextflow config, a line
  in ``tests/data/benchmarks/README.md``, and a planted ``.txt`` next to it. A Nextflow module comment
  is a very plausible home for exactly this kind of run-calibration note.
* trees it never looked at -- ``CHANGELOG.md``, ``docs/``, ``README.md``, ``scripts/``, ``.github/``.
* ``.ipynb``, where a notebook cell is plain JSON text and reads like prose to everyone but a glob.

So discovery walks the tree and *excludes*, rather than listing what to include: a new file kind, a
new top-level directory or a new prose format is scanned by default and has to be argued out of scope
in ``PRUNED_DIRS`` below. Binary files are recognised by content, not by an extension list, for the
same reason. ``test_scan_sees_the_file_kinds_that_defeated_the_narrow_glob`` pins all three defeats.

Cost, measured on this tree: 352 files / 23.2 MB, 86-91 ms warm and ~5 s on a cold page cache. That
is fast-tier work, so nothing is excluded for size -- including the ``tests/unit/data/baseline_0_7_1``
``.tsv`` baselines (252 KB together) and the four multi-megabyte benchmark CSVs under ``tests/data``.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Split so the file does not itself contain the identifier it forbids.
PRIVATE_TARGET = "MSH" + "3"

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


def _is_binary(raw: bytes) -> bool:
    """A NUL byte in the first 8 KB. Decided by content so no novel extension can hide prose."""
    return b"\x00" in raw[:8192]


def _scan(root: Path, needle: str) -> tuple[list[str], list[str]]:
    """Return ``(offenders, scanned)`` for one tree: ``path:line`` hits, and every file read.

    Symlinks are skipped: they either point back inside the tree (already walked) or outside it.
    """
    offenders: list[str] = []
    scanned: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in PRUNED_DIRS and not d.endswith(".egg-info")]
        for name in sorted(filenames):
            path = Path(dirpath) / name
            if path.is_symlink() or not path.is_file():
                continue
            raw = path.read_bytes()
            if _is_binary(raw):
                continue
            rel = path.relative_to(root).as_posix()
            scanned.append(rel)
            for n, line in enumerate(raw.decode("utf-8", errors="replace").splitlines(), start=1):
                if needle in line.lower():
                    offenders.append(f"{rel}:{n}")
    return offenders, sorted(scanned)


_OFFENDERS, _SCANNED = _scan(REPO_ROOT, PRIVATE_TARGET.lower())


@pytest.mark.unit
def test_no_private_programme_target_anywhere_in_the_repository() -> None:
    """The programme target appears in no text file in the repository, in any case.

    Case-insensitive because the species convention writes the mouse symbol in title case, and both
    spellings were present: the prose carried the upper-case human form and a conservation fixture's
    symbol column carried the title-case rodent one. A public reference gene used deliberately as a
    public baseline (``TP53``) is a different thing and is untouched.
    """
    assert not _OFFENDERS, (
        f"{len(_OFFENDERS)} site(s) name the internal programme target. Keep the measurement, drop the "
        f"identifier -- say 'one internal run' instead: {_OFFENDERS}"
    )


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
def test_scan_sees_the_file_kinds_that_defeated_the_narrow_glob(tmp_path: Path) -> None:
    """The three audited defeats, plus the exclusions, asserted against a planted tree.

    ``.config``, ``CHANGELOG.md`` and ``.ipynb`` each hid a reference from the ``*.py`` glob. The
    token is assembled here too, so this file still does not spell it.
    """
    leak = f"calibrated on the {PRIVATE_TARGET} run"
    (tmp_path / "src" / "workflows").mkdir(parents=True)
    (tmp_path / "src" / "workflows" / "nextflow.config").write_text(f"// {leak}\n", encoding="utf-8")
    (tmp_path / "CHANGELOG.md").write_text(f"- {leak}\n", encoding="utf-8")
    (tmp_path / "notebook.ipynb").write_text(f'{{"cells": [{{"source": ["# {leak}"]}}]}}\n', encoding="utf-8")
    (tmp_path / "leak_note.txt").write_text(f"{leak}\n", encoding="utf-8")
    (tmp_path / "clean.py").write_text("# one internal run\n", encoding="utf-8")

    # Exclusions must hold even when the token is present: a pruned dir and a binary file.
    (tmp_path / "_build").mkdir()
    (tmp_path / "_build" / "index.html").write_text(leak, encoding="utf-8")
    (tmp_path / "index.bin").write_bytes(b"\x00\x01" + leak.encode())

    offenders, scanned = _scan(tmp_path, PRIVATE_TARGET.lower())

    assert sorted(offenders) == [
        "CHANGELOG.md:1",
        "leak_note.txt:1",
        "notebook.ipynb:1",
        "src/workflows/nextflow.config:1",
    ]
    assert "clean.py" in scanned
    assert "index.bin" not in scanned, "binary file was decoded as prose"
    assert not [p for p in scanned if p.startswith("_build/")], "pruned directory was walked"
