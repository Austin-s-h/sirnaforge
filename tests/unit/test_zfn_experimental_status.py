"""Regression tests for the EXPERIMENTAL status of the ZFN arm.

The ZFN module ships with known unfixed defects, tracked in
https://github.com/Austin-s-h/sirnaforge/issues/82, so every entry point and every ZFN
documentation surface has to say so before a user can act on ZFN output.

The doc-page set is **discovered**, not listed. An earlier version of this file enumerated
exactly the files that had been edited, which structurally cannot fail -- and in fact missed
`docs/ccr5_zfn_benchmark.md`, a user-guide ZFN page in the toctree. Discovery means a new ZFN
page fails this file until it is marked.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from typer.testing import CliRunner

import sirnaforge.cli as cli_module
from sirnaforge.models.sirna import DesignMode, DesignParameters
from sirnaforge.models.zfn import ZFNDesignParameters
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig, ZFNWorkflowConfig
from sirnaforge.zfn import (
    ZFN_CCR5_WORKING_RIGHT_HALF_SITE,
    ZFN_EXPERIMENTAL_ISSUE_URL,
    ZFN_EXPERIMENTAL_WARNING,
    ZFNDesigner,
    emit_zfn_experimental_warning,
)
from sirnaforge.zfn.experimental import reset_zfn_experimental_warning

# A file with at least this many case-insensitive "zfn" mentions is a ZFN surface, not a page
# that happens to name the feature once. `docs/index.rst` sits at 3 (three toctree entries),
# which is the level this is calibrated to exclude.
ZFN_MENTION_THRESHOLD = 4

# Surfaces that mention ZFN often but must not carry a usage banner. Each needs a reason:
# a blanket exemption list is how the transcribed page list went wrong in the first place.
BANNER_EXEMPT_PAGES: dict[str, str] = {
    "CHANGELOG.md": "release history; describes what shipped, and pinning wording here would fight every release",
    "release_notes.md": "release history, same reason as CHANGELOG.md",
    "docs/changelog.md": "rendered release history",
    "docs/developer/index.rst": "toctree plus one-line link descriptions; every page it links carries the banner",
}

# Globs covering every prose/notebook surface a user could learn ZFN usage from.
DOC_SURFACE_GLOBS = (
    "*.md",
    "docs/**/*.md",
    "docs/**/*.rst",
    "notebooks/**/*.md",
    "notebooks/**/*.ipynb",
)

# Anchors proving discovery actually walked the tree. Without these a broken glob would make
# the parametrized banner test collect nothing and pass silently.
REQUIRED_DISCOVERED_PAGES = frozenset(
    {
        "README.md",
        "docs/zfn_module.md",
        "docs/ccr5_zfn_benchmark.md",
        "docs/cli_reference.md",
        "docs/developer/testing_guide.md",
        "notebooks/README.md",
        "notebooks/zfn_backend_runtime_comparison.ipynb",
    }
)
MIN_DISCOVERED_PAGES = 10

HALF_SITE_ARGS = [
    "--zfn-left-half-site",
    "GCGTACGTA",
    "--zfn-right-half-site",
    "TACGGCATA",
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _discover_zfn_doc_pages() -> list[str]:
    """Every prose/notebook surface with a real ZFN presence, minus the reasoned exemptions."""
    root = _repo_root()
    discovered: set[str] = set()
    for glob in DOC_SURFACE_GLOBS:
        for path in root.glob(glob):
            rel = path.relative_to(root).as_posix()
            if rel in BANNER_EXEMPT_PAGES:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            if len(re.findall("zfn", text, re.IGNORECASE)) >= ZFN_MENTION_THRESHOLD:
                discovered.add(rel)
    return sorted(discovered)


ZFN_DOC_PAGES = _discover_zfn_doc_pages()


def _unwrapped(text: str) -> str:
    """Drop wrap artefacts (rich panel borders, markdown blockquote markers) and collapse whitespace.

    Both renderings break the warning across lines, so a raw substring search would miss it.
    """
    return re.sub(r"\s+", " ", re.sub(r"[>│|]", " ", text))


def _run_cli(args: list[str], capsys: pytest.CaptureFixture[str]) -> str:
    """Invoke the CLI and return its unwrapped output.

    Capture is disabled during the call: the CLI reconfigures root logging onto whichever
    stdout is live, which does not survive nested pytest/click capture.
    """
    with capsys.disabled():
        result = CliRunner().invoke(cli_module.app, args)
    assert result.exit_code == 0, result.output
    return _unwrapped(result.output)


def _tiny_zfn_genome(tmp_path: Path, left: str, right: str) -> Path:
    """A FASTA containing exactly one matching site for the pair, so a real run is fast."""
    target = f"{left}{'A' * 5}{Seq(right).reverse_complement()}"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\nAAA{target}CCC\n", encoding="utf-8")
    return fasta


def _warning_records(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    return [record for record in caplog.records if ZFN_EXPERIMENTAL_WARNING in record.getMessage()]


@pytest.fixture(autouse=True)
def _restore_root_log_handlers() -> Iterator[None]:
    """Undo the root-logger handlers the CLI installs, after every test in this file.

    `configure_logging` attaches a console StreamHandler plus a RotatingFileHandler for the
    run's log file and never removes them, which is correct for a process that exits after
    one run. Under pytest it leaks: this is the only file that invokes the real CLI, and
    without this fixture two RotatingFileHandlers pointing into `tmp_path` survive on the
    root logger for the rest of the xdist worker's session. Later tests in that worker then
    write into directories pytest has already deleted.
    """
    root = logging.getLogger()
    preserved = list(root.handlers)
    preserved_level = root.level
    try:
        yield
    finally:
        for handler in list(root.handlers):
            if handler not in preserved:
                root.removeHandler(handler)
                handler.close()
        root.setLevel(preserved_level)


@pytest.fixture(autouse=True)
def _rearm_experimental_notice() -> Iterator[None]:
    """Re-arm the once-per-process notice latch so each test sees a fresh run."""
    reset_zfn_experimental_warning()
    yield
    reset_zfn_experimental_warning()


@pytest.mark.unit
def test_experimental_warning_text_states_scope_and_required_action() -> None:
    """The shared warning must name the status, the tracking issue and the required action."""
    assert "EXPERIMENTAL" in ZFN_EXPERIMENTAL_WARNING
    assert "independent validation" in ZFN_EXPERIMENTAL_WARNING
    assert ZFN_EXPERIMENTAL_ISSUE_URL in ZFN_EXPERIMENTAL_WARNING
    assert ZFN_EXPERIMENTAL_ISSUE_URL.endswith("/82")


@pytest.mark.unit
def test_experimental_warning_names_the_defects_that_bite_silently() -> None:
    """The two defects with no visible symptom must be named, not just alluded to.

    A user can notice a crash or an empty result set. They cannot notice that
    `best_offtarget_score` means the opposite of its name, so the warning has to say so.
    """
    assert "worst_site_score" in ZFN_EXPERIMENTAL_WARNING
    assert "best_offtarget_score" in ZFN_EXPERIMENTAL_WARNING
    assert "inverted" in ZFN_EXPERIMENTAL_WARNING

    # The hard failure on the *default* backend at the documented CCR5 mismatch budget.
    assert "pyahocorasick" in ZFN_EXPERIMENTAL_WARNING
    assert "ValueError" in ZFN_EXPERIMENTAL_WARNING
    assert "exhaustive_python" in ZFN_EXPERIMENTAL_WARNING


@pytest.mark.unit
def test_experimental_warning_gives_a_working_invocation() -> None:
    """Telling a user their example is broken without telling them what to type is not enough."""
    assert str(Seq("AAACTGCAAAAG").reverse_complement()) == ZFN_CCR5_WORKING_RIGHT_HALF_SITE
    assert f"--zfn-right-half-site {ZFN_CCR5_WORKING_RIGHT_HALF_SITE}" in ZFN_EXPERIMENTAL_WARNING


@pytest.mark.unit
def test_zfn_doc_discovery_walked_the_tree() -> None:
    """Guard the discovery itself: a broken glob would make the banner test vacuous."""
    discovered = set(ZFN_DOC_PAGES)
    missing = REQUIRED_DISCOVERED_PAGES - discovered
    assert not missing, f"ZFN doc discovery missed known ZFN pages: {sorted(missing)}"
    assert len(discovered) >= MIN_DISCOVERED_PAGES, f"only {len(discovered)} ZFN pages discovered: {sorted(discovered)}"


@pytest.mark.unit
@pytest.mark.parametrize("exempt_page", sorted(BANNER_EXEMPT_PAGES))
def test_banner_exemptions_still_exist(exempt_page: str) -> None:
    """An exemption for a deleted file silently widens the contract; make it fail instead."""
    assert (_repo_root() / exempt_page).exists(), f"exempt page {exempt_page} no longer exists; drop the exemption"


@pytest.mark.unit
@pytest.mark.parametrize("doc_page", ZFN_DOC_PAGES)
def test_zfn_docs_carry_experimental_banner(doc_page: str) -> None:
    """Each ZFN surface must warn about experimental status and point at the tracking issue."""
    text = _unwrapped((_repo_root() / doc_page).read_text(encoding="utf-8"))

    assert "EXPERIMENTAL" in text, f"{doc_page} has no EXPERIMENTAL banner"
    assert "independent validation" in text, f"{doc_page} does not require independent validation"
    assert ZFN_EXPERIMENTAL_ISSUE_URL in text or "#82" in text, f"{doc_page} does not reference issue #82"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_zfn_workflow_warns_experimental(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Entering the in-process ZFN workflow must emit the experimental warning."""
    left = "GCGTACGTA"
    right = "TACGGCATA"
    fasta = _tiny_zfn_genome(tmp_path, left, right)

    cfg = WorkflowConfig(
        gene_query="ZFN_EXPERIMENTAL",
        output_dir=tmp_path / "out",
        design_params=DesignParameters(design_mode=DesignMode.ZFN),
        zfn_config=ZFNWorkflowConfig(
            zfn_params=ZFNDesignParameters(
                left_half_site=left,
                right_half_site=right,
                search_space_fasta=str(fasta),
            )
        ),
    )

    with caplog.at_level(logging.WARNING):
        await SiRNAWorkflow(cfg).run_complete_workflow()

    assert _warning_records(caplog)


@pytest.mark.unit
def test_python_api_evaluate_pair_warns_experimental(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """`ZFNDesigner().evaluate_pair()` is documented Python API and must carry the notice itself.

    A caller who never touches the CLI or SiRNAWorkflow previously got ZFN artifacts with no
    console output and no log record saying the module is experimental.
    """
    left = "GCGTACGTA"
    right = "TACGGCATA"
    params = ZFNDesignParameters(
        left_half_site=left,
        right_half_site=right,
        search_space_fasta=str(_tiny_zfn_genome(tmp_path, left, right)),
    )

    with caplog.at_level(logging.WARNING):
        ZFNDesigner().evaluate_pair(params=params)

    assert _warning_records(caplog), "evaluate_pair() produced no experimental-status log record"


@pytest.mark.unit
def test_zfn_command_warns_experimental(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """`sirnaforge zfn` must print the experimental warning before doing any work."""

    async def _fake_workflow(**_kwargs: Any) -> dict[str, Any]:
        return {}

    monkeypatch.setattr(cli_module, "run_sirna_workflow", _fake_workflow)

    output = _run_cli(
        ["zfn", *HALF_SITE_ARGS, "--output-dir", str(tmp_path / "zfn_out")],
        capsys,
    )

    assert "EXPERIMENTAL" in output
    assert "independent validation" in output
    assert ZFN_EXPERIMENTAL_ISSUE_URL in output.replace(" ", "")


@pytest.mark.unit
def test_workflow_command_design_mode_zfn_warns_experimental(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """`workflow --design-mode zfn` must print the experimental warning too."""

    async def _fake_workflow(**_kwargs: Any) -> dict[str, Any]:
        return {}

    monkeypatch.setattr(cli_module, "run_sirna_workflow", _fake_workflow)

    output = _run_cli(
        [
            "workflow",
            "CCR5",
            "--design-mode",
            "zfn",
            *HALF_SITE_ARGS,
            "--output-dir",
            str(tmp_path / "wf_out"),
        ],
        capsys,
    )

    assert "EXPERIMENTAL" in output
    assert "independent validation" in output


@pytest.mark.unit
def test_emit_is_latched_so_layered_entry_points_do_not_repeat() -> None:
    """The notice is once per process: CLI -> workflow -> designer must not stack three copies."""
    assert emit_zfn_experimental_warning() is True
    assert emit_zfn_experimental_warning() is False
    reset_zfn_experimental_warning()
    assert emit_zfn_experimental_warning() is True


@pytest.mark.unit
@pytest.mark.asyncio
async def test_full_zfn_run_logs_the_warning_exactly_once(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """One run, one notice.

    This path stacks two emitters -- `_run_zfn_workflow` then `ZFNDesigner.evaluate_pair` --
    and previously the CLI added a third. Log records are counted rather than console text
    because the rich panel wraps the message across lines.
    """
    left = "GCGTACGTA"
    right = "TACGGCATA"
    cfg = WorkflowConfig(
        gene_query="ZFN_ONCE",
        output_dir=tmp_path / "out",
        design_params=DesignParameters(design_mode=DesignMode.ZFN),
        zfn_config=ZFNWorkflowConfig(
            zfn_params=ZFNDesignParameters(
                left_half_site=left,
                right_half_site=right,
                search_space_fasta=str(_tiny_zfn_genome(tmp_path, left, right)),
            )
        ),
    )

    with caplog.at_level(logging.WARNING):
        await SiRNAWorkflow(cfg).run_complete_workflow()

    records = _warning_records(caplog)
    assert len(records) == 1, f"expected one experimental notice per run, got {len(records)}"


@pytest.mark.unit
def test_zfn_cli_run_logs_the_warning_exactly_once(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Same contract through the real CLI, which adds a third emitter above the other two."""
    left = "GCGTACGTA"
    right = "TACGGCATA"
    fasta = _tiny_zfn_genome(tmp_path, left, right)

    with caplog.at_level(logging.WARNING), capsys.disabled():
        result = CliRunner().invoke(
            cli_module.app,
            [
                "zfn",
                "--zfn-left-half-site",
                left,
                "--zfn-right-half-site",
                right,
                "--zfn-search-space",
                str(fasta),
                "--output-dir",
                str(tmp_path / "cli_out"),
            ],
        )
    assert result.exit_code == 0, result.output

    records = _warning_records(caplog)
    assert len(records) == 1, f"expected one experimental notice per CLI run, got {len(records)}"
