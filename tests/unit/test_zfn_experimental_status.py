"""Regression tests for the EXPERIMENTAL status of the ZFN arm.

The ZFN module ships with known unfixed defects, tracked in
https://github.com/Austin-s-h/sirnaforge/issues/82, so every entry point and every ZFN
documentation surface has to say so before a user can act on ZFN output.

The doc-page set is **discovered**, not listed. An earlier version of this file enumerated
exactly the files that had been edited, which structurally cannot fail -- and in fact missed
`docs/ccr5_zfn_benchmark.md`, a user-guide ZFN page in the toctree. Discovery means a new ZFN
page fails this file until it is marked.

The runtime notice is asserted two ways, and both are needed. The `_logs_the_warning_exactly_once`
tests count *log records*, which pins the once-per-process latch across layered entry points. The
`_visible_notice_count` tests count *rendered notices on stdout/stderr*, which is what a user
actually reads: a single log record was being printed twice (raw log line plus rich panel) while
every record-count test stayed green.
"""

from __future__ import annotations

import io
import logging
import os
import re
import subprocess
import sys
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

# A file with more than one case-insensitive "zfn" mention is a ZFN surface, not a page that
# happens to name the feature in passing. This was 4, which excluded exactly one page that
# would otherwise have failed (`docs/index.rst`, at 3) -- a threshold tuned to hide the one
# failing case. That page is now an explicit, reasoned exemption below, so the threshold can
# sit at the lowest defensible value instead of carrying a hidden exemption.
ZFN_MENTION_THRESHOLD = 2

# Surfaces that mention ZFN often but must not carry a usage banner. Each needs a reason:
# a blanket exemption list is how the transcribed page list went wrong in the first place.
BANNER_EXEMPT_PAGES: dict[str, str] = {
    "CHANGELOG.md": "release history; describes what shipped, and pinning wording here would fight every release",
    "release_notes.md": "release history, same reason as CHANGELOG.md",
    "docs/changelog.md": "rendered release history",
    "docs/developer/index.rst": "toctree plus one-line link descriptions; every page it links carries the banner",
    "docs/index.rst": (
        "documentation landing page: its three ZFN mentions are bare toctree entries "
        "(zfn_module, zfn_ranking, ccr5_zfn_benchmark) with no prose and no usage, and all "
        "three linked pages carry the banner. Same rationale as docs/developer/index.rst"
    ),
}

# Globs covering every prose/notebook surface a user could learn ZFN usage from.
#
# Scope choice: prose and notebooks only, wherever they live -- Markdown, reStructuredText and
# Jupyter. Source trees are deliberately out of scope, because `src/` states the status at
# runtime instead (`emit_zfn_experimental_warning`, asserted by the entry-point tests below)
# and `tests/` prose is covered by `tests/README.md`, which these globs do reach. The
# directories below beyond `docs/`/`notebooks/` currently contain no ZFN mentions at all; they
# are listed so that a ZFN usage page added to any of them fails this file rather than
# escaping a glob set that only knew about `docs/`.
DOC_SURFACE_GLOBS = (
    "*.md",
    "docs/**/*.md",
    "docs/**/*.rst",
    "notebooks/**/*.md",
    "notebooks/**/*.ipynb",
    "examples/**/*.md",
    "examples/**/*.ipynb",
    "docker/**/*.md",
    "tests/**/*.md",
    "scripts/**/*.md",
    ".github/**/*.md",
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


# One sentence that occurs exactly once inside ZFN_EXPERIMENTAL_WARNING, so counting it counts
# whole renderings of the notice. Counting "EXPERIMENTAL" instead would be ambiguous: the rich
# panel's own title line adds a second occurrence of that word per render.
NOTICE_SENTINEL = "Do not use ZFN results for decisions without independent validation."


def _visible_notice_count(text: str) -> int:
    """How many times a user reading ``text`` is shown the experimental notice.

    Counts rendered notices, not log records. A record-count assertion cannot see the failure
    this guards: with a rich console the notice used to be *printed twice* (log line plus
    panel) from a single log record.
    """
    assert NOTICE_SENTINEL in ZFN_EXPERIMENTAL_WARNING, "sentinel drifted out of the warning text"
    return _unwrapped(text).count(_unwrapped(NOTICE_SENTINEL).strip())


def _run_in_child(args: list[str]) -> subprocess.CompletedProcess[str]:
    """Run ``args`` under this interpreter so stdout and stderr are the real, separate streams.

    In-process capture cannot answer "what does the user see": ``CliRunner`` merges the two
    streams by default, and pytest's capture replaces them after the CLI has already bound a
    log handler to one of them.
    """
    env = {k: v for k, v in os.environ.items() if k != "SIRNAFORGE_LOG_FILE"}
    env["COLUMNS"] = "100"  # deterministic rich wrapping
    return subprocess.run(  # noqa: S603 - fixed argv, same interpreter as the test session
        [sys.executable, *args],
        capture_output=True,
        text=True,
        env=env,
        timeout=300,
        check=False,
    )


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


@pytest.mark.unit
def test_console_notice_is_rendered_once_and_still_logged(caplog: pytest.LogCaptureFixture) -> None:
    """One console render *and* one log record, when both share a stream.

    This is the production topology: `logging_utils.configure_logging` puts a StreamHandler on
    the same stream rich writes to, so emitting the log record and the panel used to print the
    whole notice twice in a row. Muting the record for that stream must not delete the record
    itself -- log files and library callers still need it.
    """
    from rich.console import Console  # noqa: PLC0415 - local to keep the module import light

    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    root = logging.getLogger()
    root.addHandler(handler)
    try:
        with caplog.at_level(logging.WARNING):
            assert emit_zfn_experimental_warning(Console(file=stream, width=100)) is True
    finally:
        root.removeHandler(handler)
        handler.close()

    visible = _visible_notice_count(stream.getvalue())
    assert visible == 1, f"user saw the notice {visible} times on one stream, expected 1"
    assert len(_warning_records(caplog)) == 1, "the log record was dropped, not just kept off the console"


@pytest.mark.unit
def test_zfn_cli_notice_is_visible_exactly_once_to_a_user(tmp_path: Path) -> None:
    """`sirnaforge zfn`: exactly one notice across the real stdout and stderr of a real run.

    Counts what lands on the terminal, in a child process, because the once-per-run contract is
    about what a user reads. The sibling record-count tests above pass even when a single log
    record is rendered twice, which is exactly what shipped.
    """
    left = "GCGTACGTA"
    right = "TACGGCATA"
    fasta = _tiny_zfn_genome(tmp_path, left, right)
    output_dir = tmp_path / "cli_out"

    result = _run_in_child(
        [
            "-m",
            "sirnaforge.cli",
            "zfn",
            "--zfn-left-half-site",
            left,
            "--zfn-right-half-site",
            right,
            "--zfn-search-space",
            str(fasta),
            "--output-dir",
            str(output_dir),
        ]
    )
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"

    on_stdout = _visible_notice_count(result.stdout)
    on_stderr = _visible_notice_count(result.stderr)
    assert on_stdout + on_stderr == 1, (
        f"user saw the notice {on_stdout + on_stderr} times "
        f"(stdout={on_stdout}, stderr={on_stderr}); expected exactly once"
    )
    # Pinned to stdout deliberately: that is where the rich panel and sirnaforge's log handler
    # both write today. Moving the notice to stderr is a legitimate change -- it must update
    # this assertion and the docs, not happen silently.
    assert (on_stdout, on_stderr) == (1, 0)

    # The run's log file keeps the record even though it was muted on the console stream.
    log_file = output_dir / "logs" / "sirnaforge.log"
    assert log_file.exists(), f"no run log written to {log_file}"
    assert _visible_notice_count(log_file.read_text(encoding="utf-8")) == 1


@pytest.mark.unit
def test_library_notice_lands_on_stdout_not_stderr() -> None:
    """The no-console path: one notice, on stdout.

    An earlier revision of this module claimed the notice "reaches stderr via logging's
    last-resort handler". It does not, and cannot: `get_logger` installs a
    `StreamHandler(sys.stdout)` on the root logger on first use, so `logging.lastResort` never
    fires. This test pins the real stream so the claim cannot be reinstated untested.
    """
    result = _run_in_child(
        [
            "-c",
            "from sirnaforge.zfn import emit_zfn_experimental_warning as e; assert e() is True",
        ]
    )
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"

    assert _visible_notice_count(result.stdout) == 1, f"expected one notice on stdout, got:\n{result.stdout}"
    assert _visible_notice_count(result.stderr) == 0, f"unexpected notice on stderr:\n{result.stderr}"
