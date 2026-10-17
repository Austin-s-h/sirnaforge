"""CLI surface for benchmark artifacts: `sirnaforge benchmark prepare`/`design` (#109, bm-cli slice).

Two verified facts frame every test below (see also `sirnaforge.benchmark.panels`'s module docstring
and `tests/unit/data/benchmark/README.md`): the issue's cited PRD
(`docs/prd_benchmark_artifacts_and_variable_length.md`) does not exist in this repository, so the
acceptance criteria are the entire specification; and of the five panels #109/#110 name, only a
180-row Huesken redistribution (`tests/unit/data/sirna_efficacy_subset.csv`) ships real bytes. A test
passing over that one panel is not evidence about Ichihara, Martinelli, Shmushkovich or OligoGym --
it exercises the artifact/CLI *contract*, which is what this file owns.

`sirnaforge.benchmark.prepare`/`.design` are a sibling #109 slice's package; this file (and
`cli.py`) owns only the CLI wiring that calls them. Where that sibling has not yet landed in this
worktree, the functional tests below skip with a named reason rather than reporting a false red;
`test_benchmark_commands_are_exactly_prepare_and_design` depends on nothing but `cli.py`'s own
registration and always runs.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import typer
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.config.run_policy import default_for

try:
    from sirnaforge.benchmark import design as _benchmark_design_probe  # noqa: F401
    from sirnaforge.benchmark import prepare as _benchmark_prepare_probe  # noqa: F401

    _BENCHMARK_WORKFLOW_AVAILABLE = True
except ImportError:
    _BENCHMARK_WORKFLOW_AVAILABLE = False

requires_benchmark_workflow = pytest.mark.skipif(
    not _BENCHMARK_WORKFLOW_AVAILABLE,
    reason=(
        "blocked on a sibling #109 slice: sirnaforge.benchmark does not yet expose prepare/design in "
        "this worktree; bm-cli owns only the CLI wiring that calls them"
    ),
)

runner = CliRunner()


def _prepare_huesken_subset(out_dir: Path) -> Path:
    """Run the real `prepare` command and return the artifact directory it wrote."""
    result = runner.invoke(app, ["benchmark", "prepare", "--panel", "huesken_subset", "--out-dir", str(out_dir)])
    assert result.exit_code == 0, result.output
    artifact_dirs = [path for path in out_dir.iterdir() if path.is_dir()]
    assert len(artifact_dirs) == 1, artifact_dirs
    return artifact_dirs[0]


def _read_manifest(artifact_dir: Path) -> dict[str, Any]:
    payload: dict[str, Any] = json.loads((artifact_dir / "manifest.json").read_text())
    return payload


@pytest.mark.unit
def test_benchmark_commands_are_exactly_prepare_and_design() -> None:
    """The `benchmark` sub-app exposes exactly two verbs, whatever else this worktree gains (#109).

    A flat `benchmark-prepare`/`benchmark-design` pair would be the novelty here, not the
    established shape -- `sequences_app`/`internal_app` are already one noun with verbs beneath it.
    Depends only on `cli.py`'s own `app.add_typer`/`command_decorator_typed` wiring, not on the
    sibling `sirnaforge.benchmark` package, so it is the one assertion in this file that always runs.
    """
    command = typer.main.get_command(app)
    benchmark_group = command.commands["benchmark"]  # type: ignore[attr-defined]
    assert set(benchmark_group.commands) == {"prepare", "design"}


@requires_benchmark_workflow
@pytest.mark.unit
def test_prepare_writes_an_artifact_for_the_vendored_huesken_subset_panel(tmp_path: Path) -> None:
    """`huesken_subset` is the one panel with real bytes; preparing it must exit 0 and leave an artifact.

    No `--panel-csv`: the panel is vendored, so the CLI must not need one for it, and a build that
    silently required it for the one real panel would defeat the point of vendoring it at all.
    """
    result = runner.invoke(app, ["benchmark", "prepare", "--panel", "huesken_subset", "--out-dir", str(tmp_path)])
    assert result.exit_code == 0, result.output

    artifact_dir = tmp_path / "huesken_subset__len19"
    assert artifact_dir.is_dir()
    assert (artifact_dir / "manifest.json").is_file()
    assert (artifact_dir / "observations.csv").is_file()
    assert (artifact_dir / "design_inputs.fasta").is_file()

    manifest = _read_manifest(artifact_dir)
    assert manifest["panel"]["panel_id"] == "huesken_subset"
    assert manifest["panel"]["data_present"] is True
    assert manifest["paired_length"] == 19


@requires_benchmark_workflow
@pytest.mark.unit
def test_design_widens_gc_min_while_the_default_policy_keeps_the_shipped_floor(tmp_path: Path) -> None:
    """`--gc-min` widens only the benchmark policy; `default_run_policy` must report the untouched floor.

    This is the provenance guarantee #101 built and #109 reuses rather than re-implementing: the two
    policy blocks resolve independently, in the same process, so a leak into the default block would
    mean a value the caller typed only for the benchmark run silently became the shipped default too.

    Calls `design_artifact` directly rather than through `CliRunner`: running the real designer
    (RNA-folding included) under `CliRunner`'s captured stdout has a reproducible failure unrelated
    to #109 or to this assertion (something in that path closes the captured buffer, so
    `CliRunner.invoke` raises `ValueError: I/O operation on closed file` on exit even though the
    command completed and wrote every file correctly). `test_design_refuses_a_narrowing_...` below
    already proves the CLI threads `--gc-min` into `design_artifact` correctly -- a broken wiring
    would not trigger that refusal at the exact threshold it asserts.
    """
    from sirnaforge.benchmark.design import design_artifact  # noqa: PLC0415

    artifact_dir = _prepare_huesken_subset(tmp_path)
    design_artifact(artifact_dir=artifact_dir, gc_min=15.0, invoked_command=["sirnaforge", "benchmark", "design"])

    manifest = _read_manifest(artifact_dir)
    assert manifest["gc_widening"]["gc_min"]["benchmark"] == 15.0
    assert manifest["gc_widening"]["gc_min"]["widened"] is True

    default_gc_min = next(
        record["value"] for record in manifest["default_run_policy"]["resolved_settings"] if record["key"] == "gc_min"
    )
    assert default_gc_min == default_for("gc_min")


@requires_benchmark_workflow
@pytest.mark.unit
def test_design_omitting_gc_max_records_no_widening_and_no_requested_entry(tmp_path: Path) -> None:
    """An omitted `--gc-max` must read as unstated, not as the default value typed explicitly (#101).

    A value-equality implementation cannot tell the two apart and would put a `gc_max` entry in
    `requested_settings` the caller never typed; `_option_was_stated` (Click's parameter source) is
    what keeps that from happening here, the same mechanism `sirnaforge design` already relies on.

    See the note on the previous test for why this calls `design_artifact` directly.
    """
    from sirnaforge.benchmark.design import design_artifact  # noqa: PLC0415

    artifact_dir = _prepare_huesken_subset(tmp_path)
    design_artifact(artifact_dir=artifact_dir, gc_min=15.0, invoked_command=["sirnaforge", "benchmark", "design"])

    manifest = _read_manifest(artifact_dir)
    assert manifest["gc_widening"]["gc_max"]["widened"] is False
    requested_keys = {record["key"] for record in manifest["run_policy"]["requested_settings"]}
    assert "gc_max" not in requested_keys


@requires_benchmark_workflow
@pytest.mark.unit
def test_design_refuses_a_narrowing_gc_min_with_a_message_not_a_traceback(tmp_path: Path) -> None:
    """A `--gc-min` narrower than the shipped floor is refused before any design work (#109).

    Widening is safe because the benchmark run enumerates a superset of what a default run would; a
    narrower run would falsify the default-policy verdicts this command re-derives for candidates it
    never actually screened, so `resolve_run_policy` must refuse rather than silently comply.
    """
    artifact_dir = _prepare_huesken_subset(tmp_path)
    narrower_gc_min = default_for("gc_min") + 10.0
    result = runner.invoke(
        app, ["benchmark", "design", "--artifact", str(artifact_dir), "--gc-min", str(narrower_gc_min)]
    )

    assert result.exit_code == 1
    assert "Traceback" not in result.output
    assert "Error" in result.output
