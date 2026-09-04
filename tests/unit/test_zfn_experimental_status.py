"""Regression tests for the EXPERIMENTAL status of the ZFN arm.

The ZFN module ships with known unfixed defects (half-site orientation handling, FokI
seed-region weighting, off-target region classification), so every entry point and every
ZFN doc page has to say so before a user can act on ZFN output.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from typer.testing import CliRunner

import sirnaforge.cli as cli_module
from sirnaforge.models.sirna import DesignMode, DesignParameters
from sirnaforge.models.zfn import ZFNDesignParameters
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig, ZFNWorkflowConfig
from sirnaforge.zfn import ZFN_EXPERIMENTAL_WARNING

# Every ZFN doc page must carry the banner; the docs previously presented the recorded
# ZFN validation runs as authoritative.
ZFN_DOC_PAGES = [
    "docs/zfn_module.md",
    "docs/zfn_ranking.md",
    "docs/workflows.md",
    "docs/cli_reference.md",
    "docs/usage_examples.md",
    "docs/developer/zfn_backend_tuning.md",
    "docs/developer/zfn_hg38_primary_test_commands.md",
    "docs/developer/zfn_nextflow_bridge_validation.md",
]

HALF_SITE_ARGS = [
    "--zfn-left-half-site",
    "GCGTACGTA",
    "--zfn-right-half-site",
    "TACGGCATA",
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


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


@pytest.mark.unit
def test_experimental_warning_text_states_scope_and_required_action() -> None:
    """The shared warning must name the experimental status and the validation requirement."""
    assert "EXPERIMENTAL" in ZFN_EXPERIMENTAL_WARNING
    assert "independent validation" in ZFN_EXPERIMENTAL_WARNING


@pytest.mark.unit
@pytest.mark.asyncio
async def test_zfn_workflow_warns_experimental(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Entering the in-process ZFN workflow must emit the experimental warning."""
    left = "GCGTACGTA"
    right = "TACGGCATA"
    target = f"{left}{'A' * 5}{str(Seq(right).reverse_complement())}"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\nAAA{target}CCC\n", encoding="utf-8")

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

    with caplog.at_level(logging.WARNING, logger="sirnaforge.workflow"):
        await SiRNAWorkflow(cfg).run_complete_workflow()

    assert any(ZFN_EXPERIMENTAL_WARNING in record.getMessage() for record in caplog.records)


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
@pytest.mark.parametrize("doc_page", ZFN_DOC_PAGES)
def test_zfn_docs_carry_experimental_banner(doc_page: str) -> None:
    """Each ZFN doc page must warn about experimental status at its ZFN entry point."""
    text = _unwrapped((_repo_root() / doc_page).read_text(encoding="utf-8"))

    assert "EXPERIMENTAL" in text, f"{doc_page} has no EXPERIMENTAL banner"
    assert "independent validation" in text, f"{doc_page} does not require independent validation"
