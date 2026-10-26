"""Execute embedded shell fragments so version artifacts remain usable provenance."""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

MODULES = Path(__file__).resolve().parents[2] / "src/sirnaforge/pipeline/nextflow/workflows/modules/local"
BLOCKS = [
    (path.name, index, block)
    for path in sorted(MODULES.glob("*.nf"))
    for index, block in enumerate(
        re.findall(r"cat <<-END_VERSIONS > versions.yml\n.*?\n[ \t]*END_VERSIONS", path.read_text(), re.S)
    )
]


@pytest.fixture
def tool_environment(tmp_path: Path) -> dict[str, str]:
    """Mimic BWA-MEM2's noisy launcher without requiring host bioinformatics tools."""
    tools = tmp_path / "bin"
    tools.mkdir()
    for name, body in {
        "python": 'if [ "$1" = --version ]; then echo "Python 3.12.0"; else echo "1.0.0"; fi',
        "bwa-mem3": 'echo "3.0.0"',
        "bwa-mem2": 'echo "Launching executable bwa-mem2.avx2" >&2; echo "2.2.1"',
    }.items():
        executable = tools / name
        executable.write_text("#!/bin/sh\n" + body + "\n")
        executable.chmod(0o755)
    return {**os.environ, "PATH": str(tools) + os.pathsep + os.environ["PATH"]}


def render(text: str, aligner: str = "bwa-mem2") -> str:
    """Apply the substitutions relevant to these fragments in Groovy GStrings."""
    return (
        text.replace("${task.process}", "TEST_PROCESS")
        .replace("${species}", "human")
        .replace("${params.aligner ?: 'bwa-mem2'}", aligner)
        .replace("${params.aligner}", aligner)
        .replace(r"\$", "$")
    )


@pytest.mark.parametrize("aligner", ["bwa-mem2", "bwa-mem3"])
@pytest.mark.parametrize(("module", "index", "block"), BLOCKS, ids=[f"{m}-{i}" for m, i, _ in BLOCKS])
def test_version_fragment_has_no_terminator_or_launcher_diagnostic(
    tmp_path: Path, tool_environment: dict[str, str], module: str, index: int, block: str, aligner: str
) -> None:
    """The real shell must terminate its heredoc and capture the numeric BWA version."""
    completed = subprocess.run(
        ["bash", "--noprofile", "--norc", "-euo", "pipefail", "-c", render(block, aligner)],
        cwd=tmp_path,
        env=tool_environment,
        capture_output=True,
        text=True,
        check=True,
    )
    assert not completed.stderr, (module, index, completed.stderr)
    versions = (tmp_path / "versions.yml").read_text()
    assert versions.startswith('"TEST_PROCESS":\n')
    assert "END_VERSIONS" not in versions
    assert "python: 3.12.0" in versions
    if "params.aligner" in block:
        version = "2.2.1" if aligner == "bwa-mem2" else "3.0.0"
        assert f"{aligner}: {version}\n" in versions
        assert "Launching" not in versions


@pytest.mark.parametrize("module", ["aggregate_results.nf", "mirna_seed_analysis.nf", "offtarget_analysis.nf"])
def test_stub_evidence_heredoc_still_terminates(tmp_path: Path, tool_environment: dict[str, str], module: str) -> None:
    """Changing string indentation must not swallow the following version command."""
    match = re.search(r'stub:\s*"""(.*?)"""', (MODULES / module).read_text(), re.S)
    assert match is not None
    completed = subprocess.run(
        ["bash", "--noprofile", "--norc", "-euo", "pipefail", "-c", render(match.group(1))],
        cwd=tmp_path,
        env=tool_environment,
        capture_output=True,
        text=True,
        check=True,
    )
    assert not completed.stderr
    evidence_files = list(tmp_path.glob("*evidence.json"))
    assert evidence_files
    for path in evidence_files:
        if path.name == "mirna_reference_evidence.json":
            assert json.loads(path.read_text())["mirna"]["screening"]["status"] == "failed"
            continue
        assert json.loads(path.read_text())["schema_version"] == "2"
    assert (tmp_path / "versions.yml").exists()
