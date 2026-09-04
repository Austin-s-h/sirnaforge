"""Hold ``docs/getting_started.md``'s output inventory to what a real run writes.

Moving the ``check_off_targets`` guard to the top of ``step5_offtarget_analysis`` (so a skipped
run resolves no reference and stages no FASTA) made one documented promise false: the
"Key files in every workflow run" table listed ``off_target/input_candidates.fasta``, which a
``--skip-off-targets`` run never writes. This test runs the skip path for real -- no network, no
Nextflow, no Docker -- and asserts that

* every path in the unconditional table exists, and
* every path in the "only when off-target screening runs" table does not.

So the table cannot silently drift from the code again in either direction.
"""

import asyncio
import re
from pathlib import Path

import pytest

from sirnaforge.workflow import run_sirna_workflow

REPO_ROOT = Path(__file__).resolve().parents[2]
GETTING_STARTED = REPO_ROOT / "docs" / "getting_started.md"

EVERY_RUN_MARKER = "Key files in every workflow run:"
SCREENING_ONLY_MARKER = "Written only when off-target screening actually runs"


def _table_paths(text: str, marker: str) -> list[str]:
    """Return the backticked paths in the first markdown table after ``marker``.

    Tolerant of cell padding (a formatter may align the pipes) and skips the header and
    separator rows, which have no backticked first cell.
    """
    assert marker in text, f"marker not found in {GETTING_STARTED}: {marker!r}"
    paths: list[str] = []
    for line in text.split(marker, 1)[1].splitlines():
        stripped = line.strip()
        if not stripped.startswith("|"):
            if paths:
                break  # table finished
            continue
        first_cell = stripped.strip("|").split("|")[0].strip()
        match = re.fullmatch(r"`(.+)`", first_cell)
        if match:
            paths.append(match.group(1))
    assert paths, f"no table rows found after {marker!r} in {GETTING_STARTED}"
    return paths


@pytest.fixture(scope="module")
def skip_offtargets_output_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Run ``--skip-off-targets`` equivalent end to end and return its output directory."""
    work = tmp_path_factory.mktemp("skip_offtargets")
    fasta = work / "toy.fa"
    fasta.write_text(">toy1 toy transcript one\nATG" + "ACGTTGCAAGGCTTACGGATCCA" * 8 + "TAA\n")
    output_dir = work / "out"

    asyncio.run(
        run_sirna_workflow(
            gene_query="TOY",
            output_dir=str(output_dir),
            input_fasta=str(fasta),
            check_off_targets=False,
            top_n_candidates=5,
        )
    )
    return output_dir


@pytest.mark.integration
def test_documented_every_run_files_exist_without_off_target_screening(
    skip_offtargets_output_dir: Path,
) -> None:
    """Everything the docs promise in *every* run must survive --skip-off-targets."""
    documented = _table_paths(GETTING_STARTED.read_text(), EVERY_RUN_MARKER)

    missing = [path for path in documented if not list(skip_offtargets_output_dir.glob(path))]
    assert not missing, (
        f"docs/getting_started.md promises these in every workflow run, but a "
        f"--skip-off-targets run wrote none of them: {missing}"
    )


@pytest.mark.integration
def test_screening_only_files_absent_without_off_target_screening(
    skip_offtargets_output_dir: Path,
) -> None:
    """The screening-only table must list exactly what the skip path does not produce."""
    doc_text = GETTING_STARTED.read_text()
    screening_only = _table_paths(doc_text, SCREENING_ONLY_MARKER)

    assert "off_target/input_candidates.fasta" in screening_only, (
        "off_target/input_candidates.fasta is only staged when screening runs; it must be "
        "documented as conditional, not as an every-run file"
    )

    unexpected = [path for path in screening_only if list(skip_offtargets_output_dir.glob(path.rstrip("/")))]
    assert not unexpected, f"documented as screening-only, but a --skip-off-targets run wrote: {unexpected}"
