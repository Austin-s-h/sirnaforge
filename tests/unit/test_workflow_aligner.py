"""Backend selection reaches workflow configuration and reproducible search settings."""

from pathlib import Path

import pytest

from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("aligner", ["bwa-mem2", "bwa-mem3"])
def test_backend_selection_is_explicit_and_preserved(tmp_path: Path, aligner: str) -> None:
    """Backend selection is explicit and preserved."""
    config = WorkflowConfig(output_dir=tmp_path / "run", gene_query="test", aligner=aligner)
    assert config.aligner == aligner
    settings = SiRNAWorkflow._plan_search_settings({"aligner": config.aligner, "bwa_k": 12})
    assert settings["aligner"] == aligner
    assert settings["bwa_k"] == 12


def test_bad_backend_fails_before_output_creation(tmp_path: Path) -> None:
    """Bad backend fails before output creation."""
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="Unsupported aligner"):
        WorkflowConfig(output_dir=output, gene_query="test", aligner="bwa")
    assert not output.exists()


def test_conflicting_backend_sources_fail_before_output_creation(tmp_path: Path) -> None:
    """Conflicting backend sources fail before output creation."""
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="Conflicting aligner"):
        WorkflowConfig(
            output_dir=output,
            gene_query="test",
            aligner="bwa-mem3",
            nextflow_config={"aligner": "bwa-mem2"},
        )
    assert not output.exists()
