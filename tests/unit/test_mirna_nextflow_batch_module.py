"""Unit checks for Nextflow miRNA batch module wiring."""

from pathlib import Path


def test_nextflow_mirna_batch_module_uses_python_entrypoint_defaults() -> None:
    """Batch module should invoke the shared Python entrypoint without backend-specific knobs."""
    module_path = (
        Path(__file__).parents[2]
        / "src"
        / "sirnaforge"
        / "pipeline"
        / "nextflow"
        / "workflows"
        / "modules"
        / "local"
        / "mirna_seed_analysis.nf"
    )
    module_text = module_path.read_text()

    assert "from sirnaforge.core.off_target import run_mirna_seed_analysis" in module_text
    assert "output_path = run_mirna_seed_analysis(" in module_text
    assert "backend=" not in module_text
