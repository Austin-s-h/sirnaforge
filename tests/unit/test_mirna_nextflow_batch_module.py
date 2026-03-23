"""Unit checks for Nextflow miRNA batch module wiring."""

from pathlib import Path


def _repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").exists():
            return parent
    raise AssertionError("Could not locate repository root containing pyproject.toml")


def test_nextflow_mirna_module_uses_default_backend() -> None:
    """Batch module should invoke the shared Python entrypoint without backend-specific knobs."""
    module_path = (
        _repo_root()
        / "src"
        / "sirnaforge"
        / "pipeline"
        / "nextflow"
        / "workflows"
        / "modules"
        / "local"
        / "mirna_seed_analysis.nf"
    )
    assert module_path.exists(), (
        f"Expected Nextflow miRNA module at {module_path}. "
        "Check path under src/sirnaforge/pipeline/nextflow/workflows/modules/local/."
    )
    module_text = module_path.read_text()

    assert "from sirnaforge.core.off_target import run_mirna_seed_analysis" in module_text
    assert "output_path = run_mirna_seed_analysis(" in module_text
    assert "backend=" not in module_text
