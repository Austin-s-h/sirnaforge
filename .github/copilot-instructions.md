# GitHub Copilot Instructions for siRNAforge

## Big Picture (where to start)
- Primary entry points: CLI in `src/sirnaforge/cli.py` (Typer/Rich) → orchestration in `src/sirnaforge/workflow.py`.
- Typical data flow: gene/transcripts (`src/sirnaforge/data/gene_search.py`) → ORF analysis (`src/sirnaforge/data/orf_analysis.py`) → design (`src/sirnaforge/core/design.py`) → thermodynamics (`src/sirnaforge/core/thermodynamics.py`) → optional off-target / Nextflow (`src/sirnaforge/core/off_target.py`, `src/sirnaforge/pipeline/nextflow/`).
- Prefer typed models + validation: Pydantic models in `src/sirnaforge/models/` and Pandera schemas in `src/sirnaforge/models/schemas.py` + `src/sirnaforge/validation/`.

## Repo Conventions
- Strict `src/` layout; add new packages under `src/sirnaforge/`.
- Keep the CLI thin: new subcommands go in `src/sirnaforge/cli.py` but should delegate to orchestration/helpers (usually `workflow.py`).
- Reuse the cache/reference resolution helpers instead of ad-hoc paths: see `src/sirnaforge/config/reference_policy.py` and `src/sirnaforge/utils/unified_cache.py`. Workflow runs record resolved reference choices in `logs/workflow_summary.json`.

## Dev Workflow (use these commands)
- Setup once: `make dev` (runs `uv sync` + installs pre-commit).
- Run things via uv: `uv run sirnaforge …`.
- Lint/format/typecheck: `make lint`, `make format` (Ruff), `make check` (format + `make test-dev`).
- Quick manual sanity check: `uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test.csv`.

## Testing (how this repo behaves)
- Pytest defaults in `pyproject.toml` include xdist (`-n 2`). For debugging hangs or ordering issues, run `uv run pytest -n 0 …`.
- Tiered suites are Make targets: `make test-dev` (fast), `make test-ci`, `make test-release`, `make test`.
- Deterministic fixtures live in `tests/unit/data/` and `examples/`. Avoid network in unit tests unless explicitly marked (`requires_network`).
- Container-only integration lives in `tests/container/` and expects Docker + bundled bio tools.

## Pipeline / External Tooling Boundaries
- Nextflow assets are embedded under `src/sirnaforge/pipeline/nextflow/workflows/` and are executed via `src/sirnaforge/pipeline/nextflow/runner.py` (don’t reimplement pipeline orchestration).
- Heavy tools (BWA-MEM2, SAMtools, ViennaRNA) are pinned in `docker/Dockerfile`; on bare metal they’re optional and should be capability-gated.

## Packaging Constraints
- Packaging is Hatchling-only (`pyproject.toml`, `uv.lock`); don’t add `setup.py` or ad-hoc requirements files.
- Supported Python is 3.10–3.12 (see README); avoid relying on Python 3.13-only features.
