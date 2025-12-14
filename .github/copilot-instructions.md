# GitHub Copilot Instructions for siRNAforge

## Architecture Snapshot
- CLI (`src/sirnaforge/cli.py`) uses Typer/Rich to call the façade helpers in `workflow.py`, which orchestrate design/thermo/off-target engines under `core/`, data providers in `data/`, and Nextflow/Docker hooks in `pipeline/`. Skim `docs/developer/architecture.md` before deep rewrites.
- Outputs flow: gene search → transcript validation → siRNA design (`core/design.py`) → thermodynamic scoring (`core/thermodynamics.py`) → optional off-target + Nextflow (`core/off_target.py`, `pipeline/`). Reuse the Pydantic models in `models/` and Pandera schemas in `validation/` instead of rolling new dict structures.
- Gene/ORF retrieval (`data/gene_search.py`, `data/orf_analysis.py`) talks to Ensembl/NCBI; inject fixtures or cached FASTA files (`examples/`, `tests/data/`) when you need deterministic tests.

## Source Layout & Conventions
- Strict `src/` layout; new packages go under `src/sirnaforge/` and must be wired through `workflow.py` so the CLI stays thin.
- CLI subcommands live in a single Typer app; expose new flows by adding `@app.command` functions that delegate to orchestration functions in `workflow.py`.
- Examples and documentation leverage `examples/*.py` and `examples/modification_patterns/*.json`; prefer reusing these assets (also mirrored in `docs/usage_examples.md`).

## Developer Environment
- Bootstrap once with `make dev` (wraps `uv sync --dev` and installs pre-commit). Do not bypass this—Ruff/mypy hooks rely on the managed `.venv`.
- Daily commands run through `uv`: `uv run sirnaforge workflow TP53 --output-dir results`, `uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test.csv`.
- Packaging is Hatchling-only (`pyproject.toml`, `uv.lock`); avoid introducing `setup.py` or ad-hoc requirements files.

## CI & Make Commands
- `make lint` runs Ruff check/format (check-only) plus mypy; `make format` applies Ruff format/fixes. `make check` chains `format` then `test-dev`, so expect unstaged formatting edits.
- Test tiers map to pytest markers defined in `pyproject.toml`: `make test-dev` (`-m dev`), `make test-ci` (smoke + coverage, host-only), `make test-release` (host + Docker, merges coverage), `make test` (all, tolerant of skips).
- Requirement slices: `make test-requires-docker`, `make test-requires-nextflow`, `make test-requires-network`; use these before touching external tooling.
- Docker workflow: `make docker-build` (tags version + latest), `make docker-test` (runs `tests/container` inside image), `make docker-build-test` (clean debug dirs → build → test), `make docker-run GENE=TP53` (smoke workflow in container).
- Workflow CLI: prefer `--input-fasta` when using pre-made transcripts (container tests); override genome indices with `--offtarget-indices human:/abs/GRCh38,mouse:/abs/GRCm39` to replace cached defaults.
- Utilities worth knowing: `make docs` / `make docs-serve`, `make security` (bandit + safety JSON reports), `make cache-info` (shows transcriptome/miRNA cache mounts), `make example` (designs bundled FASTA).

## Testing Matrix
- Fast iteration: `make lint`, `make test-dev` (pytest `-m dev`, ~15s). `make check` auto-formats as part of the run, so expect files to change.
- Broader gates: `make test-ci`, `make test-release`, or marker slices (`make test-requires-docker`, `make test-requires-nextflow`, `make test-requires-network`).
- Container validation: `make docker-build`, `make docker-test`, `make docker-run GENE=TP53`; anything in `tests/container/` assumes BWA-MEM2/SAMtools/ViennaRNA and only passes inside Docker.

## Pipelines & External Integrations
- Nextflow pipeline assets live under `src/sirnaforge/pipeline/nextflow/workflows/` and are typically executed via `NextflowRunner` (Python) or via `nextflow run <embedded main.nf>`.
- Off-target steps require BWA-MEM2, SAMtools, ViennaRNA; these are pinned in `docker/Dockerfile`. On bare metal they are optional, so gate calls behind capability checks.
- Workflow outputs go under `workflow_output/` or debugging dirs like `tp53_workflow_debug/`. Preserve this structure when adding new reporters so downstream tools keep working.

## Documentation & Support Files
- Sphinx docs are under `docs/`; run `make docs` or `make docs-serve` whenever you touch `.rst` files (notably `docs/cli_reference.md`, `docs/developer/testing_guide.md`).
- `examples/complete_workflow_example.py` and `docs/tutorials/python_api.md` are the canonical references for end-to-end usage—keep them updated when signatures change.

## Manual Verification
- Host sanity check: `uv run sirnaforge design examples/sample_transcripts.fasta -o /tmp/test.csv` (completes in ~2s and produces a CSV).
- Container/pipeline check: `make docker-run GENE=TP53` or `sirnaforge workflow TP53 --output-dir tp53_workflow_debug` and inspect `tp53_workflow_debug/orf_reports/*.txt` plus the summary CSVs.
