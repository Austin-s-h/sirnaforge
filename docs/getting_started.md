# Getting Started & Quick Reference

Launch your first siRNAforge workflow in minutes, then dive deeper through the focused docs sections.

## 1. Install

- **Python/uv users:** Follow the step-by-step guide in [Installation](installation.md) for pip/uv instructions plus development setup (`make dev`) that matches the project’s CI environment.
- **Container users:** Pull or build the full bioinformatics image as described in {ref}`Installation → Docker <docker-full-bioinformatics-stack>`. The container bundles Nextflow, BWA-MEM2, SAMtools, ViennaRNA, and detects `local` mode automatically.

> Need to automate lab servers or run air-gapped? The installation guide also documents offline caching, verification commands, and how to run smoke tests with the bundled FASTA files.

## 2. Run Your First Workflow

```bash
# End-to-end design + scoring for TP53
uv run sirnaforge workflow TP53 --output-dir tp53_run
```

- Use `docker run --rm ghcr.io/austin-s-h/sirnaforge:latest sirnaforge workflow …` for the same behavior inside the prebuilt image.
- The command performs gene search → siRNA design → thermodynamic scoring → filtering and writes CSVs plus logs under `tp53_run/`.

See the [Workflows overview](workflows.md) for the full output tree, common flags (GC bounds, top-N, verbosity), and Nextflow/off-target hand-offs.

## 3. Customize Inputs & References

Pass alternative sequences or transcriptomes to match your experimental data:

```bash
uv run sirnaforge workflow TP53 \
  --input-fasta examples/sample_transcripts.fasta \
  --transcriptome-fasta ensembl_mouse_cdna \
  --output-dir tp53_custom
```

- `--input-fasta` accepts local paths or URLs, bypassing the transcript lookup while retaining familiar file naming via the positional gene argument.
- `--transcriptome-fasta` selects the reference used for transcriptome off-target analysis (local/remote/preset). Provide it whenever you need species other than the bundled Ensembl set.
- `--offtarget-indices` and `--species` feed the BWA-MEM2/Nextflow pipeline; details live in {ref}`Workflows → Nextflow Pipeline <nextflow-pipeline>`.

## 4. Inspect Results

Key files in every workflow run:

| File | Why it matters |
|------|----------------|
| `sirnaforge/*_pass.csv` | High-quality candidates ready for experiments |
| `sirnaforge/*_all.csv` | Full candidate list with all metrics for custom filtering |
| `logs/workflow_summary.json` | Summary of search/design stats, reference decisions, and QC flags |
| `off_target/input_candidates.fasta` | FASTA handed to the Nextflow off-target pipeline |

Use standard CLI tools to browse:

```bash
head -6 tp53_run/sirnaforge/TP53_pass.csv
jq '.' tp53_run/logs/workflow_summary.json | less
```

For a complete walkthrough of directory contents, refer to {ref}`Workflows → Output Structure <output-structure>`.

## 5. Understand the Scores

siRNAforge surfaces the same thermodynamic metrics described in the code under `sirnaforge/core/thermodynamics.py`. Refer to:

- [Scoring Overview](scoring.md) — optimal ranges, filtering presets, and column descriptions.
- [Data Models & Scoring Reference](models_and_scoring.md) — Pydantic models, algorithms, and literature citations.

Typical “green zone” checkpoints:

- `composite_score` ≥ 8 for top-tier candidates
- `asymmetry_score` ≥ 0.65 to ensure the guide strand loads into RISC
- `gc_content` between 35–60% (optimal 40–55%)

## 6. Next Steps

- [Usage Examples](usage_examples.md) for batch automation, miRNA design mode, and modification workflows.
- [CLI Reference](cli_reference.md) for auto-generated `--help` output from every Typer command.
- [Thermodynamic Metrics Guide](thermodynamic_guide.md) for deeper interpretation.
- [Developer Documentation](developer/index.rst) if you plan to extend the codebase.
