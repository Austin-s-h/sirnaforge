# Docker deployment

The full image contains siRNAforge, Nextflow and its matching bundled framework JAR, Java,
BWA-MEM2, BWA-MEM3, samtools and ViennaRNA. Native tools are installed from
the root `pixi.lock`. The application wheel is installed without dependency resolution.
Pin a published image digest when exact
runtime identity matters.

## Build and validate

```bash
make docker-build
make docker-test
make docker-shell
```

The build verifies `nextflow -version` with networking disabled. `make docker-test` exercises the
installed package and embedded pipeline; the chromosome 21 fixture bounds CI reference indexing.
Live Ensembl tests are excluded from the default gate; opt in with
`make docker-test DOCKER_TEST_MARKERS=runs_in_container DOCKER_TEST_PYTEST_ARGS=--run-network`
on a local machine. CI rejects this opt-in. `make docker-ensure` checks both the version
and source fingerprint before running container tests.

## Choose the deployment mode

- **Run the full image:** Nextflow detects the container and uses `local` tasks. Index construction
  and alignment both use the BWA-MEM2 already installed in the image. No nested Docker daemon or
  Docker socket mount is required.
- **Run the Python application on the host with Docker tasks:** install Java/Nextflow on the host
  and select the image with `--nextflow-docker-image` or `SIRNAFORGE_NEXTFLOW_IMAGE`. Missing indices
  are built by Nextflow in that same image; BWA-MEM2 need not be installed on the host.
- **Run natively or in an activated Conda environment:** install the complete toolchain and use
  `SIRNAFORGE_USE_LOCAL_EXECUTION=1`. The image includes the locked toolchain and the verified wheel.

Keep input, output, siRNAforge cache and Nextflow work paths accessible to the selected runtime.
Map the host UID/GID to preserve ownership. The image preserves its toolchain `PATH` in login
shells, and its Nextflow home permits an arbitrary invoking UID.

See the [installation and execution-mode guide](../docs/installation.md) for full commands,
Conda/Singularity requirements, mounted cache layout and prebuilt-index selection. Explicit
`--transcriptome-indices` replaces automatic transcriptome downloads; a neighboring cDNA FASTA
is still required for classification. `design --genome-index` has been removed.
