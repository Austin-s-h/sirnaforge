# Disposable #121 comparison image. Build from the generated context documented
# in docs/developer/performance_review.md so the project .dockerignore stays intact.
FROM ghcr.io/prefix-dev/pixi:0.76.1@sha256:540414f50b96439a6ba1475b79b2e2dd2208e86ebef33e2fbf2a176a8f037d31 AS builder
WORKDIR /workspace

COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/
COPY scripts/performance/pixi/pixi.toml scripts/performance/pixi/pixi.lock scripts/performance/pixi/

RUN pixi install --locked --manifest-path scripts/performance/pixi/pixi.toml --environment application && \
    pixi run --manifest-path scripts/performance/pixi/pixi.toml --environment application install-app

FROM debian:bookworm-slim@sha256:7ff6ba31527b9e336162fa061746c366c4b2032f4d2cff2078bea43f78efe3f8 AS runtime
ENV PIXI_PREFIX="/workspace/scripts/performance/pixi/.pixi/envs/application" \
    CONDA_PREFIX="/workspace/scripts/performance/pixi/.pixi/envs/application" \
    PATH="/workspace/scripts/performance/pixi/.pixi/envs/application/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    LD_LIBRARY_PATH="/workspace/scripts/performance/pixi/.pixi/envs/application/lib" \
    HOME="/tmp/nextflow/home" \
    JAVA_HOME="/workspace/scripts/performance/pixi/.pixi/envs/application/lib/jvm" \
    NXF_DIST="/workspace/scripts/performance/pixi/.pixi/envs/application/share/nextflow/dist" \
    NXF_HOME="/tmp/nextflow" \
    NXF_CONDA_CLI="micromamba" \
    MAMBA_ROOT_PREFIX="/tmp/nextflow/mamba"

RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /tmp/nextflow/home && chmod -R a+rwX /tmp/nextflow
WORKDIR /workspace
COPY --from=builder /workspace/src/ /workspace/src/
COPY --from=builder /workspace/scripts/performance/pixi/.pixi/envs/application/ /workspace/scripts/performance/pixi/.pixi/envs/application/
COPY scripts/performance/pixi/conda-shim.sh /usr/local/bin/conda
COPY scripts/performance/pixi/micromamba-activate.sh /workspace/scripts/performance/pixi/.pixi/envs/application/bin/activate
RUN chmod 755 /usr/local/bin/conda /workspace/scripts/performance/pixi/.pixi/envs/application/bin/activate

ENTRYPOINT ["sirnaforge"]
CMD ["--help"]
