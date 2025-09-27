# Comprehensive Docker image for siRNAforge
# Includes all dependencies: Python, bioinformatics tools, and workflow management
FROM ghcr.io/astral-sh/uv:python3.11-bookworm-slim AS builder
ARG VERSION=0.0.0+local
ENV BUILD_VERSION=${VERSION}

# Ensure we start as root to install system packages and create a non-root user
USER root

# Install system dependencies for building
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    git \
    curl \
    unzip \
    # openjdk removed from builder to avoid unnecessary download; runtime will
    # provide Java via conda (openjdk) in the runtime environment.
    \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user early and ensure /app is owned by them so build
# artifacts are not owned by root.
RUN useradd -m -s /bin/bash -u 1000 sirnauser && mkdir -p /app && chown -R sirnauser:sirnauser /app /home/sirnauser

# Switch to the non-root user for subsequent build steps so outputs are
# created with the correct ownership.
USER sirnauser
WORKDIR /app

# Copy uv configuration and dependencies as the non-root user
COPY --chown=sirnauser:sirnauser pyproject.toml uv.lock README.md LICENSE ./

# Install Python dependencies with uv (runs as sirnauser)
RUN uv sync --frozen --no-dev

# Copy source code as the non-root user
COPY --chown=sirnauser:sirnauser src/ ./src/

# Build the package (as sirnauser)
RUN uv build

# Conda builder stage: create the conda environment in a separate stage so we
# can copy it into the runtime image without an expensive recursive chown.
FROM mambaorg/micromamba:bookworm-slim AS conda-builder
# Copy environment and install into /opt/conda in this stage
COPY docker/environment-nextflow.yml /tmp/environment.yml
RUN micromamba install --prefix /opt/conda --file /tmp/environment.yml -y && \
    micromamba clean --all --yes

# Unit test stage (fast) - runs unit tests using the builder image
FROM ghcr.io/astral-sh/uv:python3.11-bookworm-slim AS unit-test
ARG VERSION=0.0.0+local
ENV BUILD_VERSION=${VERSION}
USER root

# Create non-root user and prepare home; we must create the user before COPY
# so COPY --chown can resolve the user/group name.
RUN useradd -m -s /bin/bash -u 1000 sirnauser && mkdir -p /home/sirnauser

# Copy the pre-built conda prefix from the conda-builder stage and set ownership
# to the non-root user at layer creation time (fast, avoids recursive chown).
COPY --from=conda-builder --chown=sirnauser:sirnauser /opt/conda /opt/conda

# Switch to non-root user and install the application wheel into the conda
# prefix (pip will satisfy pure-Python deps from wheels).
USER sirnauser
WORKDIR /home/sirnauser

# Copy and install the application wheel
COPY --from=builder --chown=sirnauser:sirnauser /app/dist/*.whl /tmp/
RUN /opt/conda/bin/pip install --no-cache-dir /tmp/*.whl && rm -f /tmp/*.whl

# Switch to the non-root user for installing the application wheel so files
# created during installation are owned by sirnauser.
USER sirnauser
WORKDIR /home/sirnauser

# Set environment variables
# Ensure conda/micromamba binaries are on PATH so installed tools (nextflow,
# bwa-mem2, samtools, ViennaRNA, and the installed sirnaforge entrypoint) are
# available to healthchecks and the default CMD which run via /bin/sh.
ENV PATH="/opt/conda/bin:/opt/conda/condabin:/usr/local/bin:$PATH"
ENV JAVA_HOME="/opt/conda/lib/jvm"
ENV NXF_HOME="/home/sirnauser/.nextflow"
# Configure uv to use the existing conda environment instead of creating its own venv
ENV UV_SYSTEM_PYTHON=1
ENV UV_PYTHON="/opt/conda/bin/python"
ENV VIRTUAL_ENV="/opt/conda"

# Ensure Nextflow home is created and writable by the non-root user so Nextflow
# does not attempt to create runtime directories under /usr/local/bin.
RUN mkdir -p /home/sirnauser/.nextflow && chown -R sirnauser:sirnauser /home/sirnauser/.nextflow

# Copy and setup entrypoint script
COPY --chown=sirnauser:sirnauser docker/entrypoint.sh /usr/local/bin/entrypoint.sh
USER root
RUN chmod +x /usr/local/bin/entrypoint.sh
USER sirnauser

# Health check - verify all major components
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD sirnaforge version && \
        nextflow -version && \
        bwa-mem2 version && \
        samtools --version && \
        RNAfold --version || exit 1

# Set entrypoint and default command
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []

# Metadata
LABEL org.opencontainers.image.title="siRNAforge"
LABEL org.opencontainers.image.description="siRNAforge - Multi-species gene to siRNA design, off-target prediction, and ranking. Comprehensive siRNA design toolkit for gene silencing."
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.authors="Austin S. Hovland"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.source="https://github.com/austin-s-h/sirnaforge"
