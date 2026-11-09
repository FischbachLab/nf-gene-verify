# Build for the same Linux x86_64 target as ec2-linux-64.explicit.txt:
# docker build --platform linux/amd64 -t fig1-29:1.0.0 .
FROM mambaorg/micromamba:2.9.0@sha256:7b93def90e48d6d7744938e31c29172e03b19bc9f56f0db512039f18efd21529

LABEL org.opencontainers.image.title="Fig. 1-style bacterial genome pipeline" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.description="Prokka, Panaroo/MAFFT, IQ-TREE and marker plotting for Linux x86_64"

# Install exact package artifacts, including Prokka's bundled/indexed databases.
# Preserve the base image's non-root user and environment-activation shell.
COPY --chown=$MAMBA_USER:$MAMBA_USER ec2-linux-64.explicit.txt /tmp/fig1-environment.lock
RUN test "$(uname -m)" = x86_64 && \
    micromamba install --yes --name base --file /tmp/fig1-environment.lock && \
    micromamba clean --all --yes

# Copy code only. Genome assemblies, tables and results are runtime mounts.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    run_pipeline.py fig1lib.py plot_fig1.py ortholog_pca.py /opt/fig1/

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLCONFIGDIR=/tmp/fig1-matplotlib \
    XDG_CACHE_HOME=/tmp/fig1-cache \
    LC_ALL=C

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN prokka --listdb && python /opt/fig1/run_pipeline.py --doctor

WORKDIR /tmp
# Keep the official entrypoint so the Conda environment is activated for all tools.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/fig1/run_pipeline.py"]
CMD ["--help"]
