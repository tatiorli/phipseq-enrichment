# ---------------------------------------------------------------------------
# PhIP-seq enrichment pipeline — reproducible image
#
#   docker build -t phip-pipeline:0.1 -f docker/Dockerfile docker/
#
# TWO conda environments, not one. MEME pins icu 72.x and R 4.3 wants icu 73+,
# so a single environment containing both is unsolvable — this image's first
# version tried it and failed the same way the cluster conda did. Nothing in
# the pipeline needs MEME and R in the same process, so they are kept apart
# and the runner switches PATH between stages.
#
#   /opt/conda/envs/phip     python, fastp, seqkit, bowtie2, samtools, MEME
#   /opt/conda/envs/phip-r   R, edgeR, universalmotif, ggplot2
#
# These two YAML files are the exact ones that solved on the cluster.
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="phip-pipeline" \
      org.opencontainers.image.version="0.1" \
      org.opencontainers.image.description="PhIP-seq enrichment + motif pipeline (PRJNA1438256)" \
      org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates procps gzip less bsdmainutils \
    && rm -rf /var/lib/apt/lists/*
# bsdmainutils provides `column`. The scripts no longer depend on it, but
# the manifest is far easier to read when it is there.
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER env-py.yaml env-r.yaml /tmp/
RUN micromamba env create -y -f /tmp/env-py.yaml \
 && micromamba env create -y -f /tmp/env-r.yaml \
 && micromamba clean --all --yes

ENV PY_BIN=/opt/conda/envs/phip/bin \
    R_BIN=/opt/conda/envs/phip-r/bin \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    HOME=/tmp

# Fail the BUILD if anything is missing, rather than three hours into a run.
# On a miss, print what the solver actually chose -- "MISSING streme" alone
# does not say whether the package is absent or just an old version without
# that binary. streme ships with MEME >= 5.3, which is why env-py.yaml pins it.
RUN set -eux; \
    echo "--- MEME package as resolved ---"; \
    micromamba list -n phip 2>/dev/null | grep -iE '^\s*(meme|perl)\b' || true; \
    echo "--- meme binaries present ---"; \
    ls "$PY_BIN" | grep -iE 'meme|streme|dreme|fimo' || echo "(none)"; \
    for t in python3 fastp seqkit bowtie2 samtools pigz streme mafft; do \
      test -x "$PY_BIN/$t" || { echo "MISSING $t in phip"; ls "$PY_BIN" | head -60; exit 1; }; \
    done; \
    test -x "$R_BIN/Rscript" || { echo "MISSING Rscript in phip-r"; exit 1; }; \
    "$PY_BIN/python3" -c "import pandas, pyarrow, Bio; print('python ok')"; \
    "$R_BIN/Rscript" -e 'suppressPackageStartupMessages({library(edgeR);library(data.table);library(universalmotif);library(patchwork);library(ggseqlogo)}); cat("R ok\n")'; \
    echo "image ok"

WORKDIR /work
CMD ["bash", "docker/run_in_docker.sh"]
