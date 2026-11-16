#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Run the whole pipeline inside the container, on whatever is mounted at /work.
#
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work \
#       -e PROTOCOLS="EL" phip-pipeline:0.1
#
# Env vars:
#   PROTOCOLS   space-separated list; default EL (4 runs, ~0.35 GB)
#               "EL S2 DL DLPW" = 16 runs, ~4.6 GB, the full comparison
#
# The two conda envs cannot be active at once (MEME pins icu 72, R wants 73+),
# so each stage gets the PATH it needs. No `conda activate` -- the env bin
# directories are prepended directly, which is one less thing to go wrong.
# ---------------------------------------------------------------------------
set -euo pipefail

PY_BIN="${PY_BIN:-/opt/conda/envs/phip/bin}"
R_BIN="${R_BIN:-/opt/conda/envs/phip-r/bin}"
BASE_PATH="$PATH"
use_py () { export PATH="$PY_BIN:$BASE_PATH"; }
use_r  () { export PATH="$R_BIN:$BASE_PATH"; }

export PROJECT="$PWD"
export CONTAINER_MODE="none"     # we ARE the container; do not nest

# The per-sample stages are LSF job scripts: they start with `cd "$LS_SUBCWD"`,
# the submission directory LSF exports into every job. Outside LSF that
# variable does not exist, and `set -u` turns it into a fatal error. Provide
# it here so the same scripts run unmodified under LSF and under Docker.
export LS_SUBCWD="$PWD"

mkdir -p logs fastq counts results

echo ">>> phip-pipeline $(date -u +%FT%TZ)"
echo ">>> workdir: $PROJECT"
echo ">>> protocols: ${PROTOCOLS:-EL}"

# --- select protocols ------------------------------------------------------
[ -f config/samples.full.tsv ] || cp config/samples.tsv config/samples.full.tsv
awk -F'\t' -v want="${PROTOCOLS:-EL}" \
    'BEGIN{OFS="\t"; n=split(want,a," "); for(i=1;i<=n;i++) keep[a[i]]=1}
     NR==1{print;next} {$7=($3 in keep ? "yes" : "no"); print}' \
    config/samples.full.tsv > config/samples.tsv.new
mv config/samples.tsv.new config/samples.tsv

use_py
sed -i 's/^export CONTAINER_MODE=.*/export CONTAINER_MODE="none"/' config/params.sh
bash scripts/00_manifest.sh

source config/params.sh
n=$(included_runs | wc -l)
echo ">>> $n runs selected"
[ "$n" -gt 0 ] || { echo "no runs selected"; exit 1; }

echo ">>> [1/6] download"
for i in $(seq 1 "$n"); do LSB_JOBINDEX=$i bash scripts/01_download.lsf; done

echo ">>> [2/6] count"
for i in $(seq 1 "$n"); do LSB_JOBINDEX=$i bash scripts/02_count.lsf; done

echo ">>> [3/6] matrix"
python3 scripts/03_matrix.py

echo ">>> [4/6] enrichment"
use_r;  Rscript scripts/04_enrich.R

echo ">>> [5/6] motifs"
use_py; source config/params.sh
bash scripts/05_motifs.sh || echo "!! motif stage failed, continuing"

echo ">>> [6/6] figure"
use_r;  Rscript scripts/06_figure.R || echo "!! figure stage failed"

echo ">>> DONE $(date -u +%FT%TZ)"
ls -la results results/figures 2>/dev/null || true
