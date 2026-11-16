#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Central configuration. Every script sources this file.
# Edit the top block; leave the derived paths alone unless you move the project.
# ---------------------------------------------------------------------------

# --- Wexac / LSF -----------------------------------------------------------
# Queue names taken from the Wexac container wiki examples. Confirm with
# `bqueues` and adjust if your lab has dedicated queues.
export QUEUE_SHORT="short"
export QUEUE_MEDIUM="short"

# Wexac reserves 1024 MB per job by default and KILLS anything that exceeds it.
# Every bsub in this project sets rusage[mem=] explicitly. Never remove those.
export MEM_DOWNLOAD=4000      # MB
export MEM_COUNT=16000        # MB
export MEM_MATRIX=32000       # MB
export MEM_STATS=16000        # MB
export THREADS=8

# --- Container -------------------------------------------------------------
# Wexac supports BOTH container routes (see hpcwiki /en/home/lsf/container):
#
#   docker-lsf  LSF launches the whole job inside the Docker image, via the
#               docker-cpu / nvidia-gpu application profiles. The image lives
#               in the internal registry ops:5000 (or Docker Hub). Nothing to
#               convert, nothing to copy. This is the default.
#
#   singularity You convert the image to a .sif on docker-gpu, scp it to your
#               home, and wrap each command. Useful if you want one portable
#               file, or need --fakeroot editing.
#
#   none        No container: conda env or system modules.
#
export CONTAINER_MODE="none"

# docker-lsf: the image tag LSF should pull.
#   internal registry : ops:5000/phip-pipeline:0.1   (after you docker push)
#   docker hub        : <youruser>/phip-pipeline:0.1
export IMAGE="ops:5000/phip-pipeline:0.1"

# singularity: path to the .sif in your Wexac home
export SIF="${PROJECT:-$PWD}/phip.sif"
# Wexac needs the module loaded before `singularity` is on PATH
export SINGULARITY_MODULE="Singularity"

# Lab/data paths to expose inside the container (singularity mode).
# In docker-lsf mode the bind mounts are configured by the LSF admins in the
# application profile, not by you -- which filesystems are visible is a site
# decision. Run scripts/00_smoke_test.sh once to confirm this project
# directory is readable AND writable from inside the image before you submit
# the real pipeline.
export BIND_PATHS=""   # e.g. "/home/labs/yourlab:/home/labs/yourlab"

# --- Analysis --------------------------------------------------------------
# denovo   : count unique trimmed inserts. Works with no reference. Default.
# reference: map to a VirScan peptide FASTA (set REF_FASTA below).
export COUNT_MODE="denovo"
export REF_FASTA="${PROJECT:-$PWD}/ref/virscan_peptides.fasta"

# All included runs must be compared at a common read length.
# The EL/S2/DL/DLPW arms are 50 bp; the RD arm is 151 bp (see README).
export TRUNC_LEN=50

# Adapter trimming. fastp auto-detects for single-end, but if you know the
# constant 5'/3' flanks of the oligo, put them here and they are used instead.
export ADAPTER_3P=""          # e.g. "AGATCGGAAGAGC"

# Minimum count for a sequence to enter the matrix (per sample).
export MIN_COUNT=5
# A sequence must clear MIN_COUNT in at least this many samples.
export MIN_SAMPLES=2

# Enrichment calling
export FDR_CUTOFF=0.05
export LOGFC_CUTOFF=1         # log2 fold change over the matched bead control

# Motif discovery
export N_TOP_PEPTIDES=500     # top enriched peptides per protocol fed to STREME
export MOTIF_WIDTH_MIN=5
export MOTIF_WIDTH_MAX=12

# --- Derived paths (leave alone) -------------------------------------------
export PROJECT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export SAMPLES="$PROJECT/config/samples.tsv"
export FASTQ_DIR="$PROJECT/fastq"
export COUNT_DIR="$PROJECT/counts"
export RESULT_DIR="$PROJECT/results"
export LOG_DIR="$PROJECT/logs"
mkdir -p "$FASTQ_DIR" "$COUNT_DIR" "$RESULT_DIR" "$LOG_DIR"

# Load the Singularity module if Wexac's module system is available.
load_singularity () {
  command -v singularity >/dev/null 2>&1 && return 0
  if command -v module >/dev/null 2>&1; then
    module load "$SINGULARITY_MODULE" 2>/dev/null || true
  elif [ -f /etc/profile.d/lmod.sh ]; then
    . /etc/profile.d/lmod.sh && module load "$SINGULARITY_MODULE" 2>/dev/null || true
  fi
  command -v singularity >/dev/null 2>&1
}
export -f load_singularity

# Run a command inside the container.
#   docker-lsf : the job is ALREADY inside the image (LSF put it there), so
#                this is a pass-through. Wrapping again would be wrong.
#   singularity: wrap each invocation.
#   none       : pass through.
run_tool () {
  if [ "$CONTAINER_MODE" = "singularity" ]; then
    load_singularity || { echo "singularity not available (module $SINGULARITY_MODULE)" >&2; return 127; }
    local binds=(-B "$PROJECT":"$PROJECT")
    [ -n "$BIND_PATHS" ] && binds+=(--bind "$BIND_PATHS")
    singularity exec --cleanenv "${binds[@]}" "$SIF" "$@"
  else
    "$@"
  fi
}
export -f run_tool

# The extra bsub flags needed to put a job inside the image.
# Empty for singularity/none, because there the wrapping happens inside.
container_bsub_flags () {
  if [ "$CONTAINER_MODE" = "docker-lsf" ]; then
    printf -- '-app docker-cpu -env LSB_CONTAINER_IMAGE=%s' "$IMAGE"
  fi
}
export -f container_bsub_flags

# Helper: the list of runs with include=yes
included_runs () { awk -F'\t' 'NR>1 && $7=="yes" {print $1}' "$SAMPLES"; }
export -f included_runs
