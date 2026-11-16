#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Build the bowtie2 index once, before submitting the counting array.
# Only needed when COUNT_MODE=reference.
#
#   bsub -q <queue> -R rusage[mem=16000] -n 8 bash scripts/00_build_index.sh
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/params.sh"

if [ "$COUNT_MODE" != "reference" ]; then
  echo "COUNT_MODE=$COUNT_MODE - no index needed. Nothing to do."
  exit 0
fi

if [ ! -f "$REF_FASTA" ]; then
  cat >&2 <<'EOF'
ERROR: no peptide reference found.

The VirScan oligo/peptide library (~115,000 56-aa tiles) is NOT part of
BioProject PRJNA1438256. You have to supply it yourself. Options:

  1. The supplementary tables of the VirScan library publications.
  2. Ask the corresponding author of the protocol paper
     (dev.majumdar@uvm.edu) which library build was used - the paper says
     "over 115,000 epitopes", which does not identify the version.
  3. Stay in COUNT_MODE=denovo. The pipeline runs end to end without a
     reference; you simply get sequences instead of named epitopes, and you
     can annotate the enriched ones afterwards with blastp.

Put a nucleotide FASTA of the displayed inserts at:
  REF_FASTA=$REF_FASTA
EOF
  exit 1
fi

idx="$PROJECT/ref/bt2_index"
echo "Building bowtie2 index from $REF_FASTA"
run_tool bowtie2-build --threads "$THREADS" "$REF_FASTA" "$idx"
run_tool seqkit stats "$REF_FASTA"
echo "Index written to ${idx}.*.bt2"
