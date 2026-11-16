#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Submit the whole pipeline to Wexac with LSF dependencies, so you can fire it
# once and walk away.
#
#   bash run_all.sh
#
# Each stage waits for the previous one to finish cleanly ("done(JOBID)").
# Every stage is idempotent: re-running skips work that already has a .done
# marker, so if a stage dies you resubmit the same command and only the
# missing pieces re-run.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
export PROJECT="$PWD"
source config/params.sh

N=$(included_runs | wc -l)
echo "runs to process: $N"
[ "$N" -gt 0 ] || { echo "nothing to do"; exit 1; }

case "$CONTAINER_MODE" in
  singularity)
    [ -f "$SIF" ] || { echo "ERROR: $SIF not found. See docker/README-wexac.md." >&2; exit 1; }
    echo "container: singularity $SIF" ;;
  docker-lsf)
    echo "container: LSF docker profile, image $IMAGE" ;;
  none)
    echo "container: none (using whatever is on PATH)" ;;
  *) echo "Unknown CONTAINER_MODE=$CONTAINER_MODE" >&2; exit 1 ;;
esac

# shellcheck disable=SC2207
CFLAGS=($(container_bsub_flags))

# 0. manifest (login node, one API call)
bash scripts/00_manifest.sh

# Rewrite the array ranges to match the manifest
sed -i -E "s/^#BSUB -J phip_dl\[1-[0-9]+\]/#BSUB -J phip_dl[1-${N}]/"       scripts/01_download.lsf
sed -i -E "s/^#BSUB -J phip_count\[1-[0-9]+\]/#BSUB -J phip_count[1-${N}]/" scripts/02_count.lsf

submit () { bsub "$@" | sed -nE 's/^Job <([0-9]+)>.*/\1/p'; }

# Stages 3,4,6 call a single tool, so in singularity mode they are wrapped
# here; in docker-lsf mode LSF puts the whole job inside the image and the
# command is run bare. wrap() hides the difference.
wrap () {
  if [ "$CONTAINER_MODE" = "singularity" ]; then
    printf 'ml %s; singularity exec --cleanenv -B %s:%s %s %s' \
      "$SINGULARITY_MODULE" "$PROJECT" "$PROJECT" "$SIF" "$*"
  else
    printf '%s' "$*"
  fi
}

J1=$(submit -q "$QUEUE_SHORT" "${CFLAGS[@]}" < scripts/01_download.lsf)
echo "download     : $J1"

J2=$(submit -q "$QUEUE_MEDIUM" "${CFLAGS[@]}" -w "done($J1)" < scripts/02_count.lsf)
echo "count        : $J2  (after $J1)"

J3=$(submit -q "$QUEUE_SHORT" "${CFLAGS[@]}" -w "done($J2)" \
      -J phip_matrix -o logs/matrix.%J.out -e logs/matrix.%J.err \
      -R "rusage[mem=${MEM_MATRIX}]" -n 2 \
      "$(wrap python3 scripts/03_matrix.py)")
echo "matrix       : $J3  (after $J2)"

J4=$(submit -q "$QUEUE_SHORT" "${CFLAGS[@]}" -w "done($J3)" \
      -J phip_edger -o logs/edger.%J.out -e logs/edger.%J.err \
      -R "rusage[mem=${MEM_STATS}]" -n 4 \
      "$(wrap Rscript scripts/04_enrich.R)")
echo "edgeR        : $J4  (after $J3)"

J5=$(submit -q "$QUEUE_SHORT" "${CFLAGS[@]}" -w "done($J4)" \
      -J phip_motif -o logs/motif.%J.out -e logs/motif.%J.err \
      -R "rusage[mem=${MEM_STATS}]" -n 4 \
      "bash scripts/05_motifs.sh")
echo "motifs       : $J5  (after $J4)"

J6=$(submit -q "$QUEUE_SHORT" "${CFLAGS[@]}" -w "ended($J5)" \
      -J phip_fig -o logs/fig.%J.out -e logs/fig.%J.err \
      -R "rusage[mem=8000]" -n 2 \
      "$(wrap Rscript scripts/06_figure.R)")
echo "figure       : $J6  (after $J5)"

cat <<EOF

Submitted. Watch with:
  bjobs -w
  bjobs -l $J6

Final output when it lands:
  $RESULT_DIR/figures/protocol_comparison.pdf
  $RESULT_DIR/enriched_peptides.tsv
  $RESULT_DIR/enrichment_summary.tsv

If a stage fails, fix it and just run 'bash run_all.sh' again - completed work
is skipped via the .done markers.
EOF
