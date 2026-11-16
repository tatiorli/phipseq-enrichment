#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# One cheap job that answers: will the real pipeline run at all?
#
#   bash scripts/00_smoke_test.sh            # submits itself to LSF
#   bash scripts/00_smoke_test.sh --local    # run here, no bsub
#
# Checks, in order:
#   1. are we actually inside the container
#   2. is the project directory visible AND writable from inside it
#   3. is every tool the pipeline calls on PATH
#   4. do the Python and R scripts parse
#   5. does the sample sheet produce the expected number of runs
#
# Any of these failing now costs seconds. Failing later costs a download,
# a queue wait, and a confusing log.
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/.."
export PROJECT="$PWD"
source config/params.sh

if [ "${1:-}" != "--local" ] && [ -z "${LSB_JOBID:-}" ]; then
  # shellcheck disable=SC2207
  CFLAGS=($(container_bsub_flags))
  echo "submitting smoke test (mode=$CONTAINER_MODE)"
  exec bsub -q "$QUEUE_SHORT" "${CFLAGS[@]}" \
      -J phip_smoke -o logs/smoke.%J.out -e logs/smoke.%J.err \
      -R "rusage[mem=4000]" -n 1 -I \
      "bash $PROJECT/scripts/00_smoke_test.sh --local"
fi

fail=0
ok ()   { printf '  [ ok ] %s\n' "$*"; }
bad ()  { printf '  [FAIL] %s\n' "$*"; fail=1; }

echo "=== 1. container ==="
if [ -f /.dockerenv ]; then
  ok "inside a Docker container"
elif [ -n "${SINGULARITY_CONTAINER:-}${APPTAINER_CONTAINER:-}" ]; then
  ok "inside a Singularity container"
else
  if [ "$CONTAINER_MODE" = "none" ]; then
    ok "CONTAINER_MODE=none, running on the host as intended"
  else
    bad "NOT in a container although CONTAINER_MODE=$CONTAINER_MODE"
    bad "  -> the -app/-env flags did not take effect; check the image name"
  fi
fi
echo "  host: $(hostname)   user: $(id -un)   uid: $(id -u)"
[ -r /etc/os-release ] && echo "  os:   $(. /etc/os-release; echo "$PRETTY_NAME")"

echo "=== 2. filesystem ==="
if [ -d "$PROJECT" ]; then ok "project visible: $PROJECT"; else bad "project NOT visible: $PROJECT"; fi
probe="$PROJECT/logs/.smoke.$$"
if touch "$probe" 2>/dev/null; then ok "project is writable"; rm -f "$probe"
else bad "project is NOT writable from inside the container"; fi
if [ -s config/samples.tsv ]; then ok "sample sheet readable"; else bad "cannot read config/samples.tsv"; fi
echo "  free space: $(df -h "$PROJECT" | awk 'NR==2 {print $4}') on $(df -h "$PROJECT" | awk 'NR==2 {print $1}')"

echo "=== 3. tools ==="
for t in curl pigz fastp seqkit bowtie2 samtools python3 Rscript streme mafft; do
  if command -v "$t" >/dev/null 2>&1; then ok "$t  -> $(command -v "$t")"; else bad "$t missing"; fi
done

echo "=== 4. scripts parse ==="
for p in scripts/03_matrix.py scripts/05a_peptides.py; do
  if python3 -m py_compile "$p" 2>/dev/null; then ok "$p"; else bad "$p does not compile"; fi
done
for r in scripts/04_enrich.R scripts/06_figure.R; do
  if Rscript -e "invisible(parse('$r'))" 2>/dev/null; then ok "$r"; else bad "$r does not parse"; fi
done
if Rscript -e 'suppressPackageStartupMessages({library(edgeR);library(data.table);library(universalmotif);library(patchwork)})' 2>/dev/null
then ok "R packages load"; else bad "an R package is missing from the image"; fi

echo "=== 5. sample sheet ==="
n=$(included_runs | wc -l)
if [ "$n" -gt 0 ]; then ok "$n runs marked include=yes"; else bad "no runs selected"; fi
included_runs | paste -sd' ' -

echo
if [ "$fail" -eq 0 ]; then
  echo "SMOKE TEST PASSED - safe to run: bash run_all.sh"
else
  echo "SMOKE TEST FAILED - fix the [FAIL] lines above before submitting."
fi
exit "$fail"
