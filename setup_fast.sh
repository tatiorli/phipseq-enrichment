#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Fastest route to a real result on Wexac. No Docker image needed.
#
#   cd ~/docker-testing && bash setup_fast.sh
#
# Two conda envs, not one:
#   phip    python, read tools, MEME   (MEME pins icu 72)
#   phip-r  R, edgeR, plotting          (R 4.3 wants icu 73+)
# They cannot coexist; nothing needs both at once.
#
# Then one LSF job takes the EL protocol (4 runs, ~350 MB) from download to
# figure. Resumable: rerun this script and finished steps are skipped.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
export PROJECT="$PWD"
mkdir -p logs

# --- Wexac specifics discovered by wexac_probe.sh / recorded in CLAUDE.md ---
QUEUE="${QUEUE:-short}"          # 24 h cap; use 'medium' if that is not enough
NSLOTS="${NSLOTS:-4}"

# On this cluster rusage[mem=] is PER SLOT, not per job. Asking for
# mem=24000 with -n 8 reserved 192 GB, which is why the submit output said
# "total amount of memory reservation for this job is (MB): 192000".
# Specify what the JOB needs; the per-slot figure is derived.
MEM_TOTAL_MB="${MEM_TOTAL_MB:-16000}"
MEM_MB=$(( MEM_TOTAL_MB / NSLOTS ))
[ "$MEM_MB" -lt 1024 ] && MEM_MB=1024   # never fall under the 1 GB default
# Some nodes raise "Illegal instruction" on conda numpy/pandas, and cn240 /
# cn788 have a history of killing long jobs. 03_matrix.py uses pandas, so the
# job is pinned to hosts known to work.
HOSTS="${HOSTS:-cn074 cn111 cn114 cn116 cn120 cn123 cn125 cn127 cn129 cn132 cn473 cn650 cn653 cn776 cn789 cn801 cn804 cn809 cn810 cn812 cn815 cn817 cn818 cn821 cn823 cn831 cn833 cn834 cn835 cn836 cn838 cn839 cn840 hgn14 hgn15 hgn16 hgn22 hgn23 hgn24 hgn25 hgn27 hgn28 hgn29 hgn30 hgn31 hgn32 hgn37}"

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

make_env () {           # $1 = env name, $2 = yaml
  if conda env list | awk '{print $1}' | grep -qx "$1"; then
    echo "  env '$1' exists, skipping"
  else
    echo "  creating '$1' from $2"
    conda env create -n "$1" -f "$2"
  fi
}

echo "=== 1. conda environments ==="
make_env phip   env-py.yaml
make_env phip-r env-r.yaml

echo "=== 2. verify both envs ==="
fail=0
# Resolve each env's real prefix by ASKING conda, never by guessing.
# Two wrong guesses already:
#   1. `conda run -n ENV bash -lc "command -v TOOL"` - `bash -l` is a LOGIN
#      shell, which re-sources /etc/profile and resets PATH back to base, so
#      only tools that exist in base reported "ok".
#   2. "$CONDA_BASE/envs/$env" - base conda lives under /apps and is
#      read-only, so conda writes user envs to ~/.conda/envs instead.
env_prefix () {
  conda env list 2>/dev/null | awk -v e="$1" '$1==e {print $NF; exit}'
}

check () {              # $1 = env, rest = tools
  local env="$1" p; shift
  p="$(env_prefix "$env")"
  if [ -z "$p" ] || [ ! -d "$p" ]; then
    echo "  MISS env '$env' not found by 'conda env list'"
    conda env list | sed 's/^/        /'
    fail=1; return
  fi
  printf '  env  %-8s %s\n' "$env" "$p"
  for t in "$@"; do
    if [ -x "$p/bin/$t" ]; then
      printf '  ok   %-10s\n' "$t"
    else
      printf '  MISS %-10s (no %s/bin/%s)\n' "$t" "$p" "$t"; fail=1
    fi
  done
}
check phip   python3 fastp seqkit bowtie2 samtools pigz streme mafft
check phip-r Rscript
if ! conda run -n phip-r Rscript -e \
      'suppressPackageStartupMessages({library(edgeR);library(data.table);library(universalmotif);library(patchwork);library(ggseqlogo)})' 2>/dev/null; then
  echo "  MISS one or more R packages"; fail=1
else
  echo "  ok   R packages load"
fi
if [ "$fail" -ne 0 ]; then
  echo
  echo "Environment incomplete. What actually got installed:"
  for e in phip phip-r; do
    p="$(env_prefix "$e")"
    echo "  --- $e  ($p) ---"
    ls -1 "$p/bin" 2>/dev/null | head -40 | sed 's/^/    /'
  done
  echo
  echo "If a tool is genuinely missing, add it:  conda install -n <env> -c bioconda <tool>"
  exit 1
fi

echo "=== 3. select protocols: ${PROTOCOLS:-EL} ==="
# PROTOCOLS="EL"                 4 runs,  ~0.35 GB   smoke test
# PROTOCOLS="EL S2 DL DLPW"     16 runs,  ~4.6 GB    the real comparison
[ -f config/samples.full.tsv ] || cp config/samples.tsv config/samples.full.tsv
awk -F'\t' -v want="${PROTOCOLS:-EL}" 'BEGIN{OFS="\t"; n=split(want,a," "); for(i=1;i<=n;i++) keep[a[i]]=1}
     NR==1{print;next} {$7=($3 in keep ? "yes" : "no"); print}' \
    config/samples.full.tsv > config/samples.tsv
awk -F'\t' 'NR==1||$7=="yes"' config/samples.tsv | column -t -s$'\t'

echo "=== 4. manifest ==="
sed -i 's/\r$//' scripts/*.sh scripts/*.lsf config/params.sh run_all.sh 2>/dev/null || true
sed -i 's/^export CONTAINER_MODE=.*/export CONTAINER_MODE="none"/' config/params.sh
bash scripts/00_manifest.sh

echo "=== 5. build the job script ==="
cat > quickstart.lsf <<EOF
#!/usr/bin/env bash
#BSUB -J phip_quick
#BSUB -o logs/quick.%J.out
#BSUB -e logs/quick.%J.err
#BSUB -q $QUEUE
#BSUB -n $NSLOTS
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=$MEM_MB]
#BSUB -M $MEM_MB
#BSUB -m "$HOSTS"
#BSUB -W 12:00
EOF
cat >> quickstart.lsf <<'EOF'
set -euo pipefail
cd "$LS_SUBCWD"
source "$(conda info --base)/etc/profile.d/conda.sh"
export PROJECT="$PWD"

echo ">>> host: $(hostname)  started: $(date)"

conda activate phip
source config/params.sh
n=$(included_runs | wc -l)
echo ">>> $n runs selected"

echo ">>> downloading"
for i in $(seq 1 "$n"); do LSB_JOBINDEX=$i bash scripts/01_download.lsf; done

echo ">>> counting"
for i in $(seq 1 "$n"); do LSB_JOBINDEX=$i bash scripts/02_count.lsf; done

echo ">>> matrix"
python3 scripts/03_matrix.py

echo ">>> enrichment (R env)"
conda deactivate; conda activate phip-r
Rscript scripts/04_enrich.R

echo ">>> motifs (python/MEME env)"
conda deactivate; conda activate phip
source config/params.sh
bash scripts/05_motifs.sh || echo "!! motif stage failed, continuing"

echo ">>> figure (R env)"
conda deactivate; conda activate phip-r
Rscript scripts/06_figure.R || echo "!! figure stage failed"

echo ">>> DONE: $(date)"
ls -la results results/figures 2>/dev/null || true
EOF

echo "=== 6. submit ==="
echo "  queue=$QUEUE  slots=$NSLOTS  mem=${MEM_MB}MB/slot  total=${MEM_TOTAL_MB}MB"
echo "  (check the 'total amount of memory reservation' line below matches)"
bsub < quickstart.lsf

cat <<'EOT'

Watch it:
  bjobs -w
  bpeek -J phip_quick
  tail -f logs/quick.*.out

Results when it lands:
  results/enrichment_summary.tsv
  results/enriched_peptides.tsv
  results/figures/protocol_comparison.pdf

Widen to all four protocols afterwards:
  cp config/samples.full.tsv config/samples.tsv
  bash scripts/00_manifest.sh && bash run_all.sh
EOT
