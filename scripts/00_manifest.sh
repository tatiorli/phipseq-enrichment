#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Build the download manifest by asking ENA for the current FASTQ URLs.
# Never hard-code ENA paths: the directory shard depends on the accession and
# changes between submissions.
#
# Run on the login node (it is one API call), not via bsub.
#   bash scripts/00_manifest.sh
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/params.sh"

OUT="$PROJECT/config/manifest.tsv"
API="https://www.ebi.ac.uk/ena/portal/api/filereport"

echo "Querying ENA for PRJNA1438256 ..."
curl -sf --retry 3 --retry-delay 5 \
  "${API}?accession=PRJNA1438256&result=read_run&fields=run_accession,sample_accession,read_count,base_count,fastq_bytes,fastq_ftp,fastq_md5&format=tsv" \
  -o "$PROJECT/config/ena_raw.tsv"

# Join ENA's URL list onto our curated sample sheet, keeping only include=yes.
awk -F'\t' '
  NR==FNR {
    if (FNR>1 && $7=="yes") { keep[$1]=1; name[$1]=$2; prot[$1]=$3; arm[$1]=$4; rep[$1]=$5; rl[$1]=$6 }
    next
  }
  FNR==1 { next }
  ($1 in keep) {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
      $1, name[$1], prot[$1], arm[$1], rep[$1], rl[$1], $3, $5, $6
  }
' "$SAMPLES" "$PROJECT/config/ena_raw.tsv" \
| sort -k3,3 -k4,4 -k5,5n \
| cat <(printf "run\tsample_name\tprotocol\tarm\treplicate\tread_len\treads\tbytes\turl\n") - > "$OUT"

N=$(($(wc -l < "$OUT") - 1))
TOTAL=$(awk -F'\t' 'NR>1 {s+=$8} END {printf "%.2f", s/1e9}' "$OUT")

echo "Manifest: $OUT"
echo "  runs to download : $N"
echo "  total size       : ${TOTAL} GB"
echo
column -t -s$'\t' "$OUT" | cut -c1-120

if [ "$N" -eq 0 ]; then
  echo "ERROR: no runs matched. Check config/samples.tsv." >&2
  exit 1
fi
