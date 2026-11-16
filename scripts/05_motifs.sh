#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Translate the enriched inserts and find shared sequence motifs per protocol.
#
#   bsub -q <queue> -R rusage[mem=16000] -n 4 bash scripts/05_motifs.sh
#
# WHY NOT A LOGO OF THE PEPTIDES THEMSELVES
# -----------------------------------------
# A sequence logo only means something when the sequences are aligned and the
# same length. VirScan peptides are 56-aa tiles taken from hundreds of
# different viral proteins; stacking them produces a uniform smear, not a
# motif. The honest route, and the one used in the epitope literature, is to
# ask what SHORT motif is over-represented in the enriched set relative to a
# background drawn from the same library. That is what STREME does, and its
# output PWM is what the logo in the final figure is built from.
#
# Background: peptides that were tested but NOT enriched in the same protocol.
# This controls for the library's own amino-acid composition, which a shuffled
# background would not.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/params.sh"

HITS="$RESULT_DIR/enriched_peptides.tsv"
ALL="$RESULT_DIR/enrichment_all.tsv.gz"
[ -f "$HITS" ] || { echo "Run 04_enrich.R first" >&2; exit 1; }

MOTIF_DIR="$RESULT_DIR/motifs"
mkdir -p "$MOTIF_DIR"

echo "translating enriched inserts"
run_tool python3 "$PROJECT/scripts/05a_peptides.py" \
  --hits "$HITS" --all "$ALL" --outdir "$MOTIF_DIR" \
  --top "$N_TOP_PEPTIDES" --mode "$COUNT_MODE"

for fa in "$MOTIF_DIR"/*.enriched.faa; do
  [ -e "$fa" ] || { echo "no enriched FASTA produced"; exit 1; }
  prot=$(basename "$fa" .enriched.faa)
  bg="$MOTIF_DIR/${prot}.background.faa"
  out="$MOTIF_DIR/${prot}.streme"

  n_fg=$(grep -c '^>' "$fa" || true)
  n_bg=$(grep -c '^>' "$bg" 2>/dev/null || echo 0)
  echo "[$prot] foreground=$n_fg background=$n_bg"

  if [ "$n_fg" -lt 20 ]; then
    echo "[$prot] fewer than 20 enriched peptides - skipping motif search"
    continue
  fi

  run_tool streme \
    --protein \
    --p "$fa" \
    ${n_bg:+--n "$bg"} \
    --oc "$out" \
    --minw "$MOTIF_WIDTH_MIN" --maxw "$MOTIF_WIDTH_MAX" \
    --nmotifs 5 \
    --seed 42 \
    || echo "[$prot] streme returned non-zero; continuing"
done

echo
echo "motif output under $MOTIF_DIR/<protocol>.streme/streme.xml"
ls -1 "$MOTIF_DIR"/*/streme.txt 2>/dev/null || true
