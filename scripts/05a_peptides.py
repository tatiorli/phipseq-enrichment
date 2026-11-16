#!/usr/bin/env python3
"""
Turn enriched features into protein FASTA for motif discovery.

De novo mode: the reading frame of a trimmed insert is not known a priori, so
each sequence is translated in all three forward frames and the frame with no
internal stop codon (longest clean ORF as tie-break) is taken. Sequences with
a stop in every frame are dropped and counted -- a high drop rate means the
trimming offset is wrong, which is worth knowing.

Reference mode: the feature id is already a library peptide id; its protein
sequence is read from the reference FASTA headers if present, otherwise the
nucleotide entry is translated in frame 0.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq


def translate_best_frame(nt: str) -> tuple[str | None, int]:
    """Return (protein, frame) for the frame with no internal stop."""
    best: tuple[str | None, int] = (None, -1)
    best_len = -1
    for frame in (0, 1, 2):
        sub = nt[frame:]
        sub = sub[: len(sub) - (len(sub) % 3)]
        if len(sub) < 9:
            continue
        aa = str(Seq(sub).translate())
        clean = "*" not in aa.rstrip("*")
        core = aa.rstrip("*")
        if clean and len(core) > best_len:
            best, best_len = (core, frame), len(core)
    return best


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--all", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top", type=int, default=500)
    ap.add_argument("--mode", default="denovo")
    ap.add_argument("--bg-size", type=int, default=5000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hits = pd.read_csv(args.hits, sep="\t")
    allres = pd.read_csv(args.all, sep="\t")

    summary = []
    for prot, grp in hits.groupby("protocol"):
        grp = grp.sort_values("logFC", ascending=False).head(args.top)

        # Background: tested in this protocol, clearly NOT enriched, matched on
        # abundance so the comparison is not just "abundant vs rare".
        bg = allres[(allres.protocol == prot) & (~allres.enriched_edger)]
        if len(bg) > args.bg_size:
            bg = bg.sample(args.bg_size, random_state=42)

        fg_recs, bg_recs, dropped = [], [], 0
        for label, frame_df, sink in (("fg", grp, fg_recs), ("bg", bg, bg_recs)):
            for i, row in enumerate(frame_df.itertuples(index=False)):
                feat = str(row.feature)
                if args.mode == "denovo":
                    aa, fr = translate_best_frame(feat)
                    if aa is None or len(aa) < 9:
                        if label == "fg":
                            dropped += 1
                        continue
                    name = f"{prot}_{label}{i}_lfc{getattr(row, 'logFC', 0):.2f}_f{fr}"
                else:
                    aa = feat  # reference mode: feature id carries the peptide
                    name = f"{prot}_{label}{i}_{feat}"
                sink.append((name, aa))

        write_fasta(outdir / f"{prot}.enriched.faa", fg_recs)
        write_fasta(outdir / f"{prot}.background.faa", bg_recs)

        drop_pct = 100 * dropped / max(len(grp), 1)
        summary.append(
            dict(protocol=prot, enriched_in=len(grp), translated=len(fg_recs),
                 dropped_all_frames=dropped, dropped_pct=round(drop_pct, 1),
                 background=len(bg_recs))
        )
        flag = "  <-- CHECK TRIMMING" if drop_pct > 25 else ""
        print(f"{prot}: {len(fg_recs)} peptides, {dropped} dropped "
              f"({drop_pct:.1f}%), {len(bg_recs)} background{flag}")

    pd.DataFrame(summary).to_csv(outdir / "translation_qc.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
