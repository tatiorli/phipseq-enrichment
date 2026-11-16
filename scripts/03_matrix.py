#!/usr/bin/env python3
"""
Build the feature x sample count matrix from the per-sample count files.

    python scripts/03_matrix.py

Writes:
    results/counts_matrix.parquet   features x samples, raw counts
    results/counts_matrix.tsv.gz    same, for eyeballing
    results/library_qc.tsv          per-sample QC used in the report

Memory note: de novo counting produces millions of distinct sequences. The
matrix is built in two passes so we never hold every sequence from every
sample at once -- pass 1 decides which features survive the prevalence filter,
pass 2 fills only those rows. On Wexac, give this job rusage[mem=32000].
"""
from __future__ import annotations

import gzip
import os
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
COUNT_DIR = PROJECT / "counts"
RESULT_DIR = PROJECT / "results"
MANIFEST = PROJECT / "config" / "manifest.tsv"

MIN_COUNT = int(os.environ.get("MIN_COUNT", 5))
MIN_SAMPLES = int(os.environ.get("MIN_SAMPLES", 2))


def read_counts(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            feat, _, cnt = line.rstrip("\n").partition("\t")
            if cnt:
                out[feat] = int(cnt)
    return out


def main() -> int:
    if not MANIFEST.exists():
        sys.exit("Run scripts/00_manifest.sh first")

    man = pd.read_csv(MANIFEST, sep="\t")
    samples = man["sample_name"].tolist()

    paths = {s: COUNT_DIR / f"{s}.counts.tsv.gz" for s in samples}
    missing = [s for s, p in paths.items() if not p.exists()]
    if missing:
        sys.exit(f"Missing count files for: {', '.join(missing)}")

    # ---- pass 1: prevalence -------------------------------------------------
    print(f"pass 1: prevalence filter (>= {MIN_COUNT} reads in >= {MIN_SAMPLES} samples)")
    prevalence: Counter[str] = Counter()
    per_sample_total: dict[str, int] = {}
    per_sample_feats: dict[str, int] = {}
    for s, p in paths.items():
        c = read_counts(p)
        per_sample_total[s] = sum(c.values())
        per_sample_feats[s] = len(c)
        prevalence.update(f for f, n in c.items() if n >= MIN_COUNT)
        print(f"  {s:<16} {len(c):>10,} features  {per_sample_total[s]:>12,} reads")
        del c

    keep = {f for f, k in prevalence.items() if k >= MIN_SAMPLES}
    print(f"  features kept: {len(keep):,} of {len(prevalence):,}")
    if not keep:
        sys.exit("No features survived the prevalence filter. Lower MIN_COUNT/MIN_SAMPLES.")
    keep_index = pd.Index(sorted(keep), name="feature")

    # ---- pass 2: fill -------------------------------------------------------
    print("pass 2: building matrix")
    mat = pd.DataFrame(0, index=keep_index, columns=samples, dtype="int32")
    for s, p in paths.items():
        c = read_counts(p)
        col = pd.Series(c, dtype="int64").reindex(keep_index).fillna(0).astype("int32")
        mat[s] = col.to_numpy()
        del c, col
        print(f"  filled {s}")

    RESULT_DIR.mkdir(exist_ok=True)
    mat.to_parquet(RESULT_DIR / "counts_matrix.parquet")
    mat.to_csv(RESULT_DIR / "counts_matrix.tsv.gz", sep="\t", compression="gzip")

    # ---- QC -----------------------------------------------------------------
    qc = man.set_index("sample_name").copy()
    qc["counted_reads"] = pd.Series(per_sample_total)
    qc["unique_features"] = pd.Series(per_sample_feats)
    qc["in_matrix"] = mat.sum(axis=0)
    qc["frac_in_matrix"] = (qc["in_matrix"] / qc["counted_reads"]).round(4)
    # Diversity: how concentrated is the sample? A successful IP collapses it.
    cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
    qc["top10_pct"] = (
        mat.apply(lambda c: c.nlargest(10).sum() / max(c.sum(), 1) * 100).round(2)
    )
    import numpy as np

    p = cpm / 1e6
    qc["shannon"] = (-(p * np.log(p.where(p > 0))).sum(axis=0)).round(3)
    qc.to_csv(RESULT_DIR / "library_qc.tsv", sep="\t")

    print(f"\nmatrix: {mat.shape[0]:,} features x {mat.shape[1]} samples")
    print(qc[["protocol", "arm", "counted_reads", "top10_pct", "shannon"]].to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
