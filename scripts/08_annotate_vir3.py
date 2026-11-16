#!/usr/bin/env python3
"""
Annotate every enriched feature against the VirScan Vir3 library.

    python3 scripts/08_annotate_vir3.py --vir3 vir3_annotations.csv.gz

Why not BLAST
-------------
BLAST infers what a peptide resembles. The library says what it IS. Every
peptide in this experiment either is a member of the designed oligo pool or it
is not, and that is an exact-match question, not a similarity one. BLAST also
only scaled to the 16 family representatives; this annotates all of them.

Why protein level and not nucleotide
------------------------------------
The library's `oligo` column and the sequenced inserts do not match verbatim --
the deposited oligos are codon-optimised differently from what was sequenced.
The encoded peptides do match. So the join key is the 16-aa translation of the
50 bp insert, looked up in an index of every 16-mer in every 56-aa library
tile.

A 16-mer can occur in more than one tile: adjacent tiles overlap by 28
residues, and homologous strains share sequence. All matching tiles are kept
and summarised, rather than silently taking the first.

Outputs
-------
results/vir3_annotation.tsv.gz   one row per tested feature
results/vir3_summary.tsv         per protocol: how much of the hit list is
                                 library, and which species
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
K = 16


def translate_best_frame(nt: str) -> str | None:
    """Same rule as 05a_peptides.py: the frame with no internal stop."""
    best, best_len = None, -1
    for frame in (0, 1, 2):
        sub = nt[frame:]
        sub = sub[: len(sub) - (len(sub) % 3)]
        if len(sub) < 9:
            continue
        aa = str(Seq(sub).translate())
        core = aa.rstrip("*")
        if "*" not in core and len(core) > best_len:
            best, best_len = core, len(core)
    return best


def build_index(vir3_path: Path):
    """k-mer -> list of library row indices, plus the row metadata."""
    print(f"indexing {vir3_path} ...", file=sys.stderr)
    index: dict[str, list[int]] = defaultdict(list)
    meta: list[dict] = []
    opener = gzip.open if str(vir3_path).endswith(".gz") else open
    with opener(vir3_path, "rt", newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            pep = (row.get("peptide") or "").strip()
            if not pep:
                continue
            meta.append({
                "lib_id":   row.get("id", ""),
                "species":  row.get("Species", ""),
                "organism": row.get("Organism", ""),
                "protein":  row.get("Protein names", ""),
                "gene":     row.get("Gene names", ""),
                "entry":    row.get("Entry", ""),
                "source":   row.get("source", ""),
            })
            j = len(meta) - 1
            for s in range(0, max(1, len(pep) - K + 1)):
                index[pep[s:s + K]].append(j)
    print(f"  {len(meta):,} tiles, {len(index):,} distinct {K}-mers", file=sys.stderr)
    return index, meta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vir3", default=str(PROJECT / "vir3_annotations.csv.gz"))
    ap.add_argument("--enrichment", default=str(RESULTS / "enrichment_all.tsv.gz"))
    ap.add_argument("--hits-only", action="store_true",
                    help="annotate only enriched features (faster)")
    args = ap.parse_args()

    vir3 = Path(args.vir3)
    if not vir3.exists():
        sys.exit(
            f"VirScan annotation not found: {vir3}\n"
            "Download it with:\n"
            "  curl -L -o data.tar.gz "
            "'https://zenodo.org/records/20033894/files/data.tar.gz?download=1'\n"
            "  tar xzf data.tar.gz"
        )

    index, meta = build_index(vir3)

    print("reading enrichment table ...", file=sys.stderr)
    df = pd.read_csv(args.enrichment, sep="\t")
    if args.hits_only:
        df = df[df["enriched_edger"].astype(str) == "True"]
    print(f"  {len(df):,} features", file=sys.stderr)

    # A peptide that does not match exactly is not necessarily absent from the
    # library. A single sequencing or PCR error anywhere in the 16-mer breaks
    # an exact lookup, and the hit lists are full of single-residue variant
    # families -- the parent matches, its error children do not. Without this
    # pass, "not in the library" conflates two completely different things:
    # an error variant of a real library member, and a sequence that was never
    # in the pool. Only the second one is a finding.
    AA = "ACDEFGHIKLMNPQRSTVWY"

    def near_match(kmer: str):
        """All library tiles within one substitution of this 16-mer."""
        for pos in range(len(kmer)):
            for aa_sub in AA:
                if aa_sub == kmer[pos]:
                    continue
                alt = kmer[:pos] + aa_sub + kmer[pos + 1:]
                got = index.get(alt)
                if got:
                    return got, alt
        return [], None

    rows = []
    for r in df.itertuples(index=False):
        aa = translate_best_frame(str(r.feature))
        kmer = aa[:K] if aa and len(aa) >= K else ""
        hits = index.get(kmer, []) if kmer else []
        match_type = "exact" if hits else "none"
        if not hits and kmer:
            hits, _alt = near_match(kmer)
            if hits:
                match_type = "1_mismatch"
        species = sorted({meta[j]["species"] for j in hits if meta[j]["species"]})
        proteins = sorted({meta[j]["protein"] for j in hits if meta[j]["protein"]})
        rows.append({
            "feature":        r.feature,
            "peptide":        aa or "",
            "protocol":       r.protocol,
            "logFC":          round(float(r.logFC), 3),
            "FDR":            float(r.FDR),
            "enriched_edger": r.enriched_edger,
            "enriched_paper": r.enriched_paper,
            "in_library":     bool(hits),
            "match_type":     match_type,
            "n_tiles":        len(hits),
            "n_species":      len(species),
            "species":        "; ".join(species[:5]),
            "protein":        "; ".join(proteins[:3]),
            "lib_source":     "; ".join(sorted({meta[j]["source"] for j in hits})),
        })

    out = pd.DataFrame(rows)
    RESULTS.mkdir(exist_ok=True)
    out.to_csv(RESULTS / "vir3_annotation.tsv.gz", sep="\t", index=False,
               compression="gzip")

    # ---- summary -----------------------------------------------------------
    hits = out[out.enriched_edger.astype(str) == "True"]
    summ = (
        hits.groupby("protocol")
            .agg(enriched=("feature", "size"),
                 exact=("match_type", lambda s: (s == "exact").sum()),
                 one_mismatch=("match_type", lambda s: (s == "1_mismatch").sum()),
                 absent=("match_type", lambda s: (s == "none").sum()))
            .assign(pct_library=lambda d:
                    (100 * (d.exact + d.one_mismatch) / d.enriched).round(1))
            .reset_index()
    )
    summ.to_csv(RESULTS / "vir3_summary.tsv", sep="\t", index=False)

    print("\n=== enriched peptides that are genuine library members ===")
    print(summ.to_string(index=False))

    print("\n=== species of the enriched library members ===")
    for prot, g in hits[hits.in_library].groupby("protocol"):
        top = (g.species.str.split("; ").str[0]
                 .value_counts().head(8))
        print(f"\n{prot}:")
        for sp, n in top.items():
            print(f"  {n:>4}  {sp}")

    n_abs = (hits.match_type == "none").sum()
    n_mm = (hits.match_type == "1_mismatch").sum()
    print(f"\n{n_mm} enriched peptides match a library tile with ONE substitution"
          " - sequencing or PCR error variants of a real member, which de novo"
          " counting split off as separate features.")
    if n_abs:
        print(f"{n_abs} match nothing even allowing one substitution. "
              "Those cannot be antigen hits whatever their fold change.")
        worst = (hits[hits.match_type == "none"]
                 .nlargest(5, "logFC")[["protocol", "peptide", "logFC"]])
        print("\nhighest fold change among them:")
        print(worst.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
