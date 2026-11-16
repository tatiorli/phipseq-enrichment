# BLAST identification of the enriched peptide families

Raw NCBI web BLAST reports, kept so the annotation in
[`../blast_annotation.tsv`](../blast_annotation.tsv) can be checked by hand
rather than taken on trust.

| file | database | RID | run |
|---|---|---|---|
| `blastp_swissprot_FNMV906Z016.txt` | UniProtKB/Swiss-Prot | FNMV906Z016 | 2026-11-16 |
| `blastp_nr_FNMUSCJ7016.txt` | nr (non-redundant protein) | FNMUSCJ7016 | 2026-11-16 |

**Query file:** [`../../blast_check/representatives.faa`](../../blast_check/representatives.faa)
,  16 sequences, one per variant family, taken from
`results/motifs/{DL,S2}.enriched.faa`. Each is the 16-aa translation of an
enriched 50 bp insert.

**Parameters.** blastp with short-query adjustment: word size 2, PAM30, no
compositional adjustment, expect threshold 1000, 100 target sequences.

The expect threshold matters. A perfect 16-residue match scores E ≈ 1–100
purely because the query is short, so at the default 0.05 every one of these
queries returns "No significant similarity found", including the ones that
turn out to be exact matches to respiratory syncytial virus and rhinovirus.
Short peptides are the classic case where the default cutoff hides the answer.

**Reading them.** Ignore the E-value ranking in isolation; look at alignment
length and percent identity (you want 14–16/16), then at the organism. The
calls made from these reports are in `../blast_annotation.tsv`, including the
two families that match nothing in either database and the streptococcal
IgG-binding protein that is an assay artifact rather than an epitope.

