# PhIP-seq protocol comparison — PRJNA1438256

Which peptides does each VirScan protocol variant actually enrich, how strongly,
and do the enriched peptides share a sequence motif?

Data: Kjendal et al., *Front. Virol.* 6:1761741 (2026),
[doi:10.3389/fviro.2026.1761741](https://doi.org/10.3389/fviro.2026.1761741),
BioProject **PRJNA1438256**.

---

## The question, and why only 16 of the 48 runs

The BioProject holds 48 runs (24 conditions × 2 replicates, 24.7 GB). We do not
need all of them.

The same serum — **patient 1** — was run through every protocol variant, and
each protocol has **its own bead-only control**. That gives a clean comparison:
serum vs matched control, within protocol, in duplicate.

| protocol code | `experiment` tag  | serum arm      | bead control    | runs |
|---------------|-------------------|----------------|-----------------|------|
| `EL`          | `Elledge_long`    | `EL_P1_1/2`    | `EL_BC_1/2`     | 4    |
| `S2`          | `Sergio2`         | `S2_P1_1/2`    | `S2_BC_1/2`     | 4    |
| `DL`          | `D_long_r1/r3`    | `DL3_P1_1/2`   | `DL1_BC_1/2`    | 4    |
| `DLPW`        | `D_long_*_PW`     | `DL3_P1_PW_1/2`| `DL1_BC_PW_1/2` | 4    |
| `RD`          | `round1` / `round3` | `r3_P1_1/2`  | `r1_BC_1/2`     | 4    |

**16 runs, ~4.6 GB** with `RD` excluded (the default). `DL` vs `DLPW` is the
pre-wash contrast on its own — the paper's central claim, in 8 files.

### Why `RD` is excluded by default

`EL`, `S2`, `DL` and `DLPW` were sequenced at **50 bp**. The `RD` arm was
sequenced at **151 bp**. Protocol is therefore confounded with read length and
sequencing run, and in de novo mode the longer reads produce an entirely
disjoint set of exact sequences. The pipeline truncates everything to
`TRUNC_LEN=50` precisely so this comparison is possible — but it is still a
different sequencing batch, so treat any `RD` difference as protocol *plus*
batch. To include it: set `include=yes` for those four rows in
`config/samples.tsv`.

---

## Where the sample labels came from

They are **not in the paper**, and they are **not in ENA's run report** — both
come back blank. They are in the NCBI **BioSample** records, as custom
attributes (`sample name`, `info`, `experiment`, `fastq_file`):

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=$(seq -s, 56533536 56533583)&rettype=full&retmode=xml"
```

`config/samples.tsv` is that mapping, already joined to run accessions.

---

## Running it

Build and push the image once (see
[`docker/README-wexac.md`](docker/README-wexac.md)), then:

```bash
bash scripts/00_smoke_test.sh
```

That submits one small interactive job and checks the things that otherwise
fail three hours into a run: that the job really landed inside the image, that
this directory is visible *and writable* from in there, that every tool is on
PATH, and that the Python and R scripts parse. It costs seconds.

Once it passes:

```bash
bash run_all.sh
```

Stage by stage:

| stage | script | what it does |
|---|---|---|
| 0 | `scripts/00_manifest.sh` | asks ENA for current FASTQ URLs, builds `config/manifest.tsv` |
| 1 | `scripts/01_download.lsf` | LSF array, one run per index, md5/line-count verified, `.done` markers |
| 2 | `scripts/02_count.lsf` | fastp → truncate to 50 bp → count inserts (de novo or bowtie2) |
| 3 | `scripts/03_matrix.py` | two-pass feature × sample matrix + library QC |
| 4 | `scripts/04_enrich.R` | edgeR QL F-test, serum vs matched bead control, per protocol |
| 5 | `scripts/05_motifs.sh` | translate hits, STREME against a non-enriched background |
| 6 | `scripts/06_figure.R` | the figure |

Everything is idempotent. Re-run `run_all.sh` after any failure; finished work
is skipped.

### Wexac notes

- **Containers: Docker works here.** Wexac runs IBM LSF with its Docker
  application profiles, so a job can be launched *inside* a Docker image with
  `-app docker-cpu -env LSB_CONTAINER_IMAGE=<image>`. You build on the
  dedicated server `docker-gpu.wexac.weizmann.ac.il` (the only place you have
  the privileges) and push to the internal registry `ops:5000`. Singularity is
  the alternative, not the requirement. Full instructions, both routes, in
  [`docker/README-wexac.md`](docker/README-wexac.md).
  Set `CONTAINER_MODE` in `config/params.sh` to `docker-lsf` (default),
  `singularity`, or `none`.
- **If you use Singularity, load the module first.** `singularity` is not on
  PATH until `ml Singularity`. The pipeline does this for you, but remember it
  when running commands by hand.
- **Memory is always declared.** Wexac reserves 1024 MB by default and kills
  anything that exceeds it, usually without an obvious error. Every `bsub` here
  carries `rusage[mem=]`. Do not remove them.
- **Queue names** are set in `config/params.sh` (`short`, with
  `interactive`/`short-gpu` used in the container docs). Confirm with
  `bqueues` — your lab may have dedicated queues.

---

## The analysis, and the argument behind it

The published analysis calls a peptide enriched when CV < 100, fold change over
bead control ≥ 2, and normalised count ≥ 10. Three hard thresholds, no error
model, no multiple-testing control across ~10⁵–10⁶ features.

`04_enrich.R` runs **edgeR** (TMM normalisation, negative-binomial dispersion,
quasi-likelihood F-test) against each protocol's own bead control — *and also
reimplements the published threshold rule*, so `enrichment_summary.tsv` reports
`only_edger` and `only_paper` side by side. The disagreement between the two is
the finding, not a footnote.

Each protocol is modelled separately. A single model across protocols would
confound protocol with sequencing batch and depth.

### About the sequence logo

A logo only means something for aligned, equal-length sequences. VirScan
peptides are 56-aa tiles from hundreds of different viral proteins — stacking
them yields a smear. So `05_motifs.sh` asks the right question instead: what
short motif is **over-represented in the enriched set relative to peptides from
the same library that were tested and not enriched**. STREME answers that, and
panel C renders its PWM as the logo. The background is drawn from the library
itself, so the library's own amino-acid composition is controlled for.

If no motif is found, that is a real result and the panel says so: the enriched
peptides may be genuine binders with no shared *linear* motif.

---

## Known unknowns

Two things to settle, ideally with one email to the corresponding author
(dev.majumdar@uvm.edu):

1. **Which `experiment` tag is which named protocol.** `Elledge_long`,
   `Sergio2`, `D_long` are internal lab names; the paper names Original,
   Phage-Ab Optimized, Multi-IP, Multi-IP + Pre-wash, Optimized. `_PW` =
   pre-wash and `round1`/`round3` = IP round are unambiguous; the rest is
   inference. Edit the `lab` vector in `scripts/06_figure.R` once confirmed.

2. **The peptide reference.** The ~115,000-peptide VirScan library is *not*
   deposited with the BioProject. Without it you cannot name epitopes or
   viruses. `COUNT_MODE=denovo` (the default) sidesteps this: it counts exact
   inserts, so the pipeline runs today and you annotate the handful of enriched
   sequences afterwards with blastp. Set `COUNT_MODE=reference` and point
   `REF_FASTA` at the library when you have it.

Also worth verifying: 16 of the 48 BioSamples list a `.count.csv.gz` as their
source file rather than a FASTQ. ENA reports non-zero read counts for all of
them, so the reads exist, but check those runs look like reads before trusting
them. All four `RD` runs are in that group.

---

## Output

```
results/
├── counts_matrix.parquet / .tsv.gz
├── library_qc.tsv                  depth, diversity, Shannon entropy per sample
├── enrichment_all.tsv.gz           every tested feature, both calling methods
├── enriched_peptides.tsv           the hits
├── enrichment_summary.tsv          per protocol: edgeR vs threshold rule
├── motifs/<protocol>.streme/       STREME output
└── figures/
    ├── protocol_comparison.pdf     panels A–D
    └── motif_logos.pdf             logos alone
```
