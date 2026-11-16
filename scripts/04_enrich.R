#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Differential enrichment: for each protocol, serum (P1) vs its own bead
# control (BC), using edgeR.
#
#   Rscript scripts/04_enrich.R
#
# Why edgeR and not the published thresholds
# ------------------------------------------
# Kjendal et al. call a peptide enriched when CV < 100, fold change over bead
# control >= 2, and normalised count >= 10. That treats a peptide seen 3 vs 6
# times the same as one seen 3,000 vs 6,000, and it has no multiple-testing
# control over ~10^5-10^6 features. edgeR models the count noise explicitly
# (negative binomial, dispersion estimated across features) and returns an FDR.
# We compute BOTH so the figure can show what the thresholds miss.
#
# Each protocol is tested separately against its OWN bead control. Protocols
# were run on different days and arms differ in depth, so a single model across
# all of them would confound protocol with batch.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR); library(readr); library(dplyr); library(tidyr); library(data.table)
})

# Project root: $PROJECT if exported by config/params.sh, else the parent of
# this script's directory, else the working directory.
proj <- Sys.getenv("PROJECT", unset = "")
if (!nzchar(proj)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  proj <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
}
res_dir <- file.path(proj, "results")
man     <- read_tsv(file.path(proj, "config", "manifest.tsv"), show_col_types = FALSE)

FDR_CUTOFF   <- as.numeric(Sys.getenv("FDR_CUTOFF",   "0.05"))
LOGFC_CUTOFF <- as.numeric(Sys.getenv("LOGFC_CUTOFF", "1"))

message("reading count matrix")
# fread() cannot open .gz directly without the R.utils package. Decompress
# through a pipe instead of adding a dependency just to read one file.
mat_file <- file.path(res_dir, "counts_matrix.tsv.gz")
dt <- data.table::fread(cmd = paste("gzip -dc", shQuote(mat_file)),
                        sep = "\t", header = TRUE)
rn <- dt[[1]]
counts <- as.matrix(dt[, -1, with = FALSE])
rownames(counts) <- rn
storage.mode(counts) <- "integer"
message(sprintf("  %d features x %d samples", nrow(counts), ncol(counts)))

meta <- man %>% filter(sample_name %in% colnames(counts))
counts <- counts[, meta$sample_name, drop = FALSE]

all_res <- list()

for (prot in unique(meta$protocol)) {
  m <- meta %>% filter(protocol == prot)
  if (!all(c("P1", "BC") %in% m$arm)) {
    message(sprintf("skip %s: needs both P1 and BC arms", prot)); next
  }
  message(sprintf("\n=== %s : %d samples ===", prot, nrow(m)))

  sub <- counts[, m$sample_name, drop = FALSE]
  grp <- factor(m$arm, levels = c("BC", "P1"))   # BC is the reference

  # Drop features absent from this protocol entirely
  keep <- filterByExpr(sub, group = grp, min.count = 5, min.total.count = 15)
  message(sprintf("  features tested: %d of %d", sum(keep), length(keep)))
  sub <- sub[keep, , drop = FALSE]

  y <- DGEList(counts = sub, group = grp)
  y <- calcNormFactors(y, method = "TMM")      # depth differs 2-3x between arms
  design <- model.matrix(~ grp)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)             # P1 vs BC

  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  cpms <- cpm(y, normalized.lib.sizes = TRUE)

  df <- tibble(
    feature    = rownames(tt),
    protocol   = prot,
    logFC      = tt$logFC,
    logCPM     = tt$logCPM,
    PValue     = tt$PValue,
    FDR        = tt$FDR,
    cpm_P1     = rowMeans(cpms[, m$arm == "P1", drop = FALSE]),
    cpm_BC     = rowMeans(cpms[, m$arm == "BC", drop = FALSE])
  ) %>%
    mutate(
      enriched_edger = FDR < FDR_CUTOFF & logFC > LOGFC_CUTOFF,
      # The published rule, reimplemented for side-by-side comparison:
      fc_raw         = (cpm_P1 + 1) / (cpm_BC + 1),
      enriched_paper = fc_raw >= 2 & cpm_P1 >= 10
    )

  n_e <- sum(df$enriched_edger); n_p <- sum(df$enriched_paper)
  message(sprintf("  enriched  edgeR: %d   threshold-rule: %d   both: %d",
                  n_e, n_p, sum(df$enriched_edger & df$enriched_paper)))

  # Per-protocol dispersion is itself a QC number worth reporting
  message(sprintf("  common BCV: %.3f", sqrt(y$common.dispersion)))
  all_res[[prot]] <- df
}

if (!length(all_res)) stop("No protocol had both arms; check config/samples.tsv")

res <- bind_rows(all_res)
write_tsv(res, file.path(res_dir, "enrichment_all.tsv.gz"))

hits <- res %>% filter(enriched_edger) %>% arrange(protocol, desc(logFC))
write_tsv(hits, file.path(res_dir, "enriched_peptides.tsv"))

# NOTE the column names. They must NOT match the logical columns they are
# computed from. summarise() makes each new column visible to the expressions
# after it, so `enriched_edger = sum(enriched_edger)` replaces the logical
# vector with a scalar, and the next line then evaluates
# `sum(!0 & 622)` = sum(TRUE) = 1 instead of 622. That is exactly what the
# first run reported: only_paper = 1 when it should have been 622.
safe_med <- function(x) if (length(x)) round(stats::median(x, na.rm = TRUE), 2) else NA_real_
safe_max <- function(x) if (length(x)) round(max(x, na.rm = TRUE), 2) else NA_real_

summary_tbl <- res %>%
  group_by(protocol) %>%
  summarise(
    tested       = n(),
    n_edger      = sum(enriched_edger, na.rm = TRUE),
    n_paper      = sum(enriched_paper, na.rm = TRUE),
    n_only_edger = sum(enriched_edger & !enriched_paper, na.rm = TRUE),
    n_only_paper = sum(!enriched_edger & enriched_paper, na.rm = TRUE),
    n_both       = sum(enriched_edger & enriched_paper, na.rm = TRUE),
    median_logFC = safe_med(logFC[enriched_edger]),
    max_logFC    = safe_max(logFC[enriched_edger]),
    min_FDR      = signif(min(FDR, na.rm = TRUE), 3),
    .groups = "drop"
  )
write_tsv(summary_tbl, file.path(res_dir, "enrichment_summary.tsv"))

message("\n=== summary ===")
print(as.data.frame(summary_tbl))
message(sprintf("\nwrote %s", file.path(res_dir, "enriched_peptides.tsv")))
