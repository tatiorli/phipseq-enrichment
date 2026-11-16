#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Final figure: what each protocol enriched, and what the enriched peptides
# look like.
#
#   Rscript scripts/06_figure.R
#
# Panels
#   A  how many peptides each protocol enriched (edgeR vs the published rule)
#   B  how strongly  -- distribution of log2 fold change over bead control
#   C  sequence logo of the top motif STREME found in each protocol's hits
#   D  the top peptides, and which protocols recover them
#
# Output: results/figures/protocol_comparison.{pdf,png}
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

proj <- Sys.getenv("PROJECT", unset = "")
if (!nzchar(proj)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  proj <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
}
res_dir <- file.path(proj, "results")
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)

res     <- read_tsv(file.path(res_dir, "enrichment_all.tsv.gz"), show_col_types = FALSE)
hits    <- read_tsv(file.path(res_dir, "enriched_peptides.tsv"), show_col_types = FALSE)
summ    <- read_tsv(file.path(res_dir, "enrichment_summary.tsv"), show_col_types = FALSE)

# Readable protocol labels. Edit once the author confirms the mapping to the
# five protocols named in the paper (see README, "Known unknowns").
lab <- c(EL = "Elledge_long", S2 = "Sergio2", DL = "D_long", DLPW = "D_long + pre-wash",
         RD = "round1/3")
pal <- c(EL = "#4C72B0", S2 = "#DD8452", DL = "#55A868", DLPW = "#C44E52", RD = "#8172B3")
relabel <- function(x) ifelse(x %in% names(lab), lab[x], x)

ord <- summ %>% arrange(desc(n_edger)) %>% pull(protocol)
res$protocol  <- factor(res$protocol,  levels = ord)
hits$protocol <- factor(hits$protocol, levels = ord)
summ$protocol <- factor(summ$protocol, levels = ord)

# ---- A: counts -------------------------------------------------------------
cnt <- summ %>%
  select(protocol, edgeR = n_edger, `threshold rule` = n_paper) %>%
  pivot_longer(-protocol, names_to = "method", values_to = "n")

pA <- ggplot(cnt, aes(protocol, n, fill = method)) +
  geom_col(position = position_dodge(width = .75), width = .7) +
  geom_text(aes(label = comma(n)), position = position_dodge(width = .75),
            vjust = -.35, size = 3) +
  scale_x_discrete(labels = relabel) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) +
  scale_fill_manual(values = c(edgeR = "#2B5D8C", `threshold rule` = "#BDBDBD")) +
  labs(title = "A  Peptides called enriched over the matched bead control",
       subtitle = sprintf("edgeR: FDR < %s and log2FC > %s   |   threshold rule: FC >= 2 and CPM >= 10",
                          Sys.getenv("FDR_CUTOFF", "0.05"), Sys.getenv("LOGFC_CUTOFF", "1")),
       x = NULL, y = "peptides", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.x = element_blank())

# ---- B: magnitude ----------------------------------------------------------
pB <- ggplot(hits, aes(protocol, logFC, fill = protocol)) +
  geom_violin(scale = "width", colour = NA, alpha = .55, trim = TRUE) +
  geom_boxplot(width = .13, outlier.shape = NA, fill = "white", colour = "grey25") +
  scale_x_discrete(labels = relabel) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "B  How strongly they were enriched",
       x = NULL, y = expression(log[2]~"fold change vs bead control")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

# ---- C: motif logos --------------------------------------------------------
logo_panel <- function() {
  ok <- requireNamespace("universalmotif", quietly = TRUE)
  if (!ok) return(NULL)
  plots <- list()
  for (p in levels(hits$protocol)) {
    f <- file.path(res_dir, "motifs", paste0(p, ".streme"), "streme.txt")
    if (!file.exists(f)) next
    m <- try(universalmotif::read_meme(f), silent = TRUE)
    if (inherits(m, "try-error") || !length(m)) next
    if (!is.list(m)) m <- list(m)
    top <- m[[1]]
    g <- universalmotif::view_motifs(top, use.type = "ICM") +
      ggtitle(sprintf("%s  (E = %s)", relabel(p),
                      format(top@eval %||% NA, digits = 2))) +
      theme(plot.title = element_text(size = 9))
    plots[[p]] <- g
  }
  if (!length(plots)) return(NULL)
  wrap_plots(plots, ncol = min(length(plots), 2))
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

pC <- logo_panel()
if (is.null(pC)) {
  pC <- ggplot() + theme_void() +
    annotate("text", 0, 0, size = 3.4, colour = "grey30",
             label = paste("No STREME motif found.",
                           "Either too few enriched peptides, or the enriched set",
                           "shares no short linear motif - which is itself a result.",
                           sep = "\n"))
}
pC <- pC + plot_annotation(title = NULL)

# ---- D: shared top peptides -----------------------------------------------
top_feats <- hits %>% group_by(protocol) %>% slice_max(logFC, n = 15) %>%
  ungroup() %>% distinct(feature) %>% pull(feature)

mat <- res %>% filter(feature %in% top_feats) %>%
  mutate(short = paste0(substr(feature, 1, 18), "..")) %>%
  select(short, protocol, logFC)

pD <- ggplot(mat, aes(protocol, reorder(short, logFC), fill = logFC)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_x_discrete(labels = relabel) +
  scale_fill_gradient2(low = "#3B6EA5", mid = "grey95", high = "#B2182B",
                       midpoint = 0, name = expression(log[2]~FC)) +
  labs(title = "D  Top enriched peptides across protocols",
       subtitle = "a peptide recovered by several protocols is a stronger candidate than a one-protocol hit",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(family = "mono", size = 6.5),
        panel.grid = element_blank())

# ---- assemble --------------------------------------------------------------
design <- "
AB
CC
DD
"
fig <- (pA + pB + pC + pD) +
  plot_layout(design = design, heights = c(1, 1.1, 1.5)) +
  plot_annotation(
    title = "PhIP-seq protocol comparison: enrichment over matched bead controls",
    subtitle = "BioProject PRJNA1438256 | patient-1 serum run through each protocol, duplicates, edgeR quasi-likelihood F-test",
    caption = sprintf("generated %s", format(Sys.Date())),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "protocol_comparison.pdf"), fig, width = 12, height = 14)
ggsave(file.path(fig_dir, "protocol_comparison.png"), fig, width = 12, height = 14, dpi = 200)

message("wrote ", file.path(fig_dir, "protocol_comparison.pdf"))

# Standalone logo figure, since that is the panel people ask for separately
if (!is.null(logo_panel())) {
  ggsave(file.path(fig_dir, "motif_logos.pdf"), logo_panel(), width = 9, height = 6)
  message("wrote ", file.path(fig_dir, "motif_logos.pdf"))
}
