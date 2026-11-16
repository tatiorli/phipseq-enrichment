#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Where the library annotation and the BLAST inference agree, and where they
# do not.
#
#   Rscript scripts/09_vir3_figure.R
#
# BLAST asks what a peptide resembles; the library says what it is. Running
# both and comparing them is the honest way to show how much of a BLAST-based
# annotation you should trust when the reference exists.
#
# Panel A  how much of each protocol's hit list is a genuine library member.
#          A peptide that is not in the designed pool cannot be an antigen hit
#          whatever its fold change.
# Panel B  which viruses the library-confirmed hits come from, per protocol.
# Panel C  the 16 BLAST-annotated families, side by side with what the library
#          says they are: agree / disagree / BLAST found nothing / not a
#          library member.
#
# Inputs : results/vir3_annotation.tsv.gz  (from 08_annotate_vir3.py)
#          results/blast_annotation.tsv
# Output : results/figures/vir3_vs_blast.{pdf,png}
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(patchwork); library(scales)
})

proj <- Sys.getenv("PROJECT", unset = "")
if (!nzchar(proj)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  proj <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
}
res_dir <- file.path(proj, "results")
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)

ann_file <- file.path(res_dir, "vir3_annotation.tsv.gz")
if (!file.exists(ann_file))
  stop("run scripts/08_annotate_vir3.py first - missing ", ann_file)

vir <- read_tsv(ann_file, show_col_types = FALSE) %>%
  mutate(protocol_lab = recode(protocol,
                               DL = "D_long", DLPW = "D_long + pre-wash",
                               S2 = "Sergio2", EL = "Elledge_long",
                               .default = protocol))
hits <- vir %>% filter(enriched_edger %in% c(TRUE, "True", "TRUE"))
if (!nrow(hits)) stop("no enriched features in ", ann_file)

# --- Panel A: library membership of each hit list --------------------------
memb <- hits %>%
  mutate(cls = if_else(in_library %in% c(TRUE, "True", "TRUE"),
                       "in the VirScan library", "not a library member")) %>%
  count(protocol_lab, cls)

pA <- ggplot(memb, aes(n, protocol_lab, fill = cls)) +
  geom_col(width = .55) +
  geom_text(aes(label = n), position = position_stack(vjust = .5),
            colour = "white", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c("in the VirScan library" = "#2E7D32",
                               "not a library member"   = "#B71C1C"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .05))) +
  labs(title = "A  Are the enriched peptides actually in the designed library?",
       subtitle = "exact 16-mer match against all 128,257 Vir3 tiles",
       x = "enriched peptides", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

# --- Panel B: which viruses ------------------------------------------------
sp <- hits %>%
  filter(in_library %in% c(TRUE, "True", "TRUE"), !is.na(species), species != "") %>%
  mutate(sp1 = str_trim(str_split_fixed(species, ";", 2)[, 1])) %>%
  count(protocol_lab, sp1, sort = TRUE) %>%
  group_by(protocol_lab) %>% slice_max(n, n = 10) %>% ungroup()

pB <- ggplot(sp, aes(n, reorder(sp1, n), fill = protocol_lab)) +
  geom_col(width = .7, show.legend = FALSE) +
  facet_wrap(~ protocol_lab, scales = "free", ncol = 2) +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "B  Which organisms the library-confirmed hits come from",
       x = "enriched peptides", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"))

# --- Panel C: BLAST vs library, per annotated family -----------------------
pC <- NULL
bl_file <- file.path(res_dir, "blast_annotation.tsv")
if (file.exists(bl_file)) {
  bl <- read_tsv(bl_file, comment = "#", show_col_types = FALSE)

  lib_for_pep <- vir %>%
    select(peptide, in_library, species) %>%
    distinct(peptide, .keep_all = TRUE)

  cmp <- bl %>%
    left_join(lib_for_pep, by = c("peptide" = "peptide")) %>%
    mutate(
      in_library = in_library %in% c(TRUE, "True", "TRUE"),
      blast_found = !str_detect(class, "UNASSIGNED"),
      lib_sp = str_trim(str_split_fixed(coalesce(species, ""), ";", 2)[, 1]),
      verdict = case_when(
        !in_library &  blast_found ~ "BLAST hit, not a library member",
        !in_library & !blast_found ~ "no BLAST hit, not a library member",
        in_library  & !blast_found ~ "library member BLAST missed",
        TRUE                       ~ "both agree"
      ),
      label = paste0(str_trunc(coalesce(na_if(organism, "-"), "no BLAST hit"), 36),
                     "  [", family, "]")
    )

  pC <- ggplot(cmp, aes(log2FC, reorder(label, log2FC), colour = verdict)) +
    geom_segment(aes(x = 0, xend = log2FC, yend = reorder(label, log2FC)),
                 linewidth = .4, alpha = .4) +
    geom_point(size = 3.4) +
    scale_colour_manual(values = c(
      "both agree"                         = "#2E7D32",
      "library member BLAST missed"        = "#1565C0",
      "BLAST hit, not a library member"    = "#EF6C00",
      "no BLAST hit, not a library member" = "#B71C1C"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, .08))) +
    labs(title = "C  BLAST inference vs library membership, per family",
         subtitle = "BLAST says what a peptide resembles; the library says what it is",
         x = expression(log[2]~"fold change"), y = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right", panel.grid.major.y = element_blank())

  write_tsv(cmp %>% select(family, peptide, log2FC, organism, lib_sp,
                           in_library, blast_found, verdict),
            file.path(res_dir, "blast_vs_vir3.tsv"))
}

fig <- if (is.null(pC)) pA / pB + plot_layout(heights = c(1, 3)) else
  pA / pB / pC + plot_layout(heights = c(1, 2.6, 2.2))

fig <- fig + plot_annotation(
  title = "VirScan library annotation of the enriched peptides",
  subtitle = "Vir3 annotations, Zenodo 10.5281/zenodo.20033894 (CC BY 4.0)",
  caption = paste("generated", Sys.Date()),
  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(fig_dir, "vir3_vs_blast.pdf"), fig, width = 12, height = 14)
ggsave(file.path(fig_dir, "vir3_vs_blast.png"), fig, width = 12, height = 14, dpi = 200)
message("wrote ", file.path(fig_dir, "vir3_vs_blast.pdf"))

print(as.data.frame(memb))
