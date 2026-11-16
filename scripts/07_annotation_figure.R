#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Visualise what the enriched peptides actually are.
#
#   Rscript scripts/07_annotation_figure.R
#
# The statistics produce a ranked list of sequences. This figure answers the
# question that matters next: which organism does each one come from, and how
# much of the hit list is real virome signal versus assay artifact.
#
# Panel A  every enriched peptide family, ranked by fold change, labelled with
#          the organism BLAST assigned it to, coloured by call, sized by how
#          many near-identical variants collapse into it.
# Panel B  the composition of each protocol's hit list by call. This is the
#          panel that shows one streptococcal IgG-binding family accounting
#          for most of D_long's top hits.
#
# Inputs : results/blast_annotation.tsv   (curated from the two BLAST runs)
#          results/motifs/{DL,S2}.enriched.faa
# Output : results/figures/blast_annotation.{pdf,png}
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

ann_file <- file.path(res_dir, "blast_annotation.tsv")
if (!file.exists(ann_file)) stop("missing ", ann_file)
ann <- read_tsv(ann_file, comment = "#", show_col_types = FALSE)

# --- read the enriched peptides, with their fold changes from the headers ---
read_faa <- function(path, protocol) {
  if (!file.exists(path)) return(tibble())
  ln  <- readLines(path, warn = FALSE)
  hdr <- grep("^>", ln)
  tibble(
    protocol = protocol,
    name     = sub("^>", "", ln[hdr]),
    peptide  = ln[hdr + 1]
  ) %>%
    mutate(logFC = as.numeric(str_match(name, "lfc(-?[0-9.]+)")[, 2])) %>%
    filter(!is.na(peptide), peptide != "")
}

pep <- bind_rows(
  read_faa(file.path(res_dir, "motifs", "DL.enriched.faa"), "DL"),
  read_faa(file.path(res_dir, "motifs", "S2.enriched.faa"), "S2")
)
if (!nrow(pep)) stop("no enriched peptides found under results/motifs/")

# --- assign each peptide to its nearest annotated representative ------------
# The hit lists are dominated by families of single-residue variants, which de
# novo counting splits into separate features. Collapsing them is the whole
# point: 251 "hits" is not 251 discoveries. Edit distance is enough here --
# within a family the peptides differ at one or two positions.
MAX_DIST <- 6          # of 16 residues

assign_family <- function(p, prot) {
  cand <- ann %>% filter(str_starts(family, prot))
  if (!nrow(cand)) return(NA_character_)
  d <- as.vector(adist(p, cand$peptide))
  if (min(d) > MAX_DIST) return(NA_character_)
  cand$family[which.min(d)]
}

pep$family <- mapply(assign_family, pep$peptide, pep$protocol)

fam_size <- pep %>% count(family, name = "n_variants") %>% filter(!is.na(family))

dat <- ann %>%
  left_join(fam_size, by = "family") %>%
  mutate(
    n_variants = coalesce(n_variants, 1L),
    protocol   = if_else(str_starts(family, "DL"), "D_long", "Sergio2"),
    # Shorten the call into a plotting category
    call_grp = case_when(
      str_detect(class, "ASSAY ARTIFACT")  ~ "assay artifact (binds IgG)",
      str_detect(class, "BACTERIAL")       ~ "bacterial, not viral",
      str_detect(class, "UNASSIGNED")      ~ "no database match",
      str_detect(class, "weak")            ~ "viral epitope (weak)",
      TRUE                                 ~ "viral epitope"
    ),
    label = if_else(organism == "-" | is.na(organism),
                    paste0(family, "  (no match)"),
                    paste0(str_trunc(organism, 42), "  [", family, "]"))
  )

pal <- c(
  "viral epitope"               = "#2E7D32",
  "viral epitope (weak)"        = "#8BC34A",
  "assay artifact (binds IgG)"  = "#C62828",
  "bacterial, not viral"        = "#EF6C00",
  "no database match"           = "#9E9E9E",
  "not in a BLASTed family"     = "#CFD8DC"
)

# --- Panel A: what each family is ------------------------------------------
pA <- dat %>%
  mutate(label = reorder(label, log2FC)) %>%
  ggplot(aes(log2FC, label, colour = call_grp)) +
  geom_segment(aes(x = 0, xend = log2FC, yend = label), linewidth = .5, alpha = .45) +
  geom_point(aes(size = n_variants)) +
  geom_text(aes(label = ifelse(n_variants > 1, n_variants, "")),
            colour = "white", size = 2.6, fontface = "bold", show.legend = FALSE) +
  facet_grid(protocol ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_size_continuous(range = c(3, 9), name = "variants\nin family") +
  scale_x_continuous(expand = expansion(mult = c(0, .08))) +
  labs(
    title = "A  What the enriched peptides are",
    subtitle = "each point is a family of near-identical variants; number inside = how many sequences collapse into it",
    x = expression(log[2]~"fold change over matched bead control"), y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        panel.grid.major.y = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"))

# --- Panel B: composition of each hit list ---------------------------------
# A peptide too far from every annotated representative was never BLASTed at
# all -- only the 16 family representatives were. Calling that "no database
# match" would claim a negative result that was never tested, which is a
# different and much stronger statement. Label it honestly.
comp <- pep %>%
  left_join(dat %>% select(family, call_grp), by = "family") %>%
  mutate(call_grp = coalesce(call_grp, "not in a BLASTed family"),
         protocol = if_else(protocol == "DL", "D_long", "Sergio2")) %>%
  count(protocol, call_grp)

pB <- ggplot(comp, aes(n, protocol, fill = call_grp)) +
  geom_col(width = .55) +
  geom_text(aes(label = n), position = position_stack(vjust = .5),
            colour = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = pal, name = NULL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .05))) +
  labs(title = "B  Composition of each protocol's hit list",
       subtitle = "peptides assigned to an annotated family by edit distance",
       x = "enriched peptides", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

fig <- pA / pB +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "Identity of the enriched peptides, PRJNA1438256 reanalysis",
    subtitle = "blastp against UniProtKB/Swiss-Prot and nr, short-query parameters (word size 2, PAM30, expect 1000)",
    caption = paste("generated", Sys.Date()),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "blast_annotation.pdf"), fig, width = 12, height = 9)
ggsave(file.path(fig_dir, "blast_annotation.png"), fig, width = 12, height = 9, dpi = 200)
message("wrote ", file.path(fig_dir, "blast_annotation.pdf"))

# A small text summary, so the numbers quoted anywhere can be checked
comp %>%
  group_by(protocol) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(protocol, desc(n)) %>%
  write_tsv(file.path(res_dir, "hitlist_composition.tsv"))
print(as.data.frame(comp))
