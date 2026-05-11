#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure 3 – chrX promoter/housekeeping genes - MCCS analysis
#
# Non-exclusive macro category definitions:
#   All_active       : overlaps any active gene promoter (CPM ≥ 5 in S2 RNA-seq, all WT replicates)
#   Housekeeping     : overlaps a housekeeping gene promoter (from Jayakrishnan et al., 2025)
#   Active_regulated : overlaps an active (non-HK) promoter only (setdiff)
#   Neither          : no overlap with active, HK, or active-regulated promoters
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
out_dir  <- "../../results/figure3g_chrX_polished"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msl2_bed  <- "../../results/csaw/merged_top_sites_QLF_msl2.bed"
mre_bed   <- "../../results/motifs/mre_like_motifs_from_fimo.bed"
beaf_bed  <- "../../results/motifs/beaf_like_motifs_from_fimo.bed"

active_gene_list_file <- "../../data/input/drosophila_active_genes.txt"
hk_gene_list_file     <- "../../data/input/drosophila_housekeeping_genes_jayakrishnan2025.txt"

prom_up   <- 1000L
prom_down <- 0L

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Colour palettes
motif_cols <- c(MRE = "#4393c3", BEAF = "#d6604d", Both = "#762a83")

macro_cols <- c(
  All_active       = "#276419",
  Housekeeping     = "#7fbc41",
  Active_regulated = "#d9ef8b",
  Neither          = "#fee08b"
)
macro_levels <- names(macro_cols)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
read_bed_gr <- function(path) {
  if (!file.exists(path)) stop("BED file not found: ", path)
  x <- read_tsv(path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(x) < 3) stop("BED file has fewer than 3 columns: ", path)
  GRanges(seqnames = x$X1,
          ranges   = IRanges(start = x$X2 + 1L, end = x$X3))
}

safe_read_bed_gr <- function(path) {
  if (!file.exists(path)) {
    warning("Optional BED not found, skipping: ", path)
    return(GRanges())
  }
  read_bed_gr(path)
}

load_ids <- function(path) {
  if (!file.exists(path)) stop("Gene list not found: ", path)
  x <- read_tsv(path, col_names = FALSE, show_col_types = FALSE)[[1]]
  unique(as.character(x[!is.na(x) & x != "" & !grepl("^#", x)]))
}

# Check that a sufficient fraction of gene IDs matched the txdb.
# Warns if <50% matched; stops if 0 matched.
check_id_match <- function(raw_ids, matched_ids, label) {
  n_raw     <- length(raw_ids)
  n_matched <- length(matched_ids)
  pct       <- 100 * n_matched / max(n_raw, 1)
  cat(sprintf("  %-22s %d IDs in list, %d matched txdb (%.0f%%)\n",
              paste0(label, ":"), n_raw, n_matched, pct))
  if (n_matched == 0)
    stop(label, ": no IDs matched txdb gene names — check ID format (FlyBase vs Ensembl?)")
  if (pct < 50)
    warning(label, ": only ", round(pct), "% of IDs matched txdb — verify gene ID format")
}

save_plot <- function(p, stem, w, h) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, device = cairo_pdf, width = w, height = h)
}

snippet_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 11),
      plot.title      = element_text(size = 12, face = "bold"),
      plot.subtitle   = element_text(size = 10, colour = "grey40")
    )
}

# =============================================================================
# 1. Load & annotate MSL2 peaks
# =============================================================================
cat("Loading genomic ranges...\n")
msl2 <- read_bed_gr(msl2_bed)
mre  <- safe_read_bed_gr(mre_bed)
beaf <- safe_read_bed_gr(beaf_bed)

cat(sprintf("  MSL2 sites: %d | MRE motifs: %d | BEAF motifs: %d\n",
            length(msl2), length(mre), length(beaf)))

ann    <- suppressWarnings(
  annotatePeak(msl2, TxDb = txdb, tssRegion = c(-prom_up, prom_down))
)
ann_df <- as.data.frame(ann)

site_tbl <- tibble(
  chrom             = as.character(seqnames(msl2)),
  start             = start(msl2),
  end               = end(msl2),
  is_chrX           = chrom == "chrX",
  promoter_proximal = str_detect(ann_df$annotation, "Promoter"),
  mre_overlap       = countOverlaps(msl2, mre)  > 0,
  beaf_overlap      = countOverlaps(msl2, beaf) > 0,
  both_motifs       = mre_overlap & beaf_overlap
)

# =============================================================================
# 2. Gene-set loading, ID matching, and promoter construction
# =============================================================================
cat("Loading gene lists...\n")
active_raw <- load_ids(active_gene_list_file)
hk_raw     <- load_ids(hk_gene_list_file)

genes_gr <- genes(txdb)
txdb_ids <- names(genes_gr)

# active_regulated is derived at runtime for internal consistency
active_raw_regulated <- setdiff(active_raw, hk_raw)

active           <- intersect(active_raw,           txdb_ids)
hk               <- intersect(hk_raw,               txdb_ids)
active_regulated <- intersect(active_raw_regulated, txdb_ids)

check_id_match(active_raw,           active,           "Active")
check_id_match(hk_raw,               hk,               "Housekeeping")
check_id_match(active_raw_regulated, active_regulated, "Active_regulated")

active_prom <- promoters(genes_gr[names(genes_gr) %in% active],           upstream = prom_up, downstream = prom_down)
hk_prom     <- promoters(genes_gr[names(genes_gr) %in% hk],               upstream = prom_up, downstream = prom_down)
ar_prom     <- promoters(genes_gr[names(genes_gr) %in% active_regulated], upstream = prom_up, downstream = prom_down)

# =============================================================================
# 3. Attach promoter-set overlaps
# =============================================================================
site_gr <- GRanges(seqnames = site_tbl$chrom,
                   ranges   = IRanges(site_tbl$start, site_tbl$end))

site_tbl <- site_tbl %>%
  mutate(
    active_overlap     = countOverlaps(site_gr, active_prom) > 0,
    hk_overlap         = countOverlaps(site_gr, hk_prom)     > 0,
    active_reg_overlap = countOverlaps(site_gr, ar_prom)      > 0
  )

# Convenience subsets
chrX_tbl      <- filter(site_tbl, is_chrX)
chrX_prom_tbl <- filter(site_tbl, is_chrX, promoter_proximal)

cat(sprintf("  All MCCS: %d | chrX MCCS: %d | chrX promoter MCCS: %d\n",
            nrow(site_tbl), nrow(chrX_tbl), nrow(chrX_prom_tbl)))

# =============================================================================
# 4. Motif prevalence at three filtering levels (Plot 1)
# =============================================================================
compute_motif_fractions <- function(df, label) {
  df %>%
    summarise(
      MRE  = mean(mre_overlap,  na.rm = TRUE),
      BEAF = mean(beaf_overlap, na.rm = TRUE),
      Both = mean(both_motifs,  na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "motif", values_to = "fraction") %>%
    mutate(level = label, n_sites = nrow(df))
}

motif_prev <- bind_rows(
  compute_motif_fractions(site_tbl,      "All MCCS"),
  compute_motif_fractions(chrX_tbl,      "chrX MCCS"),
  compute_motif_fractions(chrX_prom_tbl, "chrX promoter MCCS")
) %>%
  mutate(
    motif = factor(motif, levels = c("MRE", "BEAF", "Both")),
    level = factor(level, levels = c("All MCCS", "chrX MCCS", "chrX promoter MCCS"))
  )

p_motif_prev <- ggplot(motif_prev,
                       aes(x = motif, y = fraction, fill = motif)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = scales::percent(fraction, accuracy = 0.1)),
            vjust = -0.35, size = 3) +
  scale_fill_manual(values = motif_cols) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.14))) +
  facet_wrap(~ level, nrow = 1) +
  labs(x = NULL, y = "Fraction of sites",
       title = "Motif prevalence across site subsets") +
  snippet_theme() +
  theme(strip.text = element_text(size = 10, face = "bold"))

# =============================================================================
# 5. Non-exclusive macro-category fractions (Plot 2)
# =============================================================================
# Neither is defined as no overlap with any of the three gene sets.
# Because active_regulated = setdiff(active, hk), any active_reg_overlap site
# should also be active_overlap; the explicit triple-negative guard is kept
# for correctness and clarity.
macro_counts <- chrX_prom_tbl %>%
  summarise(
    total_sites      = n(),
    All_active       = sum(active_overlap,                                        na.rm = TRUE),
    Housekeeping     = sum(hk_overlap,                                            na.rm = TRUE),
    Active_regulated = sum(active_reg_overlap & !hk_overlap,                      na.rm = TRUE),
    Neither          = sum(!active_overlap & !hk_overlap & !active_reg_overlap,   na.rm = TRUE)
  ) %>%
  pivot_longer(-total_sites, names_to = "macro_group", values_to = "n_sites") %>%
  mutate(
    fraction    = n_sites / total_sites,
    macro_group = factor(macro_group, levels = macro_levels)
  )

p_macro <- ggplot(macro_counts,
                  aes(x = macro_group, y = fraction, fill = macro_group)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = scales::percent(fraction, accuracy = 0.1)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = macro_cols) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Fraction of chrX promoter MCCS",
       title = "chrX promoter MCCS – macro category overlap (non-exclusive)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# =============================================================================
# 6. Motif distribution across macro categories (Plot 3)
# =============================================================================
# Sites are assigned to each category independently (non-exclusive), so a site
# can appear in both All_active and Housekeeping. Motif fractions are then
# computed within each category's membership separately.
motif_macro <- bind_rows(
  chrX_prom_tbl %>% filter(active_overlap)                                      %>% mutate(macro_group = "All_active"),
  chrX_prom_tbl %>% filter(hk_overlap)                                          %>% mutate(macro_group = "Housekeeping"),
  chrX_prom_tbl %>% filter(active_reg_overlap & !hk_overlap)                    %>% mutate(macro_group = "Active_regulated"),
  chrX_prom_tbl %>% filter(!active_overlap & !hk_overlap & !active_reg_overlap) %>% mutate(macro_group = "Neither")
) %>%
  group_by(macro_group) %>%
  summarise(
    n_sites   = n(),
    frac_mre  = mean(mre_overlap,  na.rm = TRUE),
    frac_beaf = mean(beaf_overlap, na.rm = TRUE),
    frac_both = mean(both_motifs,  na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  pivot_longer(c(frac_mre, frac_beaf, frac_both),
               names_to  = "motif",
               values_to = "fraction") %>%
  mutate(
    motif       = recode(motif, frac_mre = "MRE", frac_beaf = "BEAF", frac_both = "Both"),
    macro_group = factor(macro_group, levels = macro_levels)
  )

p_motif_macro <- ggplot(motif_macro,
                        aes(x = motif, y = fraction, fill = macro_group)) +
  geom_col(position = position_dodge(0.8),
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = scales::percent(fraction, accuracy = 0.1)),
            position = position_dodge(0.8),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = macro_cols, breaks = macro_levels) +
  scale_x_discrete(limits = c("MRE", "BEAF", "Both")) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Fraction within macro category",
       fill  = "Macro category",
       title = "Motif distribution across macro categories (non-exclusive)") +
  theme_classic(base_size = 13)

# =============================================================================
# 7. Save plots
# =============================================================================
cat("Saving plots...\n")
save_plot(p_motif_prev,  "motif_prev_combined",      w = 10, h = 4.5)
save_plot(p_macro,       "macro_category_fractions", w = 7,  h = 4.5)
save_plot(p_motif_macro, "macro_category_motif",     w = 9,  h = 4.8)

# =============================================================================
# 8. Provenance record
# =============================================================================
tibble(
  date_run                    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  txdb                        = "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
  promoter_upstream_bp        = prom_up,
  promoter_downstream_bp      = prom_down,
  n_msl2_sites                = length(msl2),
  n_mre_motifs                = length(mre),
  n_beaf_motifs               = length(beaf),
  n_active_ids_in_list        = length(active_raw),
  n_active_ids_matched_txdb   = length(active),
  n_hk_ids_in_list            = length(hk_raw),
  n_hk_ids_matched_txdb       = length(hk),
  n_ar_ids_derived            = length(active_raw_regulated),
  n_ar_ids_matched_txdb       = length(active_regulated),
  n_all_mccs                  = nrow(site_tbl),
  n_chrX_mccs                 = nrow(chrX_tbl),
  n_chrX_promoter_mccs        = nrow(chrX_prom_tbl)
) %>% write_csv(file.path(out_dir, "run_metadata.csv"))

cat("DONE. Outputs in: ", out_dir, "\n", sep = "")
