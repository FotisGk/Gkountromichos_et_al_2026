#!/usr/bin/env Rscript
# =============================================================================
# Chromosome-level copy number estimation from INPUT sequencing
# Reads pre-computed bin-level coverage (from spikein_normalization.sh),
# scales to absolute copy numbers using WT autosomal baseline, and produces
# per-chromosome barplots with replicate-based error bars.
#
# Figures: 5C (copy number barplot)
# Inputs:  cn_per_cellline_spike/<sample>/input_50000.cov_forR.bed
# Outputs: chr_stats_canonical.csv, copy_number_barplot.{png,pdf}
# =============================================================================

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(patchwork)

output_dir <- "results/tables/spikein"
plot_dir   <- "results/plots/spikein"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
input_samples <- tribble(
  ~condition,          ~replicate, ~bin_dir,
  "WT",                1,          "cn_per_cellline_spike/s2a_wt_rep1",
 "WT",                4,          "cn_per_cellline_spike/s2a_wt_rep4",
  "WT",                5,          "cn_per_cellline_spike/s2a_wt_rep5",
  "rox2ko17",          1,          "cn_per_cellline_spike/s2a_rox2ko17_rep1",
  "rox2ko17",          2,          "cn_per_cellline_spike/s2a_rox2ko17_rep2",
  "rox2ko17",          5,          "cn_per_cellline_spike/s2a_rox2ko17_rep5",
  "rox2ko112",         1,          "cn_per_cellline_spike/s2a_rox2ko112_rep1",
  "rox2ko112",         2,          "cn_per_cellline_spike/s2a_rox2ko112_rep2",
  "rox2ko112",         3,          "cn_per_cellline_spike/s2a_rox2ko112_rep3",
  "rox2ko17rox2fl",    1,          "cn_per_cellline_spike/s2a_rox2ko17rox2fl_rep1",
  "rox2ko17rox2fl",    2,          "cn_per_cellline_spike/s2a_rox2ko17rox2fl_rep2",
  "rox2ko17rox2fl",    3,          "cn_per_cellline_spike/s2a_rox2ko17rox2fl_rep3"
)

keep_chrs    <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
autosome_list <- c("chr2L", "chr2R", "chr3L", "chr3R")
expected_copies <- data.frame(
  chromosome = keep_chrs,
  expected   = c(4, 4, 4, 4, 4, 4),
  chr_type   = c(rep("Autosome", 5), "X chromosome")
)

# ---------------------------------------------------------------------------
# Load bin-level data
# ---------------------------------------------------------------------------
chr_medians_list <- list()
for (i in seq_len(nrow(input_samples))) {
  s <- input_samples[i, ]
  bed_file <- file.path(s$bin_dir, "input_50000.cov_forR.bed")
  if (!file.exists(bed_file)) { warning("Missing: ", bed_file); next }

  bins <- read.delim(bed_file, header = FALSE,
                     col.names = c("chr", "start", "end", "counts"))
  bins <- bins %>% filter(chr %in% keep_chrs)

  # RPKM-like within-sample normalisation: convert raw coverage to a
  # density metric (counts / bin_length_kb / total_reads_millions) so
  # samples with different sequencing depths are internally comparable
  # before cross-sample scaling.
  bin_lengths_kb   <- (bins$end - bins$start + 1) / 1000
  total_reads_M    <- max(sum(bins$counts, na.rm = TRUE) / 1e6, 1)
  bins$copy_number <- (bins$counts / bin_lengths_kb) / total_reads_M

  chr_meds <- bins %>%
    group_by(chr) %>%
    summarise(chr_median = median(copy_number, na.rm = TRUE), .groups = "drop") %>%
    mutate(condition = s$condition, replicate = s$replicate)

  chr_medians_list[[paste0(s$condition, "_rep", s$replicate)]] <- chr_meds
}

chr_medians <- bind_rows(chr_medians_list)

# ---------------------------------------------------------------------------
# Scale to absolute copy numbers (WT autosomal baseline = 4 copies)
# ---------------------------------------------------------------------------
wt_auto_baseline <- chr_medians %>%
  filter(condition == "WT", chr %in% autosome_list) %>%
  pull(chr_median) %>% median(na.rm = TRUE)

chr_medians <- chr_medians %>%
  mutate(copy_number = (chr_median / wt_auto_baseline) * 4)

# Per-condition summary
chr_stats <- chr_medians %>%
  dplyr::rename(chromosome = chr) %>%
  group_by(condition, chromosome) %>%
  summarise(
    median_copy_number = median(copy_number, na.rm = TRUE),
    mean_copy_number   = mean(copy_number, na.rm = TRUE),
    sd_copy_number     = sd(copy_number, na.rm = TRUE),
    n = n(), .groups = "drop"
  ) %>%
  left_join(expected_copies, by = "chromosome")

write.csv(chr_stats, file.path(output_dir, "chr_stats_canonical.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Barplot (Fig 5C)
# ---------------------------------------------------------------------------
chr_stats <- chr_stats %>%
  mutate(
    chr_clean = factor(gsub("chr", "", chromosome),
                       levels = c("2L", "2R", "3L", "3R", "4", "X")),
    condition = factor(condition,
                       levels = c("WT", "rox2ko112", "rox2ko17", "rox2ko17rox2fl"))
  )

genotype_colors <- c(WT = "#E69F00", rox2ko112 = "#009E73",
                     rox2ko17 = "#56B4E9", rox2ko17rox2fl = "#CC79A7")

p <- ggplot(chr_stats, aes(chr_clean, median_copy_number, fill = condition)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = median_copy_number - sd_copy_number,
                    ymax = median_copy_number + sd_copy_number),
                position = position_dodge(0.8), width = 0.25) +
  geom_hline(data = expected_copies %>%
               mutate(chr_clean = factor(gsub("chr","",chromosome),
                      levels = c("2L","2R","3L","3R","4","X"))),
             aes(yintercept = expected), linetype = "dashed", colour = "red", alpha = 0.5) +
  scale_fill_manual(values = genotype_colors) +
  scale_y_continuous(breaks = 0:10) +
  labs(x = "Chromosome", y = "Estimated copy number",
       fill = "Genotype") +
  theme_classic(base_size = 14)

ggsave(file.path(plot_dir, "copy_number_barplot.png"), p,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "copy_number_barplot.pdf"), p,
       width = 10, height = 6)

cat("✅ Copy number estimation complete.\n")
