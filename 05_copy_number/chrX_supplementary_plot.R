#!/usr/bin/env Rscript
# =============================================================================
# Smoothed chrX copy number profiles (supplementary figure)
# Reads spike-in-normalised bin coverage directly from .cov_forR.bed files,
# scales to absolute copy numbers, applies LOESS smoothing, and creates
# faceted chrX copy-number profiles per genotype.
#
# Figures: SFig 4A–B (chrX CN profiles)
# Inputs:  cn_per_cellline_spike/<sample>/input_50000.cov_forR.bed
# =============================================================================

library(tidyverse)
library(patchwork)

output_dir <- "results/plots/spikein"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

EXPECTED_X_COPIES <- 2

# ---------------------------------------------------------------------------
# Sample table (same as copy_number_estimation.R)
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

keep_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
auto_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R")

# ---------------------------------------------------------------------------
# Load bin-level data from .cov_forR.bed files
# ---------------------------------------------------------------------------
bins_all <- list()
for (i in seq_len(nrow(input_samples))) {
  s <- input_samples[i, ]
  bed_file <- file.path(s$bin_dir, "input_50000.cov_forR.bed")
  if (!file.exists(bed_file)) { warning("Missing: ", bed_file); next }
  bins <- read.delim(bed_file, header = FALSE,
                     col.names = c("chr", "start", "end", "counts"))
  bins <- bins %>%
    filter(chr %in% keep_chrs) %>%
    mutate(condition = s$condition, replicate = s$replicate)

  # RPKM-like within-sample normalisation (matching original pipeline)
  bin_lengths_kb <- (bins$end - bins$start + 1) / 1000
  total_reads_M  <- max(sum(bins$counts, na.rm = TRUE) / 1e6, 1)
  bins$cov       <- (bins$counts / bin_lengths_kb) / total_reads_M

  bins_all[[i]] <- bins
}
bins_all <- bind_rows(bins_all)

# WT autosomal baseline (same logic as copy_number_estimation.R)
wt_auto_baseline <- bins_all %>%
  filter(condition == "WT", chr %in% auto_chrs) %>%
  pull(cov) %>% median(na.rm = TRUE)

bins_all <- bins_all %>%
  mutate(copy_number = (cov / wt_auto_baseline) * 4) %>%
  filter(chr == "chrX")

# ---------------------------------------------------------------------------
# Plot: Smoothed chrX profile
# ---------------------------------------------------------------------------
bins_all$condition <- factor(bins_all$condition,
  levels = c("WT", "rox2ko112", "rox2ko17", "rox2ko17rox2fl"))

# Subsample for plotting (10 % per condition for scatter)
set.seed(42)
bins_sub <- bins_all %>% group_by(condition) %>% slice_sample(prop = 0.1)

p <- ggplot(bins_sub, aes(start / 1e6, copy_number)) +
  geom_point(alpha = 0.15, size = 0.3, colour = "grey50") +
  geom_smooth(data = bins_all, method = "loess", span = 0.1,
              aes(colour = condition), se = FALSE, linewidth = 0.8) +
  geom_hline(yintercept = EXPECTED_X_COPIES, linetype = "dashed",
             colour = "red", linewidth = 0.4) +
  facet_wrap(~ condition, ncol = 1) +
  scale_colour_manual(values = c(WT = "#E69F00", rox2ko112 = "#009E73",
                                 rox2ko17 = "#56B4E9", rox2ko17rox2fl = "#CC79A7")) +
  labs(x = "chrX position (Mb)", y = "Estimated copy number") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(output_dir, "chrX_copy_number_smoothed.png"), p,
       width = 10, height = 10, dpi = 300)
ggsave(file.path(output_dir, "chrX_copy_number_smoothed.pdf"), p,
       width = 10, height = 10)

cat("✅ chrX supplementary plot complete.\n")
