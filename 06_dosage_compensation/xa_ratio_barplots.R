#!/usr/bin/env Rscript
# =============================================================================
# Three-way X:A ratio barplots (CN, raw expression, CN-corrected)
# Tests the copy-number buffering hypothesis: do elevated X copies make
# MSL-mediated upregulation dispensable?
#
# Figures: 5D
# Inputs:  genes_counts_RAW.tsv, SampleSheet.tsv, genes.gtf,
#          chr_stats_canonical.csv, all_genes_copy_numbers.csv
# =============================================================================

library(tidyverse)
library(rtracklayer)
library(patchwork)

output_dir <- "results/copy_number_hypothesis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# Gene annotation (protein-coding)
gtf <- rtracklayer::import("data/input/genes.gtf")
keep_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
genes_pc <- gtf[!is.na(mcols(gtf)$type) & mcols(gtf)$type == "gene" &
                 as.character(seqnames(gtf)) %in% keep_chrs &
                 mcols(gtf)$gene_biotype == "protein_coding"]
gene_chr_map <- data.frame(
  gene_id    = mcols(genes_pc)$gene_id,
  chromosome = as.character(seqnames(genes_pc)),
  stringsAsFactors = FALSE
)

# Metadata
sample_meta <- read.delim("data/input/SampleSheet.tsv", stringsAsFactors = FALSE) %>%
  mutate(condition_std = recode(condition,
    wt = "WT", ko17 = "rox2ko17", ko112 = "rox2ko112",
    ko17rox2fl = "rox2ko17rox2fl"))

# Raw counts (protein-coding)
raw_counts <- read.delim("data/input/genes_counts_RAW.tsv",
                         row.names = 1, stringsAsFactors = FALSE)
genes_use <- intersect(rownames(raw_counts), gene_chr_map$gene_id)
raw_counts_pc <- raw_counts[genes_use, ]

# Chromosome-level copy numbers
cn_chr <- read.csv("results/tables/spikein/chr_stats_canonical.csv",
                   stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# X and autosome copy numbers per condition
# ---------------------------------------------------------------------------
x_cn <- cn_chr %>% filter(chromosome == "chrX") %>%
  dplyr::select(condition, X_cn = median_copy_number)
a_cn <- cn_chr %>% filter(chromosome != "chrX") %>%
  group_by(condition) %>%
  summarise(A_cn = median(median_copy_number, na.rm = TRUE), .groups = "drop")

cn_comparison <- left_join(x_cn, a_cn, by = "condition") %>%
  dplyr::rename(condition_std = condition) %>%
  mutate(xa_cn_ratio = X_cn / A_cn)

# ---------------------------------------------------------------------------
# Three-way X:A ratios
# ---------------------------------------------------------------------------
gene_n <- gene_chr_map %>%
  mutate(chr_type = ifelse(chromosome == "chrX", "X", "Autosome")) %>%
  count(chr_type, name = "n_genes")

raw_long <- raw_counts_pc %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "counts") %>%
  left_join(sample_meta, by = c("sample" = "name")) %>%
  left_join(gene_chr_map, by = "gene_id") %>%
  mutate(chr_type = ifelse(chromosome == "chrX", "X", "Autosome"))

lib_sizes <- raw_long %>% group_by(sample) %>%
  summarise(lib = sum(counts, na.rm = TRUE), .groups = "drop")

total_expr <- raw_long %>%
  group_by(condition_std, replicate, sample, chr_type) %>%
  summarise(total = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  left_join(lib_sizes, by = "sample") %>%
  mutate(rpm = (total / lib) * 1e6) %>%
  left_join(gene_n, by = "chr_type") %>%
  mutate(expr_per_gene = rpm / n_genes)

# 1. Raw X:A (Zhang method)
xa_raw <- total_expr %>%
  dplyr::select(condition_std, replicate, chr_type, expr_per_gene) %>%
  pivot_wider(names_from = chr_type, values_from = expr_per_gene) %>%
  mutate(xa_ratio = X / Autosome) %>%
  group_by(condition_std) %>%
  summarise(mean_xa = mean(xa_ratio), sd_xa = sd(xa_ratio), .groups = "drop") %>%
  mutate(method = "Raw (Zhang)")

# 2. CN-corrected X:A
expr_cn <- total_expr %>%
  left_join(cn_comparison %>% dplyr::select(condition_std, X_cn, A_cn), by = "condition_std") %>%
  mutate(cn = ifelse(chr_type == "X", X_cn, A_cn),
         expr_per_copy = expr_per_gene / cn)

xa_cn <- expr_cn %>%
  dplyr::select(condition_std, replicate, chr_type, expr_per_copy) %>%
  pivot_wider(names_from = chr_type, values_from = expr_per_copy) %>%
  mutate(xa_ratio = X / Autosome) %>%
  group_by(condition_std) %>%
  summarise(mean_xa = mean(xa_ratio), sd_xa = sd(xa_ratio), .groups = "drop") %>%
  mutate(method = "CN-corrected")

# 3. Expected (copy number ratio only)
xa_expected <- cn_comparison %>%
  dplyr::select(condition_std, xa_cn_ratio) %>%
  mutate(mean_xa = xa_cn_ratio, sd_xa = 0, method = "Expected (CN only)")

xa_all <- bind_rows(xa_raw, xa_cn, xa_expected) %>%
  mutate(
    condition_std = factor(condition_std,
      levels = c("WT", "rox2ko112", "rox2ko17", "rox2ko17rox2fl")),
    method = factor(method,
      levels = c("Expected (CN only)", "Raw (Zhang)", "CN-corrected"))
  )

write.csv(xa_all, file.path(output_dir, "three_way_xa_comparison.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Plot: three-way barplot (Fig 5D)
# ---------------------------------------------------------------------------
genotype_colors <- c(WT = "#E69F00", rox2ko112 = "#009E73",
                     rox2ko17 = "#56B4E9", rox2ko17rox2fl = "#CC79A7")

facet_labels <- c(
  "Expected (CN only)" = "DNA copy number (CN)",
  "Raw (Zhang)"        = "Raw expression",
  "CN-corrected"       = "CN-corrected expression"
)

p <- ggplot(xa_all, aes(condition_std, mean_xa, fill = condition_std)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "#2A9D8F", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#E63946", linewidth = 0.5) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_errorbar(aes(ymin = mean_xa - sd_xa, ymax = mean_xa + sd_xa),
                width = 0.25, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_xa)),
            vjust = -0.5, size = 4, fontface = "bold") +
  facet_wrap(~ method, nrow = 1, labeller = labeller(method = facet_labels)) +
  scale_fill_manual(values = genotype_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5)) +
  labs(x = NULL, y = "X:A Ratio") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "three_way_xa_ratio.png"), p,
       width = 12, height = 6, dpi = 300)
ggsave(file.path(output_dir, "three_way_xa_ratio.pdf"), p,
       width = 12, height = 6)

cat("✅ X:A ratio barplots complete.\n")
