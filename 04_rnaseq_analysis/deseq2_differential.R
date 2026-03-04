#!/usr/bin/env Rscript
# =============================================================================
# Consolidated RNA-seq differential expression analysis
# Combines: DESeq2 DE, X:A ratio calculation, volcano plots, PCA, and
# sample correlation heatmap into a single pipeline.
#
# Figures: 5A (expression boxplots — via separate script),
#          SFig 3A (PCA), SFig 3B (correlation heatmap), SFig 3C (volcanos)
# Inputs:  genes_counts_RAW.tsv, SampleSheet.tsv, genes.gtf
# Outputs: de_results_*.csv, volcano PNGs, PCA, correlation heatmap
# =============================================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

# ---------------------------------------------------------------------------
# Configuration — EDIT PATHS
# ---------------------------------------------------------------------------
counts_file   <- "data/input/genes_counts_RAW.tsv"
sample_file   <- "data/input/SampleSheet.tsv"
gtf_file      <- "data/input/genes.gtf"
output_dir    <- "results/rnaseq"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Genotype colours (colour-blind friendly)
condition_colors <- c(
  WT           = "#E69F00",
  rox2ko112    = "#009E73",
  rox2ko17     = "#56B4E9",
  rox2ko17rox2fl = "#CC79A7"
)

condition_labels <- c(
  wt           = "WT",
  ko112        = "KO-A",
  ko17         = "KO-B",
  ko17rox2fl   = "Rescue"
)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
raw_counts <- read.delim(counts_file, row.names = 1, stringsAsFactors = FALSE)
sample_meta <- read.delim(sample_file, stringsAsFactors = FALSE)
sample_meta$condition_std <- dplyr::recode(sample_meta$condition, !!!setNames(
  c("WT", "rox2ko112", "rox2ko17", "rox2ko17rox2fl"),
  c("wt", "ko112", "ko17", "ko17rox2fl")
))

# Gene annotation — ALL annotated genes on main chromosomes (for DE analysis)
gtf <- rtracklayer::import(gtf_file)
keep_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
genes_all <- gtf[!is.na(mcols(gtf)$type) & mcols(gtf)$type == "gene" &
                  as.character(seqnames(gtf)) %in% keep_chrs]
gene_chr_map <- data.frame(
  gene_id    = mcols(genes_all)$gene_id,
  chromosome = as.character(seqnames(genes_all)),
  gene_biotype = mcols(genes_all)$gene_biotype,
  stringsAsFactors = FALSE
)

# For DE: use all annotated genes present in count matrix
genes_use <- intersect(rownames(raw_counts), gene_chr_map$gene_id)
counts_filt <- raw_counts[genes_use, ]

# Protein-coding subset (for X:A ratio analysis in section 6)
gene_chr_map_pc <- gene_chr_map[gene_chr_map$gene_biotype == "protein_coding", ]

# ---------------------------------------------------------------------------
# 2. DESeq2 differential expression (batch-aware)
# ---------------------------------------------------------------------------
# Auto-detect batch column from sample metadata
batch_cols <- c("batch", "batch_id", "batchID", "lane", "prep", "plate", "run")
batch_col  <- intersect(batch_cols, names(sample_meta))[1]

if (!is.na(batch_col)) {
  sample_meta[[batch_col]] <- factor(sample_meta[[batch_col]])
  design_formula <- as.formula(paste0("~ ", batch_col, " + condition_std"))
  cat(sprintf("Detected batch column '%s' — using design: %s\n",
              batch_col, deparse(design_formula)))
} else {
  design_formula <- ~ condition_std
  cat("No batch column detected — using design: ~ condition_std\n")
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = sample_meta,
  design    = design_formula
)
dds$condition_std <- relevel(factor(dds$condition_std), ref = "WT")

# Low-count filter: ≥10 counts in at least min-group-size samples
group_sizes <- table(dds$condition_std)
min_reps    <- max(2L, as.integer(min(group_sizes)))
keep        <- rowSums(counts(dds) >= 10) >= min_reps
dds         <- dds[keep, ]
cat(sprintf("Filtering: ≥10 counts in ≥%d samples → %d genes retained.\n",
            min_reps, nrow(dds)))

dds <- DESeq(dds)

comparisons <- list(
  rox2ko17_vs_WT        = c("condition_std", "rox2ko17",        "WT"),
  rox2ko112_vs_WT       = c("condition_std", "rox2ko112",       "WT"),
  rescue_vs_rox2ko17    = c("condition_std", "rox2ko17rox2fl",  "rox2ko17"),
  rescue_vs_WT          = c("condition_std", "rox2ko17rox2fl",  "WT")
)

de_results <- list()
for (name in names(comparisons)) {
  res <- results(dds, contrast = comparisons[[name]], alpha = 0.01)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_chr_map, by = "gene_id") %>%
    arrange(padj)
  de_results[[name]] <- res_df
  write.csv(res_df, file.path(output_dir, paste0("de_results_", name, ".csv")),
            row.names = FALSE)
}

cat(sprintf("DE results written for %d comparisons.\n", length(comparisons)))

# ---------------------------------------------------------------------------
# 3. Volcano plots (SFig 3C)
# ---------------------------------------------------------------------------
make_volcano <- function(df, title, padj_cut = 0.01, lfc_cut = 1) {
  df <- df %>%
    mutate(
      sig = case_when(
        is.na(padj) ~ "NS",
        padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )

  # Label top genes
  top_padj <- df %>% filter(sig != "NS") %>% slice_min(padj, n = 20)
  top_lfc  <- df %>% filter(sig != "NS") %>% slice_max(abs(log2FoldChange), n = 20)
  to_label <- distinct(bind_rows(top_padj, top_lfc), gene_id, .keep_all = TRUE)

  ggplot(df, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(size = 0.6, alpha = 0.5) +
    geom_text_repel(data = to_label,
                    aes(label = gene_id), size = 2, max.overlaps = 30) +
    scale_colour_manual(values = c(NS = "#B0B0B0", Down = "#6BAED6", Up = "#E63946")) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
}

volcano_plots <- lapply(names(de_results), function(n) {
  make_volcano(de_results[[n]], gsub("_", " ", n))
})

combined_volcanoes <- wrap_plots(volcano_plots, ncol = 2)
ggsave(file.path(output_dir, "volcano_plots.png"),
       combined_volcanoes, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "volcano_plots.pdf"),
       combined_volcanoes, width = 12, height = 10)

# ---------------------------------------------------------------------------
# 4. PCA (SFig 3A)
# ---------------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition_std", returnData = TRUE)
pca_var  <- attr(pca_data, "percentVar")

p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = condition_std)) +
  geom_point(size = 3) +
  scale_colour_manual(values = condition_colors) +
  labs(
    x = sprintf("PC1 (%.0f%%)", pca_var[1] * 100),
    y = sprintf("PC2 (%.0f%%)", pca_var[2] * 100),
    colour = "Genotype"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(output_dir, "RNAseq_PCA.png"), p_pca,
       width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "RNAseq_PCA.pdf"), p_pca,
       width = 7, height = 5)

# ---------------------------------------------------------------------------
# 5. Sample correlation heatmap (SFig 3B)
# ---------------------------------------------------------------------------
vsd_mat <- assay(vsd)
cor_mat <- cor(vsd_mat, method = "pearson")

# Sample annotation for heatmap (keep original sample IDs as row/col names)
sample_conditions <- setNames(sample_meta$condition_std, sample_meta$name)
anno_df <- data.frame(
  Condition = sample_conditions[colnames(cor_mat)],
  row.names = colnames(cor_mat)
)
anno_colors <- list(Condition = condition_colors)

pheatmap(cor_mat,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = anno_df,
         annotation_colors = anno_colors,
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize_number = 7,
         filename = file.path(output_dir, "RNAseq_correlation_heatmap.pdf"),
         width = 8, height = 7)
pheatmap(cor_mat,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = anno_df,
         annotation_colors = anno_colors,
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize_number = 7,
         filename = file.path(output_dir, "RNAseq_correlation_heatmap.png"),
         width = 8, height = 7)

# ---------------------------------------------------------------------------
# 6. X:A expression ratio (Zhang method) — protein-coding genes only
# ---------------------------------------------------------------------------
# Restrict to protein-coding genes on main chr arms for X:A analysis
# (manuscript: "restricted to protein-coding genes located on major
#  chromosome arms")
norm_counts <- counts(dds, normalized = TRUE)
pc_genes_in_dds <- intersect(rownames(norm_counts), gene_chr_map_pc$gene_id)

xa_long <- as.data.frame(norm_counts[pc_genes_in_dds, ]) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expr") %>%
  left_join(sample_meta %>% dplyr::select(name, condition_std), by = c("sample" = "name")) %>%
  left_join(gene_chr_map_pc, by = "gene_id") %>%
  mutate(chr_type = ifelse(chromosome == "chrX", "X", "Autosome"))

# RPM per-gene (Zhang method: gene-count corrected)
gene_counts_per_type <- gene_chr_map_pc %>%
  mutate(chr_type = ifelse(chromosome == "chrX", "X", "Autosome")) %>%
  count(chr_type, name = "n_genes")

xa_ratios <- xa_long %>%
  group_by(sample, condition_std, chr_type) %>%
  summarise(total = sum(expr, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    xa_long %>% group_by(sample) %>%
      summarise(lib_size = sum(expr, na.rm = TRUE), .groups = "drop"),
    by = "sample"
  ) %>%
  mutate(rpm = (total / lib_size) * 1e6) %>%
  left_join(gene_counts_per_type, by = "chr_type") %>%
  mutate(expr_per_gene = rpm / n_genes) %>%
  dplyr::select(sample, condition_std, chr_type, expr_per_gene) %>%
  pivot_wider(names_from = chr_type, values_from = expr_per_gene) %>%
  mutate(xa_ratio = X / Autosome)

xa_summary <- xa_ratios %>%
  group_by(condition_std) %>%
  summarise(
    mean_xa = mean(xa_ratio), sd_xa = sd(xa_ratio),
    n = n(), .groups = "drop"
  )

write.csv(xa_ratios,  file.path(output_dir, "xa_ratios_per_sample.csv"),  row.names = FALSE)
write.csv(xa_summary, file.path(output_dir, "xa_ratio_summary.csv"),      row.names = FALSE)

cat("\n✅ RNA-seq analysis pipeline complete.\n")
cat("   DE results:    results/rnaseq/de_results_*.csv\n")
cat("   Volcanos:      results/rnaseq/volcano_plots.{png,pdf}\n")
cat("   PCA:           results/rnaseq/RNAseq_PCA.{png,pdf}\n")
cat("   Cor. heatmap:  results/rnaseq/RNAseq_correlation_heatmap.{png,pdf}\n")
cat("   X:A ratios:    results/rnaseq/xa_ratio*.csv\n")
