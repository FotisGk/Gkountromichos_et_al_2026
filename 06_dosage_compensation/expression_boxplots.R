#!/usr/bin/env Rscript
# =============================================================================
# RNA-seq expression boxplots by chromosome arm
# Per-gene log2FC (KO vs WT) boxplots coloured by X vs autosome.
#
# Figures: 5A
# Inputs:  genes_counts_RAW.tsv, SampleSheet.tsv, genes.gtf
# =============================================================================

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(rtracklayer)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
counts_file <- "data/input/genes_counts_RAW.tsv"
sample_file <- "data/input/SampleSheet.tsv"
gtf_file    <- "data/input/genes.gtf"
output_dir  <- "results/plots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_counts <- read.delim(counts_file, row.names = 1, stringsAsFactors = FALSE)
sample_meta <- read.delim(sample_file, stringsAsFactors = FALSE)
sample_meta$condition_std <- dplyr::recode(sample_meta$condition,
  wt = "WT", ko112 = "rox2ko112", ko17 = "rox2ko17",
  ko17rox2fl = "rox2ko17rox2fl")

# Protein-coding genes on main chromosomes
gtf <- rtracklayer::import(gtf_file)
keep_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
genes_pc <- gtf[!is.na(mcols(gtf)$type) & mcols(gtf)$type == "gene" &
                 as.character(seqnames(gtf)) %in% keep_chrs &
                 mcols(gtf)$gene_biotype == "protein_coding"]
gene_chr_map <- data.frame(
  gene_id    = mcols(genes_pc)$gene_id,
  chromosome = as.character(seqnames(genes_pc)),
  chr_type   = ifelse(as.character(seqnames(genes_pc)) == "chrX",
                      "X chromosome", "Autosomes"),
  stringsAsFactors = FALSE
)

genes_use <- intersect(rownames(raw_counts), gene_chr_map$gene_id)
counts_pc <- raw_counts[genes_use, ]

# ---------------------------------------------------------------------------
# DESeq2 normalisation (size factors only — no DE testing needed)
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts_pc, sample_meta, ~ condition_std)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# ---------------------------------------------------------------------------
# Compute per-gene log2FC (condition mean / WT mean)
# ---------------------------------------------------------------------------
expr_long <- as.data.frame(norm_counts) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  left_join(sample_meta %>% dplyr::select(name, condition_std),
            by = c("sample" = "name"))

condition_means <- expr_long %>%
  group_by(gene_id, condition_std) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition_std, values_from = mean_expr) %>%
  left_join(gene_chr_map, by = "gene_id")

log2fc_data <- condition_means %>%
  pivot_longer(cols = c(rox2ko17, rox2ko112, rox2ko17rox2fl),
               names_to = "condition_std", values_to = "ko_expr") %>%
  mutate(
    log2fc = log2((ko_expr + 1) / (WT + 1)),
    chromosome = factor(chromosome, levels = keep_chrs)
  ) %>%
  filter(!is.na(log2fc) & is.finite(log2fc))

# ---------------------------------------------------------------------------
# Panel: log2FC boxplots faceted by condition
# ---------------------------------------------------------------------------
p <- ggplot(log2fc_data, aes(x = chromosome, y = log2fc, fill = chr_type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.6, outlier.alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~ condition_std, nrow = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_fill_manual(values = c("X chromosome" = "#fc8d62", "Autosomes" = "#8da0cb")) +
  scale_x_discrete(labels = function(x) gsub("^chr", "", x)) +
  labs(x = "Chromosome", y = "RNA log2FC (KO vs WT)") +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(output_dir, "expression_boxplots_by_chr.png"), p,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(output_dir, "expression_boxplots_by_chr.pdf"), p,
       width = 12, height = 5)

cat("✅ Expression boxplots complete.\n")
