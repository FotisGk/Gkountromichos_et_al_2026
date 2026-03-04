#!/usr/bin/env Rscript
# =============================================================================
# Gene-level copy number analysis
# Counts reads per gene from INPUT BAMs, normalises to WT autosomal baseline,
# and exports per-gene copy numbers used by xa_ratio_barplots.R.
#
# Figures: 5D dependency (all_genes_copy_numbers.csv)
# Inputs:  INPUT BAM files, genes.gtf, chipseq_sample_metadata.csv
# Outputs: all_genes_copy_numbers.csv, gene_level_copy_number_summary.csv
# =============================================================================

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(dplyr)
library(tidyr)

output_dir <- "results/gene_level_cnv"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load gene annotations
# ---------------------------------------------------------------------------
genes <- rtracklayer::import("data/input/genes.gtf")
main_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
genes <- genes[genes$type == "gene" & as.character(seqnames(genes)) %in% main_chrs]

chr_info <- data.frame(
  chr = main_chrs,
  expected_dna_copies = c(4, 4, 4, 4, 4, 2),
  chr_type = c(rep("autosome", 5), "sex")
)

genes_df <- data.frame(
  gene_id   = genes$gene_id,
  gene_name = ifelse(is.null(genes$gene_name), genes$gene_id, genes$gene_name),
  chr = as.character(seqnames(genes)),
  start = start(genes), end = end(genes), width = width(genes)
) %>% left_join(chr_info, by = "chr")

# ---------------------------------------------------------------------------
# Load sample metadata (INPUT samples)
# ---------------------------------------------------------------------------
metadata <- read.csv("data/chipseq_sample_metadata.csv", stringsAsFactors = FALSE)
input_samples <- metadata %>%
  dplyr::select(condition, input_control, replicate) %>%
  filter(!is.na(input_control), input_control != "") %>% distinct()

# ---------------------------------------------------------------------------
# Count reads per gene
# ---------------------------------------------------------------------------
count_reads_per_gene <- function(bam_file, gene_ranges, quality_threshold = 10) {
  if (!file.exists(bam_file)) return(rep(0, length(gene_ranges)))
  param <- ScanBamParam(
    what = c("pos", "mapq"),
    flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
    mapqFilter = quality_threshold
  )
  alns <- readGAlignments(bam_file, param = param)
  if (length(alns) == 0) return(rep(0, length(gene_ranges)))
  countOverlaps(gene_ranges, alns)
}

gene_counts <- genes_df
for (i in seq_len(nrow(input_samples))) {
  s <- input_samples[i, ]
  col_name <- paste0(s$condition, "_rep", s$replicate)
  gene_counts[[col_name]] <- count_reads_per_gene(s$input_control, genes)
}

# Integrate replicates (median)
conditions <- unique(gsub("_rep\\d+", "", grep("_rep", names(gene_counts), value = TRUE)))
for (cond in conditions) {
  reps <- grep(paste0("^", cond, "_rep"), names(gene_counts), value = TRUE)
  gene_counts[[paste0(cond, "_integrated")]] <- apply(gene_counts[reps], 1, median, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Normalise to WT autosomal baseline
# ---------------------------------------------------------------------------
auto_mask <- gene_counts$chr_type == "autosome"
wt_auto_median <- median(gene_counts$WT_integrated[auto_mask &
                         gene_counts$WT_integrated > 0], na.rm = TRUE)

int_cols <- grep("_integrated$", names(gene_counts), value = TRUE)
for (col in int_cols) {
  cn_col <- gsub("_integrated", "_copy_number", col)
  gene_counts[[cn_col]] <- (gene_counts[[col]] / wt_auto_median) * 4
}

# ---------------------------------------------------------------------------
# Export long-format table
# ---------------------------------------------------------------------------
cn_cols <- grep("_copy_number$", names(gene_counts), value = TRUE)
gene_cn_long <- gene_counts %>%
  dplyr::select(gene_id, gene_name, chr, start, end, width, expected_dna_copies,
         chr_type, all_of(cn_cols)) %>%
  pivot_longer(all_of(cn_cols), names_to = "condition",
               values_to = "observed_copy_number") %>%
  mutate(condition = gsub("_copy_number", "", condition)) %>%
  filter(observed_copy_number > 0)

write.csv(gene_cn_long, file.path(output_dir, "all_genes_copy_numbers.csv"),
          row.names = FALSE)

# Summary
gene_summary <- gene_cn_long %>%
  group_by(chr, condition) %>%
  summarise(
    n_genes = n(),
    median_cn = median(observed_copy_number, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(gene_summary, file.path(output_dir, "gene_level_copy_number_summary.csv"),
          row.names = FALSE)

cat("✅ Gene-level CNV analysis complete.\n")
