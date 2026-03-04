#!/usr/bin/env Rscript
# =============================================================================
# Differential binding site discovery using csaw
# Identifies WT-enriched MSL2 and H4K16ac CUT&RUN sites via sliding windows,
# edgeR quasi-likelihood / TREAT framework, and merging into non-redundant sites.
#
# Figures: Upstream of Fig 3B–F (produces BED sites used by downstream scripts)
# Outputs: merged_top_sites_QLF_<target>.bed  — BED files used by deeptools
# =============================================================================

library(csaw)
library(edgeR)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(ChIPseeker)

# ---------------------------------------------------------------------------
# Configuration — EDIT PATHS before running
# ---------------------------------------------------------------------------
metadata_file  <- "data/unified_sample_metadata.csv"
# NOTE: The unified metadata CSV has columns: condition, target, bam_path,
#   control_path. It maps each target to its correct control BAMs:
#   MSL2    → PPI (pre-immune serum)
#   H4K16ac → IgG (non-specific IgG)
output_dir     <- "results/csaw"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Per-target analysis parameters (must match original thresholds exactly)
target_params <- list(
  msl2    = list(window = 20,  ext = 100, max_frag = 250,
                 qlf_fdr = 0.01, qlf_lfc = 2, treat_lfc = 1, min_reps = 4),
  h4k16ac = list(window = 120, ext = 150, max_frag = 500,
                 qlf_fdr = 0.05, qlf_lfc = 1, treat_lfc = 0.5, min_reps = 3)
)

# ---------------------------------------------------------------------------
# Helper: locate WT BAM files from unified metadata CSV
# ---------------------------------------------------------------------------
find_wt_bams <- function(metadata_file, target_name) {
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  # Filter for WT condition and matching target (case-insensitive)
  meta_wt <- meta[tolower(meta$condition) == "wt" &
                   tolower(meta$target) == tolower(target_name), ]
  target_bams  <- unique(na.omit(meta_wt$bam_path))
  control_bams <- unique(na.omit(meta_wt$control_path))
  # Verify BAMs exist
  target_bams  <- target_bams[file.exists(target_bams)]
  control_bams <- control_bams[file.exists(control_bams)]
  stopifnot("No target BAMs found" = length(target_bams) > 0,
            "No control BAMs found" = length(control_bams) > 0)
  cat(sprintf("  Found %d target BAMs, %d control BAMs\n",
              length(target_bams), length(control_bams)))
  return(list(target = target_bams, control = control_bams))
}

# ---------------------------------------------------------------------------
# Main analysis loop — one iteration per target
# ---------------------------------------------------------------------------
for (target_name in names(target_params)) {
  params <- target_params[[target_name]]
  cat(sprintf("\n=== %s (window=%d, ext=%d, FDR≤%g, |logFC|≥%g) ===\n",
              toupper(target_name), params$window, params$ext,
              params$qlf_fdr, params$qlf_lfc))

  # 1. Read BAMs into windows (paired-end counting)
  bams <- find_wt_bams(metadata_file, target_name)
  all_bams <- c(bams$target, bams$control)
  group <- factor(c(rep("target", length(bams$target)),
                    rep("control", length(bams$control))))
  design <- model.matrix(~ group)

  param <- readParam(pe = "both", max.frag = params$max_frag)
  wins <- windowCounts(all_bams, width = params$window, ext = params$ext,
                       param = param)

  # 2. Minimal filter: ≥1 read in ≥ min_reps target replicates
  counts_mat <- assays(wins)$counts
  keep <- rowSums(counts_mat[, seq_along(bams$target), drop = FALSE] >= 1) >=
          params$min_reps
  wins_filt <- wins[keep, ]

  # 3. edgeR normalisation and modelling
  y <- asDGEList(wins_filt)
  y$samples$group <- group
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design)

  # 4a. QLF test (target vs control)
  qlf <- glmQLFTest(fit, coef = "grouptarget")
  tab_qlf <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab_qlf$window <- rownames(tab_qlf)

  # 4b. TREAT test
  treat <- glmTreat(fit, coef = "grouptarget", lfc = params$treat_lfc)

  # 5. Replicate-consistency filter
  cpm_mat <- cpm(y, log = FALSE)
  num_reps_with_signal <- rowSums(
    cpm_mat[, seq_along(bams$target), drop = FALSE] >= 1)
  consistent <- names(num_reps_with_signal[
    num_reps_with_signal >= params$min_reps])

  # 6. Select significant QLF windows
  sites_qlf <- subset(tab_qlf,
    FDR <= params$qlf_fdr &
    window %in% consistent &
    abs(logFC) >= params$qlf_lfc)

  # 7. Merge adjacent significant windows (min gap 50 bp)
  if (nrow(sites_qlf) > 0) {
    win_idx <- as.integer(as.character(sites_qlf$window))
    valid   <- which(!is.na(win_idx) & win_idx >= 1 &
                     win_idx <= length(rowRanges(wins_filt)))
    sig_gr  <- rowRanges(wins_filt)[win_idx[valid]]
    merged_gr <- reduce(sig_gr, min.gapwidth = 50)

    # Assign max |logFC| per merged region as score (absolute value)
    sig_gr$logFC <- as.numeric(sites_qlf$logFC[valid])
    hits <- findOverlaps(merged_gr, sig_gr)
    score_vec <- rep(NA_real_, length(merged_gr))
    score_raw <- tapply(sig_gr$logFC[subjectHits(hits)],
                        queryHits(hits),
                        function(x) x[which.max(abs(x))])
    score_vec[as.integer(names(score_raw))] <- abs(as.numeric(score_raw))

    merged_df <- data.frame(
      seqnames = as.character(seqnames(merged_gr)),
      start    = as.integer(start(merged_gr)),
      end      = as.integer(end(merged_gr)),
      name     = paste0("merged_qlf_", seq_along(merged_gr)),
      score    = round(score_vec, 2),
      strand   = "*",
      stringsAsFactors = FALSE
    )
  } else {
    merged_df <- data.frame()
    merged_gr <- GRanges()
  }

  # 8. Export BED
  bed_file <- file.path(output_dir,
                        sprintf("merged_top_sites_QLF_%s.bed", target_name))
  write.table(merged_df, bed_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  → %d merged sites written to %s\n",
              length(merged_gr), bed_file))
  # 8. Quick genomic annotation (ChIPseeker)
  if (file.exists("data/input/genes.gtf")) {
    txdb <- makeTxDbFromGFF("data/input/genes.gtf")
    anno <- annotatePeak(merged_gr, TxDb = txdb, level = "gene")
    write.csv(as.data.frame(anno),
              file.path(output_dir, sprintf("site_annotation_%s.csv", target_name)),
              row.names = FALSE)
  }
}

cat("\n✅ csaw site discovery complete.\n")
