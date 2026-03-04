#!/usr/bin/env Rscript
# =============================================================================
# Site dissection analysis: overlap classification, UpSet plots, genomic
# feature annotation, and signal enrichment by site class.
#
# Figures: 3B (chromosomal distribution stacked bar), 3C (UpSet),
#          3D (signal enrichment violins), 3E (genomic feature annotation),
#          3F (distance to TSS/HAS densities)
# Inputs:  merged_top_sites_QLF_msl2.bed, merged_top_sites_QLF_h4k16ac.bed,
#          has_sites.bed, pionx_sites.bed, WT ratio bigwigs
# Outputs: UpSet plot, Venn, overlap table, genomic feature barplot,
#          enrichment violins
# =============================================================================

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(UpSetR)
library(ggvenn)
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

output_dir <- "results/site_dissection"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load BED files
# ---------------------------------------------------------------------------
bedfiles <- list(
  msl2    = "results/csaw/merged_top_sites_QLF_msl2.bed",
  h4k16ac = "results/csaw/merged_top_sites_QLF_h4k16ac.bed",
  has     = "data/has_sites.bed",
  pionx   = "data/pionx_sites.bed"
)

read_bed_gr <- function(path) {
  bed <- read_tsv(path, col_names = c("chrom", "start", "end"),
                  col_types = cols_only(chrom = col_character(),
                                       start = col_integer(),
                                       end = col_integer()))
  GRanges(seqnames = bed$chrom,
          ranges = IRanges(start = bed$start + 1, end = bed$end))
}

gr_list <- lapply(bedfiles, read_bed_gr)
msl2_gr    <- gr_list$msl2
h4k16ac_gr <- gr_list$h4k16ac
has_gr     <- gr_list$has
pionx_gr   <- gr_list$pionx

cat(sprintf("Site counts — MSL2: %d, H4K16ac: %d, HAS: %d, PionX: %d\n",
            length(msl2_gr), length(h4k16ac_gr), length(has_gr), length(pionx_gr)))

# ---------------------------------------------------------------------------
# 2. Overlap analysis (MSL2 / HAS / PionX)
# ---------------------------------------------------------------------------
all_peaks <- reduce(c(msl2_gr, has_gr, pionx_gr), ignore.strand = TRUE)

olap_df <- data.frame(
  MSL2  = as.integer(countOverlaps(all_peaks, msl2_gr) > 0),
  HAS   = as.integer(countOverlaps(all_peaks, has_gr) > 0),
  PionX = as.integer(countOverlaps(all_peaks, pionx_gr) > 0)
)

# Label columns with counts
raw_counts <- c(length(msl2_gr), length(has_gr), length(pionx_gr))
set_labels <- paste0(c("MSL2", "HAS", "PionX"), " (", raw_counts, ")")
colnames(olap_df) <- set_labels

# Membership table
overlap_tbl <- tibble(
  chrom = as.character(seqnames(all_peaks)),
  start = start(all_peaks),
  end   = end(all_peaks),
  id    = paste(chrom, start, end, sep = "_"),
  MSL2  = countOverlaps(all_peaks, msl2_gr) > 0,
  HAS   = countOverlaps(all_peaks, has_gr) > 0,
  PionX = countOverlaps(all_peaks, pionx_gr) > 0
)
write_csv(overlap_tbl, file.path(output_dir, "DCC_site_overlap_table.csv"))

# ---------------------------------------------------------------------------
# 3. Venn diagram
# ---------------------------------------------------------------------------
venn_list <- list(
  MSL2  = overlap_tbl %>% filter(MSL2)  %>% pull(id),
  HAS   = overlap_tbl %>% filter(HAS)   %>% pull(id),
  PionX = overlap_tbl %>% filter(PionX) %>% pull(id)
)

png(file.path(output_dir, "venn_DCC_overlap.png"), width = 780, height = 680)
ggvenn(venn_list, fill_color = c("#1b7837", "#053061", "#b2182b"))
dev.off()

# ---------------------------------------------------------------------------
# 4. UpSet plot
# ---------------------------------------------------------------------------
png(file.path(output_dir, "upset_DCC_overlap.png"),
    width = 1200, height = 800, res = 250)
UpSetR::upset(
  olap_df,
  sets = set_labels,
  order.by = "freq",
  main.bar.color = "#2E86AB",
  sets.bar.color = "#F6C85F",
  matrix.color  = "#A23B72",
  sets.x.label  = "Peak Sets",
  mainbar.y.label = "Intersection Size",
  text.scale = c(1.5, 1.5, 1.5, 0.01, 1.5, 1.2)
)
dev.off()

# ---------------------------------------------------------------------------
# 5. Site-count stacked barplot
# ---------------------------------------------------------------------------
# Classify MSL2 sites by HAS/PionX overlap
has_overlap   <- countOverlaps(msl2_gr, has_gr) > 0
pionx_overlap <- countOverlaps(msl2_gr, pionx_gr) > 0
overlap_cat <- case_when(
  has_overlap & pionx_overlap ~ "Both",
  has_overlap                 ~ "HAS_only",
  pionx_overlap               ~ "PionX_only",
  TRUE                        ~ "Novel"
)
site_class_counts <- as.data.frame(table(Category = overlap_cat),
                                    responseName = "Count")

p_stacked <- ggplot(site_class_counts,
                    aes(x = "MSL2 Sites", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5), colour = "white", size = 5) +
  scale_fill_manual(values = c(Novel = "#377eb8", HAS_only = "#ff7f00",
                               PionX_only = "#4daf4a", Both = "#984ea3")) +
  labs(y = "Number of Sites", x = NULL, fill = "Category") +
  theme_minimal(base_size = 14)

ggsave(file.path(output_dir, "stacked_bar_site_counts.png"),
       p_stacked, width = 6, height = 5, dpi = 300)

# ---------------------------------------------------------------------------
# 6. Genomic feature annotation (ChIPseeker)
# ---------------------------------------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

all_anno_gr <- list(MSL2_Sites = msl2_gr, HAS = has_gr, PionX = pionx_gr)
all_anno <- lapply(all_anno_gr, function(gr) {
  annotatePeak(gr, TxDb = txdb, tssRegion = c(-1000, 1000))
})

# Normalise annotation labels
normalise_anno <- function(ann_str) {
  case_when(
    grepl("Promoter", ann_str)           ~ "promoter",
    grepl("5' UTR|5'UTR", ann_str)       ~ "5'UTR",
    grepl("Exon|exon", ann_str)          ~ "exon",
    grepl("Intron|intron", ann_str)      ~ "intron",
    grepl("3' UTR|3'UTR", ann_str)       ~ "3'UTR",
    grepl("Downstream", ann_str)         ~ "downstream",
    grepl("Intergenic|Distal", ann_str)  ~ "intergenic",
    TRUE                                 ~ "other"
  )
}

anno_summary <- map_dfr(names(all_anno), function(nm) {
  ann <- as.data.frame(all_anno[[nm]])
  tibble(Category = nm, Annotation = normalise_anno(ann$annotation))
})

anno_counts <- anno_summary %>%
  group_by(Category, Annotation) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(Percent = 100 * Count / sum(Count)) %>%
  ungroup()

anno_counts$Annotation <- factor(anno_counts$Annotation,
  levels = c("promoter", "5'UTR", "exon", "intron", "3'UTR",
             "downstream", "intergenic", "other"))
anno_counts$Category <- factor(anno_counts$Category,
  levels = c("MSL2_Sites", "HAS", "PionX"))

feature_colors <- c(
  promoter    = "#e63946", `5'UTR`     = "#2a9d8f",
  exon        = "#457b9d", intron      = "#f4a261",
  `3'UTR`     = "#e9c46a", downstream  = "#a8dadc",
  intergenic  = "#264653", other       = "#999999"
)

p_anno <- ggplot(anno_counts,
                 aes(x = Category, y = Percent, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6,
           colour = "black") +
  geom_text(aes(label = ifelse(Percent > 5,
                               paste0(round(Percent, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = feature_colors) +
  labs(y = "Percentage of Sites", fill = "Genomic Feature") +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(output_dir, "genomic_feature_annotation_stacked.png"),
       p_anno, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "genomic_feature_annotation_stacked.pdf"),
       p_anno, width = 10, height = 6)

write_csv(anno_counts, file.path(output_dir, "genomic_feature_summary_stats.csv"))

# ---------------------------------------------------------------------------
# 7. Signal enrichment by site class (violin plots)
# ---------------------------------------------------------------------------
# Requires WT ratio bigwigs
msl2_bw_path    <- "bigwigs/sorted_s2a_wt_msl2_vs_ppi_merged_ratio_mean.bigwig"
h4k16ac_bw_path <- "bigwigs/sorted_s2a_wt_h4k16ac_vs_igg_merged_ratio_mean.bigwig"

get_bw_signal <- function(bed_gr, bw_path) {
  if (!file.exists(bw_path)) { warning("Missing: ", bw_path); return(rep(NA, length(bed_gr))) }
  bw <- import(bw_path, which = bed_gr)
  sapply(seq_along(bed_gr), function(i) {
    ov <- subsetByOverlaps(bw, bed_gr[i])
    if (length(ov) == 0) NA_real_ else mean(score(ov))
  })
}

site_gr <- c(has_gr, pionx_gr, msl2_gr)
site_cat <- c(rep("HAS", length(has_gr)),
              rep("PionX", length(pionx_gr)),
              rep("MSL2", length(msl2_gr)))

enrich_df <- tibble(
  category           = factor(site_cat, levels = c("MSL2", "HAS", "PionX")),
  MSL2_enrichment    = get_bw_signal(site_gr, msl2_bw_path),
  H4K16ac_enrichment = get_bw_signal(site_gr, h4k16ac_bw_path)
)

enrich_long <- enrich_df %>%
  pivot_longer(cols = c(MSL2_enrichment, H4K16ac_enrichment),
               names_to = "target", values_to = "signal") %>%
  mutate(target = recode(target,
    MSL2_enrichment = "MSL2", H4K16ac_enrichment = "H4K16ac"))

p_enrich <- ggplot(enrich_long, aes(x = category, y = signal, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ target, nrow = 1) +
  scale_fill_manual(values = c(MSL2 = "#1b7837", HAS = "#053061",
                               PionX = "#b2182b")) +
  labs(x = "Site Category", y = "Signal") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "enrichment_violin_by_site_class.png"),
       p_enrich, width = 10, height = 5, dpi = 300)
ggsave(file.path(output_dir, "enrichment_violin_by_site_class.pdf"),
       p_enrich, width = 10, height = 5)

# ---------------------------------------------------------------------------
# 8. Classify MSL2 sites and enrichment by overlap status
# ---------------------------------------------------------------------------
msl2_annot <- tibble(
  overlap = factor(overlap_cat, levels = c("Both", "HAS_only", "PionX_only", "Novel")),
  MSL2_enrichment    = get_bw_signal(msl2_gr, msl2_bw_path),
  H4K16ac_enrichment = get_bw_signal(msl2_gr, h4k16ac_bw_path)
)

msl2_long <- msl2_annot %>%
  pivot_longer(cols = c(MSL2_enrichment, H4K16ac_enrichment),
               names_to = "target", values_to = "signal") %>%
  mutate(target = recode(target,
    MSL2_enrichment = "MSL2", H4K16ac_enrichment = "H4K16ac"))

p_overlap_enrich <- ggplot(msl2_long, aes(x = overlap, y = signal, fill = overlap)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ target, nrow = 1) +
  scale_fill_manual(values = c(Both = "#1b7837", HAS_only = "#053061",
                               PionX_only = "#b2182b", Novel = "#cccccc")) +
  labs(x = "Overlap with HAS / PionX", y = "Signal") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "overlap_enrichment_violin.png"),
       p_overlap_enrich, width = 12, height = 5, dpi = 300)
ggsave(file.path(output_dir, "overlap_enrichment_violin.pdf"),
       p_overlap_enrich, width = 12, height = 5)

cat("\n✅ Site dissection analysis complete.\n")
cat("   Outputs in:", output_dir, "\n")

# ===========================================================================
# 7. Prepare sequences for STREME motif discovery (motif_discovery.sh)
#    Extracts 200 bp summit-centred sequences and generates matched-
#    chromosome background (seed = 42) as described in the Methods.
# ===========================================================================
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(Biostrings)

motif_dir <- "motif_sequences"
dir.create(motif_dir, recursive = TRUE, showWarnings = FALSE)

dm6 <- BSgenome.Dmelanogaster.UCSC.dm6

export_motif_seqs <- function(gr, label, genome, outdir, seed = 42) {
  # Resize to 200 bp centred on midpoint (summit proxy)
  gr_200 <- resize(gr, width = 200, fix = "center")
  gr_200 <- trim(gr_200)  # clip to chromosome boundaries

  # Extract target sequences
  target_seq <- getSeq(genome, gr_200)
  names(target_seq) <- paste0(label, "_", seq_along(gr_200))
  writeXStringSet(target_seq, file.path(outdir, paste0("sites_", label, ".fa")))

  # Matched-chromosome background (same widths, same chromosomes)
  set.seed(seed)
  chr_pool <- as.character(seqnames(gr_200))
  widths   <- width(gr_200)
  bg_list  <- vector("list", length(gr_200))

  for (i in seq_along(gr_200)) {
    chr   <- chr_pool[i]
    chr_len <- seqlengths(genome)[chr]
    w     <- widths[i]
    # Sample random start, reject overlaps with target sites
    repeat {
      start_pos <- sample.int(chr_len - w, 1)
      bg_gr <- GRanges(chr, IRanges(start_pos, width = w))
      if (length(findOverlaps(bg_gr, gr_200, type = "any")) == 0) break
    }
    bg_list[[i]] <- bg_gr
  }
  bg_gr_all <- do.call(c, bg_list)
  bg_seq <- getSeq(genome, bg_gr_all)
  names(bg_seq) <- paste0("bg_", label, "_", seq_along(bg_gr_all))
  writeXStringSet(bg_seq, file.path(outdir, paste0("background_", label, ".fa")))

  cat(sprintf("  %s: %d target + %d background sequences\n",
              label, length(target_seq), length(bg_seq)))
}

# Export for each site category
# overlap_cat is a character vector aligned with msl2_gr (from section 5)
mcols(msl2_gr)$overlap <- overlap_cat

site_categories <- list(
  MCCS  = msl2_gr[mcols(msl2_gr)$overlap == "Both"],
  HAS   = msl2_gr[mcols(msl2_gr)$overlap == "HAS_only"],
  PionX = msl2_gr[mcols(msl2_gr)$overlap == "PionX_only"],
  Novel = msl2_gr[mcols(msl2_gr)$overlap == "Novel"],
  All   = msl2_gr
)

cat("\nExporting motif FASTA sequences...\n")
for (nm in names(site_categories)) {
  gr <- site_categories[[nm]]
  if (length(gr) < 10) {
    cat(sprintf("  Skipping %s: only %d sites\n", nm, length(gr)))
    next
  }
  export_motif_seqs(gr, nm, dm6, motif_dir)
}
cat("✅ Motif sequences written to", motif_dir, "\n")
