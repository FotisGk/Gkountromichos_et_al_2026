#!/usr/bin/env bash
# =============================================================================
# Gene body heatmaps and metaplots (scale-regions mode)
# Consolidated from 6+ individual scripts. Generates computeMatrix + plotHeatmap
# for X-linked genes, autosomal genes, all genes, and MSL2-bound genes, using
# both CUT&RUN and ChIP-MNase bigwigs.
#
# Figures: 4A–D (gene body heatmaps), 4E–F (metaplots)
# Supplementary: density comparison plots
# =============================================================================
set -euo pipefail
module load ngs/deeptools  2>/dev/null || true

# --- EDIT THESE ---------------------------------------------------------------
# CUT&RUN merged ratio bigwigs
CNR_MSL2_WT="sorted_s2a_wt_msl2_vs_ppi_merged_ratio_mean.bigwig"
CNR_MSL2_ALL="${CNR_MSL2_WT} \
  sorted_s2a_rox2ko112_msl2_vs_ppi_merged_ratio_mean.bigwig \
  sorted_s2a_rox2ko17_msl2_vs_ppi_merged_ratio_mean.bigwig \
  sorted_s2a_rox2ko17rox2fl_msl2_vs_ppi_merged_ratio_mean.bigwig"

CNR_H4K16AC_WT="sorted_s2a_wt_h4k16ac_vs_igg_merged_ratio_mean.bigwig"
CNR_H4K16AC_ALL="${CNR_H4K16AC_WT} \
  sorted_s2a_rox2ko112_h4k16ac_vs_igg_merged_ratio_mean.bigwig \
  sorted_s2a_rox2ko17_h4k16ac_vs_igg_merged_ratio_mean.bigwig \
  sorted_s2a_rox2ko17rox2fl_h4k16ac_vs_igg_merged_ratio_mean.bigwig"

# ChIP-MNase merged ratio bigwigs (WT only)
CHIP_MSL2_WT="s2a_wt_msl2_vs_input_merged_ratio_mean_chip.bigwig"
CHIP_H4K16AC_WT="s2a_wt_h4k16ac_vs_input_merged_ratio_mean_chip.bigwig"

SAMPLE_LABELS_ALL="WT rox2ko112 rox2ko17 Rescue"

# Region BED files
XGENES="genes_chrX.bed"
AGENES="genes_autosomes.bed"
ALLGENES="genes_all.bed"
MSL2GENES="genes_all_MSL2bound.bed"
# ------------------------------------------------------------------------------

BEFORE=2000; BODY=5000; AFTER=2000   # 2 kb flanks (manuscript + originals)

run_scale_regions() {
  local bws="$1" region="$2" prefix="$3" labels="$4"
  local mat="${prefix}.matrix.mat.gz"

  computeMatrix scale-regions -S ${bws} -R "${region}" \
    --beforeRegionStartLength ${BEFORE} --regionBodyLength ${BODY} \
    --afterRegionStartLength ${AFTER} \
    --skipZeros --missingDataAsZero -o "${mat}"

  plotHeatmap -m "${mat}" \
    --sortRegions descend --sortUsing mean --sortUsingSamples 1 \
    --heatmapWidth 6 --heatmapHeight 20 --colorMap YlGnBu \
    --samplesLabel ${labels} --interpolationMethod bilinear --dpi 300 \
    -out "${prefix}.png"

  plotProfile -m "${mat}" \
    --samplesLabel ${labels} --yAxisLabel "Mean signal" --dpi 300 \
    -out "${prefix}_metaplot.pdf"
}

# --- X-linked genes (all genotypes) ---
for TARGET in msl2 h4k16ac; do
  VAR="CNR_${TARGET^^}_ALL"
  eval BWS=\$$VAR
  run_scale_regions "${BWS}" "${XGENES}" \
    "Xgenes_scale-regions_${TARGET}" "${SAMPLE_LABELS_ALL}"
done

# --- Autosomal genes (all genotypes) ---
for TARGET in msl2 h4k16ac; do
  VAR="CNR_${TARGET^^}_ALL"
  eval BWS=\$$VAR
  run_scale_regions "${BWS}" "${AGENES}" \
    "Agenes_scale-regions_${TARGET}" "${SAMPLE_LABELS_ALL}"
done

# --- All genes (all genotypes) ---
for TARGET in msl2 h4k16ac; do
  VAR="CNR_${TARGET^^}_ALL"
  eval BWS=\$$VAR
  run_scale_regions "${BWS}" "${ALLGENES}" \
    "allgenes_scale-regions_${TARGET}" "${SAMPLE_LABELS_ALL}"
done

# --- MSL2-bound gene bodies: CUT&RUN (WT only) ---
run_scale_regions "${CNR_MSL2_WT}"    "${MSL2GENES}" "msl2genes_cnr_scale-regions"    "WT_CUT&RUN"
run_scale_regions "${CNR_H4K16AC_WT}" "${MSL2GENES}" "h4k16acgenes_cnr_scale-regions" "WT_CUT&RUN"

# --- MSL2-bound gene bodies: ChIP-MNase (WT only) ---
run_scale_regions "${CHIP_MSL2_WT}"    "${MSL2GENES}" "msl2genes_chip_scale-regions"    "WT_ChIP"
run_scale_regions "${CHIP_H4K16AC_WT}" "${MSL2GENES}" "h4k16acgenes_chip_scale-regions" "WT_ChIP"

echo "✅ Gene body profiles complete."
