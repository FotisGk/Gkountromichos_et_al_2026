#!/usr/bin/env bash
# =============================================================================
# Heatmaps and metaplots at csaw-identified binding sites
# Produces reference-point matrices centred on csaw merged sites, PionX sites,
# and HAS sites, then creates heatmaps and profile plots.
#
# Figures: SFig 2B–C (MSL2/H4K16ac signal at binding sites and gene bodies)
# =============================================================================
set -euo pipefail
module load ngs/deeptools  2>/dev/null || true

# --- EDIT THESE ---------------------------------------------------------------
# Merged ratio bigwigs (IP/Input, replicate-mean)
MSL2_BWS="sorted_s2a_wt_msl2_vs_ppi_merged_ratio_mean.bigwig \
          sorted_s2a_rox2ko112_msl2_vs_ppi_merged_ratio_mean.bigwig \
          sorted_s2a_rox2ko17_msl2_vs_ppi_merged_ratio_mean.bigwig \
          sorted_s2a_rox2ko17rox2fl_msl2_vs_ppi_merged_ratio_mean.bigwig"

H4K16AC_BWS="sorted_s2a_wt_h4k16ac_vs_igg_merged_ratio_mean.bigwig \
             sorted_s2a_rox2ko112_h4k16ac_vs_igg_merged_ratio_mean.bigwig \
             sorted_s2a_rox2ko17_h4k16ac_vs_igg_merged_ratio_mean.bigwig \
             sorted_s2a_rox2ko17rox2fl_h4k16ac_vs_igg_merged_ratio_mean.bigwig"

SAMPLE_LABELS="WT rox2ko112 rox2ko17 Rescue"

# BED files (outputs of 01_csaw_site_discovery)
CSAW_MSL2_BED="merged_top_sites_QLF_msl2.bed"
PIONX_BED="pionx_sites.bed"
HAS_BED="has_sites.bed"
# ------------------------------------------------------------------------------

FLANK=2000

# --- A. All csaw MSL2 sites (reference-point, centre) ---
# Target-specific colour maps & z-score ranges (from original analysis scripts)
#   MSL2:    YlGnBu, 0–4
#   H4K16ac: BuPu,   0–8
for TARGET in msl2 h4k16ac; do
  BWS_VAR="${TARGET^^}_BWS"   # MSL2_BWS or H4K16AC_BWS
  eval BWS=\$$BWS_VAR 2>/dev/null || BWS="${MSL2_BWS}"

  if [ "$TARGET" = "msl2" ]; then CMAP=YlGnBu; ZMAX=4; else CMAP=BuPu; ZMAX=8; fi

  computeMatrix reference-point \
    -S ${BWS} \
    -R "${CSAW_MSL2_BED}" \
    --referencePoint center -a ${FLANK} -b ${FLANK} \
    --skipZeros --missingDataAsZero \
    -o "csaw_center_ref_${TARGET}.matrix.mat.gz"

  plotHeatmap \
    -m "csaw_center_ref_${TARGET}.matrix.mat.gz" \
    --sortRegions descend --sortUsing mean --sortUsingSamples 1 \
    --heatmapWidth 6 --colorMap ${CMAP} \
    --zMin 0 --zMax ${ZMAX} \
    --samplesLabel ${SAMPLE_LABELS} \
    --refPointLabel "Center" --dpi 300 \
    -out "csaw_center_ref_${TARGET}.png"
done

# --- B. PionX / HAS / csaw site classes (multi-region) ---
for TARGET in msl2 h4k16ac; do
  BWS_VAR="${TARGET^^}_BWS"
  eval BWS=\$$BWS_VAR 2>/dev/null || BWS="${MSL2_BWS}"

  if [ "$TARGET" = "msl2" ]; then CMAP=YlGnBu; ZMAX=4; else CMAP=BuPu; ZMAX=8; fi

  computeMatrix reference-point \
    -S ${BWS} \
    -R "${PIONX_BED}" "${HAS_BED}" "${CSAW_MSL2_BED}" \
    --referencePoint center -a ${FLANK} -b ${FLANK} \
    --skipZeros --missingDataAsZero \
    -o "csaw_has_pionx_ref_${TARGET}.matrix.mat.gz"

  plotHeatmap \
    -m "csaw_has_pionx_ref_${TARGET}.matrix.mat.gz" \
    --sortRegions descend --sortUsing mean --sortUsingSamples 1 \
    --heatmapWidth 6 --heatmapHeight 20 --colorMap ${CMAP} \
    --zMin 0 --zMax ${ZMAX} \
    --regionsLabel "PionX sites" "HAS sites" "MSL2 csaw sites" \
    --samplesLabel ${SAMPLE_LABELS} \
    --refPointLabel "Center" --interpolationMethod bilinear --dpi 300 \
    -out "csaw_has_pionx_ref_${TARGET}.png"

  # Profile / metaplot
  plotProfile \
    -m "csaw_has_pionx_ref_${TARGET}.matrix.mat.gz" \
    --colors "#67001f" "#08306b" "#00441b" \
    --regionsLabel "PionX sites" "HAS sites" "MSL2 csaw sites" \
    --samplesLabel ${SAMPLE_LABELS} \
    --refPointLabel "Center" --yAxisLabel "Mean signal" --dpi 300 \
    -out "csaw_has_pionx_ref_${TARGET}_metaplot.pdf"
done

echo "✅ csaw site heatmaps and metaplots complete."
