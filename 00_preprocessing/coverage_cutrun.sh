#!/usr/bin/env bash
# =============================================================================
# CUT&RUN coverage tracks via deeptools (BAM → bigWig)
#
# Generates:
#   1. Per-replicate normalised coverage bigWigs (custom scale factor to 1M)
#   2. Per-replicate IP/Control ratio bigWigs (bamCompare)
#   3. Per-replicate log2(IP/Control) bigWigs (bamCompare)
#
# Uses a custom normalisation approach: each BAM is scaled to 1,000,000 reads
# (reads counted after filtering: -F 1804 -q 30) with assay-specific smoothing 
# (MSL2: 20 bp, H4K16ac: 80 bp).
#
# Input:  sorted_s2a_<condition>_<target>_rep<N>.bam  (in current directory)
# Output: deeptools_output_v3/<name>.bigwig
#
# Prerequisites: deeptools ≥ 3.5, samtools ≥ 1.9, bc
# Figures: Inputs for all CUT&RUN heatmaps/metaplots (Figs 3–4, Sup. Fig. 2A)
# =============================================================================
set -euo pipefail

module load ngs/deeptools  2>/dev/null || true
module load ngs/samtools     2>/dev/null || true

# ---- Parameters -------------------------------------------------------------
OUTDIR="deeptools_output_v3"
mkdir -p "${OUTDIR}"

THREADS=8
BIN_SIZE=10
TARGET_READS=1000000         # Scale each BAM to this many reads

# Smoothing lengths (bp) — assay-specific
MSL2_SMOOTH=20
H4K16AC_SMOOTH=80
DEFAULT_SMOOTH=0             # 0 = no smoothing

# bamCompare settings
COMPARE_PSEUDO=1             # Pseudocount for ratio/log2 comparisons
HIST_BIN_SIZE=50             # Larger bin size for histone comparisons

SCALING_LOG="${OUTDIR}/scaling_factors.txt"
echo -e "Sample\tTotal_Reads\tScaling_Factor" > "${SCALING_LOG}"

# ---- Helper: determine smoothing length from filename -----------------------
get_smooth_length() {
  local bam_name="$1"
  local target=$(echo "${bam_name}" | sed -E 's/.*_([^_]+)_rep[0-9]+\.bam$/\1/')
  case "${target}" in
    msl2|ppi)    echo "${MSL2_SMOOTH}" ;;
    h4k16ac|igg) echo "${H4K16AC_SMOOTH}" ;;
    *)           echo "${DEFAULT_SMOOTH}" ;;
  esac
}

# =============================================================================
# STEP 1: Per-replicate coverage bigWigs (normalised to TARGET_READS)
# =============================================================================
echo "=== Step 1: Per-replicate coverage bigWigs ==="

for BAM in sorted_s2a_*.bam; do
  [[ ! -f "${BAM}" ]] && continue
  OUT_BW="${OUTDIR}/${BAM%.bam}.bigwig"
  [[ -f "${OUT_BW}" ]] && { echo "  skip (exists): ${OUT_BW}"; continue; }

  # Ensure BAM index exists
  if [[ ! -f "${BAM}.bai" || "${BAM}.bai" -ot "${BAM}" ]]; then
    samtools index -f "${BAM}"
  fi

  # Count quality-filtered reads
  TOTAL_READS=$(samtools view -c -F 1804 -q 30 "${BAM}")
  if [[ "${TOTAL_READS}" -eq 0 ]]; then
    echo "  WARN: ${BAM} has 0 reads after filtering — skipping"
    continue
  fi

  SCALE_FACTOR=$(echo "scale=6; ${TARGET_READS} / ${TOTAL_READS}" | bc -l)
  echo -e "${BAM%.bam}\t${TOTAL_READS}\t${SCALE_FACTOR}" >> "${SCALING_LOG}"

  SMOOTH_LEN=$(get_smooth_length "${BAM}")
  BC_EXTRA=()
  [[ "${SMOOTH_LEN}" -gt 0 ]] && BC_EXTRA+=(--smoothLength "${SMOOTH_LEN}")
  BC_EXTRA+=(--skipNAs)

  echo "  bamCoverage: ${BAM}  (${TOTAL_READS} reads, scale=${SCALE_FACTOR}, smooth=${SMOOTH_LEN})"
  bamCoverage \
    --bam "${BAM}" \
    --outFileName "${OUT_BW}" \
    --binSize "${BIN_SIZE}" \
    --normalizeUsing None \
    --scaleFactor "${SCALE_FACTOR}" \
    --ignoreDuplicates \
    --minMappingQuality 30 \
    "${BC_EXTRA[@]}" \
    --numberOfProcessors "${THREADS}" \
    || { echo "  FAIL: bamCoverage for ${BAM}"; continue; }
done

# =============================================================================
# STEP 2: IP vs Control pairwise comparisons (bamCompare)
# =============================================================================
echo ""
echo "=== Step 2: IP vs Control comparisons (bamCompare) ==="

CONDITIONS="wt rox2ko17 rox2ko112 rox2ko17rox2fl"

for condition in ${CONDITIONS}; do

  # ----- MSL2 vs PPI (up to 6 replicates) -----
  for rep in $(seq 1 6); do
    msl2_bam="sorted_s2a_${condition}_msl2_rep${rep}.bam"
    ppi_bam="sorted_s2a_${condition}_ppi_rep${rep}.bam"
    [[ ! -f "${msl2_bam}" || ! -f "${ppi_bam}" ]] && continue

    out_ratio="${OUTDIR}/sorted_s2a_${condition}_msl2_vs_ppi_rep${rep}_ratio.bigwig"
    out_log2="${OUTDIR}/sorted_s2a_${condition}_msl2_vs_ppi_rep${rep}_log2.bigwig"

    # Count reads for per-pair scale factors
    TOTAL1=$(samtools view -c -F 1804 -q 30 "${msl2_bam}")
    TOTAL2=$(samtools view -c -F 1804 -q 30 "${ppi_bam}")
    [[ "${TOTAL1}" -eq 0 || "${TOTAL2}" -eq 0 ]] && continue
    SCALE1=$(echo "scale=6; ${TARGET_READS} / ${TOTAL1}" | bc -l)
    SCALE2=$(echo "scale=6; ${TARGET_READS} / ${TOTAL2}" | bc -l)

    echo "  MSL2 ${condition} rep${rep}: scaleFactors ${SCALE1}:${SCALE2}"

    # Ratio bigwig
    [[ ! -f "${out_ratio}" ]] && \
    bamCompare \
      --bamfile1 "${msl2_bam}" --bamfile2 "${ppi_bam}" \
      --outFileName "${out_ratio}" \
      --scaleFactors "${SCALE1}:${SCALE2}" \
      --ratio ratio \
      --binSize "${BIN_SIZE}" \
      --smoothLength "${MSL2_SMOOTH}" \
      --pseudocount "${COMPARE_PSEUDO}" \
      --skipZeroOverZero \
      --numberOfProcessors "${THREADS}" \
      || echo "  FAIL: bamCompare ratio for ${msl2_bam}"

    # Log2 ratio bigwig
    [[ ! -f "${out_log2}" ]] && \
    bamCompare \
      --bamfile1 "${msl2_bam}" --bamfile2 "${ppi_bam}" \
      --outFileName "${out_log2}" \
      --scaleFactors "${SCALE1}:${SCALE2}" \
      --ratio log2 \
      --binSize "${BIN_SIZE}" \
      --smoothLength "${MSL2_SMOOTH}" \
      --pseudocount "${COMPARE_PSEUDO}" \
      --skipZeroOverZero \
      --numberOfProcessors "${THREADS}" \
      || echo "  FAIL: bamCompare log2 for ${msl2_bam}"
  done

  # ----- H4K16ac vs IgG (up to 5 replicates) -----
  for rep in $(seq 1 5); do
    h4_bam="sorted_s2a_${condition}_h4k16ac_rep${rep}.bam"
    igg_bam="sorted_s2a_${condition}_igg_rep${rep}.bam"
    [[ ! -f "${h4_bam}" || ! -f "${igg_bam}" ]] && continue

    out_ratio="${OUTDIR}/sorted_s2a_${condition}_h4k16ac_vs_igg_rep${rep}_ratio.bigwig"
    out_log2="${OUTDIR}/sorted_s2a_${condition}_h4k16ac_vs_igg_rep${rep}_log2.bigwig"

    TOTAL1=$(samtools view -c -F 1804 -q 30 "${h4_bam}")
    TOTAL2=$(samtools view -c -F 1804 -q 30 "${igg_bam}")
    [[ "${TOTAL1}" -eq 0 || "${TOTAL2}" -eq 0 ]] && continue
    SCALE1=$(echo "scale=6; ${TARGET_READS} / ${TOTAL1}" | bc -l)
    SCALE2=$(echo "scale=6; ${TARGET_READS} / ${TOTAL2}" | bc -l)

    echo "  H4K16ac ${condition} rep${rep}: scaleFactors ${SCALE1}:${SCALE2}"

    [[ ! -f "${out_ratio}" ]] && \
    bamCompare \
      --bamfile1 "${h4_bam}" --bamfile2 "${igg_bam}" \
      --outFileName "${out_ratio}" \
      --scaleFactors "${SCALE1}:${SCALE2}" \
      --ratio ratio \
      --binSize "${HIST_BIN_SIZE}" \
      --smoothLength "${H4K16AC_SMOOTH}" \
      --pseudocount "${COMPARE_PSEUDO}" \
      --skipZeroOverZero \
      --numberOfProcessors "${THREADS}" \
      || echo "  FAIL: bamCompare ratio for ${h4_bam}"

    [[ ! -f "${out_log2}" ]] && \
    bamCompare \
      --bamfile1 "${h4_bam}" --bamfile2 "${igg_bam}" \
      --outFileName "${out_log2}" \
      --scaleFactors "${SCALE1}:${SCALE2}" \
      --ratio log2 \
      --binSize "${HIST_BIN_SIZE}" \
      --smoothLength "${H4K16AC_SMOOTH}" \
      --pseudocount "${COMPARE_PSEUDO}" \
      --skipZeroOverZero \
      --numberOfProcessors "${THREADS}" \
      || echo "  FAIL: bamCompare log2 for ${h4_bam}"
  done

done

echo ""
echo "✅ CUT&RUN coverage and comparison tracks complete. Results in ${OUTDIR}/"
