#!/usr/bin/env bash
# =============================================================================
# ChIP-seq spike-in normalised coverage tracks (BED → bedGraph → bigWig)
#
# Uses D. virilis spike-in read counts to scale ChIP-seq BED fragment files
# into normalised bigWig coverage tracks. Also produces unnormalised bigWigs.
# Processes both MSL2 and H4K16ac samples against matched Input controls.
#
#
# Input:
#   - <sample>.bed       (dm6 fragment BED from align_chip.sh)
#   - <sample>.spikein.bed (D. virilis fragment BED from align_chip.sh)
#   - <input>.bed / <input>.spikein.bed (matched Input control)
# Output:
#   - <sample>.spike.bw   (spike-in normalised bigWig)
#   - <sample>.bw         (unnormalised bigWig)
#
# Usage on SLURM:
#   bash coverage_chip.sh          # submits sbatch job
#
# Usage locally:
#   bash coverage_chip.sh --local
#
# Figures: Upstream of ChIP heatmaps/metaplots (Figs 3–4)
# =============================================================================
set -euo pipefail

# ---- Configuration ----------------------------------------------------------
CHROMSIZES="/work/project/becbec_006/ref/dm6.chrom.sizes.txt"
OUTDIR="outcoverage"

# ---- Launcher ---------------------------------------------------------------
if [[ "${1:-}" != "--local" ]]; then
  # Auto-detect Input and Sample BED files
  INPUT=$(ls *.bed | grep "Input" | grep -v "spikein")
  SAMPLES_MSL2=$(ls *.bed | grep "MSL2" | grep -v "spikein" || true)
  SAMPLES_H4K16AC=$(ls *.bed | grep -i "H4K16ac" | grep -v "spikein" || true)

  INPUTBASE=$(echo "${INPUT}" | sed -e 's/.bed//g')

  echo "Input: ${INPUTBASE}"

  mkdir -p "${OUTDIR}"

  # Submit MSL2 samples
  for S in ${SAMPLES_MSL2}; do
    SBASE=$(echo "${S}" | sed -e 's/.bed//g')
    echo "Submitting MSL2 sample: ${SBASE}"
    sbatch --export=FILEBASE="${SBASE}",INPUTBASE="${INPUTBASE}" \
      coverage_chip.sbatch
  done

  # Submit H4K16ac samples
  for S in ${SAMPLES_H4K16AC}; do
    SBASE=$(echo "${S}" | sed -e 's/.bed//g')
    echo "Submitting H4K16ac sample: ${SBASE}"
    sbatch --export=FILEBASE="${SBASE}",INPUTBASE="${INPUTBASE}" \
      coverage_chip.sbatch
  done

  exit 0
fi

# ---- Worker (called by sbatch or --local) -----------------------------------
module load ngs/samtools      2>/dev/null || true
module load ngs/Homer         2>/dev/null || true
module load ngs/UCSCutils         2>/dev/null || true
module load ngs/bedtools2  2>/dev/null || true

mkdir -p "${OUTDIR}"

# Count spike-in reads
INPUT_DEPTH=$(wc -l < "${INPUTBASE}.spikein.bed")
SAMPLE_DEPTH=$(wc -l < "${FILEBASE}.spikein.bed")

echo "Input  spike-in reads: ${INPUT_DEPTH}"
echo "Sample spike-in reads: ${SAMPLE_DEPTH}"

# Calculate spike-in scale factors (reads per million)
INPUT_SCALE=$(echo "1000000/${INPUT_DEPTH}" | bc -l)
SAMPLE_SCALE=$(echo "1000000/${SAMPLE_DEPTH}" | bc -l)

echo "Input  scale factor: ${INPUT_SCALE}"
echo "Sample scale factor: ${SAMPLE_SCALE}"

# --- Spike-in normalised bigWigs ---
bedtools genomecov -bg -scale "${INPUT_SCALE}" \
  -i "${INPUTBASE}.bed" -g "${CHROMSIZES}" \
  > "${INPUTBASE}.spike.bedgraph"
bedGraphToBigWig "${INPUTBASE}.spike.bedgraph" "${CHROMSIZES}" \
  "${INPUTBASE}.spike.bw"

bedtools genomecov -bg -scale "${SAMPLE_SCALE}" \
  -i "${FILEBASE}.bed" -g "${CHROMSIZES}" \
  > "${FILEBASE}.spike.bedgraph"
bedGraphToBigWig "${FILEBASE}.spike.bedgraph" "${CHROMSIZES}" \
  "${FILEBASE}.spike.bw"

# --- Unnormalised bigWigs (for comparison) ---
bedtools genomecov -bg -i "${INPUTBASE}.bed" -g "${CHROMSIZES}" \
  > "${INPUTBASE}.bedgraph"
bedGraphToBigWig "${INPUTBASE}.bedgraph" "${CHROMSIZES}" "${INPUTBASE}.bw"

bedtools genomecov -bg -i "${FILEBASE}.bed" -g "${CHROMSIZES}" \
  > "${FILEBASE}.bedgraph"
bedGraphToBigWig "${FILEBASE}.bedgraph" "${CHROMSIZES}" "${FILEBASE}.bw"

echo "Done: spike-in normalised + unnormalised bigWigs for ${FILEBASE}"
