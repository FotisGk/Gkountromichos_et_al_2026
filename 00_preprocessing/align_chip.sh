#!/usr/bin/env bash
# =============================================================================
# ChIP-seq paired-end alignment (bowtie2 → BAM → BED) with D. virilis spike-in
#
# Aligns paired-end ChIP-seq FASTQ files to dm6 (Drosophila melanogaster)
# AND to D. virilis (droVir3) for spike-in normalisation.
#
# Key differences from CUT&RUN alignment:
#   - Larger max insert size (-X 1000 vs -X 700)
#   - D. virilis spike-in alignment is always performed
#
# Input:  samples/read1_*.fastq.gz.fastqsanger.gz  (+ matching read2_*)
# Output: bowtie_output/<sample>.bam, <sample>.bed, <sample>.stats
#         bowtie_output/<sample>.spikein.bam, <sample>.spikein.bed
#
# Usage on SLURM:
#   bash align_chip.sh            # submits an array job (one per sample)
#
# Usage locally (single sample):
#   bash align_chip.sh --local samples/read1_SAMPLE.fastq.gz.fastqsanger.gz
#
# Figures: Upstream of ChIP signal tracks and spike-in normalised coverage
# =============================================================================
set -euo pipefail

# ---- Configuration ----------------------------------------------------------
BOWTIE_INDEX="/work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
BOWTIE_OPTS="-p 12 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 1000"

BOWTIE_INDEX_DVIR="/work/project/becbec_005/ChIP_Seq/Drosophila_virilis_UCSC_droVir3/bt2buildvir3/genome"
BOWTIE_OPTS_DVIR="-p 12 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 1000"

OUTPUT_DIR="bowtie_output"

# ---- Launcher (submit SLURM array or run locally) ---------------------------
if [[ "${1:-}" != "--local" ]]; then
  FILES=(samples/read1_*.fastq.gz.fastqsanger.gz)
  NUMFASTQ=${#FILES[@]}
  mkdir -p "${OUTPUT_DIR}"
  if [[ ${NUMFASTQ} -gt 0 ]]; then
    echo "Submitting array job for ${NUMFASTQ} ChIP FASTQ pairs."
    sbatch --array=1-${NUMFASTQ} align_chip.sbatch
  else
    echo "No FASTQ files found in samples/"
  fi
  exit 0
fi

# ---- Single-sample mode (called by sbatch or --local) -----------------------
FILENAME_R1="${2:-$(ls samples/read1_*.fastq.gz.fastqsanger.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")}"
FILEBASE=$(basename "${FILENAME_R1}" | sed 's/read1_//' | sed 's/.fastq.gz.fastqsanger.gz//')
FILENAME_R2="samples/read2_${FILEBASE}.fastq.gz.fastqsanger.gz"

for F in "${FILENAME_R1}" "${FILENAME_R2}"; do
  [[ ! -f "${F}" ]] && { echo "Error: file not found: ${F}"; exit 1; }
done

echo "Processing ChIP pair: ${FILENAME_R1}  +  ${FILENAME_R2}"

module load ngs/bowtie2      2>/dev/null || true
module load ngs/samtools     2>/dev/null || true
module load ngs/bedtools2    2>/dev/null || true
module load ngs/Homer        2>/dev/null || true

mkdir -p "${OUTPUT_DIR}"

# --- 1. Align to dm6 ---------------------------------------------------------
bowtie2 ${BOWTIE_OPTS} -x "${BOWTIE_INDEX}" \
  -1 "${FILENAME_R1}" -2 "${FILENAME_R2}" \
  -S "${OUTPUT_DIR}/${FILEBASE}.sam" \
  2> "${OUTPUT_DIR}/${FILEBASE}.stats"

samtools view -h -@ 8 -q 2 -b "${OUTPUT_DIR}/${FILEBASE}.sam" \
  -o "${OUTPUT_DIR}/${FILEBASE}.bam"
samtools sort -n -@ 8 "${OUTPUT_DIR}/${FILEBASE}.bam" \
  -o "${OUTPUT_DIR}/${FILEBASE}.f.bam"
bamToBed -i "${OUTPUT_DIR}/${FILEBASE}.f.bam" -bedpe \
  | cut -f 1,2,6 | sort -k1,1 \
  > "${OUTPUT_DIR}/${FILEBASE}.bed" 2>/dev/null

rm -f "${OUTPUT_DIR}/${FILEBASE}.sam"

# --- 2. Align to D. virilis spike-in genome ----------------------------------
bowtie2 ${BOWTIE_OPTS_DVIR} -x "${BOWTIE_INDEX_DVIR}" \
  -1 "${FILENAME_R1}" -2 "${FILENAME_R2}" \
  -S "${OUTPUT_DIR}/${FILEBASE}.spikein.sam" \
  2> "${OUTPUT_DIR}/${FILEBASE}.spikein.stats"

samtools view -h -@ 8 -q 2 -b "${OUTPUT_DIR}/${FILEBASE}.spikein.sam" \
  -o "${OUTPUT_DIR}/${FILEBASE}.spikein.bam"
samtools sort -n -@ 8 "${OUTPUT_DIR}/${FILEBASE}.spikein.bam" \
  -o "${OUTPUT_DIR}/${FILEBASE}.spikein.f.bam"
bamToBed -i "${OUTPUT_DIR}/${FILEBASE}.spikein.f.bam" -bedpe \
  | cut -f 1,2,6 \
  > "${OUTPUT_DIR}/${FILEBASE}.spikein.bed" 2>/dev/null

rm -f "${OUTPUT_DIR}/${FILEBASE}.spikein.sam"

echo "Done: ${FILEBASE}  (dm6 + D. virilis spike-in)"
