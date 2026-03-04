#!/usr/bin/env bash
# =============================================================================
# CUT&RUN paired-end alignment (bowtie2 → BAM → BED)
#
# Aligns paired-end CUT&RUN FASTQ files to dm6 (Drosophila melanogaster).
# Spike-in alignment to E. coli is optional (commented out by default).
#
# Input:  samples/read1_*.fastq.gz.fastqsanger.gz  (+ matching read2_*)
# Output: bowtie_output/<sample>.bam, <sample>.bed, <sample>.stats
#
# Usage on SLURM:
#   bash align_cutrun.sh          # submits an array job (one per sample)
#
# Usage locally (single sample):
#   bash align_cutrun.sh --local samples/read1_SAMPLE.fastq.gz.fastqsanger.gz
#
# Figures: Upstream of all CUT&RUN analyses (Figs 3–4, SFig 2)
# =============================================================================
set -euo pipefail

# ---- Configuration ----------------------------------------------------------
BOWTIE_INDEX="/work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
BOWTIE_OPTS="-p 12 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -I 10 -X 700"
OUTPUT_DIR="bowtie_output"
# Uncomment for E. coli spike-in alignment:
# BOWTIE_INDEX_SPIKEIN="/work/project/becbec_008/resources/bacteria/bowtie2_index/ecoli"

# ---- Launcher (submit SLURM array or run locally) ---------------------------
if [[ "${1:-}" != "--local" ]]; then
  FILES=(samples/read1_*.fastq.gz.fastqsanger.gz)
  NUMFASTQ=${#FILES[@]}
  mkdir -p "${OUTPUT_DIR}"
  if [[ ${NUMFASTQ} -gt 0 ]]; then
    echo "Submitting array job for ${NUMFASTQ} CUT&RUN FASTQ pairs."
    sbatch --array=1-${NUMFASTQ} align_cutrun.sbatch
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

echo "Processing CUT&RUN pair: ${FILENAME_R1}  +  ${FILENAME_R2}"

module load ngs/bowtie2      2>/dev/null || true
module load ngs/samtools     2>/dev/null || true
module load ngs/bedtools2    2>/dev/null || true
module load ngs/Homer        2>/dev/null || true

mkdir -p "${OUTPUT_DIR}"

# 1. Align to dm6
bowtie2 ${BOWTIE_OPTS} -x "${BOWTIE_INDEX}" \
  -1 "${FILENAME_R1}" -2 "${FILENAME_R2}" \
  -S "${OUTPUT_DIR}/${FILEBASE}.sam" \
  2> "${OUTPUT_DIR}/${FILEBASE}.stats"

# 2. SAM → name-sorted BAM → BEDPE → BED (fragment coordinates)
samtools view -h -@ 8 -q 2 -b "${OUTPUT_DIR}/${FILEBASE}.sam" \
  -o "${OUTPUT_DIR}/${FILEBASE}.bam"
samtools sort -n -@ 8 "${OUTPUT_DIR}/${FILEBASE}.bam" \
  -o "${OUTPUT_DIR}/${FILEBASE}.f.bam"
bamToBed -i "${OUTPUT_DIR}/${FILEBASE}.f.bam" -bedpe \
  | cut -f 1,2,6 | sort -k1,1 \
  > "${OUTPUT_DIR}/${FILEBASE}.bed" 2>/dev/null

rm -f "${OUTPUT_DIR}/${FILEBASE}.sam"

# 3. (Optional) Spike-in alignment — uncomment if using E. coli spike-in
# bowtie2 ${BOWTIE_OPTS} -x "${BOWTIE_INDEX_SPIKEIN}" \
#   -1 "${FILENAME_R1}" -2 "${FILENAME_R2}" \
#   -S "${OUTPUT_DIR}/${FILEBASE}.spikein.sam" \
#   2> "${OUTPUT_DIR}/${FILEBASE}.spikein.stats"
# samtools view -h -@ 8 -q 2 -b "${OUTPUT_DIR}/${FILEBASE}.spikein.sam" \
#   -o "${OUTPUT_DIR}/${FILEBASE}.spikein.bam"
# samtools sort -n -@ 8 "${OUTPUT_DIR}/${FILEBASE}.spikein.bam" \
#   -o "${OUTPUT_DIR}/${FILEBASE}.spikein.f.bam"
# bamToBed -i "${OUTPUT_DIR}/${FILEBASE}.spikein.f.bam" -bedpe \
#   | cut -f 1,2,6 > "${OUTPUT_DIR}/${FILEBASE}.spikein.bed" 2>/dev/null
# rm -f "${OUTPUT_DIR}/${FILEBASE}.spikein.sam"

echo "Done: ${FILEBASE}"
