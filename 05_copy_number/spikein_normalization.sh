#!/usr/bin/env bash
#SBATCH -J chrCN_spike
#SBATCH -o slurm-%x-%A_%a.out
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 12:00:00
# =============================================================================
# Spike-in normalised copy number estimation from INPUT BAMs
# For each sample: counts spike-in reads → computes scale factor →
# runs bamCoverage with --scaleFactor → bedtools bins → CBS segmentation.
#
# Submit: sbatch --array=0-11 spikein_normalization.sh
# Figures: 5C (chromosome-level CN barplots, via copy_number_estimation.R)
# =============================================================================
set -euo pipefail

module load ngs/deeptools
module load ngs/samtools
module load ngs/bedtools2
module load ngs/UCSCutils

# --- EDIT THESE ---------------------------------------------------------------
CHRSIZES="dm6.chrom.sizes.txt"
GTF="ensembl_genes.clean.gtf"
BIN=50000
BASELINE_CN=4
OUTDIR="cn_per_cellline_spike"
SPIKE_CHR_PREFIX="scaffold_"

PREFIXES=(
  "s2a_wt_rep1" "s2a_wt_rep4" "s2a_wt_rep5"
  "s2a_rox2ko17_rep1" "s2a_rox2ko17_rep2" "s2a_rox2ko17_rep5"
  "s2a_rox2ko112_rep1" "s2a_rox2ko112_rep2" "s2a_rox2ko112_rep3"
  "s2a_rox2ko17rox2fl_rep1" "s2a_rox2ko17rox2fl_rep2" "s2a_rox2ko17rox2fl_rep3"
)

SPIKEIN_BAMS=(
  "samples_bam_spike-in/s2a_wt_input_rep1.spikein.f.bam"
  "samples_bam_spike-in/s2a_wt_input_rep4.spikein.f.bam"
  "samples_bam_spike-in/s2a_wt_input_rep5.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko17_input_rep1.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko17_input_rep2.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko17_input_rep5.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko112_input_rep1.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko112_input_rep2.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko112_input_rep3.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko17rox2fl_input_rep1.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko17rox2fl_input_rep2.spikein.f.bam"
  "samples_bam_spike-in/s2a_rox2ko17rox2fl_input_rep3.spikein.f.bam"
)
# ------------------------------------------------------------------------------

count_spike_reads() {
  local bam="$1"
  [[ ! -f "${bam}.bai" ]] && samtools index "$bam"
  local spike_contigs
  spike_contigs=$(samtools view -H "$bam" | grep "^@SQ" | awk '{print $2}' \
    | sed 's/SN://' | grep -E "^${SPIKE_CHR_PREFIX}") || true
  if [[ -z "$spike_contigs" ]]; then echo 0; return; fi
  samtools view -c "$bam" ${spike_contigs} 2>/dev/null || echo 0
}

# Phase 1: count spike-in reads
declare -a spike_counts scale_factors
total=${#PREFIXES[@]}
for (( i=0; i<total; i++ )); do
  spike_counts[$i]=$(count_spike_reads "${SPIKEIN_BAMS[$i]}")
done

# Phase 2: compute scale factors (median / sample_count)
median_spike=$(printf '%s\n' "${spike_counts[@]}" | sort -n \
  | awk '{a[NR]=$1} END {if(NR%2==1) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')
for i in "${!spike_counts[@]}"; do
  sc="${spike_counts[$i]}"
  if [[ "$sc" -gt 0 ]]; then
    scale_factors[$i]=$(echo "scale=6; $median_spike / $sc" | bc)
  else
    scale_factors[$i]=1.0
  fi
done

# Phase 3: per-sample binned coverage + CBS segmentation
run_one_sample() {
  local idx="$1" sf="$2"
  local pref="${PREFIXES[$idx]}" bam="${SPIKEIN_BAMS[$idx]}"
  mkdir -p "${OUTDIR}/${pref}"; cd "${OUTDIR}/${pref}"

  if [[ ! -f "input_${BIN}.cov_forR.bed" ]]; then
    bedtools makewindows -g "${CHRSIZES}" -w ${BIN} \
      | awk '$1 ~ /^(chr2L|chr2R|chr3L|chr3R|chr4|chrX)$/' > genome_${BIN}.bed

    bamCoverage -b "${bam}" -o input_${BIN}.bw \
      --binSize ${BIN} --scaleFactor ${sf} --extendReads -p 4

    bigWigToBedGraph input_${BIN}.bw input_${BIN}.bdg
    awk '$1 ~ /^(chr2L|chr2R|chr3L|chr3R|chr4|chrX)$/' input_${BIN}.bdg \
      | sort -k1,1 -k2,2n > input_${BIN}.sorted.bdg

    bedtools sort -i genome_${BIN}.bed | \
      bedtools map -a stdin -b input_${BIN}.sorted.bdg -c 4 -o mean \
      | awk 'BEGIN{OFS="\t"}{if($4==".") $4=0; print}' > input_${BIN}.cov_forR.bed
  fi

  # CBS segmentation in R
  if [[ ! -f "segments_${pref}.csv" ]]; then
    export CHR_BIN=${BIN} BASELINE_CN=${BASELINE_CN} GTF_FILE="${GTF}" SAMPLE_PREFIX="${pref}"
    Rscript --vanilla -e '
library(data.table); library(DNAcopy); library(GenomicRanges); library(rtracklayer)
BIN <- as.integer(Sys.getenv("CHR_BIN"))
baseline_cn <- as.numeric(Sys.getenv("BASELINE_CN"))
pref <- Sys.getenv("SAMPLE_PREFIX")
dt <- fread(sprintf("input_%d.cov_forR.bed", BIN), header=FALSE)
setnames(dt, c("chr","start","end","cov"))
dt[is.na(cov), cov := 0]
cov1 <- dt$cov + 1
q99 <- quantile(cov1, 0.99, na.rm=TRUE)
trimmed_med <- median(cov1[cov1 <= q99], na.rm=TRUE)
dt[, log2r := log2((cov + 1) / trimmed_med)]
seg <- segment(smooth.CNA(CNA(dt$log2r, dt$chr, dt$start, data.type="logratio")), alpha=0.01, min.width=2)
fwrite(as.data.table(seg$output), sprintf("segments_%s.csv", pref))
cat(sprintf("[R] Segments written for %s\n", pref))
'
  fi
  cd - >/dev/null
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  run_one_sample "${SLURM_ARRAY_TASK_ID}" "${scale_factors[$SLURM_ARRAY_TASK_ID]}"
else
  for (( i=0; i<total; i++ )); do run_one_sample "$i" "${scale_factors[$i]}"; done
fi

echo "✅ Spike-in normalised CN estimation complete. Results in ${OUTDIR}/"
