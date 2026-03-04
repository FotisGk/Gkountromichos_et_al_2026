#!/usr/bin/env bash
# =============================================================================
# De novo motif discovery at csaw-identified MSL2 binding sites
# Uses MEME Suite STREME for discriminative motif discovery with matched
# chromosome background, plus Tomtom for motif annotation.
#
# Sequence extraction and background generation are performed by the companion
# R script site_dissection_analysis.R (BSgenome.Dmelanogaster.UCSC.dm6 +
# Biostrings). This script expects the FASTA files it produces:
#   - sites_<category>.fa     (200 bp midpoint-centred target sequences)
#   - background_<category>.fa (matched-chromosome random regions, seed=42)
# Categories: MCCS (Both HAS+PionX), HAS, PionX, Novel, All
#
# Figures: 3G (enriched motifs at MSL2 sites)
# Prerequisites: MEME Suite v5.5.8 (streme, tomtom)
# =============================================================================
set -euo pipefail

# --- EDIT THESE ---------------------------------------------------------------
SEQ_DIR="motif_sequences"          # output from site_dissection_analysis.R
OUTDIR="motif_results"
# Tomtom databases (adjust paths for local MEME Suite installation)
TOMTOM_DBS="--m db/fly_factor_survey.meme db/FlyReg_v2.meme db/OnTheFly_2014.meme db/dmmpmm2009.meme"
# ------------------------------------------------------------------------------
mkdir -p "${OUTDIR}"

# Run STREME for each binding-site category (MCCS, HAS, PionX, etc.)
for TARGET_FA in "${SEQ_DIR}"/sites_*.fa; do
  CATEGORY=$(basename "${TARGET_FA}" .fa | sed 's/^sites_//')
  BG_FA="${SEQ_DIR}/background_${CATEGORY}.fa"

  if [[ ! -f "${BG_FA}" ]]; then
    echo "⚠️  Skipping ${CATEGORY}: background FASTA not found (${BG_FA})"
    continue
  fi

  echo "--- STREME: ${CATEGORY} ---"
  streme --p "${TARGET_FA}" --n "${BG_FA}" \
    --oc "${OUTDIR}/streme_${CATEGORY}" \
    --dna --nmotifs 5 --minw 6 --maxw 20 --thresh 0.05

  # Tomtom — annotate discovered motifs against Drosophila databases
  if [[ -f "${OUTDIR}/streme_${CATEGORY}/streme.txt" ]]; then
    echo "--- Tomtom: ${CATEGORY} ---"
    tomtom -oc "${OUTDIR}/tomtom_${CATEGORY}" \
      -min-overlap 5 -dist pearson -evalue -thresh 10.0 \
      "${OUTDIR}/streme_${CATEGORY}/streme.txt" \
      ${TOMTOM_DBS}
  fi
done

echo "✅ Motif discovery + annotation complete. Results in ${OUTDIR}/"
