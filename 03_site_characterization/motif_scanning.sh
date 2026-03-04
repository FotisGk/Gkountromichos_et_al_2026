#!/usr/bin/env bash
# =============================================================================
# Scan csaw sites for known motifs (FIMO) — e.g. GA-repeat, MRE, CLAMP motif
# Classifies sites into PionX-like (GAGA motif), HAS-like, or unclassified.
#
# Figures: Upstream of site_dissection_analysis.R (provides site class BEDs)
# Prerequisites: MEME Suite (fimo), bedtools
# =============================================================================
set -euo pipefail

# --- EDIT THESE ---------------------------------------------------------------
GENOME_FA="dm6.fa"
CSAW_BED="merged_top_sites_QLF_msl2.bed"
MOTIF_DB="motif_database.meme"   # MEME-format file with target motifs
OUTDIR="fimo_results"
THRESH="1e-4"
# ------------------------------------------------------------------------------
mkdir -p "${OUTDIR}"

# Extract sequences
bedtools getfasta -fi "${GENOME_FA}" -bed "${CSAW_BED}" -fo "${OUTDIR}/sites.fa"

# Run FIMO scan
fimo --oc "${OUTDIR}/fimo_out" --thresh "${THRESH}" \
  "${MOTIF_DB}" "${OUTDIR}/sites.fa"

echo "✅ FIMO scan complete. Results in ${OUTDIR}/fimo_out/"
echo "   → Parse fimo.tsv to classify sites into PionX / HAS / other."
