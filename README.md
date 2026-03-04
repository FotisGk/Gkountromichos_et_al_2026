# Publication Code

Analysis scripts for Gkountromichos et al.  
*"Dosage compensation defects due to roX RNA deletion are rescued by recalibration of X/autosome stoichiometry"*

---

## Repository structure

```
publication_code/
├── install_dependencies.R
├── 00_preprocessing/
│   ├── align_cutrun.sh / .sbatch
│   ├── align_chip.sh   / .sbatch
│   ├── coverage_cutrun.sh
│   └── coverage_chip.sh / .sbatch
├── 01_csaw_site_discovery/
│   └── csaw_discovery_msl2_h4k16ac.R
├── 02_deeptools_signal/
│   ├── heatmaps_metaplots_csaw_sites.sh
│   └── gene_body_profiles.sh
├── 03_site_characterization/
│   ├── site_dissection_analysis.R
│   ├── motif_discovery.sh
│   └── motif_scanning.sh
├── 04_rnaseq_analysis/
│   └── deseq2_differential.R
├── 05_copy_number/
│   ├── spikein_normalization.sh
│   ├── copy_number_estimation.R
│   └── chrX_supplementary_plot.R
├── 06_dosage_compensation/
│   ├── expression_boxplots.R
│   ├── gene_level_cnv.R
│   └── xa_ratio_barplots.R
└── 07_microscopy/
    ├── max_projection_and_rename.ijm
    ├── segment_measure_crop.ijm
    └── image_feature_analysis.Rmd
```

---

## Figure → script mapping

### Main figures

| Panel | Script(s) |
|-------|-----------|
| Fig 1C – roX1/roX2 RNA-seq expression | `04/deseq2_differential.R` |
| Fig 2E – UMAP of MSL2 IF single-cell features | `07/image_feature_analysis.Rmd` |
| Fig 3B – Chromosomal distribution of binding sites | `03/site_dissection_analysis.R` |
| Fig 3C – UpSet plot (MCCS ∩ HAS ∩ PionX) | `03/site_dissection_analysis.R` |
| Fig 3D – Signal intensities at classified MCCS | `03/site_dissection_analysis.R` |
| Fig 3E – Genomic feature annotation | `03/site_dissection_analysis.R` |
| Fig 3F – Distance to TSS and HAS | `03/site_dissection_analysis.R` |
| Fig 3G – Enriched motifs (MRE, DREF/BEAF-32) | `03/motif_discovery.sh` |
| Fig 4A–F – Gene body heatmaps & metaplots | `02/gene_body_profiles.sh` |
| Fig 5A – RNA log₂FC boxplots by chromosome | `06/expression_boxplots.R` |
| Fig 5C – Chromosome copy number barplot | `05/copy_number_estimation.R` |
| Fig 5D – Three-way X:A ratio barplots | `06/xa_ratio_barplots.R` |

### Supplementary figures

| Panel | Script(s) |
|-------|-----------|
| SFig 1A – UMAP projections (additional replicates) | `07/image_feature_analysis.Rmd` |
| SFig 2B–C – MSL2/H4K16ac signal at binding sites | `02/heatmaps_metaplots_csaw_sites.sh` |
| SFig 3A – PCA | `04/deseq2_differential.R` |
| SFig 3B – Sample correlation heatmap | `04/deseq2_differential.R` |
| SFig 3C – Volcano plots | `04/deseq2_differential.R` |
| SFig 4 – Smoothed chrX CN profiles | `05/chrX_supplementary_plot.R` |

---

## Execution order

Run modules in numerical order; shell scripts before R scripts within each
module.

```
 0a  00_preprocessing/align_cutrun.sh
 0b  00_preprocessing/align_chip.sh
 0c  00_preprocessing/coverage_cutrun.sh
 0d  00_preprocessing/coverage_chip.sh
 1   01_csaw_site_discovery/csaw_discovery_msl2_h4k16ac.R
 2   02_deeptools_signal/heatmaps_metaplots_csaw_sites.sh
 3   02_deeptools_signal/gene_body_profiles.sh
 4   03_site_characterization/site_dissection_analysis.R   (needs step 1)
 5   03_site_characterization/motif_discovery.sh           (needs step 4)
 6   03_site_characterization/motif_scanning.sh
 7   04_rnaseq_analysis/deseq2_differential.R
 8   05_copy_number/spikein_normalization.sh
 9   05_copy_number/copy_number_estimation.R
10   05_copy_number/chrX_supplementary_plot.R
11   06_dosage_compensation/expression_boxplots.R
12   06_dosage_compensation/gene_level_cnv.R
13   06_dosage_compensation/xa_ratio_barplots.R            (needs steps 9 + 12)
14   07_microscopy/max_projection_and_rename.ijm            (Fiji)
15   07_microscopy/segment_measure_crop.ijm                 (Fiji, needs step 14)
16   07_microscopy/image_feature_analysis.Rmd               (needs step 15)
```

---

## Software dependencies

**R packages** — install all at once with `Rscript install_dependencies.R`.  
Key packages: csaw, edgeR, DESeq2, ChIPseeker, GenomicRanges, GenomicFeatures,
GenomicAlignments, Rsamtools, rtracklayer, BSgenome.Dmelanogaster.UCSC.dm6,
TxDb.Dmelanogaster.UCSC.dm6.ensGene, ggplot2, patchwork, pheatmap,
UpSetR, ggvenn, ggrepel, tidyverse.  
Microscopy analysis: EBImage, ComplexHeatmap, GLCMTextures, tiff, caret,
tidymodels, corrr, umap, cluster, ggfortify.

**Command-line tools:** bowtie2, samtools, bedtools, deeptools, wiggletools,
UCSC utilities (bigWigToBedGraph, bedGraphToBigWig, wigToBigWig),
MEME Suite (streme, fimo, tomtom), Homer (bamToBed), bc.

**Image processing:** Fiji/ImageJ (for `.ijm` macros in `07_microscopy/`).

**Reference files:** dm6 genome FASTA + bowtie2 index, droVir3 bowtie2 index
(ChIP spike-in), dm6.chrom.sizes.txt, genes.gtf (Ensembl dm6 annotation).

---

## Input data

```
data/
├── input/
│   ├── genes_counts_RAW.tsv    # RNA-seq raw counts (TEcount)
│   ├── SampleSheet.tsv         # RNA-seq sample metadata
│   └── genes.gtf
├── chipseq_sample_metadata.csv # CUT&RUN / ChIP metadata
└── (BAM/BED files)
```

Additional BED files consumed by deeptools scripts:
`merged_top_sites_QLF_msl2.bed`, `pionx_sites.bed`, `has_sites.bed`,
`genes_chrX.bed`, `genes_autosomes.bed`, `genes_all.bed`,
`genes_all_MSL2bound.bed`.

Microscopy input (module 07): per-nucleus multi-channel TIFF thumbnails
produced by the Fiji macros, placed in a `tiffs/` directory inside
`07_microscopy/`. The main-figure UMAP (Fig 2E) uses replicate 4 data.
Create a symlink: `ln -s /path/to/your/tiffs_rep4 07_microscopy/tiffs`
before knitting the Rmd. UMAP embedding is stochastic; cluster
assignments and composition barplots are deterministic given `set.seed(179)`.

---

## Notes

- **CUT&RUN controls differ per target:** MSL2 → PPI (pre-immune serum);
  H4K16ac → IgG. BigWig filenames reflect this convention.
- All R scripts use relative paths — run from the project root.
- Shell scripts contain `module load` commands for SLURM HPC environments;
  comment these out if running locally.
- **Genotype naming:** WT = S2A wild-type; KO-A = rox2ko112 (*roX2* KO clone
  112); KO-B = rox2ko17 (*roX2* KO clone 17); Rescue = rox2ko17rox2fl (KO-B
  with *roX2* re-expressed from a floxed transgene). Colour palette
  (colour-blind friendly): WT `#E69F00`, KO-A `#009E73`, KO-B `#56B4E9`,
  Rescue `#CC79A7`.
- **Stable line naming (microscopy, Fig 2):** rox2ko17_rox2fl = full-length
  roX2 rescue; rox2ko17_a0b0 = A0B0 (exon 3, 512 bp); rox2ko17_a3b0 = A3B0
  (5′ SL mutant); rox2ko17_minirox = miniroX (5′+3′ SLs only);
  rox2ko17_rox1d3f = roX1-D3 (3′ fragment of roX1).
- **Code consolidation.** These scripts were consolidated and factored from 
  the original analysis scripts with the assistance of Claude Opus 4 (Anthropic). 
  Numerical outputs and figures have been validated against the original production
  results.
