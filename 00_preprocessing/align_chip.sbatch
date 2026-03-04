#!/bin/bash
#SBATCH -J chip_bowtie_array
#SBATCH -p slim16
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 0-6:00
#SBATCH --mem 40000
#SBATCH -o bowtie_output/bowtie_%A_%a.out
#SBATCH -e bowtie_output/bowtie_%A_%a.err

# Called by align_chip.sh via sbatch --array
# Runs the single-sample worker from the same script
bash align_chip.sh --local ""
