#!/bin/bash
#SBATCH -J chip_coverage
#SBATCH -p slim16
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -t 0-6:00
#SBATCH --mem 40000
#SBATCH -o outcoverage/coverage_%A.out
#SBATCH -e outcoverage/coverage_%A.err

# Called by coverage_chip.sh via sbatch
# FILEBASE and INPUTBASE are passed via --export
bash coverage_chip.sh --local
