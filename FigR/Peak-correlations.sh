#!/bin/bash
#SBATCH --job-name=FigR-peak-correlations
#SBATCH --output=FigR-output-peak-correlations
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=21
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load RPlus

Rscript Peak-correlations.R &> Peak-correlations.txt