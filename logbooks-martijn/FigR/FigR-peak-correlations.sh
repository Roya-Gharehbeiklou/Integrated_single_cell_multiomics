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

cd /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/logbooks

module load RPlus

Rscript FigR-peak-correlations.R &> output.txt