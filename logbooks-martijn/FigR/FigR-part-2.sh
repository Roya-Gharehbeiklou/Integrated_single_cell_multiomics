#!/bin/bash
#SBATCH --job-name=FigR
#SBATCH --output=FigR
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=21
#SBATCH --mem=4gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

cd /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/logbooks/FigR

module load RPlus

Rscript FigR-part-2.R &> FigR-part-2.txt