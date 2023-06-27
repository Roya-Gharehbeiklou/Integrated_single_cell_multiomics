#!/bin/bash
#SBATCH --job-name=cisTopic_models
#SBATCH --output=cisTopic_models
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=64gb
#SBATCH --nodes=1

module load RPlus

Rscript Create_cisTopic.R &> Create_cisTopic.txt