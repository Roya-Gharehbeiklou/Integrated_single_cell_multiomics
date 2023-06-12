library(GenomicRanges)
library(cisTopic)

setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/output/')

# Or for non singularity server
library(GenomicRanges, lib.loc='../Fig_R_libs')
library(cisTopic, lib.loc='../Fig_R_libs')

# Example

bamFiles <- 'pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam'
regions <- 'pbmc_granulocyte_sorted_10k_atac_peaks.bed'

# Add paired = TRUE
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, project.name='3kPBMC', paired = TRUE)
