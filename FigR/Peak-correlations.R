setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()

# Load modules without issues
library(chromVAR, lib.loc='Fig_R_libs')
library(Matrix, lib.loc='Fig_R_libs')
library(S4Vectors, lib.loc='ArchR_libs')
library(GenomeInfoDb, lib.loc='ArchR_libs')
library(ggplot2, lib.loc='ArchR_libs')
library(motifmatchr, lib.loc='ArchR_libs')
library(FigR, lib.loc='Fig_R_libs')
library(pbmcapply, lib.loc="Fig_R_libs")
library(BSgenome, lib.loc='Fig_R_libs')
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='Fig_R_libs')

# Load both objects
ATAC.se <- readRDS('../Users/Martijn/Integrated_single_cell_multiomics/FigR/output/ATAC_se.rds')
RNAmat <- readRDS('../Users/Martijn/Integrated_single_cell_multiomics/FigR/output/RNA_mat.rds')

# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "hg38",
                           nCores = 16,
                           p.cut = NULL,
                           n_bg = 100)

write.table(cisCorr, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/ciscor.csv', quote=FALSE, sep="\t")
cisCorr <- read.csv('../Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/ciscor.csv', sep="\t")
head(cisCorr)