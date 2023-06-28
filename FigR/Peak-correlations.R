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

datadir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/'

# Load both objects
ATAC.se <- readRDS(paste0(datadir, 'ATAC_se.rds'))
RNAmat <- readRDS(paste0(datadir,'RNA_mat.rds'))

# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "hg38",
                           nCores = 16,
                           p.cut = NULL,
                           n_bg = 100)

write.table(cisCorr, paste0(datadir, 'ciscor.csv'), quote=FALSE, sep="\t")
cisCorr <- read.csv(paste0(datadir, 'ciscor.csv'), sep="\t")
head(cisCorr)
