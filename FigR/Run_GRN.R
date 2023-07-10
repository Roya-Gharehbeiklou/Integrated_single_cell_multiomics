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
library(BiocManager, lib.loc='Fig_R_libs')
library(BSgenome, lib.loc='Fig_R_libs')
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='Fig_R_libs')

datadir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/'

ATAC.se <- readRDS(paste0(datadir, 'ATAC_se.rds'))
RNAmat <- readRDS(paste0(datadir, 'RNA_mat.rds'))

# Load cis correlation matrix
cisCorr <- read.csv(paste0(datadir, 'ciscor.csv'), sep="\t")
# Filter correlations on p-value of <= 0.05
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
dorcMat <- readRDS(paste0(datadir, 'dorc_mat_unsmoothed.rds'))
dorcMat.s <- readRDS(paste0(datadir, 'dorcMat_smoothed.rds'))
RNAmat.s <- readRDS(paste0(datadir, 'RNAMat_smoothed.rds'))

# Load correct function with altered paths
source("../Users/Martijn/Integrated_single_cell_multiomics/FigR/FigR_GRN_function.R")

# Run on smoothed cisTopic dorcMat and rnaMat
figR.d.unsmoothed <- runFigRGRN(ATAC.se = ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat,
                     rnaMat = RNAmat, 
                     nCores = 21)

saveRDS(figR.d.unsmoothed, paste0(datadir, 'figRGRN_unsmoothed.rds'))

# Run on smoothed cisTopic dorcMat and rnaMat
figR.d <- runFigRGRN(ATAC.se = ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s, 
                     nCores = 4)

saveRDS(figR.d, paste0(datadir, 'figRGRN.rds'))
