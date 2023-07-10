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
library(lgr, lib.loc='Fig_R_libs')
library(cisTopic, lib.loc='Fig_R_libs')
library(FNN, lib.loc='Fig_R_libs')
library(pbmcapply, lib.loc='Fig_R_libs')
library(doParallel, lib.loc='Fig_R_libs')

datadir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/'

## Load cisTopic and fetch matrix
cisTopicObject <- readRDS(paste0(datadir, "best_cisTopicObject.rds"))
topic.mat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
topic.mat <- t(topic.mat)
topic.mat <- as.matrix(topic.mat)

# Derive cell kNN using this
set.seed(123)
cellkNN <- get.knn(topic.mat,k = 30)$nn.index
dim(cellkNN)
rownames(cellkNN) <- rownames(topic.mat)

# Load cis correlation matrix
cisCorr <- read.csv(paste0(datadir, 'ciscor.csv'), sep="\t")

# Filter correlations on p-value of <= 0.05
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

# Calculate DORCs with a cutoff of 3
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                         cutoff = 3, # No. sig peaks needed to be called a DORC
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=2)

ggsave(paste0(datadir, 'dorcs.png'))

numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

# Get ATAC SE and RNAmat data for DORC score calculation
ATAC.se <- readRDS(paste0(datadir, 'ATAC_se.rds'))
RNAmat <- readRDS(paste0(datadir, 'RNA_mat.rds'))

dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 21)

dim(dorcMat)
colnames(dorcMat) <- rownames(cellkNN)

# Smooth dorc scores using cell KNNs (k=30)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 4)

# Smooth RNA using cell KNNs
# This takes longer since it's all genes
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat,nCores = 4)

saveRDS(dorcMat.s, paste0(datadir, 'dorcMat_smoothed.rds'))
saveRDS(RNAmat.s, paste0(datadir, 'RNAMat_smoothed.rds'))
