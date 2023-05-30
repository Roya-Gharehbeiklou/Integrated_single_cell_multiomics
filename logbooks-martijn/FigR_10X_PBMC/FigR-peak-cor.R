## This script runs FigR on build-in data
## Add config file w/ paths in final pipline
setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()
## Load required packages

library(Matrix, lib.loc='Fig_R_libs')
library(chromVAR, lib.loc='Fig_R_libs')
library(motifmatchr, lib.loc='Fig_R_libs')
library(ggplot2, lib.loc='Fig_R_libs')
library(S4Vectors, lib.loc='Fig_R_libs')
library(GenomeInfoDb, lib.loc='Fig_R_libs')
library(FigR, lib.loc = 'R_libs')

# This first secstion (until UMAP) is the preprocessing/EDA part

## More dependencies
library(FNN, lib.loc='Fig_R_libs')

## Even more dependencies
library(doParallel, lib.loc='Fig_R_libs')
library(BuenColors, lib.loc='Fig_R_libs')
library(BSgenome, lib.loc='Fig_R_libs')
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='Fig_R_libs')
library(sp, lib.loc='ArchR_libs')
library(SeuratObject, lib.loc='ArchR_libs')
library(Seurat, lib.loc='ArchR_libs')
library(rhdf5, lib.loc='ArchR_libs')
library(ArchR, lib.loc='ArchR_libs')
library(SummarizedExperiment, lib.loc='Fig_R_libs')

seurat.object <- readRDS('../Users/Karina/pbmc_Seurat_Object_QCfiltered.rds')
RNAmat <- GetAssayData(object=seurat.object, slot="counts")
dim(RNAmat)

ATAC.data.example <- readRDS('logbooks/FigR/FigR_build_in_data/shareseq_skin_SE_final.rds')
ATAC.se <- ArchR::import10xFeatureMatrix('output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', names='')
ATAC.se

ATAC.data <- readRDS('../Users/Roya/Save-ArchR-Project.rds')
#ATAC.data <- ArchR::addPeakSet(ATAC.data)
#ATAC.data <- ArchR::addPeakMatrix(ATAC.data)

proj <- addGeneExpressionMatrix(input = ATAC.data, seRNA = ATAC.se, force = TRUE)

annoCols <- readRDS('../Users/Dilya/azimuth_results/pbmc_Seurat_Object_QCfiltered_Azimuth.rds')
#annoCols.correct <- levels(annoCols@active.ident)

## Get only specific celltypes
plasma <- 'Plasma'
b <- c('B', 'B intermediate','B naive','B memory')
cd4t <- c('CD4 TCM', 'CD4 Naive', 'Treg','MAIT','gdT','CD4 CTL' )
cd8t <- c('CD8+ T','CD8 TEM',"CD8 Naive", 'CD8 TCM')
dc <- c('pDC', 'mDC', 'cDC2','cDC1')
monocyte <- c('cMonocyte', 'ncMonocyte', 'CD14 Mono','CD16 Mono','NK_CD56bright' )
nk <- c('NK','CD56(dim) NK', 'CD56(bright) NK','NK Proliferating')
megakaryocyte <- c('Megakaryocyte')
# use groups to get lower resolution cell type
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% plasma, 'cell_type_lowerres'] <- 'plasma'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% b, 'cell_type_lowerres'] <- 'B'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% cd4t, 'cell_type_lowerres'] <- 'CD4T'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% cd8t, 'cell_type_lowerres'] <- 'CD8T'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% dc, 'cell_type_lowerres'] <- 'DC'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% monocyte, 'cell_type_lowerres'] <- 'monocyte'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% nk, 'cell_type_lowerres'] <- 'NK'
annoCols@meta.data[annoCols@meta.data$predicted.celltype.l2 %in% megakaryocyte, 'cell_type_lowerres'] <- 'megakaryocyte'

# Loosing some barcodes
annoCols <- annoCols[, !is.na(annoCols@meta.data$cell_type_lowerres)]

barcodes.atac <- colnames(ATAC.se)

# Loosing more barcodes
annoCols <- annoCols[, rownames(annoCols@meta.data) %in% barcodes.atac]

# Keep annotated cells
cellsToKeep <- sample(colnames(annoCols),replace = FALSE)
ATAC.se <- ATAC.se[,cellsToKeep]
RNAmat <- RNAmat[,cellsToKeep]

#ATAC.data
dim(ATAC.se) # Peaks x Cells
colnames(ATAC.se)<-gsub("#","",colnames(ATAC.se))

celltypes.barcode <- data.frame(rownames(annoCols@meta.data), annoCols@meta.data$predicted.celltype.l2)
celltypes.barcode



RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]

library(pbmcapply, lib.loc="R_libs")


# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "hg38", # One of hg19, mm10 or hg38 
                           nCores = 21,
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)

head(cisCorr)
write.table(cisCorr, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/', quote=FALSE, sep="\t")