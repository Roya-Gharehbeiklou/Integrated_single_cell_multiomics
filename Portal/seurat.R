#Author: Roya Gharehbeikou
## This script runs FigR on the 10X data
## Add config file w/ paths in final pipline
setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()

## For Singularity container
## Load required packages

# library(chromVAR)
# library(S4Vectors)
# library(GenomeInfoDb)
# library(FigR)

# # ## More dependencies
# library(FNN)

# ## Even more dependencies
# library(BuenColors)
# library(BSgenome)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(SeuratObject)
# library(Seurat)
# library(ArchR)
# library(SummarizedExperiment)

## No container 

library(Matrix, lib.loc='Fig_R_libs')
library(chromVAR, lib.loc='Fig_R_libs')
library(motifmatchr, lib.loc='Fig_R_libs')
library(ggplot2, lib.loc='Fig_R_libs')
library(S4Vectors, lib.loc='Fig_R_libs')
library(GenomeInfoDb, lib.loc='Fig_R_libs')
library(FigR, lib.loc = 'R_libs')

#This first secstion (until UMAP) is the preprocessing/EDA part

# More dependencies
library(FNN, lib.loc='Fig_R_libs')

# Even more dependencies
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

addArchRThreads(threads = 21) 

seurat.object <- readRDS('../Users/Dilya/azimuth_results/pbmc_Seurat_Azimuth_for_figR.rds')
RNAmat <- seurat.object@assays$RNA@data
RNAmat <- NormalizeData(RNAmat)
dim(RNAmat)

head(RNAmat)

# ATAC.data.example <- readRDS('logbooks/FigR/FigR_build_in_data/shareseq_skin_SE_final.rds')

ATAC.data <- readRDS('../Users/Roya/Save-ArchR-Project_subSet_QC_Frip.rds')
getAvailableMatrices(ATAC.data)

ATACgene.score.matrix <- getMatrixFromProject(
  ArchRProj = ATAC.data,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

colnames(ATACgene.score.matrix) <- gsub("pbmc_granulocyte_sorted_10k_HG38#","",colnames(ATACgene.score.matrix))

# Add feature matrix to ArchR project
#fL <- getFragmentsFromProject(ATAC.data)
#granges <- fL[[1]]

# Load summarized experiment
ATAC.se <- ArchR::import10xFeatureMatrix('output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', names='')
colnames(ATAC.se)<-gsub("#","",colnames(ATAC.se))

# Change rownames ATAC.data object
rownames(ATAC.data@cellColData) <- gsub("pbmc_granulocyte_sorted_10k_HG38#","",rownames(ATAC.data@cellColData))

# Get correct cells


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
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% plasma, 'cell_type_lowerres'] <- 'plasma'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% b, 'cell_type_lowerres'] <- 'B'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% cd4t, 'cell_type_lowerres'] <- 'CD4T'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% cd8t, 'cell_type_lowerres'] <- 'CD8T'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% dc, 'cell_type_lowerres'] <- 'DC'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% monocyte, 'cell_type_lowerres'] <- 'monocyte'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% nk, 'cell_type_lowerres'] <- 'NK'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% megakaryocyte, 'cell_type_lowerres'] <- 'megakaryocyte'


# Loosing some barcodes
seurat.object <- seurat.object[, !is.na(seurat.object@meta.data$cell_type_lowerres)]

barcodes.atac <- rownames(ATAC.data@cellColData)
barcodes.rnamat <- colnames(RNAmat)
# Loosing more barcodes
seurat.object <- seurat.object[, rownames(seurat.object@meta.data) %in% barcodes.atac]
seurat.object <- seurat.object[, rownames(seurat.object@meta.data) %in% barcodes.rnamat]

# Keep annotated cells
cellsToKeep <- sample(colnames(seurat.object),replace = FALSE)
ATAC.se <- ATAC.se[,colnames(ATAC.se) %in% cellsToKeep]
RNAmat <- RNAmat[,colnames(RNAmat) %in% cellsToKeep]

# Dimensions not matching
dim(ATAC.se) # Peaks x Cells
dim(RNAmat)

cellsToKeep <- intersect(colnames(ATAC.se), colnames(RNAmat))
ATAC.se <- ATAC.se[,colnames(ATAC.se) %in% cellsToKeep]
RNAmat <- RNAmat[,colnames(RNAmat) %in% cellsToKeep]
seurat.object <- seurat.object[,colnames(seurat.object) %in% cellsToKeep]

celltypes.barcode <- data.frame(rownames(seurat.object@meta.data), seurat.object@meta.data$cell_type_lowerres)

ATAC.se@colData <- DataFrame(celltypes.barcode)

dim(ATAC.se)
dim(RNAmat)

# UMAP
library(lgr)
library(cisTopic)
library(SingleCellExperiment)
library(scater)

annotation <- ATAC.se@colData$seurat.object.meta.data.cell_type_lowerres

ATAC.singlecell <- as(ATAC.se, "SingleCellExperiment")
ATAC.singlecell <- scater::logNormCounts(ATAC.singlecell)
u <- uwot::umap(as.matrix(t(counts(ATAC.singlecell))), n_neighbors=2)
reducedDim(ATAC.singlecell, "UMAP_uwot") <- u
celltypes.barcode$UMAP1 <- reducedDim(ATAC.singlecell, "UMAP_uwot")[,1]
celltypes.barcode$UMAP2 <- reducedDim(ATAC.singlecell, "UMAP_uwot")[,2]

ATAC.singlecell@colData <- DataFrame(celltypes.barcode)

# Subset ATACgene.score.matrix also
ATACgene.score.matrix <- ATACgene.score.matrix[,colnames(ATACgene.score.matrix) %in% cellsToKeep]
saveRDS(ATACgene.score.matrix, '../Users/Roya/Portal_input/gene_activity_ATAC.rds')

# Save ATAC gene scores as h5 file
library(HDF5Array, lib.loc='ArchR_libs/')

saveHDF5SummarizedExperiment(ATAC.se, dir="../Users/Roya/Portal_input", prefix="gene_scores_ATAC", replace=FALSE,
                             chunkdim=NULL, level=NULL, as.sparse=NA,
                             verbose=NA)
