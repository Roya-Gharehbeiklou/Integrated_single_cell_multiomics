setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()

# Load modules without version issues
library(Seurat, lib.loc='ArchR_libs')
library(data.table, lib.loc='ArchR_libs')
library(S4Vectors, lib.loc='ArchR_libs')
library(GenomeInfoDb, lib.loc='ArchR_libs')
library(Rcpp, lib.loc='ArchR_libs')
library(gtable, lib.loc='ArchR_libs')
library(ggplot2, lib.loc='ArchR_libs')
library(Matrix, lib.loc='ArchR_libs')
library(ArchR, lib.loc='ArchR_libs')
library(SingleCellExperiment, lib.loc='Fig_R_libs')
library(scuttle, lib.loc='Fig_R_libs')
library(scater, lib.loc='Fig_R_libs')

# Setting ArchR threads
addArchRThreads(threads = 16) 

# Load preprocessed Seurat object with Azimuth celltype annotations
seurat.object <- readRDS('../Users/Dilya/azimuth_results/pbmc_Seurat_Azimuth_for_figR.rds')
RNAmat <- seurat.object@assays$RNA@data

# Log10 normalize RNA count matrix
RNAmat <- NormalizeData(RNAmat)
dim(RNAmat)

# Load ATAC data and get available matrices
ATAC.data <- readRDS('../Users/Roya/Save-ArchR-Project_subSet_QC_Frip.rds')

# Load summarized experiment
ATAC.se <- import10xFeatureMatrix('output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', names='')
colnames(ATAC.se) <- gsub("#","",colnames(ATAC.se))

# Change rownames ATAC.data object
rownames(ATAC.data@cellColData) <- gsub("pbmc_granulocyte_sorted_10k_HG38#","",rownames(ATAC.data@cellColData))

# Get only specific celltypes
plasma <- 'Plasma'
b <- c('B', 'B intermediate','B naive','B memory')
cd4t <- c('CD4 TCM', 'CD4 Naive', 'Treg','MAIT','gdT','CD4 CTL' )
cd8t <- c('CD8+ T','CD8 TEM',"CD8 Naive", 'CD8 TCM')
dc <- c('pDC', 'mDC', 'cDC2','cDC1')
monocyte <- c('cMonocyte', 'ncMonocyte', 'CD14 Mono','CD16 Mono','NK_CD56bright' )
nk <- c('NK','CD56(dim) NK', 'CD56(bright) NK','NK Proliferating')
megakaryocyte <- c('Megakaryocyte')

# Use groups to get lower resolution cell type
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% plasma, 'cell_type_lowerres'] <- 'plasma'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% b, 'cell_type_lowerres'] <- 'B'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% cd4t, 'cell_type_lowerres'] <- 'CD4T'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% cd8t, 'cell_type_lowerres'] <- 'CD8T'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% dc, 'cell_type_lowerres'] <- 'DC'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% monocyte, 'cell_type_lowerres'] <- 'monocyte'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% nk, 'cell_type_lowerres'] <- 'NK'
seurat.object@meta.data[seurat.object@meta.data$predicted.celltype.l2 %in% megakaryocyte, 'cell_type_lowerres'] <- 'megakaryocyte'

# Lose NA values
seurat.object <- seurat.object[, !is.na(seurat.object@meta.data$cell_type_lowerres)]

# Fetch barcodes ATAC and RNA data
barcodes.atac <- rownames(ATAC.data@cellColData)
barcodes.rnamat <- colnames(RNAmat)

# Subset seurat objects
seurat.object <- seurat.object[, rownames(seurat.object@meta.data) %in% barcodes.atac]
seurat.object <- seurat.object[, rownames(seurat.object@meta.data) %in% barcodes.rnamat]

# Keep annotated cells from seurat object
cellsToKeep <- sample(colnames(seurat.object),replace = FALSE)
ATAC.se <- ATAC.se[,colnames(ATAC.se) %in% cellsToKeep]
RNAmat <- RNAmat[,colnames(RNAmat) %in% cellsToKeep]

# Dimensions not matching
dim(ATAC.se) # Peaks x Cells
dim(RNAmat)

# Dimensions not matching, find intersect
cellsToKeep <- intersect(colnames(ATAC.se), colnames(RNAmat))
ATAC.se <- ATAC.se[,colnames(ATAC.se) %in% cellsToKeep]
RNAmat <- RNAmat[,colnames(RNAmat) %in% cellsToKeep]
seurat.object <- seurat.object[,colnames(seurat.object) %in% cellsToKeep]

# Create dataframe of celltype for each barcode
celltypes.barcode <- data.frame(rownames(seurat.object@meta.data), seurat.object@meta.data$cell_type_lowerres)

# Add annotation to ATAC data
ATAC.se@colData <- DataFrame(celltypes.barcode)

# Change celtype column name
names(colData(ATAC.se))[which(names(colData(ATAC.se))=="seurat.object.meta.data.cell_type_lowerres")]="Cell type"

# Matching dimensions
dim(ATAC.se)
dim(RNAmat)

# Getting chromosome names
ATAC.chromosomes <- rowData(ATAC.se)$interval

# Different starting name chromosomes
deviating.names <- ATAC.chromosomes[!startsWith(ATAC.chromosomes, 'c')]
normal.chroms <- ATAC.chromosomes[startsWith(ATAC.chromosomes, 'c')]

# Fetching chromosome starting with a 'c'
ATAC.se <- ATAC.se[startsWith(rowData(ATAC.se)$interval, 'c')]

####################################################
# SAVE RNA MATRIX AND ATAC SE TO H5 FOR PORTAL INPUT
####################################################

# Save ATAC gene scores as h5 file
library(HDF5Array, lib.loc='ArchR_libs/')

saveHDF5SummarizedExperiment(ATAC.se, dir="../Users/Roya/Portal_input", prefix="gene_scores_ATAC", replace=FALSE,
                             chunkdim=NULL, level=NULL, as.sparse=NA,
                             verbose=NA)

# Save h5 for Portal
library(scrattch.io, lib.loc='Fig_R_libs')
library(rhdf5)
write_dgCMatrix_h5(RNAmat, cols_are = "sample_names", '../Users/Roya/Portal_input/RNA_count.h5',
  ref_name = "10Xpbmc", gene_ids = NULL)

saveRDS(ATAC.se, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/ATAC_se.rds')
saveRDS(RNAmat, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/RNA_mat.rds')

#################################
# Creating UMAP
#################################
ATAC.sce <- logNormCounts(as(ATAC.se, 'SingleCellExperiment'))

ATAC.sce <- runUMAP(ATAC.sce)
UMAP.plot <- plotReducedDim(ATAC.sce, dimred="UMAP", colour_by="Cell type") +
    theme_classic()

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/UMAP.png', UMAP.plot)
