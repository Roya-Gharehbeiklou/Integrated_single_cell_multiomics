## First we read the h5 matrix file containing the peaks
library(Seurat)
setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023')
peaks.pbmc <- Read10X_h5('data/output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', use.names=TRUE)

##  Create Seurat object
peaks.pbmc.seurat <- CreateSeuratObject(counts=peaks.pbmc$`Gene Expression`, project="pbmc10k", assay="ATAC")

## View object
str(peaks.pbmc.seurat)

## Altering to SummarizedExperiment object for FigR input
library(SummarizedExperiment)

pbmc.seurat.summarized <- SummarizedExperiment(peaks.pbmc.seurat)
str(pbmc.seurat.summarized)

library(FigR, lib.loc="/")
