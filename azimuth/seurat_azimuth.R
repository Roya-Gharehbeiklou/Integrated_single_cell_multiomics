# Loading libraries
.libPaths('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya')
library(Azimuth)
library(Seurat)
library(patchwork)

setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/azimuth_results")

# Returns a seurat object with cell-type anotations
seurat_azimuth <- RunAzimuth('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Karina/Seurat_singlets_QCfiltered.rds', reference = "pbmcref")
seurat_azimuth
#An object of class Seurat 
#52919 features across 11939 samples within 5 assays 
#Active assay: refAssay (26412 features, 4861 variable features)

# Output visualiztion
p1 <- DimPlot(seurat_azimuth, group.by = "predicted.celltype.l2")

# Data normalization
seurat_azimuth <- NormalizeData(seurat_azimuth)

# Set identity classes
Idents(seurat_azimuth) <- "predicted.celltype.l2"

# Visualization:
#expression of CCR7 on CD4 and CD8 Naive T cells
p2 <- FeaturePlot(seurat_azimuth, features = "CCR7")
# expression of FCRL3 on regulatory T cells
p3 <- FeaturePlot(seurat_azimuth, features = "FCRL3")

sort(table(seurat_azimuth$predicted.celltype.l2), decreasing = TRUE)

write_rds(seurat_azimuth, "pbmc_Seurat_Object_QCfiltered_Azimuth.rds")

#Annotation lower resolution cell type
plasma <- 'Plasma'
b <- c('B', 'B intermediate','B naive','B memory')
cd4t <- c('CD4 TCM', 'CD4 Naive', 'Treg','MAIT','gdT','CD4 CTL' )
cd8t <- c('CD8+ T','CD8 TEM',"CD8 Naive", 'CD8 TCM')
dc <- c('pDC', 'mDC', 'cDC2','cDC1')
monocyte <- c('cMonocyte', 'ncMonocyte', 'CD14 Mono','CD16 Mono','NK_CD56bright' )
nk <- c('NK','CD56(dim) NK', 'CD56(bright) NK','NK Proliferating')
megakaryocyte <- c('Megakaryocyte')

pbmc_Seurat_Azimuth_for_figR <- seurat_azimuth

# use groups to get lower resolution cell type
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% plasma, 'cell_type_lowerres'] <- 'plasma'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% b, 'cell_type_lowerres'] <- 'B'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% cd4t, 'cell_type_lowerres'] <- 'CD4T'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% cd8t, 'cell_type_lowerres'] <- 'CD8T'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% dc, 'cell_type_lowerres'] <- 'DC'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% monocyte, 'cell_type_lowerres'] <- 'monocyte'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% nk, 'cell_type_lowerres'] <- 'NK'
pbmc_Seurat_Azimuth_for_figR@meta.data[pbmc_Seurat_Azimuth_for_figR@meta.data$predicted.celltype.l2 %in% megakaryocyte, 'cell_type_lowerres'] <- 'megakaryocyte'

pbmc_Seurat_Azimuth_for_figR <- pbmc_Seurat_Azimuth_for_figR[, !is.na(pbmc_Seurat_Azimuth_for_figR@meta.data$cell_type_lowerres)]
pbmc_Seurat_Azimuth_for_figR
#An object of class Seurat
#52919 features across 11853 samples within 5 assays
#Active assay: refAssay (26412 features, 4861 variable features)
Idents(pbmc_Seurat_Azimuth_for_figR) <- "cell_type_lowerres"
p4 <- DimPlot(seurat_azimuth, group.by = "cell_type_lowerres")
