# Loading libraries
.libPaths('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya')
library(Azimuth)
library(Seurat)
library(patchwork)

setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/azimuth_results")

# Returns a seurat object with cell-type anotations
seurat_azimuth <- RunAzimuth('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Karina/Seurat_singlets_QCfiltered.rds', reference = "pbmcref")

# Output visualiztion
p1 <- DimPlot(seurat_azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()

# Data normalization
seurat_azimuth <- NormalizeData(seurat_azimuth)

# Set identity classes
Idents(seurat_azimuth) <- "predicted.celltype.l2"

# Visualization:
#expression of CCR7 on CD4 and CD8 Naive T cells
p1 <- FeaturePlot(seurat_azimuth, features = "CCR7")
# expression of FCRL3 on regulatory T cells
p2 <- FeaturePlot(seurat_azimuth, features = "FCRL3")

sort(table(seurat_azimuth$predicted.celltype.l2), decreasing = TRUE)