# Author: Karina Diaz Barba
# Code to analyse the PBMC 10X data and create a Seurat object with the data 
# filtered based on specific QC parameters

# Remove previous variables
rm(list=ls())

# Set working directory
setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Karina")

# Libraries that are going tobe used
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)

# Load the PBMC dataset. 
# Note: Files needed: barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
pbmc.data <- Read10X(data.dir = "/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/output/raw_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
# Calculate the number of unique genes and total molecules 
# Note: there was an error displayed ' Error in CreateAssayObject(counts = counts,
# min.cells = min.cells, min.features = min.features, : No cell names (colnames) 
# names present in the input matrix'
# Was solved adding the column name ($`Gene Expression`) in the matrix created above
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, 
                           project = "PBMCs10X", 
                           min.cells = 3, 
                           min.features = 200)

# Calculate mitochondrial QC metrics with the PercentageFeatureSet()
# We use the set of all genes starting with MT- as a set of mitochondrial genes
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

#calculate the dimensions of the object
dim(pbmc@meta.data) #12391     4

# Save seurat object as an rds
#saveRDS(pbmc, file = "pbmc_Seurat_Object.rds")

# Visualize QC metrics as a violin plot and create a PDF
pdf("Violin_plot_PBMC.pdf", width=10, height=10)
par(mar=c(0,0,50,0), pin=c(3,3))
# Generate the VlnPlot for three feature
vlnplot1 <- VlnPlot(pbmc, features = "nFeature_RNA", cols='#7fc97f') + NoLegend()
vlnplot2 <- VlnPlot(pbmc, features = "nCount_RNA", cols='#beaed4') + NoLegend()
vlnplot3 <- VlnPlot(pbmc, features = "percent.mt", cols='#fdc086') + NoLegend()
# Combine the plots
grid.arrange(vlnplot1, vlnplot2, vlnplot3, 
                ncol = 3 , 
                top=grid::textGrob('All PBMC data metrics', 
                gp=grid::gpar(fontsize=25)))
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, i.e. columns in object 
# metadata, PC scores etc.

#pdf("Feature_feature_relationships.pdf", width=15, height=10) 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#dev.off()


# Define a function to identify cells that are outliers based on certain MAD 
# from the median. Based on the next link:
# https://wg1-pipeline-qc.readthedocs.io/en/latest/QC/QC_Figures.html


mad_function <- function(seurat, column, number_mad){
    mad <- mad(seurat@meta.data[,column])
    low <- median(seurat@meta.data[,column]) - number_mad*mad
    high <- median(seurat@meta.data[,column]) + number_mad*mad
    print("The lower bound is:")
    print(low)
    print("The upper bound is:")
    print(high)
    seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
    return(seurat)
    if (column=='percent.mt')
        pV3 <- VlnPlot(seurat, features = "percent.mt") & geom_hline(yintercept = c(low, high))
    if (column=='nCount_RNA')
        pV2 <- VlnPlot(seurat, features = "nCount_RNA") & geom_hline(yintercept = c(low, high))
}

# Read in the seurat object if it was no loaded before
#pbmc <- readRDS("pbmc_Seurat_Object.rds")

# Identify the cells that are outliers using the MAD function created above
# Mitocondrial % - Remove %mtRNA: >= MAD_3
pbmc <- mad_function(seurat = pbmc, column = "percent.mt", number_mad = 3)
# nCount - Remove nUMIs: =< MAD_2
pbmc <- mad_function(seurat = pbmc, column = "nCount_RNA", number_mad = 2)
dim(pbmc) #26412 12391
pbmc@meta.data[1:15,]

# Visualize QC metrics as a violin plot withe the MAD lines and create a PDF
pdf("Violin_plot_MADlines.pdf", width=15, height=15) 
pV1 <- VlnPlot(pbmc, features = "nFeature_RNA")
wrap_plots(pV1, pV2, pV3, ncol = 3)
dev.off()

# Remove the outlier cells
pbmc_filtered <- subset(pbmc, subset = percent.mt_mad == "NotOutlier" & nCount_RNA_mad == "NotOutlier")

# Check the size of the seurat object after filtering to ensure that cells have 
# been removed
print(pbmc_filtered) # An object of class Seurat  26412 features across 10474 samples within 1 assay 
dim(pbmc_filtered) #26412 10474

# Save the Seurat filtered Object
saveRDS(pbmc_filtered, "pbmc_Seurat_Object_QCfiltered.rds")

# Visualize QC metrics as a violin plot and create a PDF of the FILTERED DATA
pdf("Violin_plot_PBMC_QC_filtered.pdf", width=10, height=10)
par(mar=c(0,0,50,0), pin=c(3,3))
# Generate the VlnPlot for three feature
vlnplot4 <- VlnPlot(pbmc_filtered, features = "nFeature_RNA", cols='#7fc97f') + NoLegend()
vlnplot5 <- VlnPlot(pbmc_filtered, features = "nCount_RNA", cols='#beaed4') + NoLegend()
vlnplot6 <- VlnPlot(pbmc_filtered, features = "percent.mt", cols='#fdc086') + NoLegend()
# Combine the plots
grid.arrange(vlnplot4, vlnplot5, vlnplot6, ncol = 3 , top=grid::textGrob('QC filtered PBMC data metrics', gp=grid::gpar(fontsize=25)))
dev.off()

# combine two plots (filtered and not filtered) into one
pdf("Violin_plot_PBMC_combined.pdf", width=15, height=15)
par(mar=c(0,0,50,0), pin=c(3,3))
grid.arrange(vlnplot1, vlnplot2, vlnplot3,vlnplot4, vlnplot5, vlnplot6, 
                ncol = 3 , 
                top=grid::textGrob('PBMC data QC metrics filtered ad not filtered',
                gp=grid::gpar(fontsize=25)))
dev.off()

### Normalizing the filtered data
pbmc <- NormalizeData(pbmc_filtered, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

### Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
pdf("Variable_features.pdf", width=7, height=5) 
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 + ggtitle("A. Top ten highly variable features")
dev.off()

# With the code below you can permfomed more  exploratory analyses, but is 
# comented because we do not use this in this project, as we have other tools 
# to analyse the data. 

#Scaling the data
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)

# Perform linear dimensional reduction
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#pdf("DimPlot.pdf", width=15, height=10) 
# DimPlot(pbmc, reduction = "pca")
#dev.off()

#pdf("DimHeatmap.pdf", width=15, height=10)
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#dev.off()

#### Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for 
# expediency. More approximate techniques such as those implemented in 
# ElbowPlot() can be used to reduce computation time
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# pdf("JackStrawPlot.pdf", width=15, height=10) 
# JackStrawPlot(pbmc, dims = 1:15)
# dev.off()

# pdf("ElbowPlot.pdf", width=15, height=10) 
# ElbowPlot(pbmc)
# dev.off()

#Cluster the cells
# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# head(Idents(pbmc), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
# pbmc <- RunUMAP(pbmc, dims = 1:10)

# pdf("UMAP.pdf", width=15, height=10) 
# DimPlot(pbmc, reduction = "umap")
# dev.off()

# saveRDS(pbmc, file = "pbmc_tutorial.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only t
# he positive ones
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# pbmc.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)

#  cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, 
#                                         test.use = "roc", only.pos = TRUE)

# pdf("FeaturePlot.pdf", width=15, height=10) 
# FeaturePlot(pbmc, features = c("LEF1", "BCL11B", "VCAN", "PLCB1", "CDKN1C", 
                #"FCGR3A", "CCSER1", "HLA-DQA1","CAMK4",
#     "IL7R"))
# dev.off()

# pbmc.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10
# pdf("DoHeatmap.pdf", width=15, height=10)
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()
# dev.off()
    
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                        #"FCGR3A+ Mono", "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# #pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
