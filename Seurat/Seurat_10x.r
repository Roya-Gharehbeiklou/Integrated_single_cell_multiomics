
# Set working directory
setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Karina")
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset. 
# Note: Files needed: barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
pbmc.data <- Read10X(data.dir = "/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/output/raw_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
# Calculate the number of unique genes and total molecules 
# Note: there was an error displayed ' Error in CreateAssayObject(counts = counts, min.cells = min.cells, min.features = min.features, : No cell names (colnames) names present in the input matrix'
# Was solved adding the column name ($`Gene Expression`) in the matrix created above
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "PBMCs10X", min.cells = 3, min.features = 200)
#pbmc

# Calculate mitochondrial QC metrics with the PercentageFeatureSet()
# We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# save seurat object
#saveRDS(pbmc, file = "pbmc_Seurat_Object.rds")

# Visualize QC metrics as a violin plot and create a PDF
#pdf("Violin_plot_QC.pdf", width=15, height=10) 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#pdf("Feature_feature_relationships.pdf", width=15, height=10) 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#dev.off()

#pbmc1 <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### Define a function to identify cells that are outliers based on certain MAD from the median
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
}

### Read in the seurat object
seurat <- readRDS("pbmc_Seurat_Object.rds")

### Identify the cells that are outliers using MAD function
seurat <- mad_function(seurat = seurat, column = "percent.mt", number_mad = 3)
seurat <- mad_function(seurat = seurat, column = "nCount_RNA", number_mad = 2)
# I did not filter by nFeature_RNA
#seurat <- mad_function(seurat = seurat, column = "nFeature_RNA", number_mad = 3)

dim(seurat) #26412 12391

##### Remove the outliers ##### 
# NOTE: I have a doubt in here. Is my logic correct? percent.mt_mad == "NotOutlier" & nCount_RNA_mad == "Outlier" ?
seurat_final <- subset(seurat, subset = percent.mt_mad == "NotOutlier" & nCount_RNA_mad == "Outlier")

### Check the size of the seurat object after filtering to ensure that cells have been removed
print(seurat_final) #  An object of class Seurat  26412 features across 1671 samples within 1 assay

dim(seurat_final) #26412  1671

### Save the Seurat Object
saveRDS(seurat_final, "Seurat_singlets_QCfiltered.rds")
#seurat_final <- readRDS("Seurat_singlets_QCfiltered.rds")

# Visualize QC metrics as a violin plot and create a PDF of the FILTERED DATA
#pdf("Violin_plot_QC_filtered.pdf", width=15, height=10) 
VlnPlot(seurat_final, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

### Normalizing the data
pbmc <- NormalizeData(seurat_final, normalization.method = "LogNormalize", scale.factor = 10000)

### Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
#pdf("Variable_features.pdf", width=15, height=10) 
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#dev.off()

#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#pdf("DimPlot.pdf", width=15, height=10) 
DimPlot(pbmc, reduction = "pca")
#dev.off()

#pdf("DimHeatmap.pdf", width=15, height=10)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#dev.off()

#### Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

pdf("JackStrawPlot.pdf", width=15, height=10) 
JackStrawPlot(pbmc, dims = 1:15)
dev.off()

pdf("ElbowPlot.pdf", width=15, height=10) 
ElbowPlot(pbmc)
dev.off()

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

pdf("UMAP.pdf", width=15, height=10) 
DimPlot(pbmc, reduction = "umap")
dev.off()

saveRDS(pbmc, file = "pbmc_tutorial.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

pdf("FeaturePlot.pdf", width=15, height=10) 
FeaturePlot(pbmc, features = c("LEF1", "BCL11B", "VCAN", "PLCB1", "CDKN1C", "FCGR3A", "CCSER1", "HLA-DQA1","CAMK4",
    "IL7R"))
dev.off()

pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf("DoHeatmap.pdf", width=15, height=10)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()
    
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()