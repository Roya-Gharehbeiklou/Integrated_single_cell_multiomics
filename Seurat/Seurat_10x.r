
# module load RPlus
# R
# Data
#/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/output/raw_feature_bc_matrix

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

# Visualize QC metrics as a violin plot and create a PDF of the FILTERED DATA
#pdf("Violin_plot_QC_filtered.pdf", width=15, height=10) 
VlnPlot(seurat_final, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()