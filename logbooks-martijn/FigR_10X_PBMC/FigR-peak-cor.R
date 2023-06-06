## This script runs FigR on the 10X data
## Add config file w/ paths in final pipline
setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()

## For Singularity container
## Load required packages

library(Matrix)
library(chromVAR)
library(motifmatchr)
library(ggplot2)
library(S4Vectors)
library(GenomeInfoDb)
library(FigR)

## More dependencies
library(FNN)

## Even more dependencies
library(doParallel)
library(BuenColors)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(sp)
library(SeuratObject)
library(Seurat)
library(rhdf5)
library(ArchR)
library(SummarizedExperiment)

## No container 

#library(Matrix, lib.loc='Fig_R_libs')
#library(chromVAR, lib.loc='Fig_R_libs')
#library(motifmatchr, lib.loc='Fig_R_libs')
#library(ggplot2, lib.loc='Fig_R_libs')
#library(S4Vectors, lib.loc='Fig_R_libs')
#library(GenomeInfoDb, lib.loc='Fig_R_libs')
#library(FigR, lib.loc = 'R_libs')

# This first secstion (until UMAP) is the preprocessing/EDA part

## More dependencies
#library(FNN, lib.loc='Fig_R_libs')

## Even more dependencies
#library(doParallel, lib.loc='Fig_R_libs')
#library(BuenColors, lib.loc='Fig_R_libs')
#library(BSgenome, lib.loc='Fig_R_libs')
#library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='Fig_R_libs')
#library(sp, lib.loc='ArchR_libs')
#library(SeuratObject, lib.loc='ArchR_libs')
#library(Seurat, lib.loc='ArchR_libs')
#ibrary(rhdf5, lib.loc='ArchR_libs')
#library(ArchR, lib.loc='ArchR_libs')
#ibrary(SummarizedExperiment, lib.loc='Fig_R_libs')

addArchRThreads(threads = 21) 

seurat.object <- readRDS('../Users/Dilya/azimuth_results/pbmc_Seurat_Azimuth_for_figR.rds')
RNAmat <- seurat.object@assays$RNA@data
dim(RNAmat)

ATAC.data.example <- readRDS('logbooks/FigR/FigR_build_in_data/shareseq_skin_SE_final.rds')
#ATAC.se <- ArchR::import10xFeatureMatrix('output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', names='')
#ATAC.se

ATAC.data <- readRDS('../Users/Roya/Save-ArchR-Project_subSet_QC_Frip.rds')
getAvailableMatrices(ATAC.data)

# Add feature matrix to ArchR project
fL <- getFragmentsFromProject(ATAC.data)
granges <- fL[[1]]

addArchRThreads(threads = 21) 

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
annoCols <- annoCols[, rownames(annoCols@meta.data) %in% barcodes.atac]
annoCols <- annoCols[, rownames(annoCols@meta.data) %in% barcodes.rnamat]

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

celltypes.barcode <- data.frame(rownames(seurat.object@meta.data), seurat.object@meta.data$predicted.celltype.l2)

ATAC.se@colData <- DataFrame(celltypes.barcode)

library(pbmcapply, lib.loc="R_libs")
#library(BSgenome, lib.loc='Fig_R_libs')
#library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='Fig_R_libs')

# Get correct ATAC.se levels
#levels.atac.data.example <- seqnames(ATAC.data.example)

#library(gdata)
ATAC.se <- ATAC.se[startsWith(rowData(ATAC.se)$interval, 'c')]

# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "hg38", # One of hg19, mm10 or hg38 
                           nCores = 21,
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)

write.table(cisCorr, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/ciscor.csv', quote=FALSE, sep="\t")
cisCorr <- read.csv('../Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/ciscor.csv', sep="\t")
head(cisCorr)

cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                         cutoff = 3, # No. sig peaks needed to be called a DORC
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=2)

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/dors.png')

# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

dorcMat <- getDORCScores(ATAC.se = ATAC.se.test, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 21)

dim(dorcMat)

# Smooth dorc scores using cell KNNs (k=30)
#dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 4)

# Smooth RNA using cell KNNs
# This takes longer since it's all genes
#RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat,nCores = 4)

library(BSgenome, lib.loc='Fig_R_libs')
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='ArchR_libs')

figR.d <- runFigRGRN(ATAC.se = ATAC.se.test, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat,
                     rnaMat = RNAmat, 
                     nCores = 21)

saveRDS(cisCorr, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/figRGRN.rds', quote=FALSE, sep="\t")
