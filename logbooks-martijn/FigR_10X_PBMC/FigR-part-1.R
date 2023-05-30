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
RNAmat.example <- readRDS('logbooks/FigR/FigR_build_in_data/shareseq_skin_RNAnorm_final.rds')

dim(RNAmat)

ATAC.data.example <- readRDS('logbooks/FigR/FigR_build_in_data/shareseq_skin_SE_final.rds')

ATAC.data <- readRDS('../Users/Roya/Save-ArchR-Project.rds')
seRNA <- ArchR::import10xFeatureMatrix('output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', names='')
proj <- addGeneExpressionMatrix(input = ATAC.data, seRNA = seRNA, force = TRUE)

#ATAC.data <- ArchR::addPeakSet(ATAC.data)
#ATAC.data <- ArchR::addPeakMatrix(ATAC.data)
annoCols <- readRDS('../Users/Dilya/azimuth_results/pbmc_Seurat_Object_QCfiltered_Azimuth.rds')
annoCols.correct <- levels(annoCols@active.ident)

ATAC.data
dim(ATAC.se) # Peaks x Cells
colnames(ATAC.se)<-gsub("#","",colnames(ATAC.se))

#set.seed(123)
#cellsToKeep <- sample(colnames(RNAmat),size = 10000,replace = FALSE)
#ATAC.se <- ATAC.se[,cellsToKeep]
#RNAmat <- RNAmat[,cellsToKeep]
#colnames(RNAmat)
#RNAmat
# Remove genes with zero expression across all cells
RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]

cisAssign <- readRDS("logbooks/FigR_build_in_data/shareseq_skin_cisTopicPs.rds")
cisAssign <- cisAssign
dim(cisAssign) # Cells x Topics

all(cellsToKeep %in% rownames(cisAssign))

# Subset
cisAssign <- cisAssign[cellsToKeep,]

# Derive cell kNN using this
set.seed(123)
cellkNN <- get.knn(cisAssign,k = 3)$nn.index
dim(cellkNN)
rownames(cellkNN) <- colnames(ATAC.se)

#annoCols <- readRDS("logbooks/FigR_build_in_data/shareseq_skin_annoCols.rds")

#install.packages(
#   "ggplot2",
#   repos = c("http://rstudio.org/_packages",
#   "http://cran.rstudio.com"),
#   lib='../../R_libs')
#> You may also find it useful to restart R,
#> In RStudio, that's the menu Session >> Restart R

umap <- colData(ATAC.se) %>% as.data.frame() %>% ggplot(aes(UMAP1,UMAP2,color=cellAnnot)) + 
  geom_point(size=0.5) + scale_color_manual(values=annoCols.correct)+
  theme_classic() + guides(colour = guide_legend(override.aes = list(size=2)))

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/umap.png', umap, width=25, height=10, units="cm")

## The next section is about peak-gene association testing
## The section below takes a long time to compute.

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