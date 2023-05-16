## This script runs FigR on build-in data
## Add config file w/ paths in final pipline
setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()
## Load required packages

library(Matrix, lib.loc='R_libs')
library(chromVAR, lib.loc = 'R_libs')
library(motifmatchr, lib.loc='R_libs')
library(ggplot2, lib.loc='R_libs')
library(FigR, lib.loc = 'R_libs')

# Downloaded data from:
# "https://s3.us-east-1.amazonaws.com/vkartha/FigR/FigR_SHAREseq.zip"

# This first secstion (until UMAP) is the preprocessing/EDA part

## More dependencies
library(FNN)

## Even more dependencies
library(doParallel)
library(BuenColors, lib.loc='R_libs')
library(BSgenome, lib.loc='R_libs')
library(BSgenome.Mmusculus.UCSC.mm10, lib.loc='R_libs')

ATAC.se <- readRDS('logbooks/FigR_build_in_data/shareseq_skin_SE_final.rds')
ATAC.se <- ATAC.se
print('done')
RNAmat <- readRDS("logbooks/FigR_build_in_data/shareseq_skin_RNAnorm_final.rds") # Normalized
RNAmat <- RNAmat
dim(ATAC.se) # Peaks x Cells

set.seed(123)
cellsToKeep <- sample(colnames(ATAC.se),size = 10000,replace = FALSE)

ATAC.se <- ATAC.se[,cellsToKeep]
RNAmat <- RNAmat[,cellsToKeep]

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

annoCols <- readRDS("logbooks/FigR_build_in_data/shareseq_skin_annoCols.rds")

#install.packages(
#   "ggplot2",
#   repos = c("http://rstudio.org/_packages",
#   "http://cran.rstudio.com"),
#   lib='../../R_libs')
#> You may also find it useful to restart R,
#> In RStudio, that's the menu Session >> Restart R

umap <- colData(ATAC.se) %>% as.data.frame() %>% ggplot(aes(UMAP1,UMAP2,color=cellAnnot)) + 
  geom_point(size=0.5) + scale_color_manual(values=annoCols)+
  theme_classic() + guides(colour = guide_legend(override.aes = list(size=2)))

ggsave('logbooks/FigR_plots/umap.png', umap, width=25, height=10, units="cm")

## The next section is about peak-gene association testing
## The section below takes a long time to compute.

library(pbmcapply, lib.loc="R_libs")


# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "mm10", # One of hg19, mm10 or hg38 
                           nCores = 21,
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)

head(cisCorr)
write.table(cisCorr, 'logbooks/FigR_plots/ciscorr.csv', quote=FALSE, sep="\t")