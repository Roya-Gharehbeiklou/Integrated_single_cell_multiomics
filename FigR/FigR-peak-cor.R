## This script runs FigR on the 10X data
## Add config file w/ paths in final pipline
setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()

## For Singularity container

# library(chromVAR)
# library(S4Vectors)
# library(GenomeInfoDb)
# library(FigR)

# # ## More dependencies
# library(FNN)

# ## Even more dependencies
# library(BuenColors)
# library(BSgenome)1: package ‘spatstat.core’ is not available for this version of R
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(SeuratObject)
# library(Seurat)
# library(ArchR)
# library(SummarizedExperiment)

## No container 

library(Matrix, lib.loc='Fig_R_libs')
library(chromVAR, lib.loc='Fig_R_libs')
library(motifmatchr, lib.loc='Fig_R_libs')
library(ggplot2, lib.loc='Fig_R_libs')
library(S4Vectors, lib.loc='Fig_R_libs')
library(GenomeInfoDb, lib.loc='Fig_R_libs')
library(dplyr, lib.loc='Fig_R_libs')
library(FigR, lib.loc = 'Fig_R_libs')
library(FNN, lib.loc='Fig_R_libs')
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

# Setting ArchR threads
addArchRThreads(threads = 16) 

seurat.object <- readRDS('../Users/Dilya/azimuth_results/pbmc_Seurat_Azimuth_for_figR.rds')
RNAmat <- seurat.object@assays$RNA@data

# Log10 normalize RNA count matrix
RNAmat <- NormalizeData(RNAmat)
dim(RNAmat)

ATAC.data <- readRDS('../Users/Roya/Save-ArchR-Project_subSet_QC_Frip.rds')
getAvailableMatrices(ATAC.data)

ATACgene.score.matrix <- getMatrixFromProject(
  ArchRProj = ATAC.data,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

colnames(ATACgene.score.matrix) <- gsub("pbmc_granulocyte_sorted_10k_HG38#","",colnames(ATACgene.score.matrix))

# Load summarized experiment
ATAC.se <- ArchR::import10xFeatureMatrix('output/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5', names='')
colnames(ATAC.se)<-gsub("#","",colnames(ATAC.se))

# Change rownames ATAC.data object
rownames(ATAC.data@cellColData) <- gsub("pbmc_granulocyte_sorted_10k_HG38#","",rownames(ATAC.data@cellColData))

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
seurat.object <- seurat.object[, rownames(seurat.object@meta.data) %in% barcodes.atac]
seurat.object <- seurat.object[, rownames(seurat.object@meta.data) %in% barcodes.rnamat]

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

celltypes.barcode <- data.frame(rownames(seurat.object@meta.data), seurat.object@meta.data$cell_type_lowerres)

ATAC.se@colData <- DataFrame(celltypes.barcode)

dim(ATAC.se)
dim(RNAmat)

###############
# CREATING UMAP
###############

# UMAP
library(lgr, lib.loc='Fig_R_libs/')
library(cisTopic, lib.loc='Fig_R_libs/')
library(SingleCellExperiment, lib.loc='Fig_R_libs/')
library(scuttle, lib.loc='Fig_R_libs/')
library(scater, lib.loc='Fig_R_libs/')

annotation <- ATAC.se@colData$seurat.object.meta.data.cell_type_lowerres
write.csv(annotation)

ATAC.singlecell <- as(ATAC.se, "SingleCellExperiment")
ATAC.singlecell <- scater::logNormCounts(ATAC.singlecell)
u <- uwot::umap(as.matrix(t(counts(ATAC.singlecell))), n_neighbors=2)
reducedDim(ATAC.singlecell, "UMAP_uwot") <- u
celltypes.barcode$UMAP1 <- reducedDim(ATAC.singlecell, "UMAP_uwot")[,1]
celltypes.barcode$UMAP2 <- reducedDim(ATAC.singlecell, "UMAP_uwot")[,2]

ATAC.singlecell@colData <- DataFrame(celltypes.barcode)

library(ggplot2, lib.loc='Fig_R_libs')

colData(ATAC.singlecell) %>% as.data.frame() %>% ggplot(aes(UMAP1,UMAP2)) + 
  geom_point(size=0.5) + scale_color_manual(values=annotation)+
  theme_classic() + guides(colour = guide_legend(override.aes = list(size=2)))
ggsave('../Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/UMAP.png')

##########
# END UMAP
##########

# Getting chromosome names
ATAC.chromosomes <- rowData(ATAC.se)$interval

# Different starting name chromosomes
deviating.names <- ATAC.chromosomes[!startsWith(ATAC.chromosomes, 'c')]
normal.chroms <- ATAC.chromosomes[startsWith(ATAC.chromosomes, 'c')]
#main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)


# Fetching chromosome starting with a 'c'
ATAC.se <- ATAC.se[startsWith(rowData(ATAC.se)$interval, 'c')]

# Subset ATACgene.score.matrix also
ATACgene.score.matrix <- ATACgene.score.matrix[,colnames(ATACgene.score.matrix) %in% cellsToKeep]
saveRDS(ATACgene.score.matrix, '../Users/Roya/Portal_input/gene_activity_ATAC.rds')

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

###################
# PEAK CORRELATIONS
###################

library(pbmcapply, lib.loc="Fig_R_libs")
# Don't run interactively
cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "hg38",
                           nCores = 16,
                           p.cut = NULL,
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

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/dorcs.png')

# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 21)

dim(dorcMat)

#################################
# Create cisTopicObject
#################################

library(lgr, lib.loc='Fig_R_libs')
library(cisTopic, lib.loc='Fig_R_libs')

# Change rownames to correct format
colnames(ATAC.se@assays@data$counts)<-gsub("#","",colnames(ATAC.se@assays@data$counts))
rownames(ATAC.se@assays@data$counts) <- ATAC.se@rowRanges$interval

count.matrix <- ATAC.se@assays@data$counts

atac <- as.matrix(count.matrix)
atac <- as.data.frame(atac)

# we need to work out the names of the rownames, and replace - into : to match the chromosome:start-end required format by cisTopic
rownames(atac) <- make.unique(rowRanges(ATAC.se)$interval)

cisTopicObject <- createcisTopicObject(
  atac,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE
)

saveRDS(cisTopicObject, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/cisTopic/cisTopicObject.rds')


## Load cisTopic and fetch matrix

cisTopicObject <- readRDS("../Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/cisTopic/best_cisTopicObject.rds")
topic.mat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
topic.mat <- t(topic.mat)
topic.mat <- as.matrix(topic.mat)

# Derive cell kNN using this
set.seed(123)
cellkNN <- get.knn(topic.mat,k = 30)$nn.index
dim(cellkNN)
rownames(cellkNN) <- cellsToKeep

colnames(dorcMat) <- rownames(cellkNN)
# Smooth dorc scores using cell KNNs (k=30)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 4)

# Smooth RNA using cell KNNs
# This takes longer since it's all genes
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat,nCores = 4)

#######################################################################
# CHANGED THE runFigRGRN SO IT WORKS, PROBLEM WAS THE LIBRARIES LOADING
#######################################################################

runFigRGRN <- function(ATAC.se, # SE of scATAC peak counts. Needed for chromVAR bg peaks etc.
                       dorcK=30, # How many dorc kNNs are we using to pool peaks
                       dorcTab, # peak x DORC connections (should contain indices relative to peaks in ATAC.se)
                       n_bg=50, # No. of background peaks to use for motif enrichment Z test
                       genome, # One of mm10, hg19, hg38, with no default
                       dorcMat, # Expect smoothed
                       rnaMat, # Expect smoothed
                       dorcGenes=NULL, # If only running on a subset of genes
                       nCores=1
){
  # Must be matched data
  stopifnot(all.equal(ncol(dorcMat),ncol(rnaMat)))

  # Expects "Gene" / "Peak" in dorcTab
  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")

  if(all(grepl("chr",dorcTab$Peak,ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")

    if(!(all(grepl("chr",rownames(ATAC.se),ignore.case = TRUE))))
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")

    if(!all(dorcTab$Peak %in% rownames(ATAC.se)))
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  } else{
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
  # If using index, make sure no indices are outside range of SE
    if(max(dorcTab$Peak) > nrow(ATAC.se))
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }


  if(is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  } else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }

  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
  rownames(DORC.knn) <- rownames(dorcMat)

  if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
  }

  # Set data subfolder path
  packagePath <- find.package("FigR", 'Fig_R_libs', quiet = TRUE)

  if(grepl("hg",genome)){
    pwm <- readRDS(paste0(packagePath,"/data/cisBP_human_pfms_2021.rds"))
  } else {
    pwm <- readRDS(paste0(packagePath,"/data/cisBP_mouse_pfms_2021.rds"))
  }

  # Old motif naming convention
  if(all(grepl("_",names(pwm),fixed = TRUE)))
     names(pwm) <- FigR::extractTFNames(names(pwm))

  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]
  myGeneNames <- gsub(x = rownames(rnaMat),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
  rownames(rnaMat) <- myGeneNames

  # Only non-zero expression TFs (also found in rnaMat)
  motifsToKeep <- intersect(names(pwm),myGeneNames)

  # This has to be done on the full SE (same peakset used as input to dorc calling)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se,pwms = pwm[motifsToKeep],genome=genome)

  # Keep TFs with some peak x motif match
  motif_ix <- motif_ix[,Matrix::colSums(assay(motif_ix))!=0]

  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat,1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.mat),rowRanges = granges(ATAC.se))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  } else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }

  # For each DORC, do motif enrichment among dorc sig Peaks, and correlation of DORC accessibility (smoothed) to TF RNA levels

  cat("Testing ",length(motifsToKeep)," TFs\n")
  cat("Testing ",nrow(dorcMat)," DORCs\n")
  library(doParallel, lib.loc='Fig_R_libs')
  if(nCores > 1)
    message("Running FigR using ",nCores," cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  mZtest.list <- foreach(g=dorcGenes,
                         .options.snow = opts) %dopar%   {
                         library("motifmatchr", lib.loc='Fig_R_libs')
                         library("ggplot2", lib.loc='Fig_R_libs')
                         library("Matrix", lib.loc='Fig_R_libs')
                         library("chromVAR", lib.loc='Fig_R_libs')
                         library("S4Vectors", lib.loc='Fig_R_libs')
                         library("GenomeInfoDb", lib.loc='Fig_R_libs')
                         library("FigR", lib.loc='Fig_R_libs')

                           # Take peaks associated with gene and its k neighbors
                           # Pool and use union for motif enrichment
                           DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% c(g,rownames(dorcMat)[DORC.knn[g,]])])

                           if(usePeakNames)
                             DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks) # Convert to index relative to input

                           mZ <- FigR::motifPeakZtest(peakSet = DORCNNpeaks,
                                                bgPeaks = bg,
                                                tfMat = assay(motif_ix))

                           mZ <- mZ[,c("gene","z_test")]
                           colnames(mZ)[1] <- "Motif"
                           colnames(mZ)[2] <- "Enrichment.Z"
                           mZ$Enrichment.P <- 2*pnorm(abs(mZ$Enrichment.Z),lower.tail = FALSE) # One-tailed
                           mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
                           mZ <- cbind("DORC"=g,mZ)
                           # Correlate smoothed dorc with smoothed expression, with spearman
                           corr.r <- cor(dorcMat[g,],t(as.matrix(rnaMat[mZ$Motif,])),method = "spearman")
                           stopifnot(all.equal(colnames(corr.r),mZ$Motif))

                           mZ$Corr <- corr.r[1,] # Correlation coefficient
                           mZ$Corr.Z <- scale(mZ$Corr,center = TRUE,scale = TRUE)[,1] # Z-score among all TF correlations
                           mZ$Corr.P <- 2*pnorm(abs(mZ$Corr.Z),lower.tail = FALSE) # One-tailed
                           mZ$Corr.log10P <- sign(mZ$Corr.Z)*-log10(mZ$Corr.P)
                           return(mZ)
                         }
  cat("Finished!\n")
  cat("Merging results ..\n")
  # Merge and save table for downstream filtering and plotting (network)
  TFenrich.d <- do.call('rbind',mZtest.list)
  dim(TFenrich.d)
  rownames(TFenrich.d) <- NULL

  # Make combined score based on multiplication
  # Here, we only sign by corr
  # Since sometimes we lose digit precision (1 - v small number is 1, instead of 0.9999999..)
  # Use Rmpfr, increase precision limits above default (100 here)
  TFenrich.d <- TFenrich.d %>% dplyr::mutate("Score"=sign(Corr)*as.numeric(-log10(1-(1-Rmpfr::mpfr(Enrichment.P,100))*(1-Rmpfr::mpfr(Corr.P,100)))))
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  TFenrich.d
}

library(BiocManager, lib.loc='Fig_R_libs')
library(BSgenome, lib.loc='Fig_R_libs')
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc='Fig_R_libs')

.libPaths('Fig_R_libs')

if(!require(BSgenome.Hsapiens.UCSC.hg38)){
    install.packages("BSgenome.Hsapiens.UCSC.hg38")
    library(BSgenome.Hsapiens.UCSC.hg38)
}


figR.d <- runFigRGRN(ATAC.se = ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s, 
                     nCores = 21)

saveRDS(figR.d, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/figRGRN.rds')

library(BuenColors, lib.loc='Fig_R_libs')
require(ggrastr)

figR.scatter <- figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  xlim(-5, 10) +
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/TF-DORC-scatter.png')

figR.d.subset <- figR.d[figR.d$Score > 0,]

figR.scatter.subset <- figR.d.subset %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=1,shape=16) + 
  theme_classic() + 
  xlim(-5, 10) +
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/TF-DORC-scatter-subset.png')

rank.drivers <- rankDrivers(figR.d,rankBy = "meanScore",interactive = FALSE)
rank.drivers + geom_text_repel() + geom_label_repel() + cowplot::theme_cowplot() + labs()
ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/rank_drivers.png', rank.drivers, width=20, height=6)

plt.drivers <- plotDrivers(figR.d,score.cut = 2, marker="TSG-AG1")

ggsave('logbooks/FigR/FigR_plots/plot-drivers.png', plt.drivers)

library(ComplexHeatmap)
library(png)
png("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/FigR_10X_PBMC/output/complex-heatmap.png", width=5,height=5,units="in",res=1200)

plotfigRHeatmap(figR.d = figR.d,
                score.cut = 0.1,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
dev.off()

library(networkD3, lib.loc='Fig_R_libs')
d3 <- plotfigRNetwork(figR.d,
                score.cut = 0.1,
                weight.edges = TRUE)

save_d3_html(
  d3,
  'network_tutorial_rna_n_atac.html',
)
