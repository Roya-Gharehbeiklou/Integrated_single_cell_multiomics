setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')

library(Matrix, lib.loc='Fig_R_libs')
library(chromVAR, lib.loc = 'Fig_R_libs')
library(motifmatchr, lib.loc='Fig_R_libs')
library(ggplot2, lib.loc='Fig_R_libs')
library(S4Vectors, lib.loc='Fig_R_libs')
library(GenomeInfoDb, lib.loc='Fig_R_libs')
library(FigR, lib.loc ='Fig_R_libs')

cisCorr <- read.csv('logbooks/FigR/FigR_data/ciscorr.csv', sep='\t')

cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

# Set cutoff to 2 for gaining results
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                         cutoff = 10, # No. sig peaks needed to be called a DORC
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=2)
ggsave('logbooks/FigR/FigR_plots/DORCS.png')

## Table significant associated peaks

# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

## Visualizing DORCs

ATAC.se <- readRDS('logbooks/FigR/FigR_build_in_data/shareseq_skin_SE_final.rds')
set.seed(123)
cellsToKeep <- sample(colnames(ATAC.se),size = 10000,replace = FALSE)

RNAmat <- readRDS("logbooks/FigR/FigR_build_in_data/shareseq_skin_RNAnorm_final.rds") # Normalized

ATAC.se <- ATAC.se[,cellsToKeep]
RNAmat <- RNAmat[,cellsToKeep]

# Remove genes with zero expression across all cells
RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]

cisAssign <- readRDS("logbooks/FigR/FigR_build_in_data/shareseq_skin_cisTopicPs.rds")
dim(cisAssign) # Cells x Topics

all(cellsToKeep %in% rownames(cisAssign))

# Subset
cisAssign <- cisAssign[cellsToKeep,]

# Derive cell kNN using this
library(FNN, lib.loc='Fig_R_libs')
set.seed(123)
cellkNN <- get.knn(cisAssign,k = 30)$nn.index
dim(cellkNN)
rownames(cellkNN) <- colnames(ATAC.se)

annoCols <- readRDS("logbooks/FigR/FigR_build_in_data/shareseq_skin_annoCols.rds")

# package pbmcapply not found, so loaded below
library(pbmcapply, lib.loc='Fig_R_libs')
dorcMat <- getDORCScores(ATAC.se = ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 16)

dim(dorcMat)

## Smoothing counts and visualize with expression
library(doParallel, lib.loc='Fig_R_libs')
# Smooth dorc scores using cell KNNs (k=30)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 16)

# Smooth RNA using cell KNNs
# This takes longer since it's all genes
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat,nCores = 16)

# Visualize on pre-computed UMAP
umap.d <- as.data.frame(colData(ATAC.se)[,c("UMAP1","UMAP2")])

# DORC score for Pxn
# ggrastr error so load module
library(ggrastr, lib.loc='Fig_R_libs')
library(BuenColors, lib.loc='Fig_R_libs')

dorcg <- plotMarker2D(umap.d,dorcMat.s,markers = c("Dlx3"),maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("Dlx3 DORC")

ggsave('logbooks/FigR/FigR_plots/dorc-Dlx3.png', dorcg)

# RNA for Dlx3
rnag <- plotMarker2D(umap.d,RNAmat.s,markers = c("Dlx3"),maxCutoff = "q0.99",colorPalette = "brewer_purple") + ggtitle("Dlx3 RNA")
ggsave('logbooks/FigR/FigR_plots/rna-Dlx3.png', rnag)

# Package for plotting side to side
# Installed using: devtools::install_github("thomasp85/patchwork")
library(patchwork, lib.loc='Fig_R_libs')

Dlx3.combined <- dorcg + rnag
ggsave('logbooks/FigR/FigR_plots/Dlx3-both.png', Dlx3.combined)

## TF gene associations 

#library(doParallel, lib.loc='R_libs')
#library(FigR, lib.loc='R_libs')
#library(foreach, lib.loc='R_libs')
#library(Matrix, lib.loc='R_libs')
#library(doSNOW, lib.loc='R_libs')

library(BSgenome.Mmusculus.UCSC.mm10, lib.loc='Fig_R_libs')

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

figR.d <- runFigRGRN(ATAC.se = ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "mm10",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s,
                     nCores=21)

saveRDS(figR.d, 'logbooks/FigR/FigR_data/GRN.rds')
