# Set working directory
setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya")

# Set the library path
.libPaths('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs')

library(ArchR, lib.loc="/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs")
#ArchR::installExtraPackages()

set.seed(1)
addArchRGenome("hg38")



#Get Input Fragment Files
inputFiles <- getInputFiles("pbmc_granulocyte_sorted_10k")[1]

names(inputFiles) <- "pbmc_granulocyte_sorted_10k"

# Create Arrow Files and setting QC values 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFile,
  sampleNames = names(inputFile),
  minTSS = 4,
  minFrags = 1000,
  force = TRUE
)


ArrowFiles

# Inferring Doublets
# Had to load IRanges seperately
library(IRanges, lib.loc='/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs')
# Problem with seqnames, had to load it this way
seqnames <- GenomicRanges::seqnames
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 1, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

#Import scRNA
seRNA <- import10xFeatureMatrix(
  input = c("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"),
  names = c("pbmc_granulocyte_sorted_10k")
)

seRNA

#Add scRNA
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

#Filter Cells ###  QC parameters. # FRiP QC value is missing
proj <- proj[proj$TSSEnrichment >= 4 & proj$nFrags >= 1000 & !is.na(proj$Gex_nUMI)]


#Doublet Filtration. Currently disabled just for tutorial. If you want to filter doublets uncomment below.
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj, cutEnrich = 1) # is cutEnrich = 1 ok? 

#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

#LSI-RNA
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)

proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

#Add Clusters
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)

p2 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)

p3 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)

#Print Plots
p1
p2
p3

#Save Plot
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)

#Print Session Info
sessionInfo()



