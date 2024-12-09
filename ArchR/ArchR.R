#srun --cpus-per-task=1 --mem=10gb --nodes=1 --qos=interactive --time=05:00:00 --pty bash -i
# Activate conda environment that has the macs2 function
# cd .conda/envs
# conda activate mac2_env/
# module load RPlus
# R
# run in /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs

# Erase all existing variables
rm(list=ls())

# Set working directory
setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya")

# Set the library path
.libPaths('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs')

library(ArchR, lib.loc="/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs")
#ArchR::installExtraPackages()
library('BSgenome.Hsapiens.UCSC.hg38')

set.seed(1)
# settting the threads to one, if not there is a problem when inferring the doublets and creating the arrow files. read: https://github.com/GreenleafLab/ArchR/issues/218  for more information
#addArchRThreads(threads = 1)
addArchRGenome("hg38")

#Get Input Fragment Files
inputFile <- getInputFiles("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/output")[1]
inputFile

names(inputFile) <- "pbmc_granulocyte_sorted_10k_HG38"

# Create Arrow Files and setting QC values 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFile,
  sampleNames = names(inputFile),
  minTSS = 4,
  minFrags = 1000,
  force = TRUE
)

ArrowFiles 
#saveRDS(ArrowFiles, file='ArrowFiles.rds')
#ArrowFiles<- readRDS('ArrowFiles.rds')

# Inferring Doublets
# Had to load IRanges seperately
library(IRanges, lib.loc='/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/ArchR_libs')
# Problem with seqnames, had to load it this way
seqnames <- GenomicRanges::seqnames

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

# Creating an ArchRProject 
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage. Output: Copying ArrowFiles to Ouptut Directory! If you want to save disk space set copyArrows = FALSE
)
getAvailableMatrices(proj) #  "GeneScoreMatrix" "TileMatrix"

# filter putative doublets 
proj <- filterDoublets(ArchRProj = proj)
# Output: Filtering 1479 cells from ArchRProject! pbmc_granulocyte_sorted_10k_HG38 : 1479 of 12162 (12.2%)

#Dimensionality Reduction and Clustering
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

# Visualizing in a 2D UMAP Embedding
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

# color by 'sample'
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

# color by “Clusters”
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

# To save an editable vectorized version of this plot
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Saving ArchRProject
#saveArchRProject(ArchRProj = proj)

# loadArchRProject
#proj <- loadArchRProject(path ='/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya')
getAvailableMatrices(proj)

# Calling Peaks w/ Macs2. 
# pathToMacs2 <- findMacs2()

# projTmp <- addReproduciblePeakSet(
#     ArchRProj = proj, 
#     groupBy = "Clusters", 
#     pathToMacs2 = pathToMacs2
# )
## Error: /Roya/PeakCalls/InsertionBeds/C1._.Rep1-1_summits.bed' does not exist or is non-readable.

# Calling Peaks w/ TileMatrix
projTmp <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters",
    peakMethod = "Tiles",
    method = "p"
)

getPeakSet(projTmp)

# Save
saveArchRProject(ArchRProj = projTmp)
# loadArchRProject
# projTmp<- loadArchRProject(path ='/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya')

# Adding Peak Matrix
proj_peak_matrix <- addPeakMatrix(projTmp)
getAvailableMatrices(proj_peak_matrix)

#Motif Deviations
if("Motif" %ni% names(proj_peak_matrix@peakAnnotation)){
    proj_peak_matrix <- addMotifAnnotations(ArchRProj = proj_peak_matrix, motifSet = "cisbp", name = "Motif")
}

proj_peak_matrix <- addBgdPeaks(proj_peak_matrix)

# Creating a deviations matrix in each of our Arrow files called “MotifMatrix”
proj_peak_matrix <- addDeviationsMatrix(
  ArchRProj = proj_peak_matrix, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj_peak_matrix, name = "MotifMatrix", plot = TRUE)
plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_peak_matrix, addDOC = FALSE)


#Motif Footprinting
motifPositions <- getPositions(proj_peak_matrix)

motifPositions
proj_peak_matrix <- addGroupCoverages(ArchRProj = proj_peak_matrix, force=TRUE)
# seFoot <- getFootprints(
#   ArchRProj = proj_peak_matrix, 
#   positions = motifPositions[markerMotifs], 
#   groupBy = "Clusters"
# )
#Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'i' in selecting a method for function '[': object 'markerMotifs' not found
# plotFootprints(
#   seFoot = seFoot,
#   ArchRProj = proj_peak_matrix, 
#   normMethod = "Subtract",
#   plotName = "Footprints-Subtract-Bias",
#   addDOC = FALSE,
#   smoothWindow = 5
#)


# FRiP (fraction of reads in peaks: 
#how many of total reads aligned to regions defined as “peaks” after pseudo-bulking) > 0.5
proj_peak_matrix <- addCoAccessibility(
    ArchRProj = proj_peak_matrix,
    reducedDims = "IterativeLSI"
)

# proj_peak_matrix <- addPeak2GeneLinks(
#     ArchRProj = proj_peak_matrix,
#     reducedDims = "IterativeLSI"
# )
#Error in addPeak2GeneLinks(ArchRProj = proj_peak_matrix, reducedDims = "IterativeLSI") : 
#  GeneIntegrationMatrix not in AvailableMatrices

seGroupMotif <- getGroupSE(ArchRProj = proj_peak_matrix, useMatrix = "MotifMatrix", groupBy = "Clusters")
saveRDS(seGroupMotif, file='seGroupMotif_proj.rds')

# Save
saveArchRProject(ArchRProj = proj_peak_matrix)
#load
#proj_peak_matrix<- loadArchRProject(path ='/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya')

# Dimension 
proj_peak_matrix # 10683

# Get a subset of the project based on the Frip value QC
ids <- which(proj_peak_matrix$FRIP>=0.5) #10051
cellsPass <- proj_peak_matrix$cellNames[ids]
proj_subset = proj_peak_matrix[cellsPass, ]

# Saving the project
saveArchRProject(ArchRProj = proj_subset)
















##########################################
seGeneScore_proj <- getGroupSE(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "Clusters")
saveRDS(seGeneScore_proj, file='seGeneScore_proj.rds')
seGeneScore_proj <- readRDS('seGeneScore_proj.rds')

######

# Export Group Summarized Experiment
seGeneScore <- getGroupSE(ArchRProj = proj_peak_matrix, useMatrix = "GeneScoreMatrix", groupBy = "Clusters")

# Motif Footprinting
motifPositions <- getPositions(proj_peak_matrix)
motifPositions

seGroupMotif <- getGroupSE(ArchRProj = proj_peak_matrix, useMatrix = "MotifMatrix", groupBy = "Clusters")



#######
#how data should look like

prueba<-readRDS('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/logbooks/FigR/FigR_build_in_data/shareseq_skin_SE_final.rds')

### path







prueba@colData
proj@cellColData