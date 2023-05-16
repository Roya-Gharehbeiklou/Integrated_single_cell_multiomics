setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')

library(ArchR, lib.loc='ArchR_libs')
ArchR::installExtraPackages()

set.seed(1)

addArchRThreads(threads = 16) 

addArchRGenome("hg38")

#Get Input Fragment Files
inputFile <- getInputFiles("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/output")[1]
inputFile

names(inputFile) <- "pbmc_granulocyte_sorted_10k_HG38"

#Create Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFile,
  sampleNames = names(inputFile),
  minTSS = 4,
  minFrags = 1000,
  force = TRUE
)
