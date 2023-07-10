setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')
getwd()

library(rhdf5)
h5createFile("../Users/Roya/Portal_input/ATAC_scores.h5")

library(HDF5Array)

gene_scores = loadHDF5SummarizedExperiment(dir="../Users/Roya/Portal_input", prefix="gene_scores_ATAC")

gene_scores

cell_types <- gene_scores@colData$seurat.object.meta.data.cell_type_lowerres

cell_types

barcodes <- gene_scores@colData$rownames.seurat.object.meta.data.

features <- rownames(gene_scores)

features

data <- gene_scores@assays@data$counts

h5createGroup("../Users/Roya/Portal_input/ATAC_scores.h5","matrix")

h5write(barcodes, "../Users/Roya/Portal_input/ATAC_scores.h5","matrix/barcodes")
h5write(features, "../Users/Roya/Portal_input/ATAC_scores.h5","matrix/features")
h5write(cell_types, "../Users/Roya/Portal_input/ATAC_scores.h5","matrix/celltypes")

h5ls("../Users/Roya/Portal_input/ATAC_scores.h5")
