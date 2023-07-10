setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')

library(lgr, lib.loc='Fig_R_libs')
library(cisTopic, lib.loc='Fig_R_libs')
library(png, lib.loc='Fig_R_libs')

outdir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/'

# Load ATAC SE
ATAC.se <- readRDS(paste0(outdir,'ATAC_se.rds'))

# Change rownames to correct format
colnames(ATAC.se@assays@data$counts) <- gsub("#","",colnames(ATAC.se@assays@data$counts))
rownames(ATAC.se@assays@data$counts) <- ATAC.se@rowRanges$interval

count.matrix <- ATAC.se@assays@data$counts

atac <- as.matrix(count.matrix)
atac <- as.data.frame(atac)

# we need to work out the names of the rownames, and replace - into : to match the chromosome:start-end required format by cisTopic
rownames(atac) <- make.unique(ATAC.se@rowRanges$interval)

# Create cisTopicObject
cisTopicObject <- createcisTopicObject(
  atac,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE
)

saveRDS(cisTopicObject, paste0(outdir,'cisTopicObject.rds'))

# Topics from: https://rdrr.io/github/aertslab/cisTopic/f/vignettes/10X_workflow.Rmd
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2, 5, 10, 15, 20, 25, 30), nCores=7)

png(paste0(outdir,'/logLikelihoodByIter.png'))
logLikelihoodByIter(cisTopicObject)
dev.off()

png(paste0(outdir,'best_cistopic_model.png'), height=400, width=400, units='px')
# Select model based on maximum likelihood
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
dev.off()

saveRDS(cisTopicObject, paste0(outdir, 'best_cisTopicObject.rds'))
