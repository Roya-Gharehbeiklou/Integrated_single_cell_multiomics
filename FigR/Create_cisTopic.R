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

# save topic
saveRDS(cisTopicObject, paste0(outdir,'cisTopicObject.rds'))

# load topic again if new session
cisTopicObject <- readRDS(paste0(outdir, 'cisTopicObject.rds'))

# run several models testing different number of topics
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2, 4, 10, 16, 20, 30, 35), nCores=7)

# plot log likelihood
png(paste0(outdir,'/logLikelihoodByIter.png'))
logLikelihoodByIter(cisTopicObject)
dev.off()

# plot selected models
pdf(paste0(outdir,'select_model_cistopic.pdf'), width=8, height=15)
par(mfrow=c(2,1))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
dev.off()

png(paste0(outdir,'best_cistopic_model.png'), height=400, width=400, units='px')
# Select model based on maximum likelihood
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
dev.off()

# save best model cistopic
saveRDS(cisTopicObject, paste0(outdir, 'best_cisTopicObject.rds'))

best.cisTopic <- readRDS(paste0(outdir, 'best_cisTopicObject.rds'))

# run umap
best.cisTopic <- runUmap(best.cisTopic, target='cell', seed=123, method='Probability')

saveRDS(best.cisTopic, paste0(outdir, 'best_cisTopicObject.rds'))

# generate matrix for umap
cellassign <- modelMatSelection(best.cisTopic, 'cell', 'Probability')
topic.mat <- t(cellassign)
topic.mat <- as.matrix(topic.mat)

set.seed(123)
# Run umap on cells
DR <- umap(t(cellassign))
DRdist <- dist(DR$layout)

library(densityClust, lib.loc='Fig_R_libs')
library(ggplot2, lib.loc='Fig_R_libs')

# calculate clusters
dclust <- densityClust(DRdist,gaussian=T)

dclust <- findClusters(dclust, rho = 50, delta = 2.5)

saveRDS(dclust, paste0(outdir, 'dclust_topics.rds'))

# Check thresholds and plot
png(paste0(outdir, 'Threshold_topics.png'))
options(repr.plot.width=6, repr.plot.height=6)
plot(dclust$rho,dclust$delta,pch=20,cex=0.6,xlab='rho', ylab='delta')
points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
text(dclust$rho[dclust$peaks]-2,dclust$delta[dclust$peaks]+1.5,labels=dclust$clusters[dclust$peaks])
abline(v=50)
abline(h=2.5)
dev.off()

# Add cluster information
densityClust <- dclust$clusters
densityClust <- as.data.frame(densityClust)
rownames(densityClust) <- best.cisTopic@cell.names
colnames(densityClust) <- 'densityClust'
densityClust[,1] <- as.factor(densityClust[,1])
best.cisTopic <- addCellMetadata(best.cisTopic, densityClust)
# Add celltype annotation
colorBy <- ATAC.se@colData[2]
row.names(colorBy) <- ATAC.se@colData$rownames.seurat.object.meta.data.
colnames(colorBy) <- 'cisTopic cell type colored'
colorBy[,1] <- as.factor(colorBy[,1])
best.cisTopic <- addCellMetadata(best.cisTopic, colorBy)

colnames(best.cisTopic@cell.data)[4] <- "cisTopic colored by cell type"

# Plot all 35 topics with cell type umap
pdf(paste0(outdir,'35_topics.pdf'), width=20, height=20)
par(mfrow=c(6,6))
plotFeatures(best.cisTopic, method='Umap', target='cell', topic_contr='Probability', colorBy='cisTopic colored by cell type', cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
#plotFeatures(best.cisTopic, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

# plot non-transparant, 3 plots
pdf(paste0(outdir,'35_topics_nontransparant.pdf'), width=10, height=5)
par(mfrow=c(1,3))
plotFeatures(best.cisTopic, method='Umap', target='cell', topic_contr='Probability', colorBy='cisTopic colored by cell type', cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, topics=c(9, 23))
#plotFeatures(best.cisTopic, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

# plot png
pdf(paste0(outdir,'topics.png'), width = 15, height = 20)
par(mfrow=c(6,6))
plotFeatures(best.cisTopic, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
dev.off()

##############################################################
## Analysis of regulatory topics
## Not finished but useful in the future
#############################################################

#cisTopicObject <- getRegionsScores(best.cisTopic, method='NormTop', scale=TRUE)

#cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)

#library(org.Hs.eg.db, lib.loc='Fig_R_libs')
#library(enrichplot, lib.loc='cisTopic_libs')
#cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb='org.Hs.eg.db')
## Transcription factor motif enrichment

#cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)

#cisTopicObject <- getRegionsScores(best.cisTopic)

#cisTopicObject <- binarizecisTopics(cisTopicObject)

# hg38 not available
#cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, genome='hg19')

#pred.matrix <- predictiveDistribution(best.cisTopic)

#colnames(ATAC.se@colData)[2] <- "Cell.type"
#colnames(pred.matrix) <- ATAC.se@colData$Cell.type

# Obtain signatures
#Bulk_ATAC_signatures <- paste('../data/output/signature_data/', list.files(path_to_signatures), sep='')
#labels  <- gsub('._peaks.narrowPeak', '', list.files(path_to_signatures))
#best.cisTopic <- getSignaturesRegions(best.cisTopic, Bulk_ATAC_signatures, labels=labels, minOverlap = 0.4)

# Compute cell rankings (Reference time: 9 min)
#library(AUCell, lib.loc='Fig_R_libs')
#aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Check signature enrichment in cells (Reference time: 1 min)
#best.cisTopic <- signatureCellEnrichment(best.cisTopic, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

#cellAnnotation <- ATAC.se@colData[2]
#rownames(cellAnnotation) <- ATAC.se@colData$rownames.seurat.object.meta.data.
#best.cisTopic <- addCellMetadata(best.cisTopic, cell.data=cellAnnotation)

# Plot
#pdf(paste0(outdir,'per_cell_cistopic.pdf'), width = 15, height = 20)
#par(mfrow=c(3,3))
#plotFeatures(best.cisTopic, method='Umap', target='cell', topic_contr=NULL, colorBy='Cell.type', cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)

# To only keep unique peaks per signature
#library(plyr, lib.loc='Fig_R_libs')
##best.cisTopic@signatures <- llply(1:length(best.cisTopic@signatures), function (i) best.cisTopic@signatures[[i]][-which(best.cisTopic@signatures[[i]] %in% unlist(as.vector(best.cisTopic@signatures[-i])))]) 
#names(best.cisTopic@signatures) <- labels

#png(paste0(outdir, 'nCounts_topics.png'))
#plotFeatures(best.cisTopic, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('nCounts'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
#dev.off()

# png(paste0(outdir, 'densityClust_topics.png'))
# plotFeatures(best.cisTopic, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('densityClust'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
# dev.off()

# png(paste0(outdir, 'nAcc_topics.png'))
# plotFeatures(best.cisTopic, method='tSNE', target='cell', topic_contr=NULL, colorBy=c('nAcc'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
# dev.off()

# png(paste0(outdir, 'heatmap_topics.png'))
# cellTopicHeatmap(best.cisTopic, method='Probability', colorBy=c('densityClust'))
# dev.off()

# png(paste0(outdir, 'each_topic_tnse.png'))
# plotFeatures(best.cisTopic[1], method='tSNE', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
# dev.off()
