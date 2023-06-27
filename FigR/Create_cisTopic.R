setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')

library(lgr, lib.loc='Fig_R_libs')
library(cisTopic, lib.loc='Fig_R_libs')
library(png, lib.loc='Fig_R_libs')

cisTopicObject <- readRDS("../Users/Martijn/Integrated_single_cell_multiomics/FigR/output/cisTopicObject.rds")

# Based on: http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/10X_workflow.html
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(1:5, 7, 9, 11, 13, 15, 20, 25, 30, 35, 40), seed=987, nCores=9, burnin = 120, iterations = 150, addModels=FALSE)

png('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/best_cistopic_model.png')
plot.models <- selectModel(cisTopicObject, type='maximum')
dev.off()

saveRDS(cisTopicObject, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/best_cisTopicObject.rds')