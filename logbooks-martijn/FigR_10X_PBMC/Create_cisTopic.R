setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')

library(lgr, lib.loc='Fig_R_libs')
library(cisTopic, lib.loc='Fig_R_libs')

cisTopicObject <- readRDS("../Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/cisTopic/cisTopicObject.rds")

# Based on: http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/10X_workflow.html
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2, 5, 10, 15, 20, 25, 30, 35, 40), seed=987, nCores=9, burnin = 120, iterations = 150, addModels=FALSE)

cisTopicObject <- selectModel(cisTopicObject, type='maximum')
saveRDS(cisTopicObject, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/logbooks-martijn/logbooks/best_cisTopicObject.rds')