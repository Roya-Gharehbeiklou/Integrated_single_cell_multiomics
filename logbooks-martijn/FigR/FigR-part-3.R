setwd("/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data")

figR.d <- readRDS('logbooks/FigR/FigR_data/GRN.rds')

.libPaths('Fig_R_libs')

library("ggplot2")
library("dplyr")
library("BuenColors")
library("FigR")

figR.heatmap <- figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

ggsave('logbooks/FigR/FigR_plots/heatmap.png', figR.heatmap)

# This is interactive
ranks <- rankDrivers(figR.d,score.cut = 2,rankBy = "nTargets",interactive = TRUE)
saveRDS(ranks, 'logbooks/FigR/FigR_data/rank_drivers.rds')

plt.drivers <- plotDrivers(figR.d,score.cut = 2,marker = "Lef1")
ggsave('logbooks/FigR/FigR_plots/plot-drivers.png', plt.drivers)

library(ComplexHeatmap)
library(png)
png("logbooks/FigR/FigR_plots/complex-heatmap.png")

plotfigRHeatmap(figR.d = figR.d,
                score.cut = 2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )
dev.off()

library(networkD3)
library(magrittr)
library(htmlwidgets)
library(gplots)
networks <- plotfigRNetwork(figR.d)
save(networks, file = "logbooks/FigR/FigR_data/network-plot.rda")