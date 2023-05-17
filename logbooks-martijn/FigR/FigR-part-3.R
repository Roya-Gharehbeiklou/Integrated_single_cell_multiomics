figR.heatmap <- figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

ggsave('logbooks/FigR/FigR_plots/heatmap.png', figR.heatmap)

# Interactive plotting
#ranks <- rankDrivers(figR.d,score.cut = 2,rankBy = "nTargets",interactive = TRUE)


plt.drivers <- plotDrivers(figR.d, marker='Acvr1')
ggsave('logbooks/FigR/FigR_plots/plot-drivers.png', plt.drivers)

library(ComplexHeatmap, lib.loc='R_libs')
complex.heatmap <- plotfigRHeatmap(figR.d = figR.d,
                score.cut = 1,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                )

ggsave('logbooks/FigR/FigR_plots/complexheatmap.png', complex.heatmap)

library(networkD3, lib.loc='R_libs')
library(magrittr)
library(htmlwidgets)
networks <- plotfigRNetwork(figR.d)
save(networks, file = "logbooks/FigR_plots/network.rda")
#saveWidget('logbooks/FigR_plots/network.html')