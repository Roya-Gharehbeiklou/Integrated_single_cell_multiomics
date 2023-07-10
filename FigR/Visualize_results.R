setwd('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/data/')

# Load modules without issues
library(chromVAR, lib.loc='Fig_R_libs')
library(Matrix, lib.loc='Fig_R_libs')
library(S4Vectors, lib.loc='ArchR_libs')
library(GenomeInfoDb, lib.loc='ArchR_libs')
library(ggplot2, lib.loc='ArchR_libs')
library(motifmatchr, lib.loc='ArchR_libs')
library(FigR, lib.loc='Fig_R_libs')
library(BuenColors, lib.loc='Fig_R_libs')
library(ggrastr, lib.loc='Fig_R_libs')

sourcedir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/'

figR.d.unsmoothed <- readRDS(paste0(sourcedir, 'figRGRN_unsmoothed.rds'))
figR.d <- readRDS(paste0(sourcedir, 'figRGRN.rds'))


figR.scatter <- figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=1,shape=16) + 
  theme_classic() + 
  xlim(-5, 10) +
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

ggsave('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Martijn/Integrated_single_cell_multiomics/FigR/output/TF-DORC-scatter_30_topics_unsmoothed.png')

library(ComplexHeatmap, lib.loc="Fig_R_libs")
pdf(paste0(sourcedir,'heatmap.pdf'), width=10, height=5)

plot1 <- plotfigRHeatmap(figR.d = figR.d.unsmoothed,
                score.cut=0.2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                ) %>% 
                draw() %>% 
                grid.grabExpr()

plot2 <- plotfigRHeatmap(figR.d = figR.d,
                score.cut=0.2,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
                ) %>% 
                draw() %>% 
                grid.grabExpr()

library(ggplotify, lib.loc='Fig_R_libs')              
plot1 <- as.ggplot(plot1) + ggtitle(expression(paste(bold("A. "), "GRN without using cisTopic"))) +
  theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent'),
         plot.title = element_text(size=10)
       )

plot2 <- as.ggplot(plot2) + ggtitle(expression(paste(bold("B. "), "GRN using cisTopic"))) +
  theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent'),
         plot.title = element_text(size=10)
       )

wrap_plots(plot1, plot2) + plot_annotation(
  title = 'FigR GRN topic modeling impact',
  theme = theme(plot.title = element_text(hjust = 0.5)))
dev.off()

png(paste0(sourcedir,'heatmap_smoothed.png'))
as.ggplot(plot1) + ggtitle(expression(paste(bold("A. "), "GRN without using cisTopic"))) +
  theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent'),
         plot.title = element_text(size=10)
       )
dev.off()


figR.d.ordered <- figR.d[order(figR.d$Score,decreasing=TRUE),][c(1:2,10)]
figR.d.top.ten <- rbind(head(figR.d.ordered, 5), tail(figR.d.ordered, 5))
figR.d.unsmoothed.ordered <- figR.d.unsmoothed[order(figR.d.unsmoothed$Score,decreasing=TRUE),][c(1:2,10)]
figR.d.unsmoothed.top.ten <- rbind(head(figR.d.unsmoothed.ordered, 5), tail(figR.d.unsmoothed.ordered, 5))
figR.d.top.ten <- figR.d.top.ten[order(abs(figR.d.top.ten$Score), decreasing=TRUE),]
figR.d.unsmoothed.top.ten <- figR.d.unsmoothed.top.ten[order(abs(figR.d.unsmoothed.top.ten$Score), decreasing=TRUE),]

saveRDS(figR.d.top.ten, paste0(sourcedir, 'top_ten_smoothed.rds'))
saveRDS(figR.d.unsmoothed.top.ten, paste0(sourcedir, 'top_ten_unsmoothed.rds'))

top.ten.smoothed <- readRDS(paste0(sourcedir, 'top_ten_smoothed.rds'))
top.ten.unsmoothed <- readRDS(paste0(sourcedir, 'top_ten_unsmoothed.rds'))

write.csv(top.ten.smoothed, paste0(sourcedir, 'top_ten_smoothed.csv'), row.names=FALSE, quote=FALSE)
write.csv(top.ten.unsmoothed, paste0(sourcedir, 'top_ten_unsmoothed.csv'), row.names=FALSE, quote=FALSE)

png(paste0(sourcedir, 'top_ten_smoothed.png'), width=)
grid.table(top.ten.smoothed, rows=NULL)
dev.off()
