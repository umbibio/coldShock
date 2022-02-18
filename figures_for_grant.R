library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(cowplot)
library(gam)
library(princurve)
library(parallel)
library(sme)
library(plotly)
library(patchwork)
library(ggVennDiagram)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(cowplot)



source('./util_funcs.R')


alldata.toxo.labs <- readRDS('../Input/coldShock/RData_new/S_O_list_taxo_labs_not_anchored.RData')
alldata.toxo.labs.integrated <- readRDS('../Input/coldShock/RData_new/S_O_list_taxo_labs_anchored.RData')


getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp, phase = S.O@meta.data$phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}


## Impact of cold on data: UMAP/PCA
pcaMataData.alldata <- getPcaMetaData(alldata.toxo.labs)


p1  <- ggplot(pcaMataData.alldata, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  #scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_human' = 'darkolivegreen4')) +
  #scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                             'BDiv_human' = 'darkolivegreen4')) +
  scale_color_manual(values = c("G" = "#ff6c67","SM" ='#6bb100', 'MC' = '#d473ff', 'C' = '#00c2c6')) +
  scale_fill_manual(values = c("G" = "#ff6c67","SM" ='#6bb100', 'MC' = '#d473ff', 'C' = '#00c2c6')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_grid(.~spp) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold.italic"
    )
  ) + 
  theme(legend.position="none") + 
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 0),
    axis.title.y = element_text(size=18, face="bold")
  )


plot(p1)


ggsave(filename="../Output/coldShock/figs/umap.pdf",
       plot=p1,
       width = 20, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## Proportions
tmp <- alldata.toxo.labs.integrated@meta.data

stats <- tmp %>% group_by(spp) %>% mutate(total.cells = n()) %>% 
  ungroup() %>% group_by(spp, phase) %>% summarise(counts = n(), perc = n()/total.cells[1]) 


p2 <- ggplot(data=stats, aes(x=spp, y=perc, color = phase, group = phase)) +
  geom_point(size = 1.5) + geom_path(size = 1.4) + geom_vline(xintercept = 5, col='pink', lwd=1.3, linetype=2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(face="bold", size=16, angle=45, hjust = 1.1, color = 'black')) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0, color = 'black')) +
  xlab('') + ylab('Proportion') + 
  scale_color_manual(values = c("G" = "#ff6c67","SM" ='#6bb100', 'MC' = '#d473ff', 'C' = '#00c2c6')) +
  scale_x_discrete(labels=c('BDiv0hrN' = '0hr' , 'BDiv4hrN' = '4hr',   'BDiv12hrN' = '12hr', 
                            'BDiv36hrN' = '36hr', 'BDiv7dN' = '7d',
                            'BDiv36hrY' = '36hr(R)',  'BDiv7dY' = '7d(R)')) + 
   theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + theme(legend.position = c(0.5, 0.6), legend.direction="horizontal",
             legend.title = element_text(colour="black", size=12, 
                                         face="bold"),
             legend.text = element_text(colour="black", size=12, 
                                        face="bold"))

plot(p2)

ggsave(filename="../Output/coldShock/figs/phase_proportions_trends.pdf",
       plot=p2,
       width = 4, height = 3,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


##

## Percent G vs. S/M. vs C
## 2D histogram of densities

#S.O.list <- SplitObject(alldata.toxo.labs.integrated, split.by = 'spp')
S.O.list <- SplitObject(alldata.toxo.labs, split.by = 'spp')
BDiv7dN <- getPcaMetaData(S.O.list$`7d`)

p  <- ggplot(BDiv7dN, aes(x= UMAP_1,y= UMAP_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  #geom_density_2d(color= 'red', size = 1) + 
  theme_bw(base_size = 14) +
  #theme(legend.position = "none") +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #ggtitle(titles[i]) +
  # theme(
  #   plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
  #   axis.title.x = element_text(size=18, face="bold", hjust = 1),
  #   axis.title.y = element_text(size=18, face="bold")
  # ) + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  =element_blank(),
    axis.ticks.x =element_blank(),
    axis.text.y  =element_blank(),
    axis.ticks.y =element_blank() 
  ) + 
  theme(legend.position = c(0.21, 0.77),
        legend.title = element_text(colour="black", size=16, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=16, 
                                   face="bold"))



plot(p)

ggsave(filename="../Output/coldShock/figs/density.pdf",
       plot=p,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



BDiv0hr <- getPcaMetaData(S.O.list$`0hr`)

p  <- ggplot(BDiv0hr, aes(x= UMAP_1,y= UMAP_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  #geom_density_2d(color= 'red', size = 1) + 
  theme_bw(base_size = 14) +
  #theme(legend.position = "none") +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #ggtitle(titles[i]) +
  # theme(
  #   plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
  #   axis.title.x = element_text(size=18, face="bold", hjust = 1),
  #   axis.title.y = element_text(size=18, face="bold")
  # ) + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  =element_blank(),
    axis.ticks.x =element_blank(),
    axis.text.y  =element_blank(),
    axis.ticks.y =element_blank() 
  ) + 
  theme(legend.position = 'none',
        #legend.position = c(0.2, 0.12), 
        legend.direction="horizontal",
        #legend.key.width=unit(0.4, "cm"),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold"))



plot(p)

ggsave(filename="../Output/coldShock/figs/density_0d.pdf",
       plot=p,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


BDiv7dR <- getPcaMetaData(S.O.list$`7d(R)`)

p  <- ggplot(BDiv7dR, aes(x= UMAP_1,y= UMAP_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  #geom_density_2d(color= 'red', size = 1) + 
  theme_bw(base_size = 14) +
  #theme(legend.position = "none") +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #ggtitle(titles[i]) +
  # theme(
  #   plot.title = element_text(size=18, face = "bold.italic", color = 'red'),
  #   axis.title.x = element_text(size=18, face="bold", hjust = 1),
  #   axis.title.y = element_text(size=18, face="bold")
  # ) + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  =element_blank(),
    axis.ticks.x =element_blank(),
    axis.text.y  =element_blank(),
    axis.ticks.y =element_blank() 
  ) + 
  theme(legend.position = 'none',
        #legend.position = c(0.2, 0.12), 
        legend.direction="horizontal",
        #legend.key.width=unit(0.4, "cm"),
        legend.title = element_text(colour="black", size=12, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=12, 
                                   face="bold"))



plot(p)

ggsave(filename="../Output/coldShock/figs/density_7dRd.pdf",
       plot=p,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


S.O.BDiv0hrN <- S.O.list$`0hr`


Idents(S.O.BDiv0hrN) <- 'phase'
p2 <- DimPlot(S.O.BDiv0hrN, reduction = "pca", 
              pt.size = 0.8,
              #shape.by='spp',
              label = T, repel = T, label.size = 5) + NoLegend() + 
  scale_color_manual(values = c("G1" = "#ff6c67","S/M" ='#a2a700', 'C' = '#00c377')) +
  theme(panel.spacing = unit(0.5, "lines")) + 
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(face="bold", size=16, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0)) +
  theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) 


plot(p2)

BDiv0hrN.stats <- stats %>% dplyr::filter(spp == '0hr')

p3 <- ggplot(data=BDiv0hrN.stats, aes(x=phase, y=perc, fill = phase)) +
  geom_bar(stat="identity", position=position_dodge(width = 1.0), width=0.8) +
  geom_text(aes(label=round(perc, 2)), vjust=2,  color="black", 
            size=6, fontface="bold", position=position_dodge(1.0), angle=0)+
  scale_fill_manual(values = c("G1" = "#ff6c67","S/M" ='#a2a700', 'C' = '#00c377')) +
  theme_minimal() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=16, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0)) +
  theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) #+  coord_flip()

plot(p3)



## Absolute numbers
## Calculating total mRNA levels
copy.numbers <- lapply(S.O.list, function(S.O){
  m <- mean(colSums(S.O@assays$RNA@data))
  s <- sd(colSums(S.O@assays$RNA@data))
  return(list(m = m, s = s))
})


cs <- lapply(1:length(S.O.list), function(i){
  data.frame(time = names(S.O.list)[i], copy.num = c(colSums(S.O.list[[i]]@assays$RNA@data)))
})

cs <- do.call('rbind', cs)

means <- lapply(copy.numbers, `[[`, 1)
sds <- lapply(copy.numbers, `[[`, 2)

copy.numbers <- data.frame(time = names(means), m = unlist(means), s= unlist(sds))
copy.numbers$time <- factor(copy.numbers$time,levels = c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'))



p <- ggplot(copy.numbers, aes(time , m, fill = time, color=time)) + 
  
  plot(p)

copy.numbers.filt <- copy.numbers %>% dplyr::filter(time %in% c('0hr', '7d', '7d(R)'))
p <- ggplot(data=copy.numbers.filt, aes(x = time , y = m)) +
  geom_bar(stat = "identity", fill="steelblue", width=0.8)+
  geom_errorbar(aes(ymin=m-s, ymax=m+s), width=0.1, color="black", size=0.6) + 
  geom_text(aes(label=round(m)), vjust = 6,  color="black", size=6,  fontface="bold")+
  xlab('') + ylab('mean copy number') + 
  theme_minimal() + 
  theme(panel.spacing = unit(0.2, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=18, angle=0, color = 'black')) +
  theme(axis.text.y = element_text(face="bold", size=18, angle=0, color = 'black')) +
  theme(
    axis.title.x = element_text(size=18, face="bold", color = 'black'),
    axis.title.y = element_text(size=18, face="bold", color = 'black')
  ) #+  coord_flip()

plot(p)

ggsave(filename="../Output/coldShock/figs/mean_copy_number.pdf",
       plot=p,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

cs.filt <- cs %>% dplyr::filter(time %in% c('0hr', '7d', '7d(R)'))
ggplot(cs.filt, aes(x=copy.num, fill=time)) +
  geom_density(bw = 200, alpha=.25)

p <- ggplot(cs.filt, aes(x=time, y=copy.num, fill=time)) +
  geom_violin(trim=FALSE) + stat_summary(fun.data=mean_sdl, 
                               geom="pointrange", color="black") + 
  xlab('') + ylab('copy number') + 
  theme_minimal() + 
  theme(panel.spacing = unit(0.2, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=18, angle=0, color = 'black')) +
  theme(axis.text.y = element_text(face="bold", size=18, angle=0, color = 'black')) +
  theme(
    axis.title.x = element_text(size=18, face="bold", color = 'black'),
    axis.title.y = element_text(size=18, face="bold", color = 'black')
  )  + 
  theme(legend.position="none",
    #legend.position = c(0.13, 0.83), 
    legend.title=element_text(size=18, face="bold"), 
    legend.text=element_text(size=18, face="bold"))





plot(p)

ggsave(filename="../Output/coldShock/figs/violins.pdf",
       plot=p,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


cs.filt2 <- cs %>% dplyr::filter(time %in% c('7d', '7d(R)'))
t.test(log(copy.num) ~ time, data = cs.filt2)

### DE analysis
global.shock.markers.sig <- read.xlsx('../Output/coldShock/tables/global_shock_markers_WT_base.xlsx')
matched.shock.markers.sig <- read.xlsx('../Output/coldShock/tables/matched_shock_markers.xlsx')

tmp <- global.shock.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% group_by(ref, query) %>%
  summarise(genes = list(GeneID))

venn.list <- tmp$genes


names(venn.list) <- tmp$query
names(venn.list) <- factor(names(venn.list), levels = c('4hr', '12hr', '36hr', '7d'))
g <- ggVennDiagram(venn.list, label_size = 6, set_size = 8,
                   set_color = c("4hr" = "firebrick","12hr" ="darkorchid3", '36hr' = 'darkslateblue', 
                                 '7d' = 'darkolivegreen4')) +  
  scale_color_manual(values = c("firebrick","darkorchid3", 'darkslateblue', 'darkolivegreen4')) +
  
  theme(legend.position = "none")


plot(g)


ggsave(filename="../Output/coldShock/figs/venn_degs.pdf",
       plot=g,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## barplots of global up & down genes
tmp <- global.shock.markers.sig 


tmp$direction <- ifelse(tmp$avg_log2FC > 0 , 'up', 'down')
tmp$contrast <-  tmp$query


stats.deg <- tmp %>% group_by(contrast, direction) %>% summarise(total = n()) %>% 
  group_by(contrast) %>% mutate(ylab_pos = total - 2)
stats.deg <- rbind(stats.deg, data.frame(contrast = '4hr', direction = 'up', total = 1, ylab_pos = 0))

stats.deg$contrast <- factor(stats.deg$contrast, levels = c('4hr', '12hr', '36hr', '7d'))
stats.deg$direction <- factor(stats.deg$direction, levels = c('up', 'down'))
stats.deg$label <- stats.deg$total
stats.deg$label[stats.deg$contrast == '4hr' &   stats.deg$direction == 'up'] <- 0

p <- ggplot(stats.deg, aes(fill=direction, y=total, x=contrast)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=label, group = direction), 
            position = position_dodge(width = 0.9), vjust=0, 
            color="black", size=7, fontface = 'bold') + 
  theme_minimal() + xlab('') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) +
  theme(#legend.position="top",
    legend.position = c(0.13, 0.83), 
    legend.title=element_text(size=18, face="bold"), 
    legend.text=element_text(size=18, face="bold"))



plot(p)


ggsave(filename="../Output/coldShock/figs/barplot_degs.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## barplots of matched up & down genes
tmp <- matched.shock.markers.sig 


tmp$direction <- ifelse(tmp$avg_log2FC > 0 , 'up', 'down')
tmp$contrast <-  tmp$query

lvs <- lapply(c('G1', 'SM', 'MC', 'C'), function(x) {
  paste(c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'), x, sep = '_')
})

lvs <- unlist(lvs)

stats.deg <- tmp %>% group_by(contrast, direction) %>% summarise(total = n()) %>% 
  group_by(contrast) %>% mutate(ylab_pos = total - 2)
stats.deg <- rbind(stats.deg, data.frame(contrast = '4hr_G1', direction = 'up', total = 1, ylab_pos = 0))
stats.deg <- rbind(stats.deg, data.frame(contrast = '4hr_G1', direction = 'down', total = 1, ylab_pos = 0))
stats.deg <- rbind(stats.deg, data.frame(contrast = '12hr_G1', direction = 'down', total = 1, ylab_pos = 0))

stats.deg$contrast <- factor(stats.deg$contrast, levels = lvs)
stats.deg$direction <- factor(stats.deg$direction, levels = c('up', 'down'))
stats.deg$label <- stats.deg$total
stats.deg$label[stats.deg$contrast == '4hr_G1'] <- 0
stats.deg$label[stats.deg$contrast == '12hr_G1' &   stats.deg$direction == 'down'] <- 0

p <- ggplot(stats.deg, aes(fill=direction, y=total, x=contrast)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=label, group = direction), 
            position = position_dodge(width = 0.9), vjust=0, 
            color="black", size=7, fontface = 'bold') + 
  theme_minimal() + xlab('') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) +
  theme(#legend.position="top",
    legend.position = c(0.13, 0.83), 
    legend.title=element_text(size=18, face="bold"), 
    legend.text=element_text(size=18, face="bold"))



plot(p)


ggsave(filename="../Output/coldShock/figs/barplot_matched_degs.pdf",
       plot=p,
       width = 12, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Expression distribution plots
tmp.down <- tmp %>% slice_min(order_by = avg_log2FC, n = 50) %>% group_by(GeneID) %>%
  mutate(mean_ave_log2FC = mean(avg_log2FC)) %>% dplyr::select(GeneID, gene, Product.Description) %>% 
  distinct() %>% data.frame()
tmp.up <- tmp %>% slice_max(order_by = avg_log2FC, n = 50) %>% group_by(GeneID) %>%
  mutate(mean_ave_log2FC = mean(avg_log2FC)) %>% dplyr::select(GeneID, gene, Product.Description) %>% 
  distinct() %>% data.frame()

### Genes of interest
#Bdiv_038820: E2F (EAPP) associated phosphoprotein 
#Bdiv_031680c: s-adenosylmethionine synthetase
#Bdiv_012690c: rRNA-processing CGR1
#Bdiv_023810: retinoblastoma A associated protein


## Plots of EAPP/SMAPS from Brendan synch data
library(sme)
bd.prod.desc <- read.csv('../Input/compScBdTgPb/genes/BDvi_Prod_desc.csv')
bd.prod.desc <- bd.prod.desc %>% transmute(GeneID = Gene.ID, ProductDescription = Product.Description)
sync.tc.fits <- readRDS('../Input/compScBabesia/RData/bd_sme_fits_sync_tc_20min.RData')
sync.tc.df <- readRDS('../Input/compScBabesia/RData/bd_sync_tc_df.RData')

EAPP.id <- "Bdiv_038820"
SAMPS.id <- "Bdiv_031680c"

pdf(file = "../Output/coldShock/figs/EAPP_SAMPS_expr_sync.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

par(mfrow = c(2,1), tcl=-0.5, family="serif", mai=c(0.6,0.6,0.6,0.6))
v <- EAPP.id
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('EAPP:', v), conf = T)
v <- SAMPS.id
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('SAMPS:', v), conf = T)

dev.off()


p1 <- VlnPlot(object = alldata.toxo.labs.integrated, features = 'Bdiv-038820') + 
  xlab('') + ylab('Expr') + ggtitle('EAPP: Bdiv_038820') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold") 
  ) +
  theme(legend.position="none",
        #legend.position = c(0.13, 0.83), 
        legend.title=element_text(size=16, face="bold"), 
        legend.text=element_text(size=16, face="bold"))
plot(p1)

ggsave(filename="../Output/coldShock/figs/violin_e2f.pdf",
       plot=p1,
       width = 8.5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

p2 <- VlnPlot(object = alldata.toxo.labs.integrated, features = 'Bdiv-031680c') + 
  xlab('') + ylab('Expr') + ggtitle('SAMS: Bdiv_031680c') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(legend.position="none",
        #legend.position = c(0.13, 0.83), 
        legend.title=element_text(size=16, face="bold"), 
        legend.text=element_text(size=16, face="bold"))
plot(p2)

ggsave(filename="../Output/coldShock/figs/violin_SAMS.pdf",
       plot=p2,
       width = 8.5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

p3 <- VlnPlot(object = alldata.toxo.labs.integrated, features = 'Bdiv-023810') + 
  xlab('') + ylab('Expr') + ggtitle('Rb: Bdiv_023810') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(legend.position="none",
        #legend.position = c(0.13, 0.83), 
        legend.title=element_text(size=16, face="bold"), 
        legend.text=element_text(size=16, face="bold"))
plot(p3)

ggsave(filename="../Output/coldShock/figs/violin_Rb.pdf",
       plot=p3,
       width = 8.5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



p4 <- VlnPlot(object = alldata.toxo.labs.integrated, features = 'Bdiv-011450c') + 
  xlab('') + ylab('Expr') + ggtitle('Histon2B: Bdiv_011450c') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold") 
  ) +
  theme(legend.position="none",
        #legend.position = c(0.13, 0.83), 
        legend.title=element_text(size=16, face="bold"), 
        legend.text=element_text(size=16, face="bold"))
plot(p4)


p5 <- VlnPlot(object = alldata.toxo.labs.integrated, features = 'Bdiv-031010c') + 
  xlab('') + ylab('Expr') + ggtitle('TPI: Bdiv-031010c') + 
  theme(
    axis.text = element_text(size=18, face="bold", color = 'black'),
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold") 
  ) +
  theme(legend.position="none",
        #legend.position = c(0.13, 0.83), 
        legend.title=element_text(size=16, face="bold"), 
        legend.text=element_text(size=16, face="bold"))
plot(p5)


ggsave(filename="../Output/coldShock/figs/violin_e2f.pdf",
       plot=p1,
       width = 8.5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

# p <- DotPlot(alldata.toxo.labs.integrated, features = unique(tmp.up$gene), 
#              cols = c("blue", "red", "green", "yellow"),
#              dot.scale = 8) + RotatedAxis()
# 
# plot(p)


Idents(alldata.toxo.labs.integrated) <- 'spp'
Idents(alldata.toxo.labs.integrated) <- 'phase.cond'
DefaultAssay(alldata.toxo.labs.integrated) <- 'RNA'

my.features <- c(tmp.up$gene[18], tmp.up$gene[27], 'Bdiv-031680c')

p1 <- VlnPlot(object = alldata.toxo.labs.integrated, features = my.features[1]) + 
  theme(
    axis.text = element_text(size=16, face="bold", color = 'black'),
    plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + xlab('') + ylab('Expr') + 
  theme(legend.position="none",
    #legend.position = c(0.13, 0.83), 
    legend.title=element_text(size=16, face="bold"), 
    legend.text=element_text(size=16, face="bold"))
plot(p1)
p2 <- VlnPlot(object = alldata.toxo.labs.integrated, features = my.features[2]) + 
  theme(
    axis.text = element_text(size=16, face="bold", color = 'black'),
    plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + xlab('') + ylab('') + 
  theme(legend.position="none",
    #legend.position = c(0.13, 0.83), 
    legend.title=element_text(size=16, face="bold"), 
    legend.text=element_text(size=16, face="bold"))
plot(p2)
p3 <- VlnPlot(object = alldata.toxo.labs.integrated, features = my.features[3]) + 
  theme(
    axis.text = element_text(size=16, face="bold", color = 'black'),
    plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + xlab('') + ylab('') + 
  theme(legend.position="none",
        #legend.position = c(0.13, 0.83), 
        legend.title=element_text(size=16, face="bold"), 
        legend.text=element_text(size=16, face="bold"))
plot(p3)

p <-grid.arrange(p1, p2, p3, ncol = 2)

plot_grid(p1, p2, p3)
ggsave(filename="../Output/coldShock/figs/violin_p3.pdf",
       plot=p3,
       width = 4, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## Expression:
DefaultAssay(alldata.toxo.labs) <- 'RNA'
Idents(alldata.toxo.labs) <- 'spp' 
p <- FeaturePlot(object = alldata.toxo.labs,
                  #shape.by = 'spp',
                  split.by = 'spp',
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = my.features[2],
                  cols = c("lightgrey", "red"), reduction = "pca")

plot(p)
