library(ggplot2)
library(cowplot)
library(patchwork)
library(tidytext)
library(tidyverse)
library(doParallel)
library(Seurat)

source('./util_funcs.R')
## individually processed samples, merged, but not anchored

prod.desc <- read.csv('../Input/coldShock/genes/BDiv_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

S.O.list <- readRDS('../Input/coldShock/RData_new/S_O_list.RData')
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

num.objs <- length(S.O.list)
spps <- names(S.O.list)
print(spps)


genes <- lapply(1:num.objs, function(i){
  gene.row <- rownames(S.O.list[[i]]@assays$RNA@data)
} )

common.genes <- Reduce(intersect, genes)
length(common.genes)

## Calculating total mRNA levels
copy.numbers <- lapply(S.O.list, function(S.O){
  m <- mean(colSums(as.matrix(S.O@assays$RNA@data)))
  s <- sd(colSums(as.matrix(S.O@assays$RNA@data)))
  return(list(m = m, s = s))
})

means <- lapply(copy.numbers, `[[`, 1)
sds <- lapply(copy.numbers, `[[`, 2)

copy.numbers <- data.frame(time = names(means), m = unlist(means), s= unlist(sds))
copy.numbers$time <- gsub('Y|N', "", gsub('BDiv', '', copy.numbers$time ))
copy.numbers$time <- factor(copy.numbers$time,levels = c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'))




p <- ggplot(data=copy.numbers, aes(x = time , y = m, fill = time)) +
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin=m-s, ymax=m+s), width=0.2, size=1) + 
  geom_text(aes(label=round(m), vjust=2,  color="black", size=6, fontface="bold"))+
  theme_minimal() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=16, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0)) +
  theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) #+  coord_flip()

plot(p)

## Transfer lables from toxo-based markers.
S.O.labs   <- readRDS('../Input/compScBabesia/RData_new/S_O_anchored_list_pstime_GAM_individual_phase_bound_from_toxo.Rdata')
S.O.labs   <- S.O.labs$bdiv_human
meta.data  <- S.O.labs@meta.data
trans.count.data <- S.O.labs@assays$RNA@counts
comm.id <- rownames(trans.count.data) %in% common.genes
trans.count.data <- trans.count.data[!is.na(comm.id), ]
S.O.labs <- CreateSeuratObject(counts = trans.count.data)
S.O.labs$orig.ident <- 'ToxoLabs'
S.O.labs$spp <- 'ToxoLabs'
S.O.labs <- AddMetaData(S.O.labs, meta.data)
S.O.labs <- prep_S.O(S.O.labs)

S.Os <- lapply(S.O.list, function(S.O){
  anchors <- FindTransferAnchors(reference = S.O.labs, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.labs@meta.data$cell.cycle.phase, dims = 1:30)
  predictions$phase <- predictions$predicted.id
  S.O <- AddMetaData(object = S.O, metadata = predictions)
})

# ## Transfer Toxo lables
# S.O.tg    <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')
# GT1.bdiv  <- read.xlsx("../Input/compScBdTgPb/Orthologs/rec_GT1.vs.B.divergence.xlsx")
# GT1.BD <- GT1.bdiv %>% dplyr::select('query_id', contains('subject_id'))
# colnames(GT1.BD) <- c('TGGT1', 'BDiv')
# tg_meta.data <- S.O.tg@meta.data
# tg.data <- S.O.tg@assays$RNA@counts
# bdiv.id <- GT1.BD$BDiv[match(gsub('-','_',rownames(tg.data)), GT1.BD$TGGT1)]
# orth.ind <- !is.na(bdiv.id)
# tg.data <- tg.data[orth.ind, ]
# sum(orth.ind) ## total number of orthologs
# rownames(tg.data) <- bdiv.id[orth.ind]

# S.O.toxo <- CreateSeuratObject(counts = tg.data)
# S.O.toxo$orig.ident <- 'Toxo'
# S.O.toxo$spp <- 'Toxo'
# S.O.toxo <- AddMetaData(S.O.toxo, tg_meta.data)
# 
# S.O.toxo <- prep_S.O(S.O.toxo)

## Transfer labels
# S.Os <- lapply(S.O.list, function(S.O){
#   anchors <- FindTransferAnchors(reference = S.O.toxo, query = S.O, dims = 1:30)
#   predictions <- TransferData(anchorset = anchors, refdata = S.O.toxo@meta.data$phase, dims = 1:30)
#   predictions$phase <- predictions$predicted.id
#   S.O <- AddMetaData(object = S.O, metadata = predictions)
# })

## Merge again without anchoring
alldata.toxo.labs <- merge(S.Os[[1]], S.Os[2:num.objs], add.cell.ids=spps)
alldata.toxo.labs <- NormalizeData(alldata.toxo.labs, normalization.method = "LogNormalize", scale.factor = 10000)
alldata.toxo.labs <- FindVariableFeatures(alldata.toxo.labs, selection.method = "vst", nfeatures = 3000)
alldata.toxo.labs <- ScaleData(alldata.toxo.labs)
alldata.toxo.labs <- RunPCA(alldata.toxo.labs, features = VariableFeatures(object = alldata.toxo.labs))
alldata.toxo.labs <- FindNeighbors(alldata.toxo.labs, dims = 1:30, reduction = 'pca')
alldata.toxo.labs <- FindClusters(alldata.toxo.labs, resolution = 0.2)
alldata.toxo.labs <- RunUMAP(alldata.toxo.labs, dims = 1:30)

## Change the phase label 
# alldata.toxo.labs@meta.data$phase[alldata.toxo.labs@meta.data$phase == 'G1.a'] <- 'G1'
# alldata.toxo.labs@meta.data$phase[alldata.toxo.labs@meta.data$phase == 'G1.b'] <- 'G1'
# alldata.toxo.labs@meta.data$phase[alldata.toxo.labs@meta.data$phase == 'S'] <- 'S/M'
# alldata.toxo.labs@meta.data$phase[alldata.toxo.labs@meta.data$phase == 'M'] <- 'S/M'

# alldata.toxo.labs@meta.data$phase <- factor(alldata.toxo.labs@meta.data$phase, 
#                                             levels = c('G1', 'S/M', 'C'))

alldata.toxo.labs@meta.data$phase <- factor(alldata.toxo.labs@meta.data$phase, 
                                            levels = c('G', 'SM', 'MC', 'C'))

Idents(alldata.toxo.labs) <- "phase"

## Relabel spp
alldata.toxo.labs@meta.data$spp <- gsub('Y', '\\(R\\)', gsub('N', '',gsub('BDiv', '', alldata.toxo.labs@meta.data$spp)))

alldata.toxo.labs@meta.data$spp <- factor(alldata.toxo.labs@meta.data$spp,
                                                     levels = c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'))

p <- DimPlot(alldata.toxo.labs, reduction = "pca", 
             #group.by = "cell", 
             split.by = 'spp',
             pt.size = 0.8,
             #shape.by='spp',
             label = T, repel = T, label.size = 5) + NoLegend() + 
  # scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                               'BDiv_human' = 'darkolivegreen4')) +
  # scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_human' = 'darkolivegreen4')) +
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )




plot(p)


## Save the data
saveRDS(alldata.toxo.labs, '../Input/coldShock/RData_new/S_O_list_taxo_labs_not_anchored.RData')

## Now Integrate the objects, using 0hr as reference
print(file.info)
names(S.Os)
ref.ind <- 1

alldata.toxo.labs.integrated <- processeMergedS.O(S.Os, data.ind = NA, ref.ind = ref.ind, res = 0.1, SC = FALSE)


## Change phase labels
# alldata.toxo.labs.integrated@meta.data$phase[alldata.toxo.labs.integrated@meta.data$phase == 'G1.a'] <- 'G1'
# alldata.toxo.labs.integrated@meta.data$phase[alldata.toxo.labs.integrated@meta.data$phase == 'G1.b'] <- 'G1'
# alldata.toxo.labs.integrated@meta.data$phase[alldata.toxo.labs.integrated@meta.data$phase == 'S'] <- 'S/M'
# alldata.toxo.labs.integrated@meta.data$phase[alldata.toxo.labs.integrated@meta.data$phase == 'M'] <- 'S/M'

# alldata.toxo.labs.integrated@meta.data$phase <- factor(alldata.toxo.labs.integrated@meta.data$phase, 
#                                                        levels = c('G1', 'S/M', 'C'))

alldata.toxo.labs.integrated@meta.data$phase <- factor(alldata.toxo.labs.integrated@meta.data$phase, 
                                                        levels = c('G', 'SM', 'MC', 'C'))

Idents(alldata.toxo.labs.integrated) <- "phase"

## Relabel spp
alldata.toxo.labs.integrated@meta.data$spp <- gsub('Y', '\\(R\\)', 
                                                   gsub('N', '',gsub('BDiv', '', alldata.toxo.labs.integrated@meta.data$spp)))

alldata.toxo.labs.integrated@meta.data$spp <- factor(alldata.toxo.labs.integrated@meta.data$spp,
                                                     levels = c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'))
alldata.toxo.labs.integrated@meta.data$spp <- gsub('Y', '\\(R\\)', 
                                                   gsub('N', '', gsub('BDiv', '', alldata.toxo.labs.integrated@meta.data$spp)))

alldata.toxo.labs.integrated@meta.data$spp <- factor(alldata.toxo.labs.integrated@meta.data$spp,
                                                     levels = c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'))


## Generate matched phase/spp variable
alldata.toxo.labs.integrated@meta.data$phase.cond <- paste(alldata.toxo.labs.integrated@meta.data$spp,
                                                           alldata.toxo.labs.integrated@meta.data$phase, sep = '_')

# lvs <- lapply(c('G1', 'S/M', 'C'), function(x) {
#   paste(c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'), x, sep = '_')
# })

lvs <- lapply(c('G', 'SM', 'MC', 'C'), function(x) {
  paste(c('0hr', '4hr', '12hr', '36hr', '7d', '36hr(R)', '7d(R)'), x, sep = '_')
})

lvs <- unlist(lvs)
alldata.toxo.labs.integrated@meta.data$phase.cond <- factor(alldata.toxo.labs.integrated@meta.data$phase.cond,
                                                            levels = lvs)

Idents(alldata.toxo.labs.integrated) <- "phase"
p <- DimPlot(alldata.toxo.labs.integrated, reduction = "pca", 
              #group.by = "cells", 
              split.by = 'spp',
              pt.size = 1,
              #shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

# ggsave(filename="../Output/coldShock/cold_shock_anchored_toxo_phase_transfer.pdf",
#        plot=p,
#        width = 12, height = 8,
#        units = "in", # other options are "in", "cm", "mm"
#        dpi = 300
# )


## Save the data
saveRDS(alldata.toxo.labs.integrated, '../Input/coldShock/RData_new/S_O_list_taxo_labs_anchored.RData')


## Some plots
tmp <- alldata.toxo.labs.integrated@meta.data

stats <- tmp %>% group_by(spp) %>% mutate(total.cells = n()) %>% 
  ungroup() %>% group_by(spp, phase) %>% summarise(counts = n(), perc = n()/total.cells[1]) 

## manually add the missing phase
#stats <- rbind(stats, data.frame(spp = 'BDiv7dN', phase = 'M', counts = 0, perc = 0))

stats <- stats %>% arrange(spp, phase)

p <- ggplot(data=stats, aes(x=spp, y=perc, fill = phase)) +
  geom_bar(stat="identity", position=position_dodge(width = 1.0), width=0.8) +
  geom_text(aes(label=round(perc, 2)), vjust=-1.3,  color="black", 
            size=3.5, fontface="bold", position=position_dodge(1.0), angle=0)+
  theme_minimal() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) #+  coord_flip()

plot(p)

ggsave(filename="../Output/coldShock/phase_proportions_bar_plots.pdf",
       plot=p,
       width = 12, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

p <- ggplot(data=stats, aes(x=spp, y=perc, color = phase, group = phase)) +
  geom_point() + geom_path() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) #+  coord_flip()

plot(p)

ggsave(filename="../Output/coldShock/phase_proportions_trends.pdf",
       plot=p,
       width = 12, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## 2D histogram of densities
lvs <- c("BDiv0hrN","BDiv4hrN", "BDiv12hrN", "BDiv36hrN", "BDiv7dN", "BDiv36hrY", "BDiv7dY" )
pcaMataData.alldata$spp <- factor(pcaMataData.alldata$spp, levels = lvs)
p  <- ggplot(pcaMataData.alldata, aes(x= PC_1,y=PC_2)) +
  #geom_density_2d() +
  
  #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  #geom_density_2d() + 
  
  #geom_bin2d(bins = 60) +
  #scale_fill_continuous(type = "viridis") +
  
  # geom_point(aes(#fill = lable.prob,
  #   fill = spp,
  #   color = spp
  # ), #color = 'blue', 
  # shape=21, size = 1)+ 
  #scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_human' = 'darkolivegreen4')) +
  #scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                             'BDiv_human' = 'darkolivegreen4')) +
  
  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('UMAP2') + xlab('UMAP1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  facet_wrap(spp~.) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)

BDiv0hrN <- pcaMataData.alldata %>% dplyr::filter(spp == 'BDiv0hrN')

p  <- ggplot(BDiv0hrN, aes(x= PC_1,y= PC_2)) +
  #geom_density_2d() +
  
  #stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_density_2d(color= 'red', size = 1) + 
  
  #geom_bin2d(bins = 60) +
  #scale_fill_continuous(type = "viridis") +
  
  # geom_point(aes(#fill = lable.prob,
  #   fill = spp,
  #   color = spp
  # ), #color = 'blue', 
  # shape=21, size = 1)+ 
  #scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
#                              'BDiv_human' = 'darkolivegreen4')) +
#scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
#                             'BDiv_human' = 'darkolivegreen4')) +

  theme_bw(base_size = 14) +
  #theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  #scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)

BDiv0hrN.stats <- stats %>% dplyr::filter(spp == 'BDiv0hrN')

p <- ggplot(data=BDiv0hrN.stats, aes(x=phase, y=perc, fill = phase)) +
  geom_bar(stat="identity", position=position_dodge(width = 1.0), width=0.8) +
  geom_text(aes(label=round(perc, 2)), vjust=-1.3,  color="black", 
            size=3.5, fontface="bold", position=position_dodge(1.0), angle=0)+
  theme_minimal() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) #+  coord_flip()

plot(p)


## Individual plots
DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- "spp"
Idents(all.samples.integrated) <- "phase.cond"

my.gene <- "Bdiv-038820"

p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = my.gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = my.gene, 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

p4 <- VlnPlot(object = all.samples.integrated, features = my.gene)
plot(p4)



# Differential Expression
# Identify global markers independent of cell cycle phase. Are there any global Cold Shock regulators?
# Two class comparison: fc = ident.1/ident.2

DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- "spp"
objs <- unique(all.samples.integrated@meta.data$spp)
ref.obj <- objs[1]
quesy.obj <- objs[2:length(objs)]

global.shock.markers <- FindMarkers(all.samples.integrated, ident.1 = quesy.obj,  ident.2  = ref.obj, verbose = T)
global.shock.markers$genes <- rownames(global.shock.markers)
global.shock.markers$GeneID <- gsub('-', '_', global.shock.markers$genes)
global.shock.markers.sig <- global.shock.markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5))
## Percentages are calclulated on Variable features, not all genes
global.stats <- global.shock.markers.sig %>% 
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(up.reg) %>% summarise(num.DEGs = n(), percent = n() / nrow(global.shock.markers)) 

print(global.stats)

## Top up-regulated markers
global.shock.markers.top <- global.shock.markers.sig %>% dplyr::filter(avg_log2FC > 0) %>% top_n(4, abs(avg_log2FC))
global.shock.markers.top


p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = global.shock.markers.top$genes,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = global.shock.markers.top$genes, 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

VlnPlot(object = all.samples.integrated, features = global.shock.markers.top$genes)

tmp <- left_join(global.shock.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/global_shock_markers.xlsx')


### Differential expression analysis
## Cell cycle phase specific


all.spp.list <- SplitObject(all.samples.integrated, split.by = "spp")

## Re-normalize splitted data
for (i in 1:length(all.spp.list)) {
  Idents(all.spp.list[[i]]) <- 'phase.cond'
  all.spp.list[[i]] <- NormalizeData(all.spp.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}


GM <- getCellCyclePhaseMarkers(all.spp.list)


## For each cluster, pool global markers.
pooled_markers <- bind_rows(GM$all.markers.list.sig) %>% group_by(glob.clust) %>%
  summarise(glob.markers = list(unique(gene)), num.markers = length(unique(gene)))

objs <- unique(all.samples.integrated@meta.data$spp)
ref.obj <- objs[1]
quesy.obj <- objs[2:length(objs)]
clusters <- unique(all.samples.integrated@meta.data$seurat_clusters)
ref.clusts <- data.frame(ref = paste(ref.obj, clusters, sep = '_'))
ref.clusts$cluster <- gsub('.*_', '', ref.clusts$ref)
ref.clusts$dummy <- 1
query.clusts <- data.frame(query = paste(rep(quesy.obj, each = length(clusters)), clusters, sep = '_'))
query.clusts$dummy <- 1
query.clusts$cluster <- gsub('.*_', '', query.clusts$query)

contrasts <- full_join(ref.clusts, query.clusts, by = 'dummy') %>% 
  dplyr::filter(cluster.x == cluster.y) %>%
  transmute(ref = ref, query = query)

contrasts$cluster <- gsub('.*_', '', contrasts$ref)

#ident.1 case, ident.2 is control
Idents(all.samples.integrated) <- 'phase.cond'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(all.samples.integrated, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  ind <- rownames(tmp) %in% unlist(pooled_markers$glob.markers[which(pooled_markers$glob.clust == x$cluster)]) 
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$cluster <- x$cluster
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp[ind, ])
  #return(tmp)
})


matched.DEGs <- bind_rows(matched.DEGs)

matched.DEGs.sig <- matched.DEGs %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5)) %>% arrange(desc(abs(avg_log2FC)))

matched.DEGs.stats <- matched.DEGs.sig %>%
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(query, up.reg) %>% summarise(num.DEGs = n()) 

print(matched.DEGs.stats)

matched.DEGs.stats$cluster <- as.factor(as.numeric(gsub('.*_', '', matched.DEGs.stats$query)))


matched.DEGs.top <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% group_by(cluster) %>% top_n(1, abs(avg_log2FC))



p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = matched.DEGs.top$gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



VlnPlot(object = all.samples.integrated, features = matched.DEGs.top$gene)


tmp <- left_join(matched.DEGs.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/matched_shock_markers.xlsx')


####################

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- 'seurat_clusters'


z.cells <- subset(all.samples.integrated, idents = "2")
Idents(z.cells) <- "spp"
avg.z.cells <- as.data.frame(log1p(AverageExpression(z.cells, verbose = FALSE)$RNA))
avg.z.cells$gene <- rownames(avg.z.cells)


p1 <- ggplot(avg.z.cells, aes(BDiv0hrN, BDiv7dN)) + geom_point() + ggtitle("0hr vs 12 hr")
plot(p1)


Idents(all.samples.integrated) <- 'phase.cond'

cold.response12h <- FindMarkers(all.samples.integrated, ident.1 = "BDiv12hrN_0", ident.2 = "BDiv0hrN_0", verbose = FALSE)

cold.response12h.sig <- cold.response12h %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5)) %>% 
  arrange(desc(abs(avg_log2FC)))

matched.DEGs.stats <- matched.DEGs.sig %>%
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(query, up.reg) %>% summarise(num.DEGs = n()) 

print(matched.DEGs.stats)

matched.DEGs.stats$cluster <- as.factor(as.numeric(gsub('.*_', '', matched.DEGs.stats$query)))


matched.DEGs.top <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% top_n(4, abs(avg_log2FC))




######## Plot individual
my.gene <- "Bdiv-011450c"
FeaturePlot(object = all.samples.integrated, 
            shape.by = 'spp',
            split.by = 'spp',
            label = T, pt.size = 0.6, label.size = 3, 
            features = my.gene,
            cols = c("lightgrey", "red"), reduction = "pca") 


VlnPlot(object = all.samples.integrated, features = my.gene)
#########


DefaultAssay(all.samples.integrated) <- 'RNA'
## Find cell cycle phase markers independently and cross compare
all.spp.list <- SplitObject(all.samples.integrated, split.by = "spp")

## Re-normalize splitted data
# for (i in 1:length(all.spp.list)) {
#   DefaultAssay(all.spp.list[[i]]) <- "RNA"
#   all.spp.list[[i]] <- NormalizeData(all.spp.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
#   DefaultAssay(all.spp.list[[i]]) <- "integrated"
# }


GM <- getCellCyclePhaseMarkers(all.spp.list)


## For each cluster, pool global markers.
pooled_markers <- bind_rows(GM$all.markers.list.sig) %>% group_by(glob.clust) %>%
  summarise(glob.markers = list(unique(gene)), num.markers = length(unique(gene)))

# Identify cell cycle phase specific markers. 
# Two class comparison: fc = ident.2/ident.1

objs <- names(all.spp.list)
ref.obj <- objs[1]
quesy.obj <- objs[2:length(objs)]
clusters <- unique(all.samples.integrated@meta.data$seurat_clusters)
ref.clusts <- data.frame(ref = paste(ref.obj, clusters, sep = '_'))
ref.clusts$cluster <- gsub('.*_', '', ref.clusts$ref)
ref.clusts$dummy <- 1
query.clusts <- data.frame(query = paste(rep(quesy.obj, each = length(clusters)), clusters, sep = '_'))
query.clusts$dummy <- 1
query.clusts$cluster <- gsub('.*_', '', query.clusts$query)

contrasts <- full_join(ref.clusts, query.clusts, by = 'dummy') %>% 
  dplyr::filter(cluster.x == cluster.y) %>%
  transmute(ref = ref, query = query)

contrasts$cluster <- gsub('.*_', '', contrasts$ref)

#Idents(all.samples.integrated) <- 'phase.cond'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(all.samples.integrated, ident.1 = x$ref, ident.2 = x$query, verbose = T)
  ind <- rownames(tmp) %in% unlist(pooled_markers$glob.markers[which(pooled_markers$glob.clust == x$cluster)]) 
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$cluster <- x$cluster
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp[ind, ])
})


matched.DEGs <- bind_rows(matched.DEGs)

matched.DEGs.sig <- matched.DEGs %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5)) %>% arrange(desc(abs(avg_log2FC)))

matched.DEGs.stats <- matched.DEGs.sig %>%
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(query, up.reg) %>% summarise(num.DEGs = n()) 

print(matched.DEGs.stats)

matched.DEGs.stats$cluster <- as.factor(as.numeric(gsub('.*_', '', matched.DEGs.stats$query)))


matched.DEGs.top <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% top_n(4, abs(avg_log2FC))



p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = matched.DEGs.top$gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = unique(matched.DEGs.top$gene), 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

VlnPlot(object = all.samples.integrated, features = matched.DEGs.top$gene)



####
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- "seurat_clusters"
nk.markers <- FindConservedMarkers(all.samples.integrated, ident.1 = '0', grouping.var = "spp", verbose = FALSE)
nk.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05) %>%  top_n(4, abs(avg_log2FC))

Idents(all.samples.integrated) <- "phase.cond"
Markers <- FindMarkers(all.samples.integrated, ident.1 = "BDiv7dN_2", ident.2 = "BDiv0hrN_2", verbose = FALSE)
Markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05) %>% arrange(pct.1) %>% top_n(4, abs(avg_log2FC))
####



### Can we filter cells by expression?
my.gene <- "Bdiv_011440c"
my.gene <- "Bdiv_016320"
exprs <- getNormExpr(all.spp.list[[2]])
exprs.g <- exprs %>% dplyr::filter(GeneID == my.gene)
exprs.g <- left_join(exprs.g, all.spp.list[[2]]@meta.data, by = c('Sample' = 'NAME'))
p <- ggplot(exprs.g, aes(x=seurat_clusters, y=expr)) + 
  geom_violin(trim=FALSE)
#p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# violin plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))

hist(exprs.g$expr, nclass = 50)
DefaultAssay(all.spp.list[[2]]) <- 'RNA'
VlnPlot(object = all.spp.list[[2]], features = "Bdiv-006560")

DefaultAssay(all.samples.integrated) <- 'integrated'
VlnPlot(object = all.samples.integrated, features = "Bdiv-006560")


FeaturePlot(object = all.samples.integrated, 
            shape.by = 'spp',
            split.by = 'spp',
            label = T, pt.size = 0.6, label.size = 3, 
            features = "Bdiv-006560",
            cols = c("lightgrey", "red"), reduction = "pca") 


xx <- FetchData(all.spp.list[[2]], vars = "Bdiv-016320")
hist(xx$`Bdiv-016320`)
hist(all.spp.list[[2]][['RNA']]@data["Bdiv-016320",])



## Generate a heatmap of expression of pooled markers.
exprs <- getNormExpr(all.samples.integrated)
exprs.scale <- exprs %>% dplyr::select(GeneID, Sample, expr) %>% 
  pivot_wider(names_from = GeneID, values_from = expr)
exprs.scale[,2:ncol(exprs.scale)] <- scale(exprs.scale[,2:ncol(exprs.scale)], center = T, scale = T)
exprs.scale <- exprs.scale %>% pivot_longer(-Sample, names_to = 'GeneID', values_to = 'expr')
exprs.filt <- exprs.scale %>% dplyr::filter(GeneID %in% matched.DEGs.sig$GeneID)
exprs.filt <- left_join(exprs.filt, all.samples.integrated@meta.data, by = c('Sample' = 'NAME'))
exprs.filt.groups <- exprs.filt %>% group_by(phase.cond, GeneID) %>% 
  summarise(mean_expr = mean(expr), median_expr = median(expr), sd_expr = sd(expr))
tmp <- all.samples.integrated@meta.data %>% 
  dplyr::select(c("spp", "time", "reactivate",  "phase.cond", "celltype")) %>% distinct()

exprs.filt.groups <- left_join(exprs.filt.groups, tmp, by = c('phase.cond'))


tmp <- matched.DEGs.sig %>% dplyr::select(avg_log2FC, GeneID, gene) %>% 
  mutate(reg = ifelse(avg_log2FC > 0, 1, -1)) %>% dplyr::select(GeneID, gene, reg) %>% distinct()
exprs.filt.groups <- left_join(exprs.filt.groups, tmp, by = c('GeneID'))

exprs.filt.groups <- dplyr::rename(exprs.filt.groups, cluster = celltype)
exprs.filt.groups$time <- factor(exprs.filt.groups$time,
                                 levels = c("0hr", "4hr", "12hr", "36hr", "7d"))

exprs.filt.groups.up <- exprs.filt.groups %>% dplyr::filter(reg == 1) %>% ungroup()
exprs.filt.groups.down <- exprs.filt.groups %>% dplyr::filter(reg == -1) %>% ungroup()


reorderGenes <- function(exprs){
  clust.ord <- exprs %>% dplyr::select(GeneID, cluster) %>% distinct() %>% 
    group_by(cluster) %>% summarise(GeneSet = list(GeneID)) 
  
  clusters <- unique(exprs$cluster)
  hc_eucledian.df <- {}
  for(i in 1:length(clusters)){
    class.i <- exprs %>% dplyr::filter(cluster == clusters[i]) %>% dplyr::select(GeneID, mean_expr, time) %>%
      pivot_wider(names_from = GeneID, values_from = mean_expr) %>% as.data.frame()
    
    
    ## Hierarchical clustering with Eucledian distance
    hc_eucledian <- hclust(dist(t(as.matrix(class.i[,-1] ))), method = "ward.D")
    hc_eucledian.df <- rbind(hc_eucledian.df, 
                             data.frame(GeneID = colnames(class.i[,-1]), 
                                        hc_eucledian.order = hc_eucledian$order,
                                        hc_eucledian.cluster = cutree(hc_eucledian,k = 10),
                                        cluster = clusters[i]))
  }
  
  return(hc_eucledian.df)
}


hc_eucledian_down.df <- reorderGenes(exprs.filt.groups.down)

exprs.filt.groups.down <- left_join(exprs.filt.groups.down, hc_eucledian_up.df, by = c('GeneID', 'cluster'))

p1 <- ggplot(exprs.filt.groups.down, aes(x = time, y = reorder_within(GeneID, hc_eucledian.order, cluster), fill = mean_expr)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    #axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(cluster~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p1)

matched.DEGs.top <- matched.DEGs.sig %>% group_by(cluster) %>% dplyr::filter(avg_log2FC < 0) %>% 
  top_n(2, abs(avg_log2FC))
  #top_n(2, desc(pct.1))
print(matched.DEGs.top)



p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = matched.DEGs.top$gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = matched.DEGs.top$gene, 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

VlnPlot(object = all.samples.integrated, features = matched.DEGs.top$gene)

tmp <- left_join(matched.DEGs.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/matched_shock_markers3.xlsx')


FeaturePlot(object = all.samples.integrated, 
            shape.by = 'spp',
            split.by = 'spp',
            label = T, pt.size = 0.6, label.size = 3, 
            features = c('Bdiv-005010c'),
            cols = c("lightgrey", "red"), reduction = "pca") 





#######
######
######



cluster.sim.marker1.marker2 <- crossCompareMarkers(GM$all.markers.list.sig[[1]], GM$all.markers.list.sig[[2]], 'WT', '7dN')
cluster.sim.marker1.marker1 <- crossCompareMarkers(GM$all.markers.list.sig[[1]], GM$all.markers.list.sig[[1]], 'WT1', 'WT2')
cluster.sim.marker2.marker2 <- crossCompareMarkers(GM$all.markers.list.sig[[2]], GM$all.markers.list.sig[[2]], '7dN1', '7dN2')

cluster.cont.marker1.marker2 <- markerContrasts(GM$all.markers.list.sig[[1]], GM$all.markers.list.sig[[2]], 'WT', '7dN')

write.xlsx(cluster.cont.marker1.marker2, '../Output/scClockOut/cluster_cont_WT_7dN.xlsx')  
head(unlist(cluster.cont.marker1.marker2$`7dN.only.genes`)) 

# Identify global markers independent of cell cycle phase. Are there any global Cold Shock regulators?
# Two class comparison: fc = ident.2/ident.1

Idents(all.samples.integrated) <- "spp"
global.shock.markers <- FindMarkers(all.samples.integrated, ident.1 = "BDiv0hrN",  ident.2  = "BDiv4hrN", verbose = T)
global.shock.markers$genes <- rownames(global.shock.markers)
global.shock.markers$GeneID <- gsub('-', '_', global.shock.markers$genes)
global.shock.markers.sig <- global.shock.markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5))
## Percentages are calclulated on Variable features, not all genes
global.stats <- global.shock.markers.sig %>% 
  mutate(up.reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(up.reg) %>% summarise(num.DEGs = n(), percent = n() / nrow(global.shock.markers)) 

print(global.stats)

## Top up-regulated markers
global.shock.markers.top <- global.shock.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% top_n(4, abs(avg_log2FC))
global.shock.markers.top


p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = global.shock.markers.top$genes,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = global.shock.markers.top$genes, 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

VlnPlot(object = all.samples.integrated, features = global.shock.markers.top$genes)

tmp <- left_join(global.shock.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/global_shock_markers.xlsx')



p1 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = c('Bdiv-010350c'),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p1)

RidgePlot(pbmc, features = , ncol = 2)

FeaturePlot(all.samples.integrated, features = xx$gene,  max.cutoff = 3, shape.by = 'spp', split.by = 'spp',
            cols = c("grey", "red"), reduction = 'pca')



markers.to.plot <- gsub('_', '-', head(unlist(cluster.cont.marker1.marker2$`7dN.only.genes`))) 
markers.to.plot <- GM$all.markers.list.sig[[2]] %>% 
  dplyr::filter(GeneID %in% unlist(cluster.cont.marker1.marker2$`7dN.only.genes`)) %>% 
  group_by(cluster) %>% top_n(2, avg_log2FC)

DotPlot(all.samples.integrated, features = markers.to.plot$gene, 
        cols = c("blue", "red", "green", "yellow"),
        dot.scale = 8) + 
  RotatedAxis()

FeaturePlot(all.spp.list[[2]], features = 'Bdiv-024520',  max.cutoff = 3, cols = c("grey", "red"), reduction = 'pca')

plots <- VlnPlot(all.samples.integrated, features = markers.to.plot$gene, split.by = "spp", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

VlnPlot(object = all.spp.list[[2]], features = markers.to.plot$gene)

markers.to.plot1 <- GM$all.markers.list.sig[[1]] %>% 
  group_by(cluster) %>% top_n(2, avg_log2FC)

VlnPlot(object = all.spp.list[[1]], features = markers.to.plot1$gene)


markers.to.plot2 <- GM$all.markers.list.sig[[2]] %>% 
  group_by(cluster) %>% top_n(2, avg_log2FC)

VlnPlot(object = all.spp.list[[2]], features = markers.to.plot2$gene)

VlnPlot(object = all.spp.list[[1]], features = markers.to.plot2$gene)

p1 <- FeaturePlot(all.samples.integrated, features = GM$all.markers.list.top[[1]]$gene,  split.by = 'spp',
            max.cutoff = 3, cols = c("grey", "red"), reduction = 'pca')

ggsave(filename="../Output/scClockFigs/Integrated_0hr_7D_0h_marker_expr.pdf", 
       plot=p1,
       width =8, height = 20, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



p2 <- FeaturePlot(all.samples.integrated, features = GM$all.markers.list.top[[2]]$gene,  split.by = 'spp',
                  max.cutoff = 3, cols = c("grey", "red"), reduction = 'pca')

ggsave(filename="../Output/scClockFigs/Integrated_0hr_7D_7h_marker_expr.pdf", 
       plot=p2,
       width =8, height = 20, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p3 <- FeaturePlot(all.samples.integrated, features = markers.to.plot$gene,  split.by = 'spp',
                  max.cutoff = 3, cols = c("grey", "red"), reduction = 'pca')


ggsave(filename="../Output/scClockFigs/Integrated_0hr_7D_top_non_overlapping_marker_expr.pdf", 
       plot=p3,
       width =8, height = 20, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


plot(p3)
library(plotly)
expr.df <- getNormExpr(all.spp.list[[2]])

gene.id <- 'Bdiv_024520'
gene.heat <- expr.df  %>% dplyr::filter(GeneID == gene.id) %>%
  transmute(log2.expr = expr, GeneID = GeneID, Sample = Sample)
pc <- getPCA(all.spp.list[[2]])
pc.heat <- left_join(pc, gene.heat, by = 'Sample')
fig <- plot_ly(pc.heat, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~log2.expr, colors = colorRamp(c("lightgray", "red")))
#colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers(size=2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig







library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
#library(sctransform)


source('./util_funcs.R')

### b. divergense
input.dir.bdiv <- "../Input/scRNAseqBdiv/"
bdiv.count.file <- "bdiv.expr.csv"


## Reading b. divergense data
bdiv.count <- read.csv(paste(input.dir.bdiv, bdiv.count.file, sep = ''))


genes <- bdiv.count$X
bd.expr <- bdiv.count[,-1]
rownames(bd.expr) <- genes


bdiv.pheno <- data.frame(X = colnames(bd.expr))
bdiv.pheno$spp <- 'BDiv' 


# Set initial Seurat clusters
S.O.bd <- CreateSeuratObject(counts = bd.expr, min.cells = 10, min.features = 100)
#S.O.bd <- CreateSeuratObject(counts = bd.expr)

#VlnPlot(S.O.bd, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.bd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.bd <- subset(S.O.bd, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
S.O.bd <- prep_S.O(S.O.bd)

bd.clust.info <- data.frame(X=as.character(names(S.O.bd$orig.ident)),
                            cluster=as.character(S.O.bd$seurat_clusters))
bdiv.pheno <- inner_join(bdiv.pheno, bd.clust.info, by = 'X')
bdiv.pheno$cells <- paste('BDiv', bdiv.pheno$cluster, sep = '')

bdiv.pheno$NAME <- paste(bdiv.pheno$cells, 
                         bdiv.pheno$X, sep = '_')

colnames(bdiv.pheno) <- c('Sample', 'spp', 'cluster', 'cells', 'NAME')

bd.ind <- match(colnames(bd.expr), bdiv.pheno$X)
colnames(bd.expr) <- bdiv.pheno$NAME[bd.ind]

## down-sample the data to make it more manageable
set.seed(100)
S.O.bd.filt <- subset(x = S.O.bd, downsample = 800)

# ## Differential gene expression
BD.markers <- FindAllMarkers(object = S.O.bd.filt, only.pos = TRUE, min.pct = 0)

BD.markers$GeneID <- gsub('-', '_', BD.markers$gene)
BD.markers.top <- BD.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.bd.filt, 
            features = BD.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

BD.markers.sig <- BD.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

ss <- BD.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
print(ss)
