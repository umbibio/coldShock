library(ggplot2)
library(cowplot)
library(patchwork)
library(tidytext)

## individually processed samples, merged, but not anchored


S.O.list <- readRDS('../Input/coldShock/RData_new/S_O_loom_RNAvel_list.RData')
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

num.objs <- length(S.O.list)
spps <- names(S.O.list)
print(spps)



## First merge without anchoring. Use to generate UMAP
alldata <- merge(S.O.list[[1]], S.O.list[2:num.objs], add.cell.ids=spps)

alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000)
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 3000)
alldata <- ScaleData(alldata)
alldata <- RunPCA(alldata, features = VariableFeatures(object = alldata))
alldata <- FindNeighbors(alldata, dims = 1:30, reduction = 'pca')
alldata <- FindClusters(alldata, resolution = 0.2)
alldata <- RunUMAP(alldata, dims = 1:30)

saveRDS(alldata, '../Input/coldShock/RData_new/S_O_loom_RNAvel_list_merged_not_anchored.RData')
## Extracting the data for control over PCA directions
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}


pcaMataData.alldata <- getPcaMetaData(alldata)

lvs <- c("BDiv0hrN","BDiv4hrN", "BDiv12hrN", "BDiv36hrN", "BDiv7dN", "BDiv36hrY", "BDiv7dY" )
pcaMataData.alldata$spp <- factor(pcaMataData.alldata$spp, levels = lvs)

p  <- ggplot(pcaMataData.alldata, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = spp,
    color = spp
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  #scale_color_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                              'BDiv_human' = 'darkolivegreen4')) +
  #scale_fill_manual(values = c("BBig" = "firebrick","BBov" ="darkorchid3", 'BDiv_Cow' = 'darkslateblue', 
  #                             'BDiv_human' = 'darkolivegreen4')) +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
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

Idents(alldata) <- 'phase'
alldata@meta.data$spp <- factor(alldata@meta.data$spp, 
                                levels = lvs)
p <- DimPlot(alldata, reduction = "pca", 
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


## Now Integrate the objects, using 0hr as reference
print(file.info)
names(S.O.list)
ref.ind <- 1

all.samples.integrated <- processeMergedS.O(S.Os, data.ind = NA, ref.ind = ref.ind, res = 0.1, SC = FALSE)

all.samples.integrated@meta.data$spp <- factor(all.samples.integrated@meta.data$spp, 
                                               levels = lvs)
Idents(all.samples.integrated) <- "phase"
p <- DimPlot(all.samples.integrated, reduction = "pca", 
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

saveRDS(all.samples.integrated, '../Input/coldShock/RData_new/S_O_loom_RNAvel_list_merged_anchored.RData')

tmp <- all.samples.integrated@meta.data

tmp$phase[tmp$phase == 'G1.a'] <- 'G1'
tmp$phase[tmp$phase == 'G1.b'] <- 'G1'
tmp$phase[tmp$phase == 'S'] <- 'S/M'
tmp$phase[tmp$phase == 'M'] <- 'S/M'

stats <- tmp %>% group_by(spp) %>% mutate(total.cells = n()) %>% 
  ungroup() %>% group_by(spp, phase) %>% summarise(counts = n(), perc = n()/total.cells[1]) 

## manually add the missing phase
#stats <- rbind(stats, data.frame(spp = 'BDiv7dN', phase = 'M', counts = 0, perc = 0))

#stats$phase <- factor(stats$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
stats$phase <- factor(stats$phase, levels = c("G1",  "S/M", "C"))
#stats$phase <- factor(stats$phase, levels = c("G1",  "S", "M", "C"))
stats$spp <- factor(stats$spp, levels = lvs)

stats <- stats %>% arrange(spp, phase)

p <- ggplot(data=stats, aes(x=spp, y=perc, color = phase, group = phase)) +
  geom_point(size = 1.5) + geom_path(size = 1.2) + geom_vline(xintercept = 5, col='pink', lwd=1.2, linetype=2) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(face="bold", size=16, angle=45)) +
  theme(axis.text.y = element_text(face="bold", size=16, angle=0)) +
  scale_color_manual(values = c("G1" = "#ff6c67","S/M" ='#a2a700', 'C' = '#00c377')) +
  theme(
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  ) + theme(#legend.position = 'topright',
    legend.title = element_text(colour="black", size=16, 
                                face="bold"),
    legend.text = element_text(colour="black", size=16, 
                               face="bold"))

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
S.Os <- SplitObject(all.samples.integrated, split.by = 'spp')
lvs <- c("BDiv0hrN","BDiv4hrN", "BDiv12hrN", "BDiv36hrN", "BDiv7dN", "BDiv36hrY", "BDiv7dY" )
BDiv0hrN <- getPcaMetaData(S.Os$BDiv0hrN)
BDiv7dN <- getPcaMetaData(S.Os$BDiv7dN)
p1  <- ggplot(BDiv7dN, aes(x= PC_1,y= PC_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_density_2d(color= 'red', size = 1) + 
  theme_bw(base_size = 14) +
  #theme(legend.position = "none") +
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) + 
  theme(legend.position = c(0.1, 0.2),
        legend.title = element_text(colour="black", size=16, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=16, 
                                   face="bold"))



plot(p1)


S.O.BDiv0hrN <- S.O.list$BDiv0hrN
S.O.BDiv0hrN@meta.data$Phase <- S.O.BDiv0hrN@meta.data$phase

S.O.BDiv0hrN$Phase[S.O.BDiv0hrN$Phase == 'G1.a'] <- 'G1'
S.O.BDiv0hrN$Phase[S.O.BDiv0hrN$Phase == 'G1.b'] <- 'G1'
S.O.BDiv0hrN$Phase[S.O.BDiv0hrN$Phase == 'S'] <- 'S/M'
S.O.BDiv0hrN$Phase[S.O.BDiv0hrN$Phase == 'M'] <- 'S/M'

Idents(S.O.BDiv0hrN) <- 'Phase'
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

BDiv0hrN.stats <- stats %>% dplyr::filter(spp == 'BDiv0hrN')
BDiv7dN.stats <- stats %>% dplyr::filter(spp == 'BDiv7dN')

p3 <- ggplot(data=BDiv7dN.stats, aes(x=phase, y=perc, fill = phase)) +
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


### RNAvelocity plots on merged data

#S.Os <- SplitObject(all.samples.integrated, split.by = 'spp')
S.Os <- SplitObject(alldata, split.by = 'spp')
#Idents(S.Os$BDiv7dN) <- 'phase'
S.O <- RunVelocity(object = S.Os$BDiv0hrN, reduction = 'umap', deltaT = 1, kCells = 25, 
                   fit.quantile = 0.02,  ncores = 16)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = S.O)))
names(x = ident.colors) <- levels(x = S.O)
cell.colors <- ident.colors[Idents(object = S.O)]
names(x = cell.colors) <- colnames(x = S.O)
S.O <- AddMetaData(S.O, data.frame(samp = names(cell.colors), cell.colors = cell.colors))

p <- show.velocity.on.embedding.cor(emb = Embeddings(object = S.O, reduction = "umap"),
                                    vel = Tool(object = S.O, slot = "RunVelocity"), n = 200, scale = "log",
                                    cell.colors = ac(x = cell.colors, alpha = 0.5), nPcs = 30,
                                    cex = 0.8, arrow.scale = 1, show.grid.flow = TRUE,
                                    min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                                    do.par = FALSE, cell.border.alpha = 0.1, n.cores = 16)


DimPlot(S.O, reduction = 'umap', label = T) + NoLegend()
## Some plots
S.O.list[[1]]@meta.data$cell.colors
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = S.O.list[[1]])))
names(x = ident.colors) <- levels(x = S.O.list[[1]])
cell.colors <- ident.colors[Idents(object = S.O.list[[1]])]
names(x = cell.colors) <- colnames(x = S.O.list[[1]])

show.velocity.on.embedding.cor(emb = Embeddings(object = S.O.list[[6]], reduction = "umap")[,1:2],
                               vel = Tool(object = S.O.list[[6]], slot = "RunVelocity"), n = 300, scale = "sqrt",
                               cell.colors = ac(x = S.O.list[[6]]@meta.data$cell.colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                               do.par = FALSE, cell.border.alpha = 0.1)


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
