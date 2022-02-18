library(ggplot2)
library(cowplot)
library(patchwork)
library(tidytext)
library(ggVennDiagram)
library(UpSetR)


makeGlobalContrasts <- function(all.samples.integrated, c.mode  = c('full', 'WT')){
  
  c.mode <- match.arg(c.mode)
  
  objs <- as.character(unique(all.samples.integrated@meta.data$spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.time = as.numeric(gsub("hr.*", "", gsub("7d", "168hr", gsub("BDiv", "", ref)))),
           query.time = as.numeric(gsub("hr.*", "", gsub("7d", "168hr", gsub("BDiv", "", query)))),
           ref.reactive = ifelse(grepl('(R)', ref), 'Y', 'N'),
           query.reactive = ifelse(grepl('(R)', query), 'Y', 'N'))
  if(c.mode == 'full'){
    my.contrasts <- contrasts %>%
      dplyr::filter( 
        ( (ref.time < query.time) ) | 
          ( (ref.time == query.time) & (ref.reactive != query.reactive) ) 
      ) %>% dplyr::filter(!(ref.reactive == 'N' & query.reactive == 'Y' & ref.time != 0)) %>%
      dplyr::filter(!(ref.reactive == 'Y' & query.reactive == 'Y')) %>% 
      dplyr::filter(!(ref.reactive == 'Y' & query.reactive == 'N' & ref.time != query.time ))
  }else if(c.mode == 'WT'){
    my.contrasts <- contrasts %>% dplyr::filter(ref.time == 0 & query.reactive != 'Y' & ref.time < query.time)
  }
  
  return(my.contrasts)
  
}

makeMatchedContrasts <- function(all.samples.integrated, c.mode  = c('full', 'WT')){
  
  c.mode <- match.arg(c.mode)
  tmp <- makeGlobalContrasts(all.samples.integrated, c.mode  = c.mode)
 
  clusters <- unique(all.samples.integrated@meta.data$phase)
 
  
  my.contrasts <- data.frame(ref = paste(rep(tmp$ref, each = length(clusters)), rep(sort(clusters), nrow(tmp)), sep = '_'),
                             query = paste(rep(tmp$query, each = length(clusters)), rep(sort(clusters), nrow(tmp)), sep = '_'),
                             cluster = rep(sort(clusters), length(tmp$ref))) 
  
  return(my.contrasts)
  
}

getCellCyclePhaseMarkers <- function(all.spp.list){
  all.markers.list <- mclapply(all.spp.list, function(x) FindAllMarkers(object = x, only.pos = TRUE))
  all.markers.list <- lapply(all.markers.list, function(x) {
    x$GeneID = gsub('-', '_', x$gene)
    x$glob.clust <- gsub('.*_', '', x$cluster)
    return(x)
  })
  
  
  all.markers.list.sig <- lapply(all.markers.list, function(x) {
    sig.marker = x %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
    return(sig.marker)
  })
  
  L <- list(all.markers.list = all.markers.list, 
            all.markers.list.sig = all.markers.list.sig)
  return(L)
}

## Integrate Samples
alldata.toxo.labs.integrated <- readRDS('../Input/coldShock/RData_new/S_O_list_taxo_labs_anchored.RData')



# Differential Expression
# Identify global markers independent of cell cycle phase. Are there any global Cold Shock regulators?
# Two class comparison: fc = ident.1/ident.2

DefaultAssay(alldata.toxo.labs.integrated) <- "RNA"
Idents(alldata.toxo.labs.integrated) <- "spp"

c.mode <- 'WT'
contrasts <- makeGlobalContrasts(alldata.toxo.labs.integrated, 'WT')


global.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(alldata.toxo.labs.integrated, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})

global.shock.markers <- bind_rows(global.DEGs)
saveRDS(global.shock.markers, '../Input/coldShock/RData_new/global_shock_markers.RData')

global.shock.markers.sig <- global.shock.markers %>% 
  dplyr::filter(p_val < 0.01 & abs(avg_log2FC) > log2(1.5)) %>% arrange(desc(abs(avg_log2FC)))

global.shock.markers.overlap <- global.shock.markers.sig %>% group_by(GeneID) %>% 
  summarise(num.contrasts = length(query)) %>% arrange(desc(num.contrasts))

global.shock.markers.sig <- left_join(global.shock.markers.sig, global.shock.markers.overlap, by = 'GeneID')
global.shock.markers.stats <- global.shock.markers.sig %>%
  mutate(reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(ref, query, reg) %>% summarise(num.DEGs = n()) 

print(global.shock.markers.stats)

global.shock.markers.sig <- left_join(global.shock.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(global.shock.markers.sig, '../Output/coldShock/tables/global_shock_markers_WT_base.xlsx')




## plot
top.up.markers <- global.shock.markers.sig %>% slice_max(order_by = avg_log2FC, n= 20)
Idents(alldata.toxo.labs.integrated) <- 'spp'
DefaultAssay(alldata.toxo.labs.integrated) <- 'RNA'

VlnPlot(object = alldata.toxo.labs.integrated, features = 'Bdiv-011450c')


### Differential expression analysis
## Cell cycle phase specific


all.spp.list <- SplitObject(alldata.toxo.labs.integrated, split.by = "spp")

## Re-normalize splitted data
for (i in 1:length(all.spp.list)) {
  Idents(all.spp.list[[i]]) <- 'phase'
  all.spp.list[[i]] <- NormalizeData(all.spp.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}


contrasts <- makeMatchedContrasts(alldata.toxo.labs.integrated, c.mode  = 'WT')

GM <- getCellCyclePhaseMarkers(all.spp.list)

## For each cluster, pool global markers.
pooled_markers <- bind_rows(GM$all.markers.list.sig) %>% group_by(glob.clust) %>%
  summarise(glob.markers = list(unique(gene)), num.markers = length(unique(gene)))

#ident.1 case, ident.2 is control

Idents(alldata.toxo.labs.integrated) <- 'phase.cond'
DefaultAssay(alldata.toxo.labs.integrated) <- 'RNA'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(alldata.toxo.labs.integrated, ident.1 = x$query, ident.2 = x$ref, verbose = T)
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
saveRDS(matched.DEGs, '../Input/coldShock/RData_new/matched_DEGs.Rdata')

matched.DEGs.sig <- matched.DEGs %>% dplyr::filter(p_val < 0.01 & abs(avg_log2FC) > log2(1.5)) %>% 
  arrange(desc(abs(avg_log2FC)))

matched.DEGs.overlap <- matched.DEGs.sig %>% group_by(GeneID) %>% 
  summarise(num.contrasts = length(query)) %>% arrange(desc(num.contrasts))

matched.DEGs.sig <- left_join(matched.DEGs.sig, matched.DEGs.overlap, by = 'GeneID')

matched.DEGs.sig <- left_join(matched.DEGs.sig, prod.desc, by = 'GeneID')
write.xlsx(matched.DEGs.sig, '../Output/coldShock/tables/matched_shock_markers.xlsx')


matched.DEGs.stats <- matched.DEGs.sig %>%
  mutate(reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(query, reg) %>% summarise(num.DEGs = n()) 

print(matched.DEGs.stats)

matched.DEGs.stats$phase <- factor(gsub('.*_', '', matched.DEGs.stats$query),
                                      levels = c('G1', 'S/M', 'C'))


matched.DEGs.top <- matched.DEGs.sig %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 1)



p2 <- FeaturePlot(object = alldata.toxo.labs.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = matched.DEGs.top$gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)


tmp.down <- matched.DEGs %>% slice_min(order_by = avg_log2FC, n = 20)
tmp.up <- matched.DEGs %>% slice_max(order_by = avg_log2FC, n = 20)

my.features <- c(tmp.up$gene[1:2], 'Bdiv-031680c')


p1 <- VlnPlot(object = alldata.toxo.labs.integrated, features = my.features[2]) + 
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


tmp <- left_join(matched.DEGs.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/coldShock/tables/matched_shock_markers.xlsx')


tmp <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% group_by(ref, query) %>%
  summarise(genes = list(GeneID))

venn.list <- tmp$genes
names(venn.list) <- paste(tmp$ref, "_vs_", tmp$query, sep = '')

ggVennDiagram(venn.list)
upset(fromList(venn.list))

