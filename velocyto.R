library(velocyto.R)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), phase = S.O@meta.data$phase, 
                          predicted.id = S.O@meta.data$predicted.id, 
                          lable.prob = S.O@meta.data$prediction.score.max,
                          spp = S.O@meta.data$spp)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

## loom files generated using velocyto
## On Markov:
## source .bashrc
## conda activate RNAvel
## velocyto run10x /H3/scRNAseqBabesCS/count/bd7dOUT_count/ /H3/scRNAseqBabesCS/Genome/PiroplasmaDB-55_Bdivergens1802A.gtf
## conda deactivate


input.dir.bdiv.loom <- "../Input/coldShock/loom/"
loom.files <- list.files(input.dir.bdiv.loom)
num.total.files <- length(loom.files)

file.info <- data.frame(time = rep(NA, num.total.files), 
                        reactivate = rep(NA, num.total.files), 
                        filename = rep(NA, num.total.files))

loom.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- loom.files[i]
  file.info$time[i] <- gsub('OUT', '', gsub('bd', '', strsplit(file.info$filename[i], split = '_')[[1]][1]))
  file.info$reactivate[i] = ifelse(grepl('OUT', file.info$filename[i]), 'Y', 'N')
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  ldat <- read.loom.matrices(paste(input.dir.bdiv.loom, loom.files[i], sep = ''))
  loom.list <- c(loom.list, list(ldat))
}

names(loom.list) <- paste("Bdiv", file.info$time, file.info$reactivate, sep = '')

loom.list <- lapply(loom.list, function(ldat){
  ldat <- lapply(ldat,function(x) {
    colnames(x) <-  gsub("x","\\.1",gsub(".*:","",colnames(x)))
    x
  })
})



S.Os <- readRDS('../Input/coldShock/RData_new/S_O_list_taxo_labs_not_anchored.RData')
S.O.list <- SplitObject(S.Os, split.by = 'spp')
S.O.list$BDiv0hrN@meta.data



prep.Loom <- function(ldat){
  bm <- as.Seurat(x = ldat)
  bm <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  bm <- SCTransform(object = bm, assay = "spliced")
  bm <- RunPCA(object = bm, verbose = FALSE)
  bm <- FindNeighbors(object = bm, dims = 1:20)
  bm <- FindClusters(object = bm)
  bm <- RunUMAP(object = bm, dims = 1:20)
  bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
  return(bm)
}


bms <- lapply(loom.list, function(ldat){
  bm <- prep.Loom(ldat)
  return(bm)
})


ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bms$Bdiv36hrN)))
names(x = ident.colors) <- levels(x = bms$Bdiv36hrN)
cell.colors <- ident.colors[Idents(object = bms$Bdiv36hrN)]
names(x = cell.colors) <- colnames(x = bms$Bdiv36hrN)


show.velocity.on.embedding.cor(emb = Embeddings(object = bms$Bdiv36hrN, reduction = "umap"), 
                               vel = Tool(object = bms$Bdiv36hrN, slot = "RunVelocity"), n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


emb <- pca.dat %>% select(PC_1, PC_2)
rownames(emb) <- pca.dat$Sample

cell.colors <- pca.dat %>% select(predicted.id)

cell.colors <- cell.colors %>% 
  mutate(color = case_when(predicted.id == "G1.a" ~ "firebrick",
                           predicted.id == "G1.b" ~ "darkorange2", 
                           predicted.id == 'S' ~ 'gold3', 
                           predicted.id == 'M' ~ 'darkolivegreen4', 
                           predicted.id == 'C' ~ 'darkorchid2'))
rownames(cell.colors) <- pca.dat$Sample
cc <- cell.colors$color
names(cc) <- rownames(cell.colors)


cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))
emb <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/embedding.rds"))


hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cc,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cc,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cc,min.max.cluster.average = 0.8)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))

