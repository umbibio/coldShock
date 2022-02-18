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
library(sctransform)
library(cowplot)
library(patchwork)
library(doParallel)


source("./util_funcs.R")

set.seed(100)

Bdiv0hrN.cutoff  <- c(200, 1250)
Bdiv12hrN.cutoff <- c(170, 1000)
Bdiv36hrN.cutoff <- c(150, 900)
Bdiv36hrY.cutoff <- c(220, 1250)
Bdiv4hrN.cutoff  <- c(110, 900)
Bdiv7dN.cutoff   <- c(75, 650)
Bdiv7dY.cutoff   <- c(220, 1250)

processLoomColdShock <- function(input.dir, filename, tt, rr, S.O.toxo, down.sample = T){
  ldat <- read.loom.matrices(paste(input.dir, filename, sep = ''))
  ldat <- lapply(ldat,function(x) {
    colnames(x) <-  gsub("x","\\.1",gsub(".*:","",colnames(x)))
    x
  })
  
  S.O <- as.Seurat(x = ldat)
  #VlnPlot(S.O, features = c("nFeature_spliced", "nCount_spliced"),ncol = 2)
  #FeatureScatter(S.O, feature1 = "nCount_spliced", feature2 = "nFeature_RNA")
  #VlnPlot(S.O, features = c("nFeature_spliced"), ncol = 1) + scale_y_continuous(name="Stopping distance", limits=c(1000, 1400))
  
  #cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  #print(cutoffs)
  
  if(tt == '0hr'){
    cutoff = Bdiv0hrN.cutoff
  }else if(tt == '4hr'){
    cutoff = Bdiv4hrN.cutoff
  }else if(tt == '12hr'){
    cutoff = Bdiv12hrN.cutoff
  }else if(tt == '36hr' & rr == 'N'){
    cutoff = Bdiv36hrN.cutoff
  }else if(tt == '36hr' & rr == 'Y'){
    cutoff = Bdiv36hrY.cutoff
  }else if(tt == '7d' & rr == 'N'){
    cutoff = Bdiv7dN.cutoff
  }else if(tt == '7d' & rr == 'Y'){
    cutoff = Bdiv7dY.cutoff
  }
  
  S.O <- subset(S.O, subset = nFeature_spliced > cutoff[1] & nFeature_spliced < cutoff[2] )
  
  if(down.sample){
    set.seed(100)
    S.O <- subset(x = S.O, downsample = 4000)
  }
  
  S.O <- prep_S.O(S.O)
  
  
  cluster <- as.character(S.O$seurat_clusters)
  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  spp <- paste('BDiv', tt, rr, sep = '')
  pheno$spp <- spp
  pheno$time <- tt
  pheno$reactivate <- rr
  pheno$cluster <- cluster
  pheno$cell <- paste(pheno$spp, pheno$cluster , sep = "")
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
  rownames(pheno) <- pheno$Sample
  
  S.O <- AddMetaData(S.O, metadata = pheno)
  
  ## Transfer lables from Toxo
  anchors <- FindTransferAnchors(reference = S.O.toxo, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.toxo@meta.data$phase, dims = 1:30)
  predictions$phase <- factor(predictions$predicted.id, 
                              levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  
  Idents(S.O) <- 'phase'
  S.O <- RunVelocity(object = S.O, deltaT = 1, kCells = 25, fit.quantile = 0.02)
  
  ident.colors <- (scales::hue_pal())(n = length(x = levels(x = S.O)))
  names(x = ident.colors) <- levels(x = S.O)
  cell.colors <- ident.colors[Idents(object = S.O)]
  names(x = cell.colors) <- colnames(x = S.O)
  S.O <- AddMetaData(S.O, data.frame(samp = names(cell.colors), cell.colors = cell.colors))
  
  p <- show.velocity.on.embedding.cor(emb = Embeddings(object = S.O, reduction = "umap"),
                                 vel = Tool(object = S.O, slot = "RunVelocity"), n = 300, scale = "sqrt",
                                 cell.colors = ac(x = cell.colors, alpha = 0.5),
                                 cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE,
                                 min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                                 do.par = FALSE, cell.border.alpha = 0.1)
  
  
  return(S.O)
}


## Read Toxo Data for label transfer
S.O.tg    <- readRDS('../Input/compScBdTgPb/RData/S.O.tg.RData')
GT1.bdiv  <- read.xlsx("../Input/compScBdTgPb/Orthologs/rec_GT1.vs.B.divergence.xlsx")
GT1.BD <- GT1.bdiv %>% dplyr::select('query_id', contains('subject_id'))
colnames(GT1.BD) <- c('TGGT1', 'BDiv')
tg_meta.data <- S.O.tg@meta.data
tg.data <- S.O.tg@assays$RNA@counts
bdiv.id <- GT1.BD$BDiv[match(gsub('-','_',rownames(tg.data)), GT1.BD$TGGT1)]
orth.ind <- !is.na(bdiv.id)
tg.data <- tg.data[orth.ind, ]
rownames(tg.data) <- bdiv.id[orth.ind]

S.O.toxo <- CreateSeuratObject(counts = tg.data)
S.O.toxo$orig.ident <- 'Toxo'
S.O.toxo$spp <- 'Toxo'
S.O.toxo <- AddMetaData(S.O.toxo, tg_meta.data)
S.O.toxo <- prep_S.O(S.O.toxo)

prod.desc <- read.csv('../Input/coldShock/genes/BDiv_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)


### b. divergense
input.dir  <- "../Input/coldShock/loom/"
loom.files <- list.files(input.dir)
num.total.files <- length(loom.files)

file.info <- data.frame(time = rep(NA, num.total.files), 
                        reactivate = rep(NA, num.total.files), 
                        filename = rep(NA, num.total.files))

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- loom.files[i]
  file.info$time[i] <- gsub('OUT', '', gsub('bd', '', strsplit(file.info$filename[i], split = '_')[[1]][1]))
  file.info$reactivate[i] = ifelse(grepl('OUT', file.info$filename[i]), 'Y', 'N')
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processLoomColdShock(input.dir, file.info$filename[i], file.info$time[i], file.info$reactivate[i], S.O.toxo)
  S.O.list <- c(S.O.list, list(L))
}

file.info$spp <- paste('BDiv', file.info$time, file.info$reactivate, sep = '')

print(file.info)


spps <- unlist(lapply(S.O.list, function(S.O) S.O@meta.data$spp[1]))
names(S.O.list) <- spps

#saveRDS(file.info, '../Input/coldShock/RData/file_info.RData')
#saveRDS(S.O.list, '../Input/coldShock/RData/S_O_list.RData')

saveRDS(file.info, '../Input/coldShock/RData_new/file_info_loom.RData')
saveRDS(S.O.list, '../Input/coldShock/RData_new/S_O_loom_RNAvel_list.RData')


## Some plots
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = S.O.list[[1]])))
names(x = ident.colors) <- levels(x = S.O.list[[1]])
cell.colors <- ident.colors[Idents(object = S.O.list[[1]])]
names(x = cell.colors) <- colnames(x = S.O.list[[1]])

show.velocity.on.embedding.cor(emb = Embeddings(object = S.O.list[[6]], reduction = "pca")[,1:2],
                               vel = Tool(object = S.O.list[[6]], slot = "RunVelocity"), n = 300, scale = "sqrt",
                               cell.colors = ac(x = cell.colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                               do.par = FALSE, cell.border.alpha = 0.1)

# emat <- S.O.list[[1]]@assays$spliced@data
# nmat <- S.O.list[[1]]@assays$unspliced@data
# rvel.cd.pooled <- gene.relative.velocity.estimates(emat, nmat,
#                                                    fit.quantile = 0.05,
#                                                    min.nmat.emat.correlation = 0.2, 
#                                                    min.nmat.emat.slope = 0.2, 
#                                                    kCells = 200)
# pca.velocity.plot(rvel.cd.pooled,
#                   nPcs=2,
#                   plot.cols=1,
#                   cell.colors=cell.colors,
#                   pc.multipliers=c(1,-1), ## adjust as needed to orient pcs
#                   show.grid.flow = TRUE, 
#                   grid.n=40 ## adjust as needed
# )


# 
# genes <- lapply(1:num.objs, function(i){
#   gene.row <- rownames(S.Os[[i]]@assays$RNA@data)
# } )
# 
# common.genes <- Reduce(intersect, genes)
# 
# 
