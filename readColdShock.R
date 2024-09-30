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

processCountColdShock <- function(input.dir, filename, tt, rr, down.sample = T){
  file.counts <- read.csv(paste(input.dir, filename, sep = ''))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 5, min.features = 20)
  
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #VlnPlot(S.O, features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(name="Stopping distance", limits=c(1000, 1400))
  
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
  
  S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2] )
  
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
  
  
 return(S.O)
}



prod.desc <- read.csv('../../Input/coldShock/genes/BDiv_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

#saveRDS(prod.desc, '../Input/scClock/prod.desc.RData')

### b. divergense
#input.dir.bdiv <- "../Input/scRNAseqBdivCS/"
input.dir.bdiv <- "../Input/coldShock/counts/"
count.files <- list.files(input.dir.bdiv)
num.total.files <- length(count.files)

file.info <- data.frame(time = rep(NA, num.total.files), 
                        reactivate = rep(NA, num.total.files), 
                        filename = rep(NA, num.total.files))

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- count.files[i]
  file.info$time[i] <- gsub('OUT', '', gsub('bd', '', strsplit(file.info$filename[i], split = '\\.')[[1]][1]))
  file.info$reactivate[i] = ifelse(grepl('OUT', file.info$filename[i]), 'Y', 'N')
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCountColdShock(input.dir.bdiv, file.info$filename[i], file.info$time[i], file.info$reactivate[i])
  S.O.list <- c(S.O.list, list(L))
}

file.info$spp <- paste('BDiv', file.info$time, file.info$reactivate, sep = '')

print(file.info)


spps <- unlist(lapply(S.O.list, function(S.O) S.O@meta.data$spp[1]))
names(S.O.list) <- spps
  
#saveRDS(file.info, '../Input/coldShock/RData/file_info.RData')
#saveRDS(S.O.list, '../Input/coldShock/RData/S_O_list.RData')

saveRDS(file.info, '../../Input/coldShock/RData_new/file_info.RData')
saveRDS(S.O.list, '../../Input/coldShock/RData_new/S_O_list.RData')

num.objs <- length(spps)
genes <- lapply(1:num.objs, function(i){
  gene.row <- rownames(S.O.list[[i]]@assays$RNA@data)
} )

common.genes <- Reduce(intersect, genes)
length(common.genes)

