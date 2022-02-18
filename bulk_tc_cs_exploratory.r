library(tidytext)
library(jcolors)
library(openxlsx)
library(ggplot2)
require(gridExtra)
library(grid)
library(tidyverse)
library(tidytext)
library(ggrepel)



final.CS.tab <- read.xlsx("../Input/coldShock/bulk_tc/babesia_ColdSchock_new_logCPM_logFC_qval_edgeR_V4.xlsx")

sync.vs.bulk.fc <- final.CS.tab %>% 
  select(GeneName, Sync4C.2.vs.tc.6_log_fc, Sync4C.4.vs.tc.8_log_fc, Sync4C.6.vs.tc.10_log_fc)

sync.vs.bulk.fc.lng <-sync.vs.bulk.fc %>%
  pivot_longer(-GeneName, names_to = "contrast", values_to = "log_fc")
sync.vs.bulk.fc.lng$contrast <- gsub("_log_fc", "", sync.vs.bulk.fc.lng$contrast)


sync.vs.bulk.qval <- final.CS.tab %>% select(GeneName, Sync4C.2.vs.tc.6_qval,  Sync4C.4.vs.tc.8_qval, Sync4C.6.vs.tc.10_qval)
sync.vs.bulk.qval.lng <- sync.vs.bulk.qval %>% pivot_longer(-GeneName, names_to = "contrast", values_to = "qval")
sync.vs.bulk.qval.lng$contrast <- gsub("_qval", "", sync.vs.bulk.qval.lng$contrast)

sync.vs.bulk.df <- inner_join(sync.vs.bulk.fc.lng, sync.vs.bulk.qval.lng, by = c("GeneName", "contrast"))

sync.vs.bulk.df <- sync.vs.bulk.df %>% mutate(sig = ifelse(abs(log_fc) > 0.58 & qval < 0.05, "yes", "no")) 

prod.desc <- read.csv('../Input/coldShock/genes/BDiv_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, ProductDescriptionBdiv = Product.Description) %>% distinct()
sync.vs.bulk.prod.desc <- inner_join(sync.vs.bulk.df, prod.desc, by = c("GeneName" = "GeneID"))

#write.xlsx(sync.vs.bulk.prod.desc, "../output/ColdShock/cold_vs_hot.xlsx")

toxo.bdiv.orth <- read.xlsx('../Input/coldShock/Orthologs/GT1_BDiv.xlsx')
colnames(toxo.bdiv.orth) <- c('GeneID_Toxo', 'GeneID', 'ProductDescriptionToxo')

######### volcano plot #########

volc.df <- sync.vs.bulk.df 
volc.df$contrast <- gsub("Sync4C.2.vs.tc.6", "2hr", 
                         gsub("Sync4C.4.vs.tc.8", "4hr", gsub("Sync4C.6.vs.tc.10", "6hr", volc.df$contrast)))
volc.df$geneLable <- ''
show.lab <- c('Bdiv_000320c', 'Bdiv_023990c', 'Bdiv_036750c', 'Bdiv_031010c')
# volc.df <- volc.df %>% group_by(contrast) %>% 
#   mutate(GeneLable = ifelse(-log(qval) > 20  & abs(log_fc) > 1.7 , "True", "False"))
volc.df <- volc.df %>% group_by(contrast) %>% 
  mutate(GeneLable = ifelse(GeneName %in%  show.lab, "True", "False"))

volc.df$Name <- gsub("Bdiv_", "", volc.df$GeneName)
volc.df$Name <- volc.df$GeneName

# remooving some random genes to make some space and (nicer plot) 
# more genes to be removed? 
volc.df$Name <- gsub("009180", "", 
                     gsub("021910", "",
                          gsub("013880c", "", 
                               gsub("006490c", "", 
                                    gsub("010240c", "", volc.df$Name)))))

p <- ggplot(volc.df) + geom_point(aes(log_fc, -log(qval),color=sig))+ 
  geom_text_repel(aes(log_fc, -log(qval)), size = 5, fontface = "bold",
                  label = ifelse(volc.df$GeneLable == "True", as.character(volc.df$Name),""), 
                  box.padding = unit(0.45, "lines"),
                  hjust=1,
                  segment.angle = 180,
                  nudge_x = -1, 
                  nudge_y = 1,
                  segment.size = 0.4) + 
  theme(legend.title=element_blank(),text = element_text(size=20))+ 
  scale_color_manual(values = c("yes" = "red", "no" = "black"))+
  facet_grid(. ~ contrast,  scales = "fixed",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) + 
  theme_bw()+
  theme(
    #plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    #legend.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18, angle = 90, vjust = 0.45, color = "black"),  
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black", face = "bold"),
    strip.text = element_text(color = "black", size = 18, face = "bold"),
    strip.background = element_rect(fill = "white"))+
  theme(legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  xlab("log2 FC")

p


ggsave(filename="../Output/coldShock/figs/volcano.pdf",
       plot=p,
       width = 10, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

toxo.bdiv.orth <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/compScBabesia/Orthologs/GT1_BDiv.xlsx')
colnames(toxo.bdiv.orth) <- c('GeneID_Toxo', 'GeneID', 'ProductDescriptionToxo')

sync.vs.bulk.filt <-  sync.vs.bulk.df %>%  
  filter(abs(log_fc) > 0.58 & qval < 0.05) %>% select(GeneName, contrast, log_fc)
sync.vs.bulk.filt <- sync.vs.bulk.filt %>% mutate(direction = ifelse(log_fc > 0 , "up", "down"))
sync.vs.bulk.filt$contrast <- gsub("Sync4C.2.vs.tc.6", "2hr", 
                                   gsub("Sync4C.4.vs.tc.8", "4hr", 
                                        gsub("Sync4C.6.vs.tc.10", "6hr", sync.vs.bulk.filt$contrast)))


sync.vs.bulk.filt.orth <- inner_join(sync.vs.bulk.filt, toxo.bdiv.orth, by = c("GeneName" = "GeneID"))
sync.vs.bulk.filt.orth.stat <- sync.vs.bulk.filt.orth %>% 
  group_by(contrast,  direction) %>% summarise(genes = list(GeneID_Toxo), num.degs = n()) 

#write.xlsx(sync.vs.bulk.filt.orth.stat, "../output/ColdShock/CS_vs_TC_DEGs_orth_stat.xlsx")

## GO trm 
in.dir <- '../output/ColdShock/GO/'
all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  GF <- strsplit(nn, split = '_')[[1]][1]
  phase <- strsplit(nn, split = '_')[[1]][2]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$phase <- phase
  all.clust.items <- c(all.clust.items, list(tmp))
}

all.clust.items <- do.call(rbind, all.clust.items)
names(all.clust.items) <- gsub(" ", ".", gsub("phase", "Sample", names(all.clust.items)))

#saveRDS(all.clust.items, '../Input/compScBabesia/RData/GO_Conserved_Markers_TGGT1.RData')
filtered.Go <- all.clust.items %>% arrange(Sample, Benjamini) %>% distinct() %>%
  group_by(Sample) %>% mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank < 40 & Bgd.count > 10) %>% 
  arrange(Sample, Benjamini) 

names(filtered.Go)


GO <- filtered.Go %>% 
  mutate(Result.gene.list = strsplit(as.character(Result.gene.list), ",")) %>% 
  unnest(Result.gene.list)
toxo.bdiv.orth <- read.xlsx('~/work/ToxoPlasmaGondiiR/Input/compScBabesia/Orthologs/GT1_BDiv.xlsx')
GO <- inner_join(toxo.bdiv.orth, GO , by = c("GT1" = "Result.gene.list"))
GO.df <- GO %>% group_by(ID, Name, GF, Sample) %>% summarise(Bdiv_ID = list(BDiv), Toxo_ID = list(GT1), total = n())
GO.final <- left_join(GO.df , filtered.Go,  by = c("ID", "Name", "GF", "Sample"))

write.xlsx(GO.final, "../output/ColdShock/filtered_GO_term.xlsx")

############# bar plot ##########
sync.vs.tc.stat  <- sync.vs.bulk.filt %>% mutate(direction = ifelse(log_fc > 0 , "up", "down"))

stats.deg <- sync.vs.tc.stat %>% group_by(contrast, direction) %>% 
  summarise(Genes = list(GeneName), total = n()) %>% group_by(contrast) %>% 
  mutate(ylab_pos = total - 2)

stats.deg$contrast <- gsub("Sync4C.2.vs.tc.6", "2hr", 
                           gsub("Sync4C.4.vs.tc.8", "4hr", gsub("Sync4C.6.vs.tc.10", "6hr", stats.deg$contrast)))

stats.deg$contrast <- factor(stats.deg$contrast, levels = c('2hr', '4hr', '6hr'))

stats.deg$direction <- factor(stats.deg$direction, levels = c('up', 'down'))
stats.deg$label <- stats.deg$total


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

ggsave(filename="../output/ColdShock/bar_plots_cs.png",
       plot=p,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


sync.vs.tc  <- sync.vs.bulk.filt %>% mutate(direction = ifelse(log_fc > 0 , "up", "down")) 
sync.vs.tc$Category <- paste( sync.vs.tc$direction , sync.vs.tc$contrast, sep = "_")

## up eg
sync.vs.tc.up <- sync.vs.tc %>% filter(direction == "up") 
sync.vs.tc.up <- split(sync.vs.tc.up , f = sync.vs.tc.up$Category)


sync.vs.tc.up.list <- list()
for(i in 1:length(sync.vs.tc.up)){
  categ <- names(sync.vs.tc.up)[i] 
  genes <- sync.vs.tc.up[[categ]]$GeneName
  sync.vs.tc.up.list <- c(sync.vs.tc.up.list, list(genes))
}
names(sync.vs.tc.up.list) <- names(sync.vs.tc.up)

up.comm.genes  <- Reduce(intersect, sync.vs.tc.up.list)

p.up <- ggVennDiagram(sync.vs.tc.up.list, 
                      category.names = c("2 hr"," 4 hr","6 hr"),
                      set_size = 8, label_size = 6,
                      lwd = 0.8, lty = 1, label_alpha = 0.7, edge_size = 1) +
  #scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") 
  ggtitle("up-regulated") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 14, face = 'bold'), 
        legend.title = element_text(size = 16,  face = "bold"))


p.up


# down reg
sync.vs.tc.down <- sync.vs.tc %>% filter(direction == "down") 
sync.vs.tc.down <- split(sync.vs.tc.down , f = sync.vs.tc.down$Category)



sync.vs.tc.down.list <- list()

for(i in 1:length(sync.vs.tc.down)){
  
  categ <- names(sync.vs.tc.down)[i] 
  genes <- sync.vs.tc.down[[categ]]$GeneName
  
  sync.vs.tc.down.list <- c(sync.vs.tc.down.list, list(genes))
  
}
names(sync.vs.tc.down.list) <- names(sync.vs.tc.down)


p.down <- ggVennDiagram(sync.vs.tc.down.list, 
                        category.names = c("2 hr"," 4 hr","6 hr"),
                        set_size = 8,label_size = 6,
                        lwd = 0.8, lty = 1, label_alpha = 0.7, edge_size = 1) +
  #scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") 
  ggtitle("down-regulated") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 14, face = 'bold'), 
        legend.title = element_text(size = 16,  face = "bold"))

p.down

p <- grid.arrange(p.up, p.down, ncol = 2)

ggsave(filename="../output/ColdShock/figures/venn_up_down.png",
       plot=p,
       width = 12, height = 12,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)
