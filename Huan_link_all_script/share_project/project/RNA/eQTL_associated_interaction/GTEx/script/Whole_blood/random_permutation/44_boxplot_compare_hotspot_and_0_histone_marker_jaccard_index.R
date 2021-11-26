library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_cutoff_0.176_histone_marker_jaccard_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# random<-fread("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t")
random<-read.table("0_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t") %>% as.data.frame()


active_mark <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3")
repressed_mark <-c("H3K27me3","H3K9me3")

hotspot$Random_number <- 1
hotspot$Class <- "hotspot"
random$Class <- "random"

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                # axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

rs<- bind_rows(hotspot,random)

# figure_list <-list()
my_comparisons <- list(c("hotspot","random"))

large <-function(state){
    aa <-filter(rs,Marker==state)
    title_name<-state

    p <-ggplot(aa,aes(x=Class,y=jaacard_index,fill = Class))+ geom_violin(outlier.colour = NA)+ theme(legend.position ="none") +ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+ylab("Jaccard index")+p_theme+coord_cartesian(ylim = c(0, 1.1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
}
 
plist = lapply(active_mark,large)

pdf("./figure/0/compare_0_active_mark_jaccard_index.pdf",width=6.3, height=4.5)
CombinePlots(plist,ncol=3,nrow=2)
dev.off()

plist = lapply(repressed_mark,large)

pdf("./figure/0/compare_0_repressed_mark_jaccard_index.pdf",width=4, height=2.25)
CombinePlots(plist,ncol=2)
dev.off()
print("aaa")


#----
