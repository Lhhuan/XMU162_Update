library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/enrichment/figure/0/")
org<-read.table("compare_0_jaacard_index_fisher_test_histone_tfbs_two_side.txt",header = T,sep = "\t") %>% as.data.frame()
org$marker<-str_replace(org$marker,"CHROMATIN_Accessibility","CA")
org <-org%>%arrange(OR) #order to plot
org$FDR <-p.adjust(org$p_value,method="fdr")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
credplot.gg_main <- function(d){
  # org$marker<-str_replace(org$marker,"_"," ")
  p1<-ggplot(data=org, aes(x=marker,y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up)))
  p1<-p1+geom_pointrange(size = 0.001) 
  p1<-p1+coord_flip()# +  # flip coordinates (puts labels on y axis)
  p1<-p1+ylab("Odds ratio (log scale)")+xlab("Markers") 
  p1 <-p1+p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
  p1 <-p1+geom_text(aes(x=marker,y=1.1,label = significant)) 
}

# pdf("fisher_test_0_jaacard_index_forest_plot.pdf",height = 3.2,width = 4)
pdf("fisher_test_0_jaacard_index_forest_plot_log_histone_tfbs.pdf",height = 3.5,width = 3.8) #
org$marker <- factor(org$marker, levels=org$marker) #marker ranking
p<-credplot.gg_main(org)
p
dev.off()

# #---------------------log2
# credplot.gg_main <- function(d){
#   # org$marker<-str_replace(org$marker,"_"," ")
#   p1<-ggplot(data=org, aes(x=marker,y=log2(OR), ymin=log2(conf_int_bottom) , ymax=log2(conf_int_up)))
#   p1<-p1+geom_pointrange(size = 0.3) 
#   p1<-p1+coord_flip()# +  # flip coordinates (puts labels on y axis)
#   p1<-p1+ylab("Odds ratio (log2 scale)")+xlab("Markers") 
#   p1 <-p1+p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
#   p1 <-p1+geom_text(aes(x=marker,y=1.85,label = significant)) 
# }

# # pdf("fisher_test_0_jaacard_index_forest_plot.pdf",height = 3.2,width = 4)
# pdf("fisher_test_0_jaacard_index_forest_plot_log2_histone_tfbs.pdf",height = 3.5,width = 3.8) #
# org$marker <- factor(org$marker, levels=org$marker) #marker ranking
# p<-credplot.gg_main(org)
# p
# dev.off()
