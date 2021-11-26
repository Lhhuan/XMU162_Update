library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/figure/0_0.176/")
org<-read.table("compare_0_0.176_jaacard_index_fisher_test_two_side.txt",header = T,sep = "\t") %>% as.data.frame()

org <-org%>%arrange(OR) #order to plot
org$FDR <-p.adjust(org$p_value,method="fdr")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
credplot.gg_main <- function(d){
p1<-ggplot(data=org, aes(x=marker,y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up)))
  p1<-p1+geom_pointrange(size = 0.3) 
  p1<-p1+coord_flip()# +  # flip coordinates (puts labels on y axis)
  p1<-p1+ylab("Odds ratio (log scale)")+xlab("Markers") 
  p1 <-p1+p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
  p1 <-p1+geom_text(aes(x=marker,y=1.1,label = significant)) 
}

# pdf("fisher_test_0_0.176_jaacard_index_forest_plot.pdf",height = 3.2,width = 4)
pdf("fisher_test_0_0.176_jaacard_index_forest_plot_log.pdf",height = 3.5,width = 4) #
org$marker <- factor(org$marker, levels=org$marker) #marker ranking
p<-credplot.gg_main(org)
p
dev.off()
