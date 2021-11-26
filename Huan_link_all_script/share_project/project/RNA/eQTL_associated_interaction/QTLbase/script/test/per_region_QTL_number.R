library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

org<-read.table("../output/count_per_500bp_QTL_number_region.txt",header = T,sep = "\t") %>% as.data.frame()
org2<-read.table("../output/count_per_750bp_QTL_number_region.txt",header = T,sep = "\t") %>% as.data.frame()
org3<-read.table("../output/count_per_1000bp_QTL_number_region.txt",header = T,sep = "\t") %>% as.data.frame()

#--------------------------------------------
#-----------------------------------------------------------------------------------500bp
figure_list <-list ()
#----------------------------------------------------------

#-------------------------------存到一起

for(i in c(1:22)){
    chr_specific<-filter(org,chr == i)
    title_name<-paste(i)
    figure_list[[i]] <- print(ggplot(chr_specific,aes(x = region, y = number)) + 
        geom_bar(stat = 'identity', fill = "black")+
        xlab("position region") +
        ylab("QTL number") +
        # p.theme+
        ggtitle(title_name)+theme(plot.title=element_text(size=9,hjust = 0.5)))
}
pdf("chr_specific_bar_plot_distance_of_locus_region.pdf",height = 24,width = 16) 
marrangeGrob(figure_list,nrow=6,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# marrangeGrob(figure_list[1:13],nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
dev.off()
#--------------------------------------------------------
